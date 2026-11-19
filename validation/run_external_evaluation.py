#!/usr/bin/env python3
"""External evaluation of Locus-Recon on public isolates with a complete genome.

This script carries out section 2 of ``validation/external-evaluation/PROTOCOL.md``,
which was written before any of these analyses was run.  For each isolate of
``validation/external-evaluation/isolates.tsv`` and each of eight loci it

1. downloads the Illumina run from the European Nucleotide Archive and
   subsamples it to 100x of the genome length with a fixed seed;
2. assembles a draft with SPAdes in isolate mode;
3. takes one bait per locus from the RefSeq reference genome of the species
   and defines the truth by aligning that bait to the complete genome of the
   isolate;
4. runs Locus-Recon with default options, then ``locus-recon-depth-ratio``
   and ``locus-recon-copy-number`` on the outputs of that run;
5. runs ARIBA and KMA with default options on the same reads and baits, with
   the functions of ``run_comparator_benchmark.py``;
6. scores every unit of isolate and locus against the truth with the scoring
   code of ``run_comparator_benchmark.py``, and applies the decision rules
   that the protocol fixes for KMA and ARIBA.

The work is split in stages so that each can run in the environment that
holds its tools, and so that an interrupted run resumes where it stopped::

    # environment of Locus-Recon (environment.yml) plus seqtk
    python validation/run_external_evaluation.py --workdir W --stage prepare --seqtk /path/to/seqtk
    python validation/run_external_evaluation.py --workdir W --stage assemble
    python validation/run_external_evaluation.py --workdir W --stage reconstruct
    python validation/run_external_evaluation.py --workdir W --stage copy-number
    python validation/run_external_evaluation.py --workdir W --stage block-size
    # environment with ARIBA, KMA and MUMmer
    python validation/run_external_evaluation.py --workdir W --stage comparators
    python validation/run_external_evaluation.py --workdir W --stage concordance
    # environment of Locus-Recon again
    python validation/run_external_evaluation.py --workdir W --stage score

Nothing in ``locus_recon/`` is changed or configured here.  The only options
passed to Locus-Recon are the inputs, the number of threads and
``--noncoding-locus`` for the rRNA genes.

Choices that the protocol leaves to the implementation are fixed in this file
and were fixed before any tool was run on the reads.  They are listed in the
``README.md`` of the output directory.
"""

import argparse
import csv
import gzip
import hashlib
import http.client
import json
import math
import platform
import random
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "validation"))

from locus_recon.io import reverse_complement  # noqa: E402
from validation.common import (  # noqa: E402
    build_provenance, command_version, sha256_file, source_sha256, write_json,
)
import run_comparator_benchmark as bench  # noqa: E402
import run_tier_calibration as calibration  # noqa: E402

VALIDATION = REPO_ROOT / "validation"
EVALUATION_DIR = VALIDATION / "external-evaluation"
FROZEN_VERSION = "1.0.1"
FROZEN_DIGEST = "2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0"

SEED = 20261118                 # subsampling, and resampling of isolates
TARGET_DEPTH = 100              # x of the genome length
LOCI = ("gyrB", "rpoB", "recA", "dnaK", "atpD", "23S", "16S", "tuf")
RRNA_PRODUCT = {"23S": "23S ribosomal RNA", "16S": "16S ribosomal RNA"}
TRUTH_MIN_IDENTITY = 80.0       # protocol 2.3
TRUTH_MIN_BAIT_FRACTION = 0.90  # protocol 2.3
BLOCK_SIZES = (100, 250, 500, 1000)
BOOTSTRAP_RESAMPLES = 10000     # resampling of isolates for intervals

# Protocol 2.2 finds a protein-coding bait by its gene symbol and records the
# unit as not available when the reference annotation lacks it.  Units found
# only by the two rules below are NOT part of the protocol set.  They are run
# as a supplementary set, labelled as such in every table, and never pooled
# with the protocol set without a label.  The rules were written after the
# reference annotations were read and before any tool was run on the reads.
FALLBACK_PARALOGUE_SUFFIXES = ("A", "B")          # e.g. tufA, tufB for tuf
FALLBACK_PRODUCT = {
    "gyrB": ("DNA gyrase subunit B",),
    "rpoB": ("DNA-directed RNA polymerase subunit beta",),
    "recA": ("recombinase RecA", "recombinase A"),
    "dnaK": ("molecular chaperone DnaK",),
    "atpD": ("F0F1 ATP synthase subunit beta",),
    "tuf": ("elongation factor Tu",),
}
RULE_PROTOCOL = "protocol"
SET_PROTOCOL, SET_SUPPLEMENTARY, SET_UNAVAILABLE = (
    "protocol", "supplementary (bait by fallback rule)", "not available")

ENA_FILEREPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run}&result=read_run"
    "&fields=fastq_ftp,fastq_md5,fastq_bytes,read_count,base_count,instrument_model,"
    "library_layout&format=tsv"
)

# ARIBA and KMA with default options only (protocol section 1).  The second
# ARIBA configuration of the comparator benchmark is not part of the protocol.
bench.CONFIGURATIONS = tuple(item for item in bench.CONFIGURATIONS if item[1] == "default")

IDENTICAL = bench.RELATIONS[0]
INCOMPLETE = bench.RELATIONS[2]
ABSENT_BUT_RETURNED = bench.RELATIONS[7]
NO_SEQUENCE = "no sequence returned"
STRATA = ("one copy", "several identical copies", "several different copies",
          "no copy (negative control)",
          "no copy under the protocol rule, bait aligned in parts")
# A unit without a copy under the rule of protocol 2.3 is a negative control
# only if the genome really lacks the locus.  When alignments of at least 80%
# identity still cover this fraction of the bait, in pieces that an insertion
# or a deletion separates, the locus is present and no truth sequence is
# defined for it.  Such units form the fifth stratum and are reported apart.
ALIGNED_IN_PARTS_MIN_FRACTION = 0.50
NOT_DEFINED = "sequence returned; no truth copy under the protocol rule"
# Not in the protocol: each draft is compared with the complete genome of the same
# isolate (MUMmer dnadiff).  The truth of protocol 2.3 holds only if the reads come
# from the genome that was deposited.  An isolate below either limit is kept in every
# table and is left out of one additional, labelled set of summaries.
CONCORDANCE_MIN_IDENTITY = 99.9          # % identity of the 1-to-1 alignments
CONCORDANCE_MIN_GENOME_ALIGNED = 98.0    # % of the complete genome aligned by the draft
SET_CONCORDANT = "protocol, isolates whose draft agrees with the complete genome"
ARIBA_BIT = {name: 1 << bit for bit, name in enumerate(bench.ARIBA_FLAGS)}
ARIBA_MUST_BE_SET = ("assembled", "assembled_into_one_contig", "unique_contig")
ARIBA_MUST_BE_UNSET = ("region_assembled_twice", "scaffold_graph_bad", "assembly_fail",
                       "variants_suggest_collapsed_repeat")
TOOLS = ("Locus-Recon", "ARIBA", "KMA")


# ------------------------------------------------------------------------- helpers

def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, rows, header=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if header is None:
        header = []
        for row in rows:
            for key in row:
                if key not in header:
                    header.append(key)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _cell(row.get(key, "")) for key in header})


def _cell(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.6g}"
    return value


def load_json(path, default):
    path = Path(path)
    return json.loads(path.read_text()) if path.is_file() else default


def run_logged(command, log, cwd=None):
    """Run a command with its output in a log; return exit code and seconds."""
    started = time.time()
    with open(log, "w") as handle:
        handle.write("$ " + " ".join(shlex.quote(str(part)) for part in command) + "\n")
        handle.flush()
        completed = subprocess.run([str(part) for part in command], stdout=handle,
                                   stderr=subprocess.STDOUT, cwd=cwd, check=False)
    return completed.returncode, round(time.time() - started, 1)


def file_digests(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            md5.update(block)
            sha.update(block)
    return md5.hexdigest(), sha.hexdigest()


def gc_percent(sequence):
    sequence = sequence.upper()
    called = sum(sequence.count(letter) for letter in "ACGT")
    return round(100.0 * (sequence.count("G") + sequence.count("C")) / called, 2) if called else ""


def read_genome(path):
    """Ordered ``{name: sequence}`` of a FASTA file, upper case."""
    return {name: sequence.upper() for name, sequence in bench.read_fasta(path)}


def stage_record(work, stage, args, tools, extra=None):
    record = {
        "stage": stage, "finished_utc": now(), "platform": platform.platform(),
        "python": sys.version.split()[0], "script_sha256": sha256_file(__file__),
        "threads": args.threads,
        "isolates": args.isolates or "all",
        "tools": {name: command_version(command) for name, command in sorted(tools.items())},
    }
    record.update(extra or {})
    directory = Path(work) / "stage_records"
    directory.mkdir(parents=True, exist_ok=True)
    write_json(directory / f"{stage}.json", record)
    return record


def check_frozen_package():
    """Version and source digest of the package that is imported and of the one installed."""
    import locus_recon
    imported_root = Path(locus_recon.__file__).resolve().parent.parent
    record = {
        "imported_from": str(imported_root),
        "imported_version": locus_recon.VERSION,
        "imported_source_sha256": source_sha256(imported_root),
    }
    probe = subprocess.run(
        [sys.executable, "-c", "import locus_recon; print(locus_recon.__file__)"],
        capture_output=True, text=True, cwd="/", check=False)
    if probe.returncode == 0 and probe.stdout.strip():
        installed_root = Path(probe.stdout.strip()).resolve().parent.parent
        record["installed_from"] = str(installed_root)
        record["installed_source_sha256"] = source_sha256(installed_root)
    script = subprocess.run(["locus-recon", "--version"], capture_output=True, text=True,
                            check=False)
    record["console_script_version"] = (script.stdout + script.stderr).strip()
    record["frozen_source_sha256"] = FROZEN_DIGEST
    problems = []
    for key in ("imported_source_sha256", "installed_source_sha256"):
        if record.get(key) != FROZEN_DIGEST:
            problems.append(f"{key} is {record.get(key)}")
    if FROZEN_VERSION not in record["console_script_version"]:
        problems.append(f"locus-recon --version says {record['console_script_version']!r}")
    if problems:
        raise RuntimeError("the software is not the frozen release: " + "; ".join(problems))
    return record


def isolate_dir(work, strain):
    return Path(work) / "isolates" / strain


def load_isolates(work):
    """Records of the prepare stage, one file per isolate: ``{strain: record}``."""
    return {path.stem: json.loads(path.read_text())
            for path in sorted((Path(work) / "prepare" / "isolates").glob("*.json"))}


def load_units(work):
    """Units of the prepare stage, one file per isolate: ``{strain|locus: unit}``."""
    units = {}
    for path in sorted((Path(work) / "prepare" / "units").glob("*.json")):
        units.update(json.loads(path.read_text()))
    for unit in units.values():
        unit["truth_defined"] = unit["unit_set"] != SET_UNAVAILABLE
        if unit["unit_set"] == SET_UNAVAILABLE or unit["truth_copies"]:
            continue
        table = isolate_dir(work, unit["strain"]) / "truth" / f"{unit['locus']}.blastn.tsv"
        fraction, longest = bait_alignment(table, unit["bait_bp"])
        unit["bait_fraction_aligned"] = fraction
        unit["bait_fraction_in_longest_alignment"] = longest
        if fraction >= ALIGNED_IN_PARTS_MIN_FRACTION:
            unit["stratum"] = STRATA[4]
            unit["truth_defined"] = False
    return units


def bait_alignment(table, bait_length):
    """Fraction of the bait covered by all alignments of >=80% identity, and by the longest."""
    covered, longest = set(), 0
    with open(table) as handle:
        for line in handle:
            parts = line.split("\t")
            if float(parts[2]) < TRUTH_MIN_IDENTITY:
                continue
            start, end = sorted((int(parts[4]), int(parts[5])))
            covered.update(range(start, end + 1))
            longest = max(longest, end - start + 1)
    return round(len(covered) / bait_length, 4), round(longest / bait_length, 4)


def unit_key(strain, locus):
    return f"{strain}|{locus}"


def selected(isolates, args):
    wanted = set(filter(None, (args.isolates or "").split(",")))
    return [entry for entry in isolates if not wanted or entry["strain"] in wanted]


def selected_loci(args):
    wanted = set(filter(None, (args.loci or "").split(",")))
    return [locus for locus in LOCI if not wanted or locus in wanted]


# ------------------------------------------------------------------ stage: prepare

def refseq_file(accession, cache, suffix):
    """One file of a RefSeq assembly directory, downloaded once."""
    cached = [path for path in Path(cache).glob(f"{accession}_*_{suffix}")
              if path.stat().st_size and (
                  "_from_genomic" in suffix or "_from_genomic" not in path.name)]
    if len(cached) == 1:
        return cached[0]
    directory, name = calibration._refseq_directory(accession)
    return Path(calibration._fetch(directory + f"{name}_{suffix}",
                                   Path(cache) / f"{name}_{suffix}"))


def annotation_entries(records):
    """Fields of the records of a ``*_from_genomic.fna`` file, in coordinate order."""
    entries, order = [], {}
    for header, sequence in records:
        identifier = header.split()[0]
        match = re.match(r"lcl\|(.+?)_(?:cds|rrna|trna|ncrna|tmrna|mrna|misc_rna|"
                         r"precursor_rna|miscrna)_", identifier)
        replicon = match.group(1) if match else identifier
        order.setdefault(replicon, len(order))
        fields = dict(re.findall(r"\[(\w+)=([^\]]*)\]", header))
        location = fields.get("location", "")
        numbers = [int(number) for number in re.findall(r"\d+", location)]
        entries.append({
            "identifier": identifier, "replicon": replicon,
            "replicon_order": order[replicon],
            "start": min(numbers) if numbers else 0, "end": max(numbers) if numbers else 0,
            "strand": "-" if "complement" in location else "+",
            "location": location, "gene": fields.get("gene", ""),
            "product": fields.get("protein") or fields.get("product", ""),
            "locus_tag": fields.get("locus_tag", ""), "gbkey": fields.get("gbkey", ""),
            "pseudo": fields.get("pseudo", "") == "true",
            "partial": "<" in location or ">" in location,
            "sequence": sequence.upper(),
        })
    entries.sort(key=lambda entry: (entry["replicon_order"], entry["start"]))
    return entries


def select_bait(locus, cds, rna):
    """The bait of a locus in a reference annotation: (entry, rule, candidates).

    The protocol rule is the gene symbol for a protein-coding locus and the
    product name for an rRNA gene, first copy by coordinate.  The product name
    of an rRNA gene is compared without regard to case.  The two fallback
    rules are described at ``FALLBACK_PRODUCT``.
    """
    if locus in RRNA_PRODUCT:
        found = [entry for entry in rna if entry["gbkey"] == "rRNA"
                 and entry["product"].lower() == RRNA_PRODUCT[locus].lower()]
        return (found[0] if found else None), (RULE_PROTOCOL if found else "none"), len(found)
    found = [entry for entry in cds if entry["gene"] == locus]
    if found:
        return found[0], RULE_PROTOCOL, len(found)
    symbols = {locus + suffix for suffix in FALLBACK_PARALOGUE_SUFFIXES}
    found = [entry for entry in cds if entry["gene"] in symbols]
    if found:
        return found[0], "fallback: gene symbol with a paralogue letter", len(found)
    products = {product.lower() for product in FALLBACK_PRODUCT[locus]}
    found = [entry for entry in cds if entry["product"].lower() in products]
    if found:
        return found[0], "fallback: product name of the CDS", len(found)
    return None, "none", 0


def blast_bait_to_genome(bait, genome, output, task):
    command = ["blastn", "-task", task, "-query", str(bait), "-subject", str(genome),
               "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send bitscore evalue",
               "-evalue", "1e-10", "-out", str(output)]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode:
        raise RuntimeError(f"blastn failed: {completed.stderr[-500:]}")
    hits = []
    with open(output) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            hits.append({
                "replicon": parts[1], "identity": float(parts[2]), "length": int(parts[3]),
                "qstart": int(parts[4]), "qend": int(parts[5]),
                "sstart": int(parts[6]), "send": int(parts[7]), "bitscore": float(parts[8]),
            })
    return hits, command


def truth_copies(hits, bait_length, genome):
    """Copies of a locus in a complete genome (protocol 2.3).

    Every alignment of at least 80% identity that covers at least 90% of the
    bait is one copy, and the copy is extended to the full span of the bait,
    base for base.  Two alignments whose extended intervals overlap by more
    than half of the shorter one, on the same replicon, are the same copy, and
    the one with the higher bit score is kept.
    """
    copies = []
    replicon_order = {name: index for index, name in enumerate(genome)}
    for hit in sorted(hits, key=lambda item: -item["bitscore"]):
        covered = hit["qend"] - hit["qstart"] + 1
        if hit["identity"] < TRUTH_MIN_IDENTITY or covered < TRUTH_MIN_BAIT_FRACTION * bait_length:
            continue
        before, after = hit["qstart"] - 1, bait_length - hit["qend"]
        if hit["sstart"] <= hit["send"]:
            strand, low, high = "+", hit["sstart"] - before, hit["send"] + after
        else:
            strand, low, high = "-", hit["send"] - after, hit["sstart"] + before
        length = len(genome[hit["replicon"]])
        clipped = max(0, 1 - low) + max(0, high - length)
        low, high = max(1, low), min(length, high)
        duplicate = False
        for copy in copies:
            if copy["replicon"] != hit["replicon"]:
                continue
            overlap = min(high, copy["end"]) - max(low, copy["start"]) + 1
            if overlap > 0.5 * min(high - low + 1, copy["end"] - copy["start"] + 1):
                duplicate = True
                break
        if duplicate:
            continue
        sequence = genome[hit["replicon"]][low - 1:high]
        if strand == "-":
            sequence = reverse_complement(sequence)
        copies.append({
            "replicon": hit["replicon"], "start": low, "end": high, "strand": strand,
            "alignment_identity_pct": hit["identity"], "alignment_bp": hit["length"],
            "bait_start": hit["qstart"], "bait_end": hit["qend"],
            "extension_clipped_bp": clipped, "sequence": sequence,
        })
    copies.sort(key=lambda copy: (replicon_order[copy["replicon"]], copy["start"]))
    return copies


def stratum_of(copies):
    if not copies:
        return STRATA[3]
    if len(copies) == 1:
        return STRATA[0]
    return STRATA[1] if len({copy["sequence"] for copy in copies}) == 1 else STRATA[2]


def prepare_references(entry, work, cache):
    """Baits, truth copies and genome statistics of one isolate."""
    strain = entry["strain"]
    directory = isolate_dir(work, strain)
    (directory / "baits").mkdir(parents=True, exist_ok=True)
    (directory / "truth").mkdir(parents=True, exist_ok=True)

    reference = entry["bait_reference_assembly"]
    cds = annotation_entries(calibration._annotated_records(reference, cache, "cds"))
    rna = annotation_entries(calibration._annotated_records(reference, cache, "rna"))
    reference_files = [refseq_file(reference, cache, f"{kind}_from_genomic.fna.gz")
                       for kind in ("cds", "rna")]

    archive = refseq_file(entry["isolate_assembly"], cache, "genomic.fna.gz")
    genome_path = directory / "truth" / "complete_genome.fna"
    genome = read_genome(archive)
    bench.write_fasta(genome_path, list(genome.items()))
    genome_bp = sum(len(sequence) for sequence in genome.values())

    units = {}
    for locus in LOCI:
        bait, rule, candidates = select_bait(locus, cds, rna)
        unit = {
            "species": entry["species"], "strain": strain, "locus": locus,
            "locus_type": "rRNA gene" if locus in RRNA_PRODUCT else "protein-coding",
            "isolate_assembly": entry["isolate_assembly"],
            "illumina_run": entry["illumina_run"],
            "bait_reference_assembly": reference,
            "bait_rule": rule, "bait_candidates_in_annotation": candidates,
        }
        if bait is None:
            unit.update(unit_set=SET_UNAVAILABLE, stratum="not available", truth_copies="")
            units[unit_key(strain, locus)] = unit
            continue
        bait_path = directory / "baits" / f"{locus}.fasta"
        bench.write_fasta(bait_path, [(locus, bait["sequence"])])
        unit.update(
            unit_set=SET_PROTOCOL if rule == RULE_PROTOCOL else SET_SUPPLEMENTARY,
            bait_path=str(bait_path), bait_bp=len(bait["sequence"]),
            bait_gc_pct=gc_percent(bait["sequence"]),
            bait_source_record=bait["identifier"], bait_gene_symbol=bait["gene"],
            bait_product=bait["product"], bait_locus_tag=bait["locus_tag"],
            bait_location=f"{bait['replicon']}:{bait['location']}",
            bait_pseudo=bait["pseudo"], bait_partial=bait["partial"],
            bait_sha256=sha256_file(bait_path),
        )
        hits, command = blast_bait_to_genome(
            bait_path, genome_path, directory / "truth" / f"{locus}.blastn.tsv", "blastn")
        copies = truth_copies(hits, len(bait["sequence"]), genome)
        megablast, _ = blast_bait_to_genome(
            bait_path, genome_path, directory / "truth" / f"{locus}.megablast.tsv", "megablast")
        for index, copy in enumerate(copies, start=1):
            copy["name"] = f"{strain}|{locus}|copy{index}"
            copy["gc_pct"] = gc_percent(copy["sequence"])
            copy["identical_to_copy1"] = copy["sequence"] == copies[0]["sequence"]
        bench.write_fasta(directory / "truth" / f"{locus}.truth.fasta",
                          [(copy["name"], copy["sequence"]) for copy in copies])
        unit.update(
            truth_copies=len(copies),
            truth_distinct_sequences=len({copy["sequence"] for copy in copies}),
            truth_copies_megablast=len(truth_copies(megablast, len(bait["sequence"]), genome)),
            stratum=stratum_of(copies), truth_blast_command=" ".join(
                Path(part).name if part.startswith("/") else part for part in command),
            copies=copies,
        )
        units[unit_key(strain, locus)] = unit
    record = {
        "genome_fasta": str(genome_path), "genome_archive": archive.name,
        "genome_archive_sha256": sha256_file(archive),
        "genome_bp_observed": genome_bp, "genome_bp_table": int(entry["genome_bp"]),
        "genome_replicons": len(genome), "genome_gc_pct": gc_percent("".join(genome.values())),
        "bait_reference_files": [
            {"file": path.name, "sha256": sha256_file(path)} for path in reference_files],
    }
    return record, units


def ena_report(run):
    rows = list(csv.DictReader(
        calibration._read_url(ENA_FILEREPORT.format(run=run)).decode().splitlines(),
        delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"ENA returned {len(rows)} records for {run}, expected one")
    return rows[0]


def fetch_resumable(url, destination, expected_bytes, attempts=60):
    """Download a file of known size; an interrupted transfer is resumed.

    A connection that closes early leaves a shorter file without any error, so
    the size that ENA reports is the test of completeness.  The bytes already
    received are kept in ``<name>.part`` and the next request asks for the rest.
    """
    destination = Path(destination)
    partial = destination.with_name(destination.name + ".part")
    if destination.is_file():
        if destination.stat().st_size == expected_bytes:
            return destination
        destination.replace(partial)
    for attempt in range(attempts):
        have = partial.stat().st_size if partial.is_file() else 0
        if have == expected_bytes:
            break
        if have > expected_bytes:
            partial.unlink()
            have = 0
        request = urllib.request.Request(
            url, headers={"Range": f"bytes={have}-"} if have else {})
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                resumed = have and response.status == 206
                with open(partial, "ab" if resumed else "wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
        except (urllib.error.URLError, http.client.HTTPException, OSError) as error:
            print(f"  {destination.name}: {error!r}; resuming", flush=True)
        if partial.is_file() and partial.stat().st_size == expected_bytes:
            break
        time.sleep(min(60, 5 * (attempt + 1)))
    if not partial.is_file() or partial.stat().st_size != expected_bytes:
        raise RuntimeError(f"{destination.name} is incomplete after {attempts} requests")
    partial.replace(destination)
    return destination


def fastq_size(path, seqtk):
    """Number of records and of bases of a FASTQ file."""
    completed = subprocess.run([seqtk, "size", str(path)], capture_output=True, text=True,
                               check=False)
    fields = completed.stdout.split()
    if completed.returncode == 0 and len(fields) >= 2:
        return int(fields[0]), int(fields[1])
    records = bases = 0
    with gzip.open(path, "rt") as handle:
        for index, line in enumerate(handle):
            if index % 4 == 1:
                records += 1
                bases += len(line) - 1
    return records, bases


def pairs_in_step(first, second):
    """Number of records whose names differ between two FASTQ files, and of records."""
    name = "awk 'NR%4==1{split($1,part,\"/\"); print part[1]}'"
    command = (f"paste <(gzip -cd {shlex.quote(str(first))} | {name}) "
               f"<(gzip -cd {shlex.quote(str(second))} | {name}) "
               "| awk '$1!=$2{bad++} END{print bad+0, NR}'")
    completed = subprocess.run(["bash", "-c", command], capture_output=True, text=True,
                               check=False)
    bad, records = completed.stdout.split()[:2]
    return int(bad), int(records)


def subsample_file(source, destination, pairs, seqtk):
    """``seqtk sample`` with the fixed seed; the same seed keeps mates in step."""
    partial = Path(str(destination) + ".part")
    with open(partial, "wb") as handle:
        sampler = subprocess.Popen(
            [seqtk, "sample", "-s", str(SEED), str(source), str(pairs)],
            stdout=subprocess.PIPE)
        packer = subprocess.Popen(["gzip", "-n", "-c"], stdin=sampler.stdout, stdout=handle)
        sampler.stdout.close()
        codes = (packer.wait(), sampler.wait())
    if any(codes):
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"subsampling of {source} failed (exit codes {codes})")
    partial.replace(destination)


def prepare_reads(entry, work, cache, seqtk, previous, verify, processors):
    """Download, verify and subsample the reads of one isolate."""
    strain, run = entry["strain"], entry["illumina_run"]
    reads_dir = isolate_dir(work, strain) / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)
    if previous and not verify and all(
            Path(path).is_file() for path in (previous["r1"], previous["r2"])):
        return previous
    report = ena_report(run)
    fetched = []
    for url, md5, size in zip(report["fastq_ftp"].split(";"), report["fastq_md5"].split(";"),
                              report["fastq_bytes"].split(";")):
        mate = re.search(r"_([12])\.fastq\.gz$", url)
        if mate:
            fetched.append((mate, url, md5, size, fetch_resumable(
                "https://" + url, Path(cache) / Path(url).name, int(size))))
    # Downloads wait on the network; what follows uses processors and is rationed.
    with processors:
        return _verify_and_subsample(entry, report, fetched, reads_dir, seqtk, verify)


def _verify_and_subsample(entry, report, fetched, reads_dir, seqtk, verify):
    strain, run = entry["strain"], entry["illumina_run"]
    files = {}
    for mate, url, md5, size, path in fetched:
        local_md5, sha256 = file_digests(path)
        if local_md5 != md5:
            path.unlink()
            raise RuntimeError(f"MD5 mismatch for {path.name}; file removed, run again")
        records, bases = fastq_size(path, seqtk)
        files[mate.group(1)] = {
            "file": path.name, "path": str(path), "url": "https://" + url, "ena_md5": md5,
            "local_md5": local_md5, "sha256": sha256, "size_bytes": path.stat().st_size,
            "ena_bytes": int(size), "records": records, "bases": bases,
        }
    if sorted(files) != ["1", "2"]:
        raise RuntimeError(f"{run} does not have the two mate files in ENA")
    if files["1"]["records"] != files["2"]["records"]:
        raise RuntimeError(f"{run}: the two mate files hold different numbers of reads")
    pairs = files["1"]["records"]
    bases = files["1"]["bases"] + files["2"]["bases"]
    genome_bp = int(entry["genome_bp"])
    target = TARGET_DEPTH * genome_bp
    record = {
        "strain": strain, "illumina_run": run, "instrument": report.get("instrument_model", ""),
        "library_layout": report.get("library_layout", ""),
        "ena_read_count": int(report["read_count"]), "ena_base_count": int(report["base_count"]),
        "pairs": pairs, "bases": bases, "depth_full_x": round(bases / genome_bp, 1),
        "files": files, "seed": SEED, "target_depth_x": TARGET_DEPTH,
    }
    if bases <= target:
        record.update(subsampled=False, pairs_kept=pairs, r1=files["1"]["path"],
                      r2=files["2"]["path"], depth_used_x=record["depth_full_x"])
        return record
    keep = math.ceil(target * pairs / bases)
    outputs = {}
    for mate in ("1", "2"):
        destination = reads_dir / f"{run}_{TARGET_DEPTH}x_seed{SEED}_{mate}.fastq.gz"
        if verify and destination.is_file():
            check = Path(str(destination) + ".verify")
            subsample_file(files[mate]["path"], check, keep, seqtk)
            same = sha256_file(check) == sha256_file(destination)
            check.unlink()
            if not same:
                raise RuntimeError(f"{destination.name} is not reproduced by a second run")
        elif not destination.is_file():
            subsample_file(files[mate]["path"], destination, keep, seqtk)
        records, kept_bases = fastq_size(destination, seqtk)
        outputs[mate] = {"file": destination.name, "path": str(destination),
                         "sha256": sha256_file(destination), "records": records,
                         "bases": kept_bases, "size_bytes": destination.stat().st_size}
    out_of_step, compared = pairs_in_step(outputs["1"]["path"], outputs["2"]["path"])
    if (out_of_step or compared != keep or outputs["1"]["records"] != keep
            or outputs["2"]["records"] != keep):
        raise RuntimeError(
            f"{run}: subsampled mates are not in step ({out_of_step} of {compared} names "
            f"differ; {outputs['1']['records']} and {outputs['2']['records']} reads for "
            f"{keep} requested)")
    kept = outputs["1"]["bases"] + outputs["2"]["bases"]
    record.update(
        subsampled=True, pairs_kept=keep, subsampled_files=outputs,
        r1=outputs["1"]["path"], r2=outputs["2"]["path"],
        depth_used_x=round(kept / genome_bp, 1), subsample_reproduced=bool(verify),
        subsample_command=f"seqtk sample -s {SEED} <mate file> {keep} | gzip -n -c",
    )
    return record


def stage_prepare(args, isolates, work):
    cache = work / "downloads"
    cache.mkdir(parents=True, exist_ok=True)
    for name in ("isolates", "units"):
        (work / "prepare" / name).mkdir(parents=True, exist_ok=True)
    isolate_records = load_isolates(work)
    for entry in isolates:
        strain = entry["strain"]
        print(f"[prepare] {strain}: references, baits and truth", flush=True)
        record, isolate_units = prepare_references(entry, work, cache)
        isolate_records.setdefault(strain, {}).update(record)
        write_json(work / "prepare" / "isolates" / f"{strain}.json", isolate_records[strain])
        write_json(work / "prepare" / "units" / f"{strain}.json", isolate_units)

    # One isolate at a time per four threads verifies and subsamples (at most
    # four busy processes each); the downloads themselves wait on the network.
    processors = threading.Semaphore(max(1, args.threads // 4))

    def reads(entry):
        strain = entry["strain"]
        previous = isolate_records.get(strain, {}).get("reads")
        result = prepare_reads(entry, work, cache, args.seqtk, previous,
                               args.verify_subsample, processors)
        isolate_records[strain]["reads"] = result
        write_json(work / "prepare" / "isolates" / f"{strain}.json", isolate_records[strain])
        print(f"[prepare] {strain}: {result['pairs']} pairs, "
              f"{result['depth_full_x']}x; used {result['pairs_kept']} pairs, "
              f"{result['depth_used_x']}x", flush=True)

    with ThreadPoolExecutor(max(1, args.download_jobs)) as pool:
        list(pool.map(reads, isolates))
    stage_record(work, "prepare", args, {
        "seqtk": ["bash", "-c", f"{shlex.quote(args.seqtk)} 2>&1 | grep Version"],
        "blastn": ["blastn", "-version"], "gzip": ["gzip", "--version"],
    }, {"seed": SEED, "target_depth_x": TARGET_DEPTH,
        "verify_subsample": bool(args.verify_subsample)})


# ----------------------------------------------------------------- stage: assemble

def draft_statistics(path):
    lengths = sorted((len(sequence) for _, sequence in bench.read_fasta(path)), reverse=True)
    total, running, n50 = sum(lengths), 0, 0
    for length in lengths:
        running += length
        if running * 2 >= total:
            n50 = length
            break
    return {"contigs": len(lengths), "total_bp": total, "n50_bp": n50,
            "longest_bp": lengths[0] if lengths else 0}


def assemble(entry, work, isolate_records, threads):
    """Draft assembly with SPAdes in isolate mode, as ``run_tier_calibration._draft``."""
    strain = entry["strain"]
    directory = isolate_dir(work, strain)
    output = directory / "spades"
    record_path = directory / "spades.run.json"
    draft = output / "contigs.fasta"
    graph = output / "assembly_graph_after_simplification.gfa"
    if draft.is_file() and graph.is_file() and record_path.is_file():
        return json.loads(record_path.read_text())
    shutil.rmtree(output, ignore_errors=True)
    reads = isolate_records[strain]["reads"]
    command = ["spades.py", "--isolate", "-1", reads["r1"], "-2", reads["r2"],
               "-o", str(output), "-t", str(threads), "-m", str(calibration.SPADES_MEMORY_GB)]
    code, seconds = run_logged(command, directory / "spades.stdout.log")
    if code or not draft.is_file() or not graph.is_file():
        raise RuntimeError(f"SPAdes failed for {strain} (exit {code})")
    record = {"strain": strain, "exit_code": code, "seconds": seconds,
              "command": [Path(part).name if part.startswith("/") else part
                          for part in command],
              "draft": str(draft), "graph": str(graph),
              "draft_sha256": sha256_file(draft), "graph_sha256": sha256_file(graph)}
    record.update(draft_statistics(draft))
    write_json(record_path, record)
    print(f"[assemble] {strain}: {record['contigs']} contigs, {record['total_bp']} bp, "
          f"N50 {record['n50_bp']} bp, {seconds} s", flush=True)
    return record


def stage_assemble(args, isolates, work):
    isolate_records = load_isolates(work)
    jobs = max(1, min(args.assembly_jobs, len(isolates)))
    # One thread count for every draft, whatever the number of drafts made at a time.
    threads = args.assembly_threads or max(1, args.threads // jobs)
    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(lambda entry: assemble(entry, work, isolate_records, threads), isolates))
    stage_record(work, "assemble", args, {"spades": ["spades.py", "--version"]},
                 {"threads_per_assembly": threads, "assemblies_at_a_time": jobs,
                  "spades_memory_gb": calibration.SPADES_MEMORY_GB})


# -------------------------------------------------------------- stage: reconstruct

def unit_paths(work, unit):
    strain, locus = unit["strain"], unit["locus"]
    output = isolate_dir(work, strain) / "locus_recon" / locus
    sample = output / strain
    copy_number = isolate_dir(work, strain) / "copy_number"
    return {
        "output": output, "sample": sample,
        "report": output / f"locus_recon_report_{locus}.tsv",
        "candidate": sample / f"{strain}_{locus}_reconstructed.fasta",
        "qc_report": sample / f"{strain}_{locus}_qc_report.txt",
        "support": sample / "allele_remap.support.tsv",
        "discovery": sample / "blast_step1.tsv",
        "bam": sample / "full_map.sorted.bam",
        "run": output / "run.json", "log": output.parent / f"{locus}.stdout.log",
        "copy_number": copy_number,
        "depth_tsv": copy_number / f"{locus}.depth_ratio.tsv",
        "copy_json": copy_number / f"{locus}.copy_number.json",
        "copy_tsv": copy_number / f"{locus}.copy_number.tsv",
        "copy_run": copy_number / f"{locus}.run.json",
        "block_json": copy_number / f"{locus}.block_size.json",
    }


def reconstruct(unit, work, threads):
    """One Locus-Recon run with default options."""
    paths = unit_paths(work, unit)
    if paths["run"].is_file() and paths["report"].is_file():
        return json.loads(paths["run"].read_text())
    strain, locus = unit["strain"], unit["locus"]
    spades = json.loads((isolate_dir(work, strain) / "spades.run.json").read_text())
    reads = load_isolates(work)[strain]["reads"]
    shutil.rmtree(paths["output"], ignore_errors=True)
    paths["output"].parent.mkdir(parents=True, exist_ok=True)
    samplesheet = paths["output"].parent / f"samples_{locus}.tsv"
    samplesheet.write_text("sample_id\tassembly_path\tr1_path\tr2_path\n"
                           f"{strain}\t{spades['draft']}\t{reads['r1']}\t{reads['r2']}\n")
    command = ["locus-recon", "--samplesheet", str(samplesheet),
               "--main-output-dir", str(paths["output"]), "--bait", unit["bait_path"],
               "--locus", locus, "--threads", str(threads)]
    if locus in RRNA_PRODUCT:
        command.append("--noncoding-locus")
    code, seconds = run_logged(command, paths["log"])
    record = {"strain": strain, "locus": locus, "exit_code": code, "seconds": seconds,
              "command": [Path(part).name if part.startswith("/") else part
                          for part in command]}
    if paths["report"].is_file():
        write_json(paths["run"], record)
    row = read_table(paths["report"])[0] if paths["report"].is_file() else {}
    print(f"[reconstruct] {strain} {locus}: exit {code}, {seconds} s, "
          f"{row.get('result_disposition', 'no report')} "
          f"{row.get('sequence_confidence', '')}", flush=True)
    return record


def runnable_units(args, isolates, work):
    units = load_units(work)
    loci = selected_loci(args)
    return [units[unit_key(entry["strain"], locus)] for entry in isolates for locus in loci
            if units[unit_key(entry["strain"], locus)]["unit_set"] != SET_UNAVAILABLE]


def stage_reconstruct(args, isolates, work):
    package = check_frozen_package()
    units = runnable_units(args, isolates, work)
    jobs = max(1, args.threads // args.threads_per_run)
    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(lambda unit: reconstruct(unit, work, args.threads_per_run), units))
    stage_record(work, "reconstruct", args, {
        "locus-recon": ["locus-recon", "--version"], "spades": ["spades.py", "--version"],
        "bwa-mem2": ["bash", "-c", "bwa-mem2 version 2>/dev/null | tail -1"],
        "samtools": ["samtools", "--version"], "blastn": ["blastn", "-version"],
    }, {"package": package, "threads_per_run": args.threads_per_run, "runs_at_a_time": jobs})


# -------------------------------------------------------------- stage: copy-number

def discovery_intervals(path, bait_length):
    """Intervals of the locus on the draft, from ``blast_step1.tsv``.

    The selection is that of ``validation/workflows/run_low_copy_23S.py``:
    alignments of at least 80% identity and at least max(200 bp, 10% of the
    bait) are taken in decreasing order of bit score, and one that overlaps an
    interval already taken on the same contig is dropped.
    """
    hits = []
    if not Path(path).is_file():
        return [], []
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            _, contig, qstart, qend, sstart, send, length, identity, _, bitscore = parts[:10]
            if float(identity) < 80.0 or int(length) < max(200, int(0.10 * bait_length)):
                continue
            hits.append({
                "contig": contig,
                "subject_start": min(int(sstart), int(send)),
                "subject_end": max(int(sstart), int(send)),
                "query_start": min(int(qstart), int(qend)),
                "query_end": max(int(qstart), int(qend)),
                "bitscore": float(bitscore),
            })
    retained = []
    for hit in sorted(hits, key=lambda item: item["bitscore"], reverse=True):
        overlaps = any(
            old["contig"] == hit["contig"]
            and hit["subject_start"] <= old["subject_end"]
            and old["subject_start"] <= hit["subject_end"]
            for old in retained)
        if not overlaps:
            retained.append(hit)
    retained.sort(key=lambda item: (item["contig"], item["subject_start"]))
    regions = [f"{hit['contig']}:{hit['subject_start']}-{hit['subject_end']}"
               for hit in retained]
    spans = [(hit["query_start"], hit["query_end"]) for hit in retained]
    return regions, spans


def copy_number(unit, work):
    """``locus-recon-depth-ratio`` and ``locus-recon-copy-number`` for one unit."""
    paths = unit_paths(work, unit)
    if paths["copy_run"].is_file():
        return json.loads(paths["copy_run"].read_text())
    paths["copy_number"].mkdir(parents=True, exist_ok=True)
    strain, locus = unit["strain"], unit["locus"]
    regions, spans = discovery_intervals(paths["discovery"], unit["bait_bp"])
    record = {"strain": strain, "locus": locus, "regions": regions,
              "gene_spans": [f"{start}-{end}" for start, end in spans],
              "gene_length": unit["bait_bp"]}
    if not regions or not paths["bam"].is_file():
        record["status"] = ("no discovery interval" if not regions
                            else "no full_map.sorted.bam")
        write_json(paths["copy_run"], record)
        return record
    spades = json.loads((isolate_dir(work, strain) / "spades.run.json").read_text())
    shared = ["--bam", str(paths["bam"])]
    for region, (start, end) in zip(regions, spans):
        shared += ["--locus", region, "--gene-span", f"{start}-{end}"]
    shared += ["--gene-length", str(unit["bait_bp"])]
    depth = ["locus-recon-depth-ratio"] + shared + ["--tsv", str(paths["depth_tsv"])]
    code, seconds = run_logged(depth, paths["copy_number"] / f"{locus}.depth_ratio.log")
    record["depth_ratio"] = {"exit_code": code, "seconds": seconds, "command": [
        Path(part).name if part.startswith("/") else part for part in depth]}
    combined = ["locus-recon-copy-number"] + shared + [
        "--gfa", spades["graph"], "--bait", unit["bait_path"],
        "--output-json", str(paths["copy_json"]), "--output-tsv", str(paths["copy_tsv"])]
    code, seconds = run_logged(combined, paths["copy_number"] / f"{locus}.copy_number.log")
    record["copy_number"] = {"exit_code": code, "seconds": seconds, "command": [
        Path(part).name if part.startswith("/") else part for part in combined]}
    record["status"] = "run"
    write_json(paths["copy_run"], record)
    print(f"[copy-number] {strain} {locus}: {len(regions)} interval(s), exit codes "
          f"{record['depth_ratio']['exit_code']} and {code}", flush=True)
    return record


def stage_copy_number(args, isolates, work):
    package = check_frozen_package()
    units = runnable_units(args, isolates, work)
    with ThreadPoolExecutor(max(1, args.threads)) as pool:
        list(pool.map(lambda unit: copy_number(unit, work), units))
    stage_record(work, "copy-number", args, {
        "locus-recon-depth-ratio": ["locus-recon-depth-ratio", "--version"],
        "locus-recon-copy-number": ["locus-recon-copy-number", "--version"],
        "samtools": ["samtools", "--version"], "blastn": ["blastn", "-version"],
    }, {"package": package})


# --------------------------------------------------------------- stage: block-size

def block_size_unit(job):
    """Depth ratio of one unit at each block size, through the library function."""
    from locus_recon.depth_ratio import estimate_locus_copy_number
    bam, regions, spans, gene_length, destination = job
    rows = []
    for block in BLOCK_SIZES:
        started = time.time()
        result = estimate_locus_copy_number(
            bam, regions, gene_spans=[tuple(span) for span in spans],
            gene_length=gene_length, block=block)
        row = result.as_row()
        row["seconds"] = round(time.time() - started, 1)
        rows.append(row)
    write_json(destination, rows)
    return destination


def stage_block_size(args, isolates, work):
    package = check_frozen_package()
    jobs = []
    for unit in runnable_units(args, isolates, work):
        paths = unit_paths(work, unit)
        record = load_json(paths["copy_run"], {})
        if record.get("status") != "run" or paths["block_json"].is_file():
            continue
        spans = [tuple(int(value) for value in span.split("-")) for span in record["gene_spans"]]
        jobs.append((str(paths["bam"]), record["regions"], spans, record["gene_length"],
                     str(paths["block_json"])))
    with ProcessPoolExecutor(max(1, args.threads)) as pool:
        for destination in pool.map(block_size_unit, jobs):
            print(f"[block-size] {destination}", flush=True)
    stage_record(work, "block-size", args, {"samtools": ["samtools", "--version"]},
                 {"package": package, "block_sizes": list(BLOCK_SIZES),
                  "function": "locus_recon.depth_ratio.estimate_locus_copy_number(block=...)",
                  "other_arguments": "defaults (n_boot 1000, seed 20260817, "
                                     "max_backbone_bootstrap_blocks 512)"})


# -------------------------------------------------------------- stage: comparators

def comparator_dir(work):
    return Path(work) / "comparators"


def comparators(entry, work, units, threads):
    """ARIBA and KMA on the reads and baits of one isolate, as the comparator benchmark."""
    strain = entry["strain"]
    directory = comparator_dir(work) / strain
    directory.mkdir(parents=True, exist_ok=True)
    available = [unit for unit in (units[unit_key(strain, locus)] for locus in LOCI)
                 if unit["unit_set"] != SET_UNAVAILABLE]
    records = [(unit["locus"], bench.read_fasta(unit["bait_path"])[0][1].upper())
               for unit in available]
    names = [name for name, _ in records]
    fasta = directory / "catalogue.fasta"
    bench.write_fasta(fasta, records)
    coding = {name: name not in RRNA_PRODUCT for name in names}
    relaxed = directory / "declared_noncoding_for_ariba.txt"
    if relaxed.is_file():
        for name in relaxed.read_text().split():
            coding[name] = False
    metadata = directory / "catalogue_ariba_metadata.tsv"
    with open(metadata, "w") as handle:
        for name in names:
            handle.write(f"{name}\t{1 if coding[name] else 0}\t0\t.\t.\t{name}\n")
    prepared, database, noncoding = bench.prepare_references(
        strain, fasta, metadata, names, comparator_dir(work))
    reads = load_isolates(work)[strain]["reads"]
    sample = {"sample": strain, "r1": reads["r1"], "r2": reads["r2"]}
    bench.run_tools(strain, sample, prepared, database, comparator_dir(work), threads)
    write_json(directory / "catalogue.json", {
        "names": names, "declared_noncoding_for_ariba": noncoding,
        "catalogue_sha256": sha256_file(fasta)})
    return strain, noncoding


def stage_comparators(args, isolates, work):
    units = load_units(work)
    jobs = max(1, args.threads // args.threads_per_run)
    with ThreadPoolExecutor(jobs) as pool:
        relaxed = dict(pool.map(
            lambda entry: comparators(entry, work, units, args.threads_per_run), isolates))
    versions = bench.tool_versions()
    stage_record(work, "comparators", args, {
        "kma": ["kma", "-v"], "blastn": ["blastn", "-version"],
    }, {"tool_versions_text": versions, "threads_per_run": args.threads_per_run,
        "isolates_at_a_time": jobs,
        "configurations": [{"tool": tool, "configuration": name, "extra_arguments": extra}
                           for tool, name, extra in bench.CONFIGURATIONS],
        "catalogue_sequences_declared_noncoding_for_ariba": relaxed})


# -------------------------------------------------------------- stage: concordance

def concordance(entry, work):
    """Agreement of the draft of one isolate with its complete genome (dnadiff)."""
    strain = entry["strain"]
    directory = Path(work) / "concordance"
    directory.mkdir(parents=True, exist_ok=True)
    record_path = directory / f"{strain}.json"
    if record_path.is_file():
        return json.loads(record_path.read_text())
    genome = isolate_dir(work, strain) / "truth" / "complete_genome.fna"
    draft = isolate_dir(work, strain) / "spades" / "contigs.fasta"
    command = ["dnadiff", "-p", strain, str(genome), str(draft)]
    code, seconds = run_logged(command, directory / f"{strain}.dnadiff.log", cwd=directory)
    report = directory / f"{strain}.report"
    if code or not report.is_file():
        raise RuntimeError(f"dnadiff failed for {strain} (exit {code})")
    text = report.read_text()

    def first(name, pattern):
        return re.search(rf"^{name}\s+{pattern}", text, re.MULTILINE)

    aligned = first("AlignedBases", r"(\d+)\(([\d.]+)%\)\s+(\d+)\(([\d.]+)%\)")
    record = {
        "strain": strain, "exit_code": code, "seconds": seconds,
        "command": ["dnadiff", "-p", strain, genome.name, draft.name],
        "genome_aligned_bp": int(aligned.group(1)),
        "genome_aligned_pct": float(aligned.group(2)),
        "draft_aligned_bp": int(aligned.group(3)),
        "draft_aligned_pct": float(aligned.group(4)),
        "identity_1_to_1_pct": float(first("AvgIdentity", r"([\d.]+)").group(1)),
        "snps": int(first("TotalSNPs", r"(\d+)").group(1)),
        "indels": int(first("TotalIndels", r"(\d+)").group(1)),
    }
    record["reads_agree_with_genome"] = (
        record["identity_1_to_1_pct"] >= CONCORDANCE_MIN_IDENTITY
        and record["genome_aligned_pct"] >= CONCORDANCE_MIN_GENOME_ALIGNED)
    write_json(record_path, record)
    print(f"[concordance] {strain}: {record['genome_aligned_pct']}% of the genome aligned, "
          f"{record['identity_1_to_1_pct']}% identity, {record['snps']} SNPs", flush=True)
    return record


def stage_concordance(args, isolates, work):
    with ThreadPoolExecutor(max(1, args.threads // 2)) as pool:
        list(pool.map(lambda entry: concordance(entry, work), isolates))
    stage_record(work, "concordance", args, {
        "nucmer": ["bash", "-c", "nucmer --version 2>&1 | grep -i version"],
    }, {"min_identity_pct": CONCORDANCE_MIN_IDENTITY,
        "min_genome_aligned_pct": CONCORDANCE_MIN_GENOME_ALIGNED})


# -------------------------------------------------------------------- stage: score

def result_class(comparison, unit):
    """The outcome classes of protocol 2.4."""
    if unit["stratum"] == STRATA[4]:
        return ("no sequence; truth not defined" if comparison == NO_SEQUENCE
                else "sequence returned; truth not defined")
    if unit["stratum"] == STRATA[3]:
        return ("no sequence, locus absent" if comparison == NO_SEQUENCE
                else "sequence returned, locus absent")
    if comparison == NO_SEQUENCE:
        return "no sequence"
    if comparison == IDENTICAL:
        return "exact"
    if comparison == INCOMPLETE:
        return "incomplete, identical over its length"
    return "wrong, extra or missing bases"


def score_sequences(found, truth_records, blast_work):
    """Best relation of the returned sequences to the truth, as ``bench.score_sample``."""
    if not found:
        return {"sequences_returned": 0, "returned_bp": 0, "uncalled_bases": 0,
                "lowercase_bases": 0, "comparison": NO_SEQUENCE,
                "truth_sequences_identical": bench.recovered([], truth_records)}
    scored = []
    for name, sequence in found:
        plain = sequence.replace("-", "")
        entry = {"returned_name": name, "returned_sequence": plain,
                 "returned_bp": len(plain),
                 "uncalled_bases": plain.upper().count("N"),
                 "lowercase_bases": sum(1 for letter in plain if letter.islower())}
        if truth_records:
            entry.update(bench.compare(plain, truth_records, blast_work))
            called = plain.upper()
            entry["matching_truth_copies"] = ";".join(
                name for name, truth in truth_records
                if called in (truth, reverse_complement(truth)))
        else:
            entry["comparison"] = ABSENT_BUT_RETURNED
        scored.append(entry)
    detail = "; ".join(f"{entry['returned_bp']} bp: {entry['comparison']}" for entry in scored)
    scored.sort(key=lambda entry: (bench.RELATIONS.index(entry["comparison"]),
                                   -int(entry.get("aligned_bp") or 0)))
    best = dict(scored[0])
    best.update(sequences_returned=len(found), all_returned_sequences=detail,
                truth_sequences_identical=bench.recovered(
                    [sequence for _, sequence in found], truth_records))
    return best


def flag_present(flags, name):
    return bool(re.search(rf"(^|[;|,\s]){name}\b", flags or ""))


def locus_recon_row(unit, work, truth_records, blast_work):
    paths = unit_paths(work, unit)
    run = load_json(paths["run"], {})
    report = read_table(paths["report"])[0] if paths["report"].is_file() else {}
    found = []
    if paths["candidate"].is_file():
        found = [(name, sequence) for name, sequence in bench.read_fasta(paths["candidate"])]
    row = score_sequences(found, truth_records, blast_work)
    disposition = report.get("result_disposition", "")
    flags = report.get("qc_flags", "")
    row.update(
        tool="Locus-Recon", configuration=f"{FROZEN_VERSION}, default options",
        exit_code=run.get("exit_code", ""), seconds=run.get("seconds", ""),
        decision="accepted" if disposition == "PASS" else "not accepted",
        decision_rule="disposition PASS",
        disposition=disposition, workflow_status=report.get("workflow_status", ""),
        tier=report.get("sequence_confidence", ""), qc_flags=flags,
        catalogue_status=report.get("catalogue_status", ""),
        nearest_allele_identity_pct=report.get("nearest_allele_identity_pct", ""),
        remap_mean_depth=report.get("remap_mean_depth", ""),
        uncertain_base_count=report.get("uncertain_base_count", ""),
        interior_uncertain_base_count=report.get("interior_uncertain_base_count", ""),
        candidate_mixed_sites=report.get("candidate_mixed_sites", ""),
        mixed_site_count=report.get("mixed_site_count", ""),
        strand_biased_sites=report.get("strand_biased_sites", ""),
        mixture_detected=report.get("mixture_detected", ""),
        span_clipped_bp=report.get("span_clipped_bp", ""),
        flag_STRAND_BIASED_ALT_SUPPORT=flag_present(flags, "STRAND_BIASED_ALT_SUPPORT"),
        flag_MIXED_SITE_BELOW_MIXTURE_CALL=flag_present(flags, "MIXED_SITE_BELOW_MIXTURE_CALL"),
        flag_ALLELE_MIXTURE=flag_present(flags, "ALLELE_MIXTURE"),
        tool_report=report.get("message", "") or "no report row",
    )
    return row, report


def ariba_row(unit, raw, timings, truth_records, blast_work):
    tag = "ariba_default"
    report, sequences = bench.ariba_results(raw, tag)
    locus = unit["locus"]
    found = sequences.get(locus, [])
    row = score_sequences(found, truth_records, blast_work)
    timing = timings.get(tag, {})
    entry = report.get(locus)
    flag = int(entry["first"]["flag"]) if entry else None
    reasons = []
    if len(found) != 1:
        reasons.append(f"{len(found)} sequences returned")
    if flag is None:
        reasons.append("locus absent from the report")
    else:
        reasons += [f"{name} not set" for name in ARIBA_MUST_BE_SET
                    if not flag & ARIBA_BIT[name]]
        reasons += [f"{name} set" for name in ARIBA_MUST_BE_UNSET if flag & ARIBA_BIT[name]]
        if unit["locus_type"] == "protein-coding" and not flag & ARIBA_BIT["complete_gene"]:
            reasons.append("complete_gene not set")
    row.update(
        tool="ARIBA", configuration="default options",
        exit_code=timing.get("exit_code", ""), seconds=timing.get("seconds", ""),
        decision="not accepted" if reasons else "accepted",
        decision_rule="protocol 2.5, ARIBA rule", decision_detail="; ".join(reasons),
        ariba_flag="" if flag is None else flag,
        tool_report=bench.ariba_report_text(entry) if entry else "locus absent from the report",
    )
    return row


def kma_row(unit, raw, timings, truth_records, blast_work):
    tag = "kma_default"
    report, sequences = bench.kma_results(raw, tag)
    locus = unit["locus"]
    found = sequences.get(locus, [])
    row = score_sequences(found, truth_records, blast_work)
    timing = timings.get(tag, {})
    entry = report.get(locus)
    reasons, coverage = [], ""
    if not found:
        reasons.append("no consensus returned")
    else:
        if row["lowercase_bases"]:
            reasons.append(f"{row['lowercase_bases']} lower-case bases")
        if row["uncalled_bases"]:
            reasons.append(f"{row['uncalled_bases']} n")
    secondary = "not accepted" if reasons else "accepted"
    if entry:
        coverage = entry.get("Template_Coverage", "")
    try:
        full = float(coverage) == 100.0
    except ValueError:
        full = False
    if not full:
        reasons.append(f"template coverage {coverage or 'not reported'}")
    row.update(
        tool="KMA", configuration="default options",
        exit_code=timing.get("exit_code", ""), seconds=timing.get("seconds", ""),
        decision="not accepted" if reasons else "accepted",
        decision_rule="protocol 2.5, KMA primary rule", decision_detail="; ".join(reasons),
        decision_secondary=secondary, kma_template_coverage=coverage,
        kma_template_identity=entry.get("Template_Identity", "") if entry else "",
        kma_depth=entry.get("Depth", "") if entry else "",
        tool_report=bench.kma_report_text(entry) if entry else "locus absent from the report",
    )
    return row


def map_position(candidate, truth, position):
    """Base of a truth copy that aligns to a 1-based position of the candidate."""
    if len(candidate) == len(truth):
        return truth[position - 1]
    offset = truth.find(candidate)
    if offset >= 0:
        return truth[offset + position - 1]
    try:
        from Bio import Align
    except ImportError:
        return "?"
    try:
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score, aligner.mismatch_score = 1, -1
        aligner.open_gap_score, aligner.extend_gap_score = -5, -1
        aligner.end_gap_score = 0
        alignment = aligner.align(candidate, truth)[0]
        for (cand_start, cand_end), (truth_start, _) in zip(*alignment.aligned):
            if cand_start < position <= cand_end:
                return truth[int(truth_start) + position - 1 - int(cand_start)]
    except Exception:  # the site is then reported as not mapped
        return "?"
    return "-"


def mixed_sites(unit, row, work, truth_records):
    """Credible mixed sites of a candidate, with the truth copies at each site."""
    paths = unit_paths(work, unit)
    if not paths["support"].is_file() or not row.get("returned_sequence"):
        return []
    candidate = row["returned_sequence"].upper()
    oriented = []
    for name, truth in truth_records:
        forward_distance = bench.edit_distance(candidate, truth)
        reverse_distance = bench.edit_distance(candidate, reverse_complement(truth))
        oriented.append((name, truth if forward_distance <= reverse_distance
                         else reverse_complement(truth)))
    sites = []
    for record in read_table(paths["support"]):
        if record.get("credible_mixed_site") != "true":
            continue
        position = int(record["position"])
        bases = [(name, map_position(candidate, truth, position)) for name, truth in oriented]
        distinct = {base for _, base in bases}
        sites.append({
            "species": unit["species"], "strain": unit["strain"], "locus": unit["locus"],
            "unit_set": unit["unit_set"], "stratum": unit["stratum"],
            "truth_copies": unit["truth_copies"], "disposition": row["disposition"],
            "tier": row["tier"], "decision": row["decision"], "exact": row["exact"],
            "flag_MIXED_SITE_BELOW_MIXTURE_CALL": row["flag_MIXED_SITE_BELOW_MIXTURE_CALL"],
            "flag_ALLELE_MIXTURE": row["flag_ALLELE_MIXTURE"],
            "position": position, "candidate_base": record["reference_base"],
            "alt_base": record["alt_base"], "alt_fraction": record["alt_fraction"],
            "depth": record["depth"], "alt_forward": record["alt_forward"],
            "alt_reverse": record["alt_reverse"],
            "truth_bases_at_site": ";".join(f"{name.rsplit('|', 1)[1]}:{base}"
                                            for name, base in bases),
            "truth_copies_differ_at_site": len(distinct) > 1,
            "alt_base_is_in_a_truth_copy": record["alt_base"].upper() in distinct,
        })
    return sites


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def copy_number_row(unit, work):
    """Outcomes of protocol 2.6 for one unit."""
    paths = unit_paths(work, unit)
    run = load_json(paths["copy_run"], {})
    # No copy number is defined when the bait aligns only in parts (fifth stratum).
    truth = unit["truth_copies"] if unit["truth_defined"] else ""
    row = {
        "species": unit["species"], "strain": unit["strain"], "locus": unit["locus"],
        "unit_set": unit["unit_set"], "stratum": unit["stratum"], "truth_copies": truth,
        "truth_distinct_sequences": unit.get("truth_distinct_sequences", ""),
        "locus_bp": unit["bait_bp"],
        "locus_gc_pct": unit["copies"][0]["gc_pct"] if unit["copies"] else unit["bait_gc_pct"],
        "discovery_intervals": len(run.get("regions", [])),
        "status": run.get("status", "not run"),
        "locus_regions": ";".join(run.get("regions", [])),
        "gene_spans": ";".join(run.get("gene_spans", [])),
    }
    depth = read_table(paths["depth_tsv"])[0] if paths["depth_tsv"].is_file() else {}
    combined = read_table(paths["copy_tsv"])[0] if paths["copy_tsv"].is_file() else {}
    source = combined or depth
    for key in ("ratio", "ratio_ci_low", "ratio_ci_high", "depth_call", "dosage_estimate",
                "dosage_status", "dosage_flags", "dosage_ci_low", "dosage_ci_high",
                "dosage_lower", "dosage_upper", "geometry_expected_dosage",
                "assembly_region_count", "ambiguity_index", "reliable",
                "locus_median_depth", "backbone_median_depth", "backbone_mad",
                "gene_overlap_bp", "gene_overlap_fraction", "query_coverage_fraction",
                "region_ratios", "notes"):
        row[key] = source.get(key, "")
    for key in ("graph_context_count", "graph_context_lower_bound", "graph_context_status",
                "graph_context_flags", "left_context_count", "right_context_count",
                "graph_bait_coverage_fraction", "copy_number_call", "copy_number_kind",
                "copy_number_method", "copy_number_status", "copy_number_lower_bound",
                "copy_number_upper_bound", "consensus_flags"):
        row[key] = combined.get(key, "")
    if depth and combined:
        row["depth_ratio_command_agrees_with_copy_number_command"] = all(
            depth.get(key, "") == combined.get(key, "")
            for key in ("ratio", "ratio_ci_low", "ratio_ci_high", "depth_call"))
    ratio, low, high = (number(row[key]) for key in ("ratio", "ratio_ci_low", "ratio_ci_high"))
    dosage = number(row["dosage_estimate"])
    dosage_low, dosage_high = number(row["dosage_ci_low"]), number(row["dosage_ci_high"])
    if ratio is not None and truth != "":
        row["ratio_error"] = round(ratio - truth, 4)
        row["ratio_relative_error"] = round((ratio - truth) / truth, 4) if truth else ""
        row["ratio_ci_width"] = round(high - low, 4)
        row["ratio_ci_contains_truth"] = low <= truth <= high
        row["called_multicopy_depth"] = row["depth_call"] == "MULTICOPY_DEPTH"
    if dosage is not None and truth != "":
        row["dosage_error"] = round(dosage - truth, 4)
        row["dosage_relative_error"] = round((dosage - truth) / truth, 4) if truth else ""
        if dosage_low is not None and dosage_high is not None:
            row["dosage_ci_contains_truth"] = dosage_low <= truth <= dosage_high
    count, bound = number(row["graph_context_count"]), number(row["graph_context_lower_bound"])
    if combined and truth != "":
        if count is not None:
            row["graph_against_truth"] = (
                "equal to the truth" if count == truth
                else "wrong, below the truth" if count < truth else "wrong, above the truth")
        elif bound is not None:
            row["graph_against_truth"] = ("correct lower bound" if bound <= truth
                                          else "wrong, lower bound above the truth")
        else:
            row["graph_against_truth"] = "no count and no bound"
        call, kind = number(row["copy_number_call"]), row["copy_number_kind"]
        lower = number(row["copy_number_lower_bound"])
        if kind == "INTEGER_CONTEXT_COUNT":
            correct = call == truth
        elif kind == "MEAN_DEPTH_DOSAGE":
            # A continuous dosage is counted as correct when it rounds to the truth.
            correct = call is not None and math.floor(call + 0.5) == truth
        elif kind == "LOWER_BOUND":
            correct = lower is not None and lower <= truth
        else:
            correct = None
        row["consensus_against_truth"] = (
            "no value" if correct is None else "correct" if correct else "wrong")
        if call is not None:
            row["consensus_error"] = round(call - truth, 4)
        row["supported_but_wrong"] = (
            row["copy_number_status"] == "SUPPORTED" and correct is False)
    return row


def block_size_rows(unit, work, copy_row):
    paths = unit_paths(work, unit)
    rows = []
    truth = unit["truth_copies"] if unit["truth_defined"] else ""
    for record in load_json(paths["block_json"], []):
        low, high = record["ratio_ci_low"], record["ratio_ci_high"]
        row = {
            "species": unit["species"], "strain": unit["strain"], "locus": unit["locus"],
            "unit_set": unit["unit_set"], "stratum": unit["stratum"], "truth_copies": truth,
            "block_bp": record["bootstrap_block_size"], "ratio": record["ratio"],
            "ratio_ci_low": low, "ratio_ci_high": high,
            "ratio_ci_width": round(high - low, 4),
            "ratio_ci_contains_truth": "" if truth == "" else low <= truth <= high,
            "depth_call": record["depth_call"],
            "dosage_estimate": record["dosage_estimate"],
            "dosage_ci_low": record["dosage_ci_low"], "dosage_ci_high": record["dosage_ci_high"],
            "locus_blocks": record["bootstrap_locus_blocks"],
            "backbone_blocks_total": record["bootstrap_backbone_blocks_total"],
            "backbone_blocks_used": record["bootstrap_backbone_blocks_used"],
            "resamples": record["bootstrap_resamples"], "seed": record["bootstrap_seed"],
            "assembly_region_count": record["assembly_region_count"],
            "seconds": record["seconds"],
        }
        if record["dosage_ci_low"] is not None and record["dosage_ci_high"] is not None:
            row["dosage_ci_width"] = round(record["dosage_ci_high"] - record["dosage_ci_low"], 4)
            row["dosage_ci_contains_truth"] = "" if truth == "" else (
                record["dosage_ci_low"] <= truth <= record["dosage_ci_high"])
        if record["bootstrap_block_size"] == 500:
            command = (number(copy_row.get("ratio")), number(copy_row.get("ratio_ci_low")),
                       number(copy_row.get("ratio_ci_high")))
            row["agrees_with_the_command_at_500_bp"] = all(
                value is not None and abs(value - other) < 5e-7
                for value, other in zip(command, (record["ratio"], low, high)))
        rows.append(row)
    return rows


def isolate_interval(pairs):
    """Percentile interval of sum(numerator)/sum(denominator), resampling isolates."""
    pairs = list(pairs)
    generator = random.Random(SEED)
    values = []
    for _ in range(BOOTSTRAP_RESAMPLES):
        numerator = denominator = 0
        for _ in pairs:
            first, second = pairs[generator.randrange(len(pairs))]
            numerator += first
            denominator += second
        if denominator:
            values.append(numerator / denominator)
    if not values:
        return "", ""
    values.sort()
    return (round(values[int(0.025 * len(values))], 4),
            round(values[min(len(values) - 1, int(0.975 * len(values)))], 4))


class Summary:
    """Rows of the summary table; proportions carry an interval over isolates."""

    def __init__(self, strains):
        self.strains = strains
        self.rows = []

    def proportion(self, metric, unit_set, group, tool, rows, numerator, source, note=""):
        rows = list(rows)
        hits = [row for row in rows if numerator(row)]
        per_isolate = [(sum(1 for row in hits if row["strain"] == strain),
                        sum(1 for row in rows if row["strain"] == strain))
                       for strain in self.strains]
        low, high = isolate_interval(per_isolate) if rows else ("", "")
        self.rows.append({
            "metric": metric, "unit_set": unit_set, "group": group, "tool": tool,
            "numerator": len(hits), "denominator": len(rows),
            "value": round(len(hits) / len(rows), 4) if rows else "",
            "interval_low": low, "interval_high": high,
            "interval": f"95% percentile, {BOOTSTRAP_RESAMPLES} resamples of isolates, "
                        f"seed {SEED}" if rows else "",
            "isolates_with_units": sum(1 for _, total in per_isolate if total),
            "source_table": source, "note": note,
        })

    def value(self, metric, unit_set, group, tool, value, count, source, note=""):
        self.rows.append({
            "metric": metric, "unit_set": unit_set, "group": group, "tool": tool,
            "numerator": "", "denominator": count, "value": value,
            "source_table": source, "note": note,
        })


def unit_sets(rows):
    """The protocol set, the supplementary set, both together, and a set added afterwards.

    The fourth set holds the protocol units of the isolates whose draft agrees
    with the complete genome.  It was not foreseen by the protocol and exists
    only when an isolate fails that comparison.
    """
    rows = list(rows)
    sets = [
        ("protocol", [row for row in rows if row["unit_set"] == SET_PROTOCOL]),
        ("supplementary", [row for row in rows if row["unit_set"] == SET_SUPPLEMENTARY]),
        ("protocol and supplementary", [row for row in rows
                                        if row["unit_set"] != SET_UNAVAILABLE]),
    ]
    if any(row.get("reads_agree_with_genome") is False for row in rows):
        sets.append((SET_CONCORDANT, [
            row for row in rows if row["unit_set"] == SET_PROTOCOL
            and row.get("reads_agree_with_genome") is not False]))
    return tuple(sets)


def decisions_of(row):
    """(rule label, accepted) for every decision rule that applies to a row."""
    if row["tool"] == "Locus-Recon":
        return [("Locus-Recon, PASS", row["decision"] == "accepted")]
    if row["tool"] == "ARIBA":
        return [("ARIBA, protocol rule", row["decision"] == "accepted")]
    return [("KMA, primary rule", row["decision"] == "accepted"),
            ("KMA, secondary rule", row["decision_secondary"] == "accepted")]


def decision_tables(unit_rows):
    table = []
    for set_name, rows in unit_sets(unit_rows):
        for stratum in ("all strata with a copy",) + STRATA:
            if stratum == "all strata with a copy":
                chosen = [row for row in rows if row["stratum"] in STRATA[:3]]
            else:
                chosen = [row for row in rows if row["stratum"] == stratum]
            cells = {}
            for row in chosen:
                for rule, accepted in decisions_of(row):
                    cell = cells.setdefault(rule, {
                        "units": 0, "accepted_exact": 0, "accepted_not_exact": 0,
                        "not_accepted_exact": 0, "not_accepted_not_exact": 0,
                        "of_which_no_sequence": 0})
                    cell["units"] += 1
                    exact = row["exact"] == "yes"
                    key = (("accepted_" if accepted else "not_accepted_")
                           + ("exact" if exact else "not_exact"))
                    cell[key] += 1
                    if not accepted and row["comparison"] == NO_SEQUENCE:
                        cell["of_which_no_sequence"] += 1
            for rule in ("Locus-Recon, PASS", "KMA, primary rule", "KMA, secondary rule",
                         "ARIBA, protocol rule"):
                cell = cells.get(rule)
                if not cell:
                    continue
                accepted = cell["accepted_exact"] + cell["accepted_not_exact"]
                rejected = cell["not_accepted_exact"] + cell["not_accepted_not_exact"]
                entry = {
                    "unit_set": set_name, "stratum": stratum, "decision_rule": rule, **cell,
                    "accepted": accepted, "not_accepted": rejected,
                    "exact": cell["accepted_exact"] + cell["not_accepted_exact"],
                    "proportion_accepted": round(accepted / cell["units"], 4),
                    "proportion_exact_among_accepted":
                        round(cell["accepted_exact"] / accepted, 4) if accepted else "",
                    "proportion_exact_among_not_accepted":
                        round(cell["not_accepted_exact"] / rejected, 4) if rejected else "",
                }
                if stratum in STRATA[3:]:
                    # Without a truth sequence a result is neither exact nor inexact.
                    for key in ("accepted_exact", "accepted_not_exact", "not_accepted_exact",
                                "not_accepted_not_exact", "exact",
                                "proportion_exact_among_accepted",
                                "proportion_exact_among_not_accepted"):
                        entry[key] = ""
                    entry["note"] = "no truth sequence in this stratum; exactness is not defined"
                table.append(entry)
    return table


def summarise(unit_rows, copy_rows, block_rows, site_rows, strains):
    summary = Summary(strains)
    units_table, copy_table = "external_units.tsv", "external_copy_number.tsv"
    every = [row for row in unit_rows if row["tool"] == "Locus-Recon"]
    for set_name in (SET_PROTOCOL, SET_SUPPLEMENTARY, SET_UNAVAILABLE):
        chosen = [row for row in every if row["unit_set"] == set_name]
        summary.value("units", set_name, "all", "", len(chosen), len(every), units_table,
                      "units of isolate and locus; the protocol set is the one found by the "
                      "bait rule of protocol 2.2")
        for stratum in STRATA:
            count = sum(1 for row in chosen if row["stratum"] == stratum)
            if count:
                summary.value("units", set_name, stratum, "", count, len(chosen), units_table)
    for set_name, rows in unit_sets(unit_rows):
        with_copy = [row for row in rows if row["stratum"] in STRATA[:3]]
        for tool in TOOLS:
            mine = [row for row in with_copy if row["tool"] == tool]
            summary.proportion("sequence returned", set_name, "units with a copy", tool, mine,
                               lambda row: row["comparison"] != NO_SEQUENCE, units_table)
            summary.proportion("sequence identical to a truth copy", set_name,
                               "units with a copy", tool, mine,
                               lambda row: row["exact"] == "yes", units_table)
            for stratum in STRATA[:3]:
                summary.proportion(
                    "sequence identical to a truth copy", set_name, stratum, tool,
                    [row for row in mine if row["stratum"] == stratum],
                    lambda row: row["exact"] == "yes", units_table)
            for rule in [rule for rule, _ in decisions_of({"tool": tool, "decision": "",
                                                           "decision_secondary": ""})]:
                def accepted(row, rule=rule):
                    return dict(decisions_of(row))[rule]
                summary.proportion("accepted", set_name, "units with a copy", rule, mine,
                                   accepted, units_table)
                summary.proportion(
                    "accepted sequences identical to a truth copy", set_name,
                    "units with a copy", rule, [row for row in mine if accepted(row)],
                    lambda row: row["exact"] == "yes", units_table)
                summary.proportion(
                    "not accepted although identical to a truth copy", set_name,
                    "units with a copy", rule, [row for row in mine if not accepted(row)],
                    lambda row: row["exact"] == "yes", units_table)
                for stratum in STRATA[:3]:
                    summary.proportion(
                        "accepted", set_name, stratum, rule,
                        [row for row in mine if row["stratum"] == stratum], accepted,
                        units_table)
        recon = [row for row in with_copy if row["tool"] == "Locus-Recon"]
        rejected = [row for row in recon if row["decision"] != "accepted"]
        for label in ("exact", "incomplete, identical over its length",
                      "wrong, extra or missing bases", "no sequence"):
            summary.proportion("not accepted by Locus-Recon, by result", set_name, label,
                               "Locus-Recon", rejected,
                               lambda row, label=label: row["result_class"] == label,
                               units_table)
        for disposition in ("PASS", "REVIEW", "HOLD"):
            summary.proportion("disposition", set_name, disposition, "Locus-Recon", recon,
                               lambda row, disposition=disposition:
                               row["disposition"] == disposition, units_table)
        one_copy = [row for row in recon if row["stratum"] == STRATA[0]]
        summary.proportion("STRAND_BIASED_ALT_SUPPORT", set_name, STRATA[0], "Locus-Recon",
                           [row for row in one_copy if row["comparison"] != NO_SEQUENCE],
                           lambda row: row["flag_STRAND_BIASED_ALT_SUPPORT"] is True,
                           units_table, "one-copy units with a candidate")
        accepted_rows = [row for row in recon if row["decision"] == "accepted"]
        summary.proportion("MIXED_SITE_BELOW_MIXTURE_CALL among accepted sequences", set_name,
                           "units with a copy", "Locus-Recon", accepted_rows,
                           lambda row: row["flag_MIXED_SITE_BELOW_MIXTURE_CALL"] is True,
                           units_table)
        for stratum, label in ((STRATA[3], "sequence returned although the locus is absent"),
                               (STRATA[4], "sequence returned, truth not defined")):
            chosen = [row for row in rows if row["stratum"] == stratum]
            for tool in TOOLS:
                mine = [row for row in chosen if row["tool"] == tool]
                summary.proportion(label, set_name, stratum, tool, mine,
                                   lambda row: row["comparison"] != NO_SEQUENCE, units_table)
                for rule in [rule for rule, _ in decisions_of(
                        {"tool": tool, "decision": "", "decision_secondary": ""})]:
                    summary.proportion(
                        "accepted", set_name, stratum, rule, mine,
                        lambda row, rule=rule: dict(decisions_of(row))[rule], units_table)
                summary.proportion(
                    "returned sequence occurs in the complete genome", set_name, stratum, tool,
                    [row for row in mine if row["comparison"] != NO_SEQUENCE],
                    lambda row: row.get("returned_occurs_in_complete_genome") is True,
                    units_table, "not part of the protocol")

    for set_name, rows in unit_sets(copy_rows):
        evaluated = [row for row in rows if row.get("ratio") not in ("", None)
                     and row["truth_copies"] not in ("", 0)]
        summary.proportion("depth ratio obtained", set_name, "units with a copy",
                           "locus-recon-depth-ratio",
                           [row for row in rows if row["truth_copies"] not in ("", 0)],
                           lambda row: row.get("ratio") not in ("", None), copy_table)
        single = [row for row in evaluated if row["truth_copies"] == 1]
        multiple = [row for row in evaluated if row["truth_copies"] >= 2]
        summary.proportion("called MULTICOPY_DEPTH", set_name, "truth: one copy",
                           "locus-recon-depth-ratio", single,
                           lambda row: row["depth_call"] == "MULTICOPY_DEPTH", copy_table,
                           "false alerts of the 1.5 threshold")
        summary.proportion("called MULTICOPY_DEPTH", set_name, "truth: two or more copies",
                           "locus-recon-depth-ratio", multiple,
                           lambda row: row["depth_call"] == "MULTICOPY_DEPTH", copy_table)
        for copies in sorted({row["truth_copies"] for row in evaluated}):
            group = [row for row in evaluated if row["truth_copies"] == copies]
            label = f"truth: {copies} cop{'y' if copies == 1 else 'ies'}"
            summary.proportion("called MULTICOPY_DEPTH", set_name, label,
                               "locus-recon-depth-ratio", group,
                               lambda row: row["depth_call"] == "MULTICOPY_DEPTH", copy_table)
            for key, name in (("ratio", "depth ratio"), ("dosage_estimate", "dosage")):
                values = [number(row[key]) for row in group if number(row[key]) is not None]
                if not values:
                    continue
                summary.value(f"median {name}", set_name, label, "locus-recon-depth-ratio",
                              round(statistics.median(values), 4), len(values), copy_table,
                              f"range {min(values):.4g} to {max(values):.4g}")
                errors = [abs(value - copies) for value in values]
                summary.value(f"median absolute error of the {name}", set_name, label,
                              "locus-recon-depth-ratio", round(statistics.median(errors), 4),
                              len(values), copy_table, f"largest {max(errors):.4g}")
        for label, group in (("truth: one copy", single),
                             ("truth: two or more copies", multiple),
                             ("units with a copy", evaluated)):
            summary.proportion("interval of the depth ratio contains the truth", set_name,
                               label, "locus-recon-depth-ratio", group,
                               lambda row: row.get("ratio_ci_contains_truth") is True,
                               copy_table)
            with_dosage = [row for row in group
                           if row.get("dosage_ci_contains_truth") in (True, False)]
            summary.proportion("interval of the dosage contains the truth", set_name, label,
                               "locus-recon-depth-ratio", with_dosage,
                               lambda row: row["dosage_ci_contains_truth"] is True,
                               copy_table)
            graph = [row for row in group if row.get("graph_against_truth")]
            for outcome in ("equal to the truth", "correct lower bound",
                            "wrong, below the truth", "wrong, above the truth",
                            "wrong, lower bound above the truth", "no count and no bound"):
                summary.proportion("graph context count against the truth", set_name,
                                   f"{label}; {outcome}", "locus-recon-copy-number", graph,
                                   lambda row, outcome=outcome:
                                   row["graph_against_truth"] == outcome, copy_table)
            consensus = [row for row in group if row.get("consensus_against_truth")]
            for status in sorted({row["copy_number_status"] for row in consensus}):
                chosen = [row for row in consensus if row["copy_number_status"] == status]
                summary.proportion("consensus copy number correct", set_name,
                                   f"{label}; status {status}", "locus-recon-copy-number",
                                   chosen,
                                   lambda row: row["consensus_against_truth"] == "correct",
                                   copy_table,
                                   "a continuous dosage is correct when it rounds to the truth")
            summary.proportion("consensus SUPPORTED and wrong", set_name, label,
                               "locus-recon-copy-number", consensus,
                               lambda row: row.get("supported_but_wrong") is True, copy_table)

    for set_name, rows in unit_sets(block_rows):
        for block in BLOCK_SIZES:
            group = [row for row in rows if row["block_bp"] == block
                     and row["truth_copies"] not in ("", 0)]
            if not group:
                continue
            widths = [row["ratio_ci_width"] for row in group]
            summary.value("median width of the interval of the depth ratio", set_name,
                          f"block {block} bp", "estimate_locus_copy_number",
                          round(statistics.median(widths), 4), len(widths),
                          "block_size_sensitivity.tsv",
                          f"range {min(widths):.4g} to {max(widths):.4g}")
            for label, chosen in (
                    ("all", group),
                    ("truth: one copy", [row for row in group if row["truth_copies"] == 1]),
                    ("truth: two or more copies",
                     [row for row in group if row["truth_copies"] >= 2])):
                summary.proportion("interval of the depth ratio contains the truth", set_name,
                                   f"block {block} bp; {label}", "estimate_locus_copy_number",
                                   chosen, lambda row: row["ratio_ci_contains_truth"] is True,
                                   "block_size_sensitivity.tsv")
    for set_name, rows in unit_sets(site_rows):
        accepted_sites = [row for row in rows if row["decision"] == "accepted"]
        summary.proportion("mixed sites of accepted sequences at which the truth copies differ",
                           set_name, "sites", "Locus-Recon", accepted_sites,
                           lambda row: row["truth_copies_differ_at_site"] is True,
                           "mixed_sites.tsv")
    return summary.rows


def counts_by(unit_rows, key):
    """Counts per isolate or per locus, for each decision rule (protocol 2.4)."""
    table = []
    for set_name, rows in unit_sets(unit_rows):
        values = []
        for row in rows:
            if row[key] not in values:
                values.append(row[key])
        for value in values:
            for tool in TOOLS:
                mine = [row for row in rows if row[key] == value and row["tool"] == tool]
                with_copy = [row for row in mine if row["stratum"] in STRATA[:3]]
                for rule in [rule for rule, _ in decisions_of(
                        {"tool": tool, "decision": "", "decision_secondary": ""})]:
                    accepted = [row for row in with_copy if dict(decisions_of(row))[rule]]
                    rejected = [row for row in with_copy if not dict(decisions_of(row))[rule]]
                    table.append({
                        key: value, "unit_set": set_name, "decision_rule": rule,
                        "units": len(mine), "units_with_a_copy": len(with_copy),
                        "sequence_returned": sum(
                            1 for row in with_copy if row["comparison"] != NO_SEQUENCE),
                        "exact": sum(1 for row in with_copy if row["exact"] == "yes"),
                        "accepted": len(accepted),
                        "accepted_exact": sum(1 for row in accepted if row["exact"] == "yes"),
                        "accepted_not_exact": sum(
                            1 for row in accepted if row["exact"] != "yes"),
                        "not_accepted_exact": sum(
                            1 for row in rejected if row["exact"] == "yes"),
                        "not_accepted_not_exact": sum(
                            1 for row in rejected if row["exact"] != "yes"),
                        "units_without_a_truth_sequence": len(mine) - len(with_copy),
                    })
    return table


FACTOR_BINS = {
    "locus_gc_pct": ((0, 40, "below 40%"), (40, 55, "40 to 55%"), (55, 101, "55% or more")),
    "genome_gc_pct": ((0, 40, "below 40%"), (40, 55, "40 to 55%"), (55, 101, "55% or more")),
    "locus_bp": ((0, 1500, "below 1,500 bp"), (1500, 2500, "1,500 to 2,499 bp"),
                 (2500, 10 ** 9, "2,500 bp or more")),
    "discovery_intervals": ((1, 2, "1"), (2, 3, "2"), (3, 10 ** 9, "3 or more")),
}


def error_by_factor(copy_rows):
    """Relative error of the depth ratio and of the dosage by the factors of protocol 2.6."""
    table = []
    for set_name, rows in unit_sets(copy_rows):
        rows = [row for row in rows if row.get("ratio_relative_error") not in ("", None)]
        for copies, group in (("one copy", [row for row in rows if row["truth_copies"] == 1]),
                              ("two or more copies",
                               [row for row in rows if row["truth_copies"] >= 2])):
            for factor, bins in FACTOR_BINS.items():
                for low, high, label in bins:
                    chosen = [row for row in group if low <= float(row[factor]) < high]
                    if not chosen:
                        continue
                    entry = {"unit_set": set_name, "truth": copies, "factor": factor,
                             "bin": label, "units": len(chosen),
                             "isolates": len({row["strain"] for row in chosen})}
                    for key, name in (("ratio_relative_error", "ratio"),
                                      ("dosage_relative_error", "dosage")):
                        values = [row[key] for row in chosen if row.get(key) not in ("", None)]
                        if values:
                            entry[f"{name}_median_relative_error"] = round(
                                statistics.median(values), 4)
                            entry[f"{name}_median_absolute_relative_error"] = round(
                                statistics.median(abs(value) for value in values), 4)
                            entry[f"{name}_largest_absolute_relative_error"] = round(
                                max(abs(value) for value in values), 4)
                    entry["called_multicopy_depth"] = sum(
                        1 for row in chosen if row["depth_call"] == "MULTICOPY_DEPTH")
                    table.append(entry)
    return table


def headline(rows):
    """The rows of the summary that answer the questions of protocol 2.4 to 2.6."""
    detailed = ("graph context count against the truth", "consensus copy number correct",
                "consensus SUPPORTED and wrong", "interval of the dosage contains the truth")
    chosen = []
    for row in rows:
        group = row["group"]
        if row["unit_set"] == "supplementary" or row["denominator"] in ("", 0):
            continue
        if row["metric"] == "units" or re.match(r"truth: \d+ cop", group):
            if row["metric"] != "units":
                continue
        elif row["metric"] in detailed and not group.startswith("units with a copy"):
            continue
        elif row["metric"].startswith("median") and "block" not in group:
            continue
        elif "block" in group and not group.endswith("; all") and "median" not in row["metric"]:
            continue
        chosen.append(row)
    return chosen


def stage_score(args, isolates, work, outdir):
    package = check_frozen_package()
    outdir.mkdir(parents=True, exist_ok=True)
    units = load_units(work)
    isolate_records = load_isolates(work)
    loci = selected_loci(args)
    blast_work = work / "score_blast"
    strains = [entry["strain"] for entry in isolates]

    # A stage that did not run must not be read as a tool that returned nothing.
    missing = []
    for entry in isolates:
        raw = comparator_dir(work) / "raw" / entry["strain"] / entry["strain"]
        for tool, configuration, _ in bench.CONFIGURATIONS:
            if not (raw / f"{tool.lower()}_{configuration}.run.json").is_file():
                missing.append(f"{entry['strain']}: {tool} was not run")
        for locus in loci:
            unit = units[unit_key(entry["strain"], locus)]
            if unit["unit_set"] == SET_UNAVAILABLE:
                continue
            paths = unit_paths(work, unit)
            if locus == loci[0] and not (
                    work / "concordance" / f"{entry['strain']}.json").is_file():
                missing.append(f"{entry['strain']}: draft was not compared with the genome")
            for key, label in (("run", "Locus-Recon"), ("copy_run", "copy number")):
                if not paths[key].is_file():
                    missing.append(f"{entry['strain']} {locus}: {label} was not run")
            if (load_json(paths["copy_run"], {}).get("status") == "run"
                    and not paths["block_json"].is_file()):
                missing.append(f"{entry['strain']} {locus}: block sizes were not run")
    if missing and not args.allow_incomplete:
        raise RuntimeError("incomplete run:\n  " + "\n  ".join(missing))
    write_table(outdir / "not_run.tsv", [{"item": item} for item in missing], header=["item"])

    unit_rows, report_rows, copy_rows, block_rows, site_rows = [], [], [], [], []
    candidates, returned, truth_fasta, bait_fasta = [], [], [], []
    bait_rows, truth_rows, input_rows, draft_rows = [], [], [], []
    evidence = outdir / "per_unit"

    for entry in isolates:
        strain = entry["strain"]
        record = isolate_records[strain]
        reads = record["reads"]
        spades = load_json(isolate_dir(work, strain) / "spades.run.json", {})
        for mate in ("1", "2"):
            source = reads["files"][mate]
            kept = reads.get("subsampled_files", {}).get(mate, {})
            input_rows.append({
                "species": entry["species"], "strain": strain,
                "illumina_run": entry["illumina_run"], "mate": mate,
                "instrument": reads["instrument"], "url": source["url"],
                "ena_md5": source["ena_md5"], "local_md5": source["local_md5"],
                "sha256": source["sha256"], "size_bytes": source["size_bytes"],
                "reads": source["records"], "bases": source["bases"],
                "ena_read_count": reads["ena_read_count"],
                "ena_base_count": reads["ena_base_count"],
                "genome_bp_table": record["genome_bp_table"],
                "genome_bp_observed": record["genome_bp_observed"],
                "depth_full_x": reads["depth_full_x"], "subsampled": reads["subsampled"],
                "seed": reads["seed"], "pairs_kept": reads["pairs_kept"],
                "depth_used_x": reads["depth_used_x"],
                "subsampled_file": kept.get("file", ""),
                "subsampled_sha256": kept.get("sha256", ""),
                "subsampled_reads": kept.get("records", ""),
                "subsampled_bases": kept.get("bases", ""),
                "subsample_reproduced_by_a_second_run": reads.get("subsample_reproduced", ""),
            })
        draft_rows.append({
            "species": entry["species"], "strain": strain,
            "genome_bp": record["genome_bp_observed"],
            "genome_replicons": record["genome_replicons"],
            "genome_gc_pct": record["genome_gc_pct"],
            "truth_genome_file": record["genome_archive"],
            "truth_genome_sha256": record["genome_archive_sha256"],
            "draft_contigs": spades.get("contigs", ""), "draft_bp": spades.get("total_bp", ""),
            "draft_n50_bp": spades.get("n50_bp", ""),
            "draft_longest_bp": spades.get("longest_bp", ""),
            "draft_sha256": spades.get("draft_sha256", ""),
            "graph_sha256": spades.get("graph_sha256", ""),
            "spades_command": " ".join(spades.get("command", [])),
            "spades_seconds": spades.get("seconds", ""),
        })
        raw = comparator_dir(work) / "raw" / strain / strain
        timings = {}
        genome = read_genome(record["genome_fasta"])
        agreement = load_json(work / "concordance" / f"{strain}.json", {})
        agrees = agreement.get("reads_agree_with_genome", "")
        draft_rows[-1].update({
            "genome_aligned_by_draft_pct": agreement.get("genome_aligned_pct", ""),
            "draft_aligned_to_genome_pct": agreement.get("draft_aligned_pct", ""),
            "identity_1_to_1_pct": agreement.get("identity_1_to_1_pct", ""),
            "snps": agreement.get("snps", ""), "indels": agreement.get("indels", ""),
            "reads_agree_with_genome": agrees,
        })
        for tool, configuration, _ in bench.CONFIGURATIONS:
            tag = f"{tool.lower()}_{configuration}"
            timings[tag] = load_json(raw / f"{tag}.run.json", {})
        if raw.is_dir():
            target = outdir / "comparator_raw" / strain
            target.mkdir(parents=True, exist_ok=True)
            for path in sorted(raw.iterdir()):
                if path.is_file() and not path.name.endswith(".assemblies.fa.gz"):
                    shutil.copy2(path, target / path.name)
            for name in ("catalogue.fasta", "catalogue_ariba_metadata.tsv", "catalogue.json",
                         "declared_noncoding_for_ariba.txt"):
                source = comparator_dir(work) / strain / name
                if source.is_file():
                    shutil.copy2(source, target / name)

        for locus in loci:
            unit = units[unit_key(strain, locus)]
            base = {key: unit.get(key, "") for key in (
                "species", "strain", "isolate_assembly", "illumina_run", "locus", "locus_type",
                "unit_set", "bait_rule", "bait_bp", "truth_copies", "truth_distinct_sequences",
                "stratum", "bait_fraction_aligned", "bait_fraction_in_longest_alignment")}
            base["genome_gc_pct"] = record["genome_gc_pct"]
            base["reads_agree_with_genome"] = agrees
            bait_rows.append({key: unit.get(key, "") for key in (
                "species", "strain", "locus", "locus_type", "unit_set",
                "bait_reference_assembly", "bait_rule", "bait_candidates_in_annotation",
                "bait_source_record", "bait_gene_symbol", "bait_product", "bait_locus_tag",
                "bait_location", "bait_bp", "bait_gc_pct", "bait_pseudo", "bait_partial",
                "bait_sha256", "truth_copies", "truth_distinct_sequences",
                "truth_copies_megablast", "stratum", "bait_fraction_aligned",
                "bait_fraction_in_longest_alignment")})
            if unit["unit_set"] == SET_UNAVAILABLE:
                for tool in TOOLS:
                    unit_rows.append(dict(base, tool=tool, comparison="not available",
                                          result_class="not available", exact="",
                                          decision="not available"))
                continue
            copies = unit["copies"]
            truth_records = [(copy["name"], copy["sequence"]) for copy in copies]
            truth_fasta += truth_records
            bait_fasta.append((f"{strain}|{locus}|bait {unit['bait_source_record']}",
                               bench.read_fasta(unit["bait_path"])[0][1]))
            for copy in copies:
                truth_rows.append({
                    "species": unit["species"], "strain": strain, "locus": locus,
                    "unit_set": unit["unit_set"], "copy": copy["name"],
                    "replicon": copy["replicon"], "start": copy["start"], "end": copy["end"],
                    "strand": copy["strand"], "bp": len(copy["sequence"]),
                    "gc_pct": copy["gc_pct"],
                    "alignment_identity_pct": copy["alignment_identity_pct"],
                    "alignment_bp": copy["alignment_bp"], "bait_start": copy["bait_start"],
                    "bait_end": copy["bait_end"],
                    "extension_clipped_bp": copy["extension_clipped_bp"],
                    "identical_to_copy1": copy["identical_to_copy1"],
                    "sha256": hashlib.sha256(copy["sequence"].encode()).hexdigest(),
                })

            row, report = locus_recon_row(unit, work, truth_records, blast_work)
            rows = [row,
                    ariba_row(unit, raw, timings, truth_records, blast_work),
                    kma_row(unit, raw, timings, truth_records, blast_work)]
            for tool_row in rows:
                if unit["stratum"] == STRATA[4] and tool_row["comparison"] == ABSENT_BUT_RETURNED:
                    tool_row["comparison"] = NOT_DEFINED
                tool_row["result_class"] = result_class(tool_row["comparison"], unit)
                tool_row["exact"] = ("not defined" if unit["stratum"] in STRATA[3:]
                                     else "yes" if tool_row["comparison"] == IDENTICAL else "no")
                sequence = tool_row.pop("returned_sequence", "")
                name = tool_row.pop("returned_name", "")
                if sequence:
                    # Not part of the protocol: whether the returned sequence occurs,
                    # base for base and without a gap, in the complete genome.
                    called = sequence.upper()
                    tool_row["returned_occurs_in_complete_genome"] = any(
                        called in replicon or reverse_complement(called) in replicon
                        for replicon in genome.values())
                    header = (f"{strain}|{locus}|{tool_row['tool']} "
                              f"decision={tool_row['decision'].replace(' ', '_')} "
                              f"result={tool_row['result_class'].replace(' ', '_')} {name}")
                    (candidates if tool_row["tool"] == "Locus-Recon" else returned).append(
                        (header, sequence))
                    tool_row["returned_sha256"] = hashlib.sha256(
                        sequence.upper().encode()).hexdigest()
                unit_rows.append(dict(base, **tool_row))
            row["returned_sequence"] = (
                bench.read_fasta(unit_paths(work, unit)["candidate"])[0][1]
                if unit_paths(work, unit)["candidate"].is_file() else "")
            for site in mixed_sites(unit, row, work, truth_records):
                site["reads_agree_with_genome"] = agrees
                site_rows.append(site)
            row.pop("returned_sequence", None)
            if report:
                report_rows.append(dict(
                    {"species": unit["species"], "strain": strain,
                     "unit_set": unit["unit_set"]}, **report))
            copy_row = copy_number_row(unit, work)
            copy_row["genome_gc_pct"] = record["genome_gc_pct"]
            copy_row["reads_agree_with_genome"] = agrees
            copy_rows.append(copy_row)
            for block_row in block_size_rows(unit, work, copy_row):
                block_row["reads_agree_with_genome"] = agrees
                block_rows.append(block_row)

            paths = unit_paths(work, unit)
            target = evidence / strain / locus
            target.mkdir(parents=True, exist_ok=True)
            for key in ("candidate", "qc_report", "report", "discovery", "depth_tsv",
                        "copy_json", "copy_tsv", "block_json", "run", "copy_run"):
                if paths[key].is_file():
                    shutil.copy2(paths[key], target / paths[key].name)
            if paths["support"].is_file():
                with open(paths["support"], "rb") as source, gzip.GzipFile(
                        target / "allele_remap.support.tsv.gz", "wb", mtime=0) as packed:
                    shutil.copyfileobj(source, packed)

    write_table(outdir / "external_units.tsv", unit_rows)
    write_table(outdir / "locus_recon_report_rows.tsv", report_rows)
    write_table(outdir / "external_copy_number.tsv", copy_rows)
    write_table(outdir / "block_size_sensitivity.tsv", block_rows)
    write_table(outdir / "mixed_sites.tsv", site_rows)
    write_table(outdir / "decision_tables.tsv", decision_tables(unit_rows))
    write_table(outdir / "baits.tsv", bait_rows)
    write_table(outdir / "truth_copies.tsv", truth_rows)
    write_table(outdir / "input_reads.tsv", input_rows)
    write_table(outdir / "isolates_and_drafts.tsv", draft_rows)
    write_table(outdir / "counts_per_isolate.tsv", counts_by(unit_rows, "strain"))
    write_table(outdir / "counts_per_locus.tsv", counts_by(unit_rows, "locus"))
    write_table(outdir / "copy_number_error_by_factor.tsv", error_by_factor(copy_rows))
    summary_rows = summarise(unit_rows, copy_rows, block_rows, site_rows, strains)
    write_table(outdir / "external_summary.tsv", summary_rows)
    columns = []
    for row in summary_rows:
        columns += [key for key in row if key not in columns]
    write_table(outdir / "external_headline.tsv", headline(summary_rows), header=columns)
    listed = [row for row in unit_rows if row["tool"] == "Locus-Recon"
              and row["decision"] == "accepted" and row["exact"] != "yes"]
    write_table(outdir / "accepted_not_identical.tsv", listed,
                header=list(unit_rows[0]) if not listed else None)
    wrong = [row for row in copy_rows if row.get("supported_but_wrong") is True]
    write_table(outdir / "copy_number_supported_but_wrong.tsv", wrong,
                header=list(copy_rows[0]) if not wrong and copy_rows else None)
    bench.write_fasta(outdir / "locus_recon_candidates.fasta", candidates)
    bench.write_fasta(outdir / "comparator_returned_sequences.fasta", returned)
    bench.write_fasta(outdir / "truth_sequences.fasta", truth_fasta)
    bench.write_fasta(outdir / "baits.fasta", bait_fasta)
    for name in ("PROTOCOL.md", "isolates.tsv", "isolate_selection_log.tsv"):
        shutil.copy2(EVALUATION_DIR / name, outdir / f"protocol_{name}")

    stage_record(work, "score", args, {"blastn": ["blastn", "-version"]}, {"package": package})
    stages = {path.stem: json.loads(path.read_text())
              for path in sorted((work / "stage_records").glob("*.json"))}
    provenance = build_provenance(
        workflow="external_evaluation",
        repo_root=REPO_ROOT,
        inputs=[EVALUATION_DIR / "PROTOCOL.md", EVALUATION_DIR / "isolates.tsv",
                EVALUATION_DIR / "isolate_selection_log.tsv", Path(__file__),
                VALIDATION / "run_comparator_benchmark.py",
                VALIDATION / "run_tier_calibration.py", VALIDATION / "common.py",
                outdir / "baits.fasta", outdir / "truth_sequences.fasta"],
        parameters={
            "subsampling_seed": SEED, "target_depth_x": TARGET_DEPTH, "loci": list(LOCI),
            "isolates": strains,
            "draft_assembler": "SPAdes --isolate, contigs.fasta",
            "truth": {"blast_task": "blastn", "evalue": "1e-10",
                      "min_identity_pct": TRUTH_MIN_IDENTITY,
                      "min_fraction_of_bait": TRUTH_MIN_BAIT_FRACTION},
            "locus_recon_options": "defaults; --threads; --noncoding-locus for 23S and 16S",
            "copy_number_options": "defaults; --locus and --gene-span from blast_step1.tsv; "
                                   "--gene-length; --gfa of the draft; --bait",
            "depth_ratio_seed_default": 20260817, "bootstrap_resamples_default": 1000,
            "block_sizes_bp": list(BLOCK_SIZES),
            "comparator_configurations": [
                {"tool": tool, "configuration": name, "extra_arguments": extra}
                for tool, name, extra in bench.CONFIGURATIONS],
            "isolate_resampling": {"resamples": BOOTSTRAP_RESAMPLES, "seed": SEED},
            "draft_against_complete_genome": {
                "tool": "MUMmer dnadiff", "min_identity_pct": CONCORDANCE_MIN_IDENTITY,
                "min_genome_aligned_pct": CONCORDANCE_MIN_GENOME_ALIGNED,
                "in_protocol": False},
            "fallback_bait_rules": {
                "paralogue_suffixes": list(FALLBACK_PARALOGUE_SUFFIXES),
                "product_names": {key: list(value) for key, value in FALLBACK_PRODUCT.items()}},
        },
        tools={"python": [sys.executable, "--version"], "blastn": ["blastn", "-version"],
               "locus-recon": ["locus-recon", "--version"]},
    )
    provenance["frozen_package"] = package
    provenance["stages"] = stages
    provenance["commands"] = {
        "subsample": f"seqtk sample -s {SEED} <mate file> <pairs> | gzip -n -c",
        "draft": f"spades.py --isolate -1 R1 -2 R2 -o DIR -t N -m {calibration.SPADES_MEMORY_GB}",
        "truth": "blastn -task blastn -query BAIT -subject COMPLETE_GENOME -evalue 1e-10",
        "locus-recon": "locus-recon --samplesheet TSV --main-output-dir DIR --bait BAIT "
                       "--locus NAME --threads N [--noncoding-locus]",
        "depth-ratio": "locus-recon-depth-ratio --bam full_map.sorted.bam "
                       "(--locus CONTIG:START-END --gene-span START-END)... "
                       "--gene-length BP --tsv OUT",
        "copy-number": "locus-recon-copy-number --bam full_map.sorted.bam "
                       "(--locus ... --gene-span ...)... --gene-length BP --gfa "
                       "assembly_graph_after_simplification.gfa --bait BAIT "
                       "--output-json OUT --output-tsv OUT",
        "ariba": "ariba prepareref -f CATALOGUE -m METADATA DIR; "
                 "ariba run --threads N DIR R1 R2 OUT",
        "kma": "kma index -i CATALOGUE -o DB; kma -ipe R1 R2 -o OUT -t_db DB -t N",
        "concordance": "dnadiff -p PREFIX COMPLETE_GENOME DRAFT_CONTIGS",
    }
    provenance["read_files"] = input_rows
    provenance["outputs"] = [
        {"path": path.name, "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in sorted(outdir.iterdir()) if path.is_file() and path.suffix in (
            ".tsv", ".fasta")]
    write_json(outdir / "PROVENANCE.json", provenance)
    print(f"[score] {len(unit_rows)} rows for {len(unit_rows) // 3} units in {outdir}",
          flush=True)


# ----------------------------------------------------------------------------- main

STAGES = ("prepare", "assemble", "reconstruct", "copy-number", "block-size", "comparators",
          "concordance", "score")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--workdir", required=True,
                        help="directory for downloads, drafts and tool outputs")
    parser.add_argument("--outdir", default=str(VALIDATION / "external-evaluation-results"))
    parser.add_argument("--stage", required=True, choices=STAGES)
    parser.add_argument("--manifest", default=str(EVALUATION_DIR / "isolates.tsv"))
    parser.add_argument("--isolates", default="",
                        help="comma-separated strain names; default: all")
    parser.add_argument("--loci", default="", help="comma-separated loci; default: all")
    parser.add_argument("--threads", type=int, default=12,
                        help="threads used at any one time")
    parser.add_argument("--threads-per-run", type=int, default=4,
                        help="threads of one Locus-Recon, ARIBA or KMA run")
    parser.add_argument("--assembly-jobs", type=int, default=2,
                        help="draft assemblies made at the same time")
    parser.add_argument("--assembly-threads", type=int, default=0,
                        help="threads of each draft assembly; default: threads / assembly jobs")
    parser.add_argument("--download-jobs", type=int, default=3,
                        help="isolates whose reads are fetched and subsampled at a time")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="score: write the tables although a stage did not finish; "
                             "what is missing is listed in not_run.tsv")
    parser.add_argument("--seqtk", default="seqtk")
    parser.add_argument("--verify-subsample", action="store_true",
                        help="prepare: subsample again and require identical files")
    args = parser.parse_args()

    work, outdir = Path(args.workdir).resolve(), Path(args.outdir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    isolates = selected(read_table(args.manifest), args)
    if not isolates:
        parser.error("no isolate selected")
    print(f"[{args.stage}] {len(isolates)} isolate(s), started {now()}", flush=True)
    if args.stage == "prepare":
        stage_prepare(args, isolates, work)
    elif args.stage == "assemble":
        stage_assemble(args, isolates, work)
    elif args.stage == "reconstruct":
        stage_reconstruct(args, isolates, work)
    elif args.stage == "copy-number":
        stage_copy_number(args, isolates, work)
    elif args.stage == "block-size":
        stage_block_size(args, isolates, work)
    elif args.stage == "comparators":
        stage_comparators(args, isolates, work)
    elif args.stage == "concordance":
        stage_concordance(args, isolates, work)
    else:
        stage_score(args, isolates, work, outdir)
    print(f"[{args.stage}] finished {now()}", flush=True)


if __name__ == "__main__":
    main()
