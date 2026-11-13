#!/usr/bin/env python3
"""Sequences returned by ARIBA and KMA for loci that Locus-Recon was tested on.

Locus-Recon, ARIBA and KMA share one output: a nucleotide sequence for a target
locus.  This script runs ARIBA and KMA on two data sets that already have a
Locus-Recon result and a known truth, and compares what each tool returns with
that truth.

- ``hpylori``: the five closed *Helicobacter pylori* genomes of the tier
  calibration (stages B and C of ``run_tier_calibration.py``), five loci each.
  The reference catalogue is the same one-allele bait per locus, from strain
  26695, and the truth is the annotation of each closed genome.
- ``mock``: the 13 simulated samples of ``validation/mock_dataset``, with the
  same three-allele catalogue.
- ``liba6656``: *Clostridioides difficile* LIBA-6656 (run ERR467623) with the
  111 *tcdB* alleles of ``validation/liba6656-gfa-rerun``.  No independent
  truth exists for this isolate.  The returned sequences are compared with the
  two analyst-built candidates described in the article, which are hypotheses.

Both tools are run with their default settings.  ARIBA is run a second time
with ``--nucmer_min_id 80``, because its default of 90% is close to the
identity of one locus (*cagA*) to its bait.  The result table also carries one
row per case for Locus-Recon, read from the tables already deposited under
``validation/``; Locus-Recon is not run here.

The comparison is limited to the sequence and to what each tool reports about
it.  Read-support tiers, the measurement of bases lost at a contig end, copy
number and graph evidence have no counterpart in the other tools and are not
compared here.

Requirements on ``PATH``: ``ariba``, ``kma`` and ``blastn``.  They are not part
of the Locus-Recon environment; ``README.md`` in the output directory gives the
environment that was used.

    python validation/run_comparator_benchmark.py --workdir /scratch/comparators

``--score-only`` recomputes the tables from the raw tool outputs kept under
``<workdir>/raw`` without running the tools again.
"""

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "validation"))

from locus_recon.io import reverse_complement  # noqa: E402
from validation.common import build_provenance, write_json  # noqa: E402
import run_tier_calibration as calibration  # noqa: E402

VALIDATION = REPO_ROOT / "validation"
CALIBRATION_DIR = VALIDATION / "tier-calibration"
MOCK_DIR = VALIDATION / "mock_dataset"
LIBA_BAITS = VALIDATION / "liba6656-gfa-rerun" / "tcdB_111_ungapped_baits.fasta"
LIBA_RUN = "ERR467623"
LIBA_CANDIDATES = VALIDATION / "comparator-benchmark" / "liba6656_analyst_built_candidates.fasta"
REFERENCE_KIND = {
    "hpylori": "annotation of the closed genome",
    "mock": "simulated allele",
    "liba6656": "analyst-built candidate, not an independent truth",
}

HPYLORI_LOCI = ("gyrB", "23S", "ureB", "glmM", "cagA")
NONCODING = {"23S"}
CONFIGURATIONS = (
    ("ARIBA", "default", []),
    ("ARIBA", "nucmer_min_id_80", ["--nucmer_min_id", "80"]),
    ("KMA", "default", []),
)
ARIBA_FLAGS = (
    "assembled", "assembled_into_one_contig", "region_assembled_twice",
    "complete_gene", "unique_contig", "scaffold_graph_bad", "assembly_fail",
    "variants_suggest_collapsed_repeat", "hit_both_strands", "has_variant",
    "ref_seq_choose_fail",
)
RELATIONS = (
    "identical",
    "uncalled bases, otherwise identical",
    "incomplete, identical over its own length",
    "truth contained, extra bases",
    "identical over the overlap, different boundaries",
    "differs",
    "no alignment to the truth",
    "sequence returned although the locus is absent",
)
MOCK_GROUPS = {
    "baseline": "one allele", "fragmented": "one allele", "depth_series": "one allele",
    "mixture": "mixture of two alleles", "paralog": "paralogue control",
    "negative": "negative control",
}
RESULT_HEADER = [
    "dataset", "sample", "locus", "tool", "configuration", "exit_code", "seconds",
    "sequences_returned", "returned_bp", "uncalled_bases", "lowercase_bases",
    "reference_kind", "truth_name", "truth_bp", "truth_copies",
    "truth_sequences_identical", "comparison", "identity_pct",
    "aligned_bp", "mismatches", "gap_columns", "length_difference_bp", "edit_distance",
    "also_identical_to", "all_returned_sequences", "tool_reference", "tool_report",
]


# --------------------------------------------------------------------------- inputs

def read_fasta(path):
    records, name, parts = [], None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                name, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def write_fasta(path, records):
    with open(path, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start:start + 70] + "\n")


def build_catalogue(dataset, workdir):
    """Catalogue FASTA with plain names, and the ARIBA metadata table."""
    if dataset == "hpylori":
        records = []
        for locus in HPYLORI_LOCI:
            bait = read_fasta(CALIBRATION_DIR / f"bait_{locus}_26695.fasta")
            if len(bait) != 1:
                raise RuntimeError(f"expected one bait record for {locus}")
            records.append((locus, bait[0][1].upper()))
        coding = {locus: locus not in NONCODING for locus in HPYLORI_LOCI}
    else:
        source = LIBA_BAITS if dataset == "liba6656" else MOCK_DIR / "mock_locus_alleles.fasta"
        # A dot in a name would be read as a field separator in ARIBA output.
        records = [(name.replace(".", "_"), seq.upper()) for name, seq in read_fasta(source)]
        coding = {name: True for name, _ in records}
    directory = Path(workdir) / dataset
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / "catalogue.fasta"
    write_fasta(fasta, records)
    relaxed = directory / "declared_noncoding_for_ariba.txt"
    if relaxed.is_file():
        for name in relaxed.read_text().split():
            coding[name] = False
    metadata = directory / "catalogue_ariba_metadata.tsv"
    with open(metadata, "w") as handle:
        for name, _ in records:
            handle.write(f"{name}\t{1 if coding[name] else 0}\t0\t.\t.\t{name}\n")
    return fasta, metadata, [name for name, _ in records]


def samples_of(dataset, cache, wanted=(), reads=True):
    """Sample name, read pair and truth records per locus.

    With ``reads=False`` the public read files are not fetched.
    """
    samples = []
    if dataset == "hpylori":
        for entry in calibration.read_tsv(CALIBRATION_DIR / "closed_genome_manifest.tsv"):
            if wanted and entry["strain"] not in wanted:
                continue
            r1, r2 = calibration._reads(entry["illumina_run"], cache) if reads else ("", "")
            accession = entry["assembly_accession"]
            cds = calibration._annotated_records(accession, cache, "cds")
            rna = calibration._annotated_records(accession, cache, "rna")
            truth = {}
            for locus in HPYLORI_LOCI:
                if locus in NONCODING:
                    found = calibration._select(rna, product="23S ribosomal RNA")
                else:
                    found = calibration._select(cds, gene=locus)
                truth[locus] = [
                    (f"{entry['strain']}|{locus}|copy{index + 1}", sequence.upper())
                    for index, (_, sequence) in enumerate(found)
                ]
            samples.append({"sample": entry["strain"], "r1": r1, "r2": r2,
                            "truth": truth, "other": {}})
        return samples
    if dataset == "liba6656":
        r1, r2 = calibration._reads(LIBA_RUN, cache) if reads else ("", "")
        candidates = [(name, sequence.upper()) for name, sequence in read_fasta(LIBA_CANDIDATES)]
        return [{"sample": "LIBA6656", "r1": r1, "r2": r2,
                 "truth": {"tcdB": candidates}, "other": {}}]
    alleles = dict(read_fasta(MOCK_DIR / "truth" / "truth_alleles.fasta"))
    alleles.update(dict(read_fasta(MOCK_DIR / "mock_locus_alleles.fasta")))
    manifest = json.loads((MOCK_DIR / "dataset_manifest.json").read_text())
    for entry in manifest["samples"]:
        if wanted and entry["sample_id"] not in wanted:
            continue
        truth = []
        if entry["truth_allele"]:
            truth = [(entry["truth_allele"], alleles[entry["truth_allele"]].upper())]
        if entry["secondary_allele"]:
            # A mixture holds two alleles, and both are truth sequences.
            truth.append((entry["secondary_allele"],
                          alleles[entry["secondary_allele"]].upper()))
        other = {}
        if entry["category"] == "paralog":
            other["mockLocus_paralog"] = alleles["mockLocus_paralog"].upper()
        samples.append({
            "sample": entry["sample_id"], "r1": str(MOCK_DIR / entry["r1"]),
            "r2": str(MOCK_DIR / entry["r2"]), "truth": {"mockLocus": truth},
            "other": other, "group": MOCK_GROUPS[entry["category"]],
        })
    return samples


# ---------------------------------------------------------------------------- tools

def run(command, log):
    started = time.time()
    with open(log, "w") as handle:
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT,
                                   check=False)
    return completed.returncode, round(time.time() - started, 1)


def prepare_references(dataset, fasta, metadata, names, workdir):
    """ARIBA reference directory and KMA database for one catalogue.

    ``ariba prepareref`` drops a sequence that is declared as a gene when it
    finds no complete reading frame in it.  Such sequences are declared
    non-coding and the references are prepared again, so that every catalogue
    sequence is available to ARIBA.  Their names are returned and recorded.
    """
    directory = Path(workdir) / dataset
    prepared = directory / "ariba_prepareref"
    relaxed = directory / "declared_noncoding_for_ariba.txt"

    def build():
        shutil.rmtree(prepared, ignore_errors=True)
        code, _ = run(["ariba", "prepareref", "-f", str(fasta), "-m", str(metadata),
                       str(prepared)], directory / "ariba_prepareref.log")
        if code:
            raise RuntimeError(f"ariba prepareref failed for {dataset}; see the log")

    def kept():
        found = set()
        clusters = prepared / "02.cdhit.clusters.tsv"
        if clusters.is_file():
            for line in clusters.read_text().splitlines():
                found.update(line.split("\t")[1:])
        return found

    if not (prepared / "02.cdhit.clusters.tsv").is_file():
        build()
    dropped = sorted(set(names) - kept())
    if dropped:
        print(f"  ariba prepareref dropped {len(dropped)} sequence(s) declared as genes; "
              "declaring them non-coding", flush=True)
        relaxed.write_text("\n".join(dropped) + "\n")
        lines = []
        for line in Path(metadata).read_text().splitlines():
            fields = line.split("\t")
            if fields[0] in dropped:
                fields[1] = "0"
            lines.append("\t".join(fields))
        Path(metadata).write_text("\n".join(lines) + "\n")
        build()
        for record in (Path(workdir) / "raw" / dataset).glob("*/ariba_*.run.json"):
            record.unlink()
        if set(names) - kept():
            raise RuntimeError(f"ariba prepareref still drops sequences of {dataset}")
    noncoding = relaxed.read_text().split() if relaxed.is_file() else []
    database = directory / "kma_db"
    if not Path(str(database) + ".name").is_file():
        code, _ = run(["kma", "index", "-i", str(fasta), "-o", str(database)],
                      directory / "kma_index.log")
        if code:
            raise RuntimeError(f"kma index failed for {dataset}; see the log")
    return prepared, database, noncoding


def run_tools(dataset, sample, prepared, database, workdir, threads):
    """Run the three configurations for one sample; keep the raw outputs."""
    raw = Path(workdir) / "raw" / dataset / sample["sample"]
    raw.mkdir(parents=True, exist_ok=True)
    timings = {}
    for tool, configuration, extra in CONFIGURATIONS:
        tag = f"{tool.lower()}_{configuration}"
        record = raw / f"{tag}.run.json"
        if record.is_file() and json.loads(record.read_text()).get("exit_code") == 0:
            timings[tag] = json.loads(record.read_text())
            continue
        scratch = raw / f"{tag}.work"
        shutil.rmtree(scratch, ignore_errors=True)
        if tool == "ARIBA":
            command = ["ariba", "run", "--threads", str(threads)] + extra + [
                str(prepared), sample["r1"], sample["r2"], str(scratch)]
            code, seconds = run(command, raw / f"{tag}.log")
            for source, target in (("report.tsv", f"{tag}.report.tsv"),
                                   ("assembled_seqs.fa.gz", f"{tag}.assembled_seqs.fa.gz"),
                                   ("assembled_genes.fa.gz", f"{tag}.assembled_genes.fa.gz"),
                                   ("assemblies.fa.gz", f"{tag}.assemblies.fa.gz")):
                if (scratch / source).is_file():
                    shutil.copy2(scratch / source, raw / target)
        else:
            scratch.mkdir()
            command = ["kma", "-ipe", sample["r1"], sample["r2"], "-o",
                       str(scratch / "kma"), "-t_db", str(database), "-t", str(threads)]
            code, seconds = run(command, raw / f"{tag}.log")
            for suffix in ("res", "fsa"):
                if (scratch / f"kma.{suffix}").is_file():
                    shutil.copy2(scratch / f"kma.{suffix}", raw / f"{tag}.{suffix}")
        shutil.rmtree(scratch, ignore_errors=True)
        timings[tag] = {"exit_code": code, "seconds": seconds,
                        "command": [Path(part).name if part.startswith("/") else part
                                    for part in command]}
        record.write_text(json.dumps(timings[tag]))
        print(f"    {tag}: exit {code}, {seconds} s", flush=True)
    return raw, timings


# -------------------------------------------------------------------------- parsing

def ariba_results(raw, tag):
    """Per reference name: report fields and the sequences that ARIBA returned."""
    report = {}
    path = raw / f"{tag}.report.tsv"
    if path.is_file():
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name = row["ref_name"]
                entry = report.setdefault(name, {"rows": 0, "multi_nt_rows": 0, "first": row})
                entry["rows"] += 1
                if "," in (row.get("smtls_nts") or ""):
                    entry["multi_nt_rows"] += 1
    sequences = {}
    path = raw / f"{tag}.assembled_seqs.fa.gz"
    if path.is_file() and path.stat().st_size:
        for name, sequence in read_fasta(path):
            sequences.setdefault(name.split(".")[0], []).append((name, sequence))
    return report, sequences


def ariba_report_text(entry):
    row = entry["first"]
    flag = int(row["flag"])
    names = [name for bit, name in enumerate(ARIBA_FLAGS) if flag & (1 << bit)]
    return (
        f"flag {flag} ({', '.join(names) or 'none'}); "
        f"ref_base_assembled {row['ref_base_assembled']}/{row['ref_len']}; "
        f"pc_ident {row['pc_ident']}; ctg_cov {row['ctg_cov']}; reads {row['reads']}; "
        f"report rows {entry['rows']}, with more than one nucleotide in the reads "
        f"{entry['multi_nt_rows']}"
    )


def kma_results(raw, tag):
    """Per template: the .res fields and the consensus."""
    report = {}
    path = raw / f"{tag}.res"
    if path.is_file():
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                clean = {key.strip().lstrip("#"): (value or "").strip()
                         for key, value in row.items() if key}
                report[clean["Template"]] = clean
    sequences = {}
    path = raw / f"{tag}.fsa"
    if path.is_file() and path.stat().st_size:
        for name, sequence in read_fasta(path):
            sequences.setdefault(name, []).append((name, sequence))
    return report, sequences


def kma_report_text(row):
    keys = ("Template_Identity", "Template_Coverage", "Query_Identity",
            "Query_Coverage", "Depth", "Score")
    return "; ".join(f"{key} {row[key]}" for key in keys if key in row)


# -------------------------------------------------------------------------- scoring

def blast_pair(query, subject, work):
    """Longest HSP of the query against the subject."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    (work / "q.fasta").write_text(f">q\n{query}\n")
    (work / "s.fasta").write_text(f">s\n{subject}\n")
    completed = subprocess.run(
        ["blastn", "-query", str(work / "q.fasta"), "-subject", str(work / "s.fasta"),
         "-outfmt", "6 pident length mismatch gaps", "-max_hsps", "1"],
        capture_output=True, text=True, check=False,
    )
    best = None
    for line in completed.stdout.splitlines():
        pident, length, mismatch, gaps = line.split("\t")[:4]
        if best is None or int(length) > best[1]:
            best = (float(pident), int(length), int(mismatch), int(gaps))
    return best


def edit_distance(first, second):
    """Global Levenshtein distance (bit-vector algorithm of Myers and Hyyro)."""
    if not first or not second:
        return len(first) + len(second)
    length = len(first)
    mask, high = (1 << length) - 1, 1 << (length - 1)
    positions = {}
    for index, letter in enumerate(first):
        positions[letter] = positions.get(letter, 0) | (1 << index)
    plus, minus, score = mask, 0, length
    for letter in second:
        equal = positions.get(letter, 0)
        vertical = equal | minus
        horizontal = (((equal & plus) + plus) ^ plus) | equal
        plus_h = (minus | ~(horizontal | plus)) & mask
        minus_h = plus & horizontal
        if plus_h & high:
            score += 1
        elif minus_h & high:
            score -= 1
        plus_h = ((plus_h << 1) | 1) & mask
        minus_h = (minus_h << 1) & mask
        plus = (minus_h | ~(vertical | plus_h)) & mask
        minus = plus_h & vertical
    return score


def compare(sequence, truth_records, work):
    """Relation of one returned sequence to the truth copies of a locus.

    ``identical``: the same sequence, on either strand.
    ``uncalled bases, otherwise identical``: same length, differs only at N.
    ``incomplete, identical over its own length``: an exact part of the truth.
    ``truth contained, extra bases``: the whole truth plus bases at an end.
    ``identical over the overlap, different boundaries``: the aligned part is
    identical and at least one sequence continues beyond it.
    ``differs``: at least one aligned base or gap column differs.

    Identity, aligned length, mismatches and gap columns describe the longest
    local alignment.  ``edit_distance`` is the global Levenshtein distance
    between the whole returned sequence and the whole truth copy named in the
    row, so it also counts bases that are missing or extra.
    """
    called = sequence.upper().replace("-", "")
    result = _compare(called, truth_records, work)
    truth = dict(truth_records)[result["truth_name"]]
    result["edit_distance"] = min(edit_distance(called, truth),
                                  edit_distance(reverse_complement(called), truth))
    return result


def recovered(sequences, truth_records):
    """How many distinct truth sequences are returned exactly, as ``n of m``."""
    distinct = {truth for _, truth in truth_records}
    if not distinct:
        return ""
    called = {sequence.upper().replace("-", "") for sequence in sequences}
    hits = [truth for truth in distinct
            if truth in called or reverse_complement(truth) in called]
    return f"{len(hits)} of {len(distinct)}"


def _compare(called, truth_records, work):
    exact, aligned = None, None
    for name, truth in truth_records:
        for oriented in (called, reverse_complement(called)):
            if oriented == truth:
                relation = RELATIONS[0]
            elif len(oriented) == len(truth) and "N" in oriented and all(
                a == b or a == "N" for a, b in zip(oriented, truth)
            ):
                relation = RELATIONS[1]
            elif oriented and oriented in truth:
                relation = RELATIONS[2]
            elif truth in oriented:
                relation = RELATIONS[3]
            else:
                continue
            candidate = (RELATIONS.index(relation), name, truth, relation)
            if exact is None or candidate[0] < exact[0]:
                exact = candidate
        hsp = blast_pair(called, truth, work)
        if hsp and (aligned is None or hsp[0] * hsp[1] > aligned[0][0] * aligned[0][1]):
            aligned = (hsp, name, truth)
    if exact is not None:
        _, name, truth, relation = exact
        return {"truth_name": name, "truth_bp": len(truth), "comparison": relation,
                "identity_pct": "100.00", "aligned_bp": min(len(called), len(truth)),
                "mismatches": 0, "gap_columns": 0,
                "length_difference_bp": len(called) - len(truth)}
    if aligned is None:
        name, truth = truth_records[0]
        return {"truth_name": name, "truth_bp": len(truth), "comparison": RELATIONS[6],
                "length_difference_bp": len(called) - len(truth)}
    (pident, length, mismatch, gaps), name, truth = aligned
    relation = RELATIONS[4] if mismatch == 0 and gaps == 0 else RELATIONS[5]
    return {"truth_name": name, "truth_bp": len(truth), "comparison": relation,
            "identity_pct": f"{pident:.2f}", "aligned_bp": length,
            "mismatches": mismatch, "gap_columns": gaps,
            "length_difference_bp": len(called) - len(truth)}


def score_sample(dataset, sample, catalogue_names, raw, timings, work):
    """One row per locus and configuration, and the returned sequences."""
    rows, returned = [], []
    for tool, configuration, _ in CONFIGURATIONS:
        tag = f"{tool.lower()}_{configuration}"
        timing = timings.get(tag, {})
        if tool == "ARIBA":
            report, sequences = ariba_results(raw, tag)
        else:
            report, sequences = kma_results(raw, tag)
        for locus, truth in sample["truth"].items():
            if dataset == "hpylori":
                references = [locus]
            else:
                references = [name for name in catalogue_names
                              if name in sequences or name in report]
            if tool == "KMA" and len(references) > 1:
                references.sort(
                    key=lambda name: -float(report.get(name, {}).get("Score") or 0))
            base = {
                "dataset": dataset, "sample": sample["sample"], "locus": locus,
                "tool": tool, "configuration": configuration,
                "exit_code": timing.get("exit_code", ""), "seconds": timing.get("seconds", ""),
                "truth_copies": len(truth), "reference_kind": REFERENCE_KIND[dataset],
            }
            found = [(name, seq) for reference in references
                     for name, seq in sequences.get(reference, [])]
            reported = [reference for reference in references if reference in report]
            if tool == "ARIBA":
                text = " | ".join(f"{name}: {ariba_report_text(report[name])}"
                                  for name in reported)
            else:
                text = " | ".join(f"{name}: {kma_report_text(report[name])}"
                                  for name in reported)
            base["tool_report"] = text or "locus absent from the report"
            base["tool_reference"] = ";".join(reported)
            if truth:
                base.update(truth_name=truth[0][0], truth_bp=len(truth[0][1]))
            base["truth_sequences_identical"] = recovered([seq for _, seq in found], truth)
            if not found:
                rows.append(dict(base, sequences_returned=0, returned_bp=0,
                                 uncalled_bases=0, lowercase_bases=0,
                                 comparison="no sequence returned"))
                continue
            scored = []
            for index, (name, sequence) in enumerate(found):
                plain = sequence.replace("-", "")
                identifier = (f"{dataset}|{sample['sample']}|{locus}|{tool}|"
                              f"{configuration}|{index + 1}")
                returned.append((identifier + " " + name, plain))
                entry = {"returned_bp": len(plain),
                         "uncalled_bases": plain.upper().count("N"),
                         "lowercase_bases": sum(1 for letter in plain if letter.islower())}
                if truth:
                    entry.update(compare(plain, truth, work))
                    also = [other for other, other_seq in sample["other"].items()
                            if plain.upper() in (other_seq, reverse_complement(other_seq))]
                    entry["also_identical_to"] = ";".join(also)
                else:
                    entry["comparison"] = RELATIONS[7]
                entry["label"] = f"{name.split('.')[0]}, {len(plain)} bp: {entry['comparison']}"
                if entry.get("truth_name"):
                    entry["label"] += (
                        f" [{entry['truth_name']}; {entry.get('mismatches', '')} mismatches, "
                        f"{entry.get('gap_columns', '')} gap columns]")
                if entry.get("also_identical_to"):
                    entry["label"] += f" (identical to {entry['also_identical_to']})"
                scored.append(entry)
            detail = "; ".join(entry.pop("label") for entry in scored)
            scored.sort(key=lambda entry: (RELATIONS.index(entry["comparison"]),
                                           -int(entry.get("aligned_bp") or 0)))
            rows.append(dict(base, sequences_returned=len(found),
                             all_returned_sequences=detail, **scored[0]))
    return rows, returned


# --------------------------------------------------------- deposited Locus-Recon rows

LOCUS_RECON = ("Locus-Recon", "1.0.0, deposited result")


def _evidence_report(record, keys):
    return "; ".join(f"{label} {record.get(key) or 'none'}" for label, key in keys)


def locus_recon_rows(dataset, samples, work):
    """Results of Locus-Recon for the same cases, read from the deposited tables.

    They are not recomputed here.  For *H. pylori* and for the simulated set
    the deposit holds the comparison with the truth and not the candidate
    sequences, so the relation is derived from the deposited fields.
    """
    rows = []
    wanted = {sample["sample"] for sample in samples}
    if dataset == "hpylori":
        for table in ("tier_calibration_closed_genomes.tsv",
                      "tier_calibration_divergent_loci.tsv"):
            for record in calibration.read_tsv(CALIBRATION_DIR / table):
                strain, locus = record["case"].split()
                if strain not in wanted:
                    continue
                found = re.search(r"truth (\d+) bp x (\d+) copy/ies; candidate (\d+) bp at "
                                  r"([\d.]+)% over (\d+) bp", record["evidence"])
                truth_bp, copies, returned, aligned = (int(found.group(i)) for i in (1, 2, 3, 5))
                identity = float(found.group(4))
                if record["matches_truth"] == "yes":
                    relation = RELATIONS[0]
                elif identity == 100 and aligned == returned < truth_bp:
                    relation = RELATIONS[2]
                elif identity == 100 and aligned == truth_bp < returned:
                    relation = RELATIONS[3]
                else:
                    relation = RELATIONS[4] if identity == 100 else RELATIONS[5]
                contained = relation in RELATIONS[:4]
                truth = next(sample["truth"][locus] for sample in samples
                             if sample["sample"] == strain)
                distinct = len({sequence for _, sequence in truth})
                rows.append({
                    "truth_name": truth[0][0],
                    "truth_sequences_identical":
                        f"{int(record['matches_truth'] == 'yes')} of {distinct}",
                    "dataset": dataset, "sample": strain, "locus": locus,
                    "sequences_returned": 1, "returned_bp": returned, "truth_bp": truth_bp,
                    "truth_copies": copies, "comparison": relation,
                    "identity_pct": f"{identity:.2f}", "aligned_bp": aligned,
                    "mismatches": 0 if identity == 100 else "",
                    "gap_columns": 0 if identity == 100 else "",
                    "length_difference_bp": returned - truth_bp,
                    "edit_distance": abs(returned - truth_bp) if contained else "",
                    "tool_reference": f"tier-calibration/{table}",
                    "tool_report": _evidence_report(record, (
                        ("tier", "sequence_confidence"), ("disposition", "disposition"),
                        ("flags", "qc_flags"))),
                })
    elif dataset == "mock":
        directory = VALIDATION / "reference_results"
        checked = {row["sample_id"]: row
                   for row in calibration.read_tsv(directory / "validation_results.tsv")}
        for record in calibration.read_tsv(directory / "locus_recon_report_mockLocus.tsv"):
            name = record["sample_id"]
            if name not in wanted:
                continue
            check = checked[name]
            returned = int(record["allele_length"] or 0)
            identity = float(check["truth_identity_pct"] or 0)
            row = {"dataset": dataset, "sample": name, "locus": "mockLocus",
                   "sequences_returned": 1 if returned else 0, "returned_bp": returned,
                   "tool_reference": "reference_results/locus_recon_report_mockLocus.tsv"}
            truth = next(sample["truth"]["mockLocus"] for sample in samples
                         if sample["sample"] == name)
            row["truth_copies"] = len(truth)
            if truth:
                row["truth_name"], row["truth_bp"] = truth[0][0], len(truth[0][1])
                row["truth_sequences_identical"] = (
                    f"{int(check['exact_truth_match'] == 'True')} of {len(truth)}")
            if not returned:
                row["comparison"] = "no sequence returned"
            elif check["exact_truth_match"] == "True":
                row.update(comparison=RELATIONS[0], identity_pct="100.00", aligned_bp=returned,
                           mismatches=0, gap_columns=0, length_difference_bp=0,
                           edit_distance=0)
            else:
                row.update(comparison=RELATIONS[5], identity_pct=f"{identity:.2f}",
                           length_difference_bp=returned - len(truth[0][1]))
            row["tool_report"] = _evidence_report(record, (
                ("status", "workflow_status"), ("tier", "sequence_confidence"),
                ("disposition", "result_disposition"))) + _evidence_report(check, (
                    ("; mixture detected", "mixture_detected"),
                    ("ambiguity flagged", "ambiguity_flagged")))
            rows.append(row)
    else:
        directory = VALIDATION / "liba6656-gfa-rerun"
        record = calibration.read_tsv(directory / "standard_workflow_report_tcdB.tsv")[0]
        sequence = read_fasta(directory / "standard_workflow_reconstructed_tcdB.fasta")[0][1]
        truth = samples[0]["truth"]["tcdB"]
        row = {"dataset": dataset, "sample": samples[0]["sample"], "locus": "tcdB",
               "sequences_returned": 1, "returned_bp": len(sequence),
               "truth_copies": len(truth),
               "tool_reference": "liba6656-gfa-rerun/standard_workflow_report_tcdB.tsv",
               "tool_report": _evidence_report(record, (
                   ("tier", "sequence_confidence"), ("disposition", "result_disposition"),
                   ("flags", "qc_flags")))}
        row.update(compare(sequence, truth, work))
        row["truth_sequences_identical"] = recovered([sequence], truth)
        rows.append(row)
    for row in rows:
        row["tool"], row["configuration"] = LOCUS_RECON
        row["reference_kind"] = REFERENCE_KIND[dataset]
        row.setdefault("uncalled_bases", 0)
        row.setdefault("lowercase_bases", 0)
    return rows


def summarise(rows, groups):
    """Counts of each relation per dataset, group of cases and configuration."""
    relations = list(RELATIONS) + ["no sequence returned"]
    table = {}
    for row in rows:
        key = (row["dataset"], groups[(row["dataset"], row["sample"], row["locus"])],
               row["tool"], row["configuration"])
        entry = table.setdefault(key, {"cases": 0, "all": 0, "seconds": [],
                                       **{relation: 0 for relation in relations}})
        entry["cases"] += 1
        entry[row["comparison"]] += 1
        found, _, present = str(row.get("truth_sequences_identical") or "").partition(" of ")
        entry["all"] += int(bool(present) and found == present)
        if row.get("seconds") not in ("", None):
            entry["seconds"].append(float(row["seconds"]))
    summary = []
    for (dataset, group, tool, configuration), entry in table.items():
        seconds = sorted(entry["seconds"])
        summary.append({
            "dataset": dataset, "group": group, "tool": tool, "configuration": configuration,
            "cases": entry["cases"], **{relation: entry[relation] for relation in relations},
            "every_truth_sequence_identical": entry["all"],
            "seconds_min": f"{seconds[0]:.1f}" if seconds else "",
            "seconds_max": f"{seconds[-1]:.1f}" if seconds else "",
        })
    return summary


# ----------------------------------------------------------------------------- main

def tool_versions():
    versions = {}
    for name, command in (("ariba", ["ariba", "version"]), ("kma", ["kma", "-v"]),
                          ("blastn", ["blastn", "-version"])):
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        versions[name] = (completed.stdout + completed.stderr).strip()
    return versions


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--workdir", required=True,
                        help="directory for downloads and raw tool outputs")
    parser.add_argument("--outdir", default=str(VALIDATION / "comparator-benchmark"))
    parser.add_argument("--datasets", default="mock,hpylori,liba6656")
    parser.add_argument("--samples", default="",
                        help="comma-separated sample names; default: all")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--score-only", action="store_true",
                        help="use the raw outputs already under <workdir>/raw")
    args = parser.parse_args()

    workdir, outdir = Path(args.workdir).resolve(), Path(args.outdir).resolve()
    cache = workdir / "downloads"
    workdir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    rows, returned, truth_records, catalogues, removed, groups = [], [], [], [], {}, {}
    for dataset in args.datasets.split(","):
        print(f"[{dataset}] catalogue and samples", flush=True)
        fasta, metadata, names = build_catalogue(dataset, workdir)
        catalogues += [fasta, metadata]
        samples = samples_of(dataset, cache, set(filter(None, args.samples.split(","))),
                             reads=not args.score_only)
        rows += locus_recon_rows(dataset, samples, workdir / "blast")
        if not args.score_only:
            prepared, database, removed[dataset] = prepare_references(
                dataset, fasta, metadata, names, workdir)
        for sample in samples:
            for locus in sample["truth"]:
                groups[(dataset, sample["sample"], locus)] = sample.get("group") or locus
            print(f"  [{dataset}] {sample['sample']}", flush=True)
            raw = workdir / "raw" / dataset / sample["sample"]
            timings = {}
            if args.score_only:
                for tool, configuration, _ in CONFIGURATIONS:
                    record = raw / f"{tool.lower()}_{configuration}.run.json"
                    if record.is_file():
                        timings[f"{tool.lower()}_{configuration}"] = json.loads(
                            record.read_text())
            else:
                raw, timings = run_tools(dataset, sample, prepared, database, workdir,
                                         args.threads)
            sample_rows, sample_sequences = score_sample(
                dataset, sample, names, raw, timings, workdir / "blast")
            rows += sample_rows
            returned += sample_sequences
            for records in sample["truth"].values():
                truth_records += [(f"{dataset}|{name}", seq) for name, seq in records]

    with open(outdir / "comparator_benchmark_results.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_HEADER, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        tools = [LOCUS_RECON] + [(tool, name) for tool, name, _ in CONFIGURATIONS]
        position = {}
        for row in rows:
            position.setdefault((row["dataset"], row["sample"], row["locus"]), len(position))
        rows.sort(key=lambda row: (position[(row["dataset"], row["sample"], row["locus"])],
                                   tools.index((row["tool"], row["configuration"]))))
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in RESULT_HEADER})
    summary = summarise(rows, groups)
    with open(outdir / "comparator_benchmark_summary.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    write_fasta(outdir / "comparator_returned_sequences.fasta", returned)
    seen, unique_truth = set(), []
    for name, sequence in truth_records:
        if name not in seen:
            seen.add(name)
            unique_truth.append((name, sequence))
    write_fasta(outdir / "comparator_truth_sequences.fasta", unique_truth)
    for path in catalogues:
        shutil.copy2(path, outdir / f"{path.parent.name}_{path.name}")

    if not args.score_only:
        versions = tool_versions()
        (outdir / "TOOL_VERSIONS.txt").write_text(
            "\n\n".join(f"== {name}\n{text}" for name, text in versions.items()) + "\n")
        provenance = build_provenance(
            workflow="comparator_benchmark",
            repo_root=REPO_ROOT,
            inputs=[CALIBRATION_DIR / f"bait_{locus}_26695.fasta" for locus in HPYLORI_LOCI]
            + [CALIBRATION_DIR / "closed_genome_manifest.tsv",
               MOCK_DIR / "mock_locus_alleles.fasta", MOCK_DIR / "dataset_manifest.json",
               LIBA_BAITS, LIBA_CANDIDATES],
            parameters={
                "datasets": args.datasets.split(","), "threads": args.threads,
                "configurations": [
                    {"tool": tool, "configuration": name, "extra_arguments": extra}
                    for tool, name, extra in CONFIGURATIONS],
                "catalogue_sequences_declared_noncoding_for_ariba": removed,
            },
            tools={"python": [sys.executable, "--version"],
                   "kma": ["kma", "-v"], "blastn": ["blastn", "-version"]},
        )
        provenance["tool_versions_text"] = versions
        write_json(outdir / "PROVENANCE.json", provenance)
    print(f"{len(rows)} rows written to {outdir}", flush=True)


if __name__ == "__main__":
    main()
