#!/usr/bin/env python3
"""Calibration of the confidence tiers against known truth.

Every other validation stage asks whether a reconstruction is correct. This one
asks a different question: when the tool reports a tier, what does that tier
buy the reader?  Two quantities answer it.

- Precision of the top tier: of the results reported as deposit-ready, how many
  match truth exactly.  A single incorrect result in that tier would make the
  tier worthless, so the count of false accepts is the primary number.
- Acceptance rate: what fraction of results reach the top tier at all.  This is
  a property of the case mix, not of the tool: a benchmark built mostly from
  degraded inputs is expected to accept few of them.  The rate is therefore
  reported per case class, and the class definition is prespecified below
  rather than read off the results.

Case classes:
  supported       the locus is single copy and intact in the draft, depth is
                  adequate, and the sample is a pure culture.  A correct tool
                  accepts these.
  divergent       a single-copy locus whose allele sits further than
                  DIVERGENT_IDENTITY_PCT from the one-allele bait.  Selection
                  ranks candidates by measured distance, so a locus that ranks
                  high but still lands within a few points of the housekeeping
                  locus is not treated as a divergence stress: the class
                  follows the measured identity, not which arm ran it.
  low-depth       the same single-copy locus reconstructed from reads
                  subsampled to a fraction of their original depth.  The tier
                  is expected to fall as depth falls, and a degraded
                  reconstruction must be withheld rather than accepted.
  multi-copy      the locus is annotated in more than one copy in the closed
                  genome, so near-identical copies collapse in a short-read
                  draft.  A correct tool withholds the top tier even when the
                  consensus it reports is exact, because the reads cannot
                  establish which copy was reconstructed.
  degraded        depth below the stated floor, a controlled allele mixture, a
                  paralogue, or a target-negative genome.  A correct tool
                  withholds these, and withholding is not an error.
  truncated       the locus runs off a contig end by construction.  The tier is
                  expected to fall with the size of the loss.

Stage A uses the deposited benchmark results: the 13-sample mock truth set
(``reference_results/validation_results.tsv``) and the 14-case completeness
benchmark (``completeness-benchmark/completeness_benchmark_results.tsv``).
Both carry a per-case truth verdict alongside the tier, so the calibration is
recomputed from them rather than asserted.

Stage B repeats the closed-genome audit end to end on public data, which is
where the top tier has to be earned on real reads: five closed *H. pylori*
chromosomes, their Illumina runs, and a fixed external bait from strain 26695.
Sample-specific truth is never used as a bait.  For each genome and locus the
draft is assembled, the locus reconstructed, and the candidate compared to the
annotated gene sequence of that same closed genome.

Usage:
    # Stage A only, deterministic, no network or external tools
    python validation/run_tier_calibration.py --outdir validation/tier-calibration

    # Stage A and B
    python validation/run_tier_calibration.py \
        --outdir validation/tier-calibration \
        --manifest validation/tier-calibration/closed_genome_manifest.tsv \
        --workdir /scratch/tier-calibration --threads 24
"""

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from locus_recon.io import read_fasta_sequences, reverse_complement  # noqa: E402
from validation.common import (  # noqa: E402
    build_provenance,
    run_checked,
    sha256_file,
    write_json,
)

NCBI_GENOMES = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
ENA_FILEREPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run}"
    "&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
)
BAIT_ASSEMBLY = "GCF_000008525.1"          # H. pylori 26695, the fixed bait source

# Prespecified shortlist for the divergence arm: single-copy H. pylori genes
# with no known paralogue family, spanning housekeeping conservation to the
# mosaic virulence genes.  Which of them are used is decided by measurement in
# select_divergent_loci(), not here.
CANDIDATE_LOCI = {
    "glmM": "phosphoglucosamine mutase, housekeeping",
    "recA": "recombinational repair, housekeeping",
    "ftsZ": "cell division protein, housekeeping",
    "ureB": "urease beta subunit",
    "katA": "catalase",
    "flaA": "flagellin A",
    "rocF": "arginase",
    "vacA": "vacuolating cytotoxin, mosaic and highly divergent",
    "cagA": "cytotoxin-associated antigen, divergent and variably present",
}
TOP_TIER = "HIGH"
# A single-copy locus within this distance of the bait is as close to it as the
# housekeeping loci are, so it belongs in `supported` rather than `divergent`.
DIVERGENT_IDENTITY_PCT = 95.0
TOP_DISPOSITION = "PASS"

# Prespecified class of every mock sample.  "supported" means a correct tool is
# expected to accept; "degraded" means it is expected to withhold.
MOCK_CLASSES = {
    "baseline": "supported",
    "fragmented": "supported",
    "depth_series": "degraded",
    "mixture": "degraded",
    "paralog": "degraded",
    "negative": "degraded",
}


def write_tsv(path, header, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _truthy(value):
    return str(value).strip().lower() in {"true", "yes", "1"}


def stage_a(repo_root, outdir):
    """Recompute the tier calibration from the deposited benchmark results."""
    validation = Path(repo_root) / "validation"
    rows = []

    dispositions = {
        record["sample_id"]: record.get("result_disposition", "")
        for record in read_tsv(
            validation / "reference_results" / "locus_recon_report_mockLocus.tsv"
        )
    }
    for record in read_tsv(validation / "reference_results" / "validation_results.tsv"):
        if record["category"] not in MOCK_CLASSES:
            raise RuntimeError(
                f"mock category {record['category']!r} has no prespecified class"
            )
        rows.append([
            "mock truth set",
            record["sample_id"],
            MOCK_CLASSES[record["category"]],
            record["status"],
            dispositions.get(record["sample_id"], ""),
            record["qc_confidence"],
            "yes" if _truthy(record["exact_truth_match"]) else "no",
            record["expected_outcome"],
        ])

    for record in read_tsv(
        validation / "completeness-benchmark" / "completeness_benchmark_results.tsv"
    ):
        withheld = int(record["truth_clipped_bp"])
        rows.append([
            "completeness benchmark",
            record["case"],
            "truncated" if withheld else "supported",
            "SUCCESS",
            record["disposition_with_span"],
            record["sequence_confidence_with_span"],
            "yes" if _truthy(record["exact"]) else "no",
            f"{withheld} bp withheld by construction",
        ])

    write_tsv(
        Path(outdir) / "tier_calibration_benchmarks.tsv",
        ["stage", "case", "case_class", "workflow_status", "disposition",
         "sequence_confidence", "matches_truth", "case_design"],
        rows,
    )
    return summarise(rows, stage="benchmarks")


def summarise(rows, stage):
    """Acceptance rate and top-tier precision, overall and per case class."""
    def bucket(subset):
        top = [row for row in subset if row[5] == TOP_TIER]
        exact_in_top = [row for row in top if row[6] == "yes"]
        withheld_exact = [
            row for row in subset if row[5] != TOP_TIER and row[6] == "yes"
        ]
        return {
            "cases": len(subset),
            "reached_top_tier": len(top),
            "acceptance_rate": round(len(top) / len(subset), 3) if subset else None,
            "top_tier_exact": len(exact_in_top),
            "top_tier_precision": (
                round(len(exact_in_top) / len(top), 3) if top else None
            ),
            "false_accepts": len(top) - len(exact_in_top),
            "withheld_but_exact": len(withheld_exact),
        }

    classes = sorted({row[2] for row in rows})
    return {
        "stage": stage,
        "top_tier": TOP_TIER,
        "overall": bucket(rows),
        "by_case_class": {
            name: bucket([row for row in rows if row[2] == name]) for name in classes
        },
    }


def _fetch(url, destination):
    destination = Path(destination)
    if destination.is_file() and destination.stat().st_size:
        return str(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=600) as response, open(
        destination, "wb"
    ) as handle:
        shutil.copyfileobj(response, handle)
    return str(destination)


def _refseq_directory(accession):
    """Resolve the RefSeq download directory for one assembly accession."""
    prefix, digits = accession.split("_")
    body = digits.split(".")[0]
    parent = f"{NCBI_GENOMES}/{prefix}/{body[0:3]}/{body[3:6]}/{body[6:9]}/"
    with urllib.request.urlopen(parent, timeout=300) as response:
        listing = response.read().decode("utf-8", "replace")
    matches = sorted(set(re.findall(rf'href="({re.escape(accession)}_[^"/]+)/"', listing)))
    if not matches:
        raise RuntimeError(f"no RefSeq directory found for {accession}")
    return parent + matches[0] + "/", matches[0]


def _annotated_records(accession, cache, kind):
    """Download and read one RefSeq annotation FASTA (``cds`` or ``rna``)."""
    directory, name = _refseq_directory(accession)
    filename = f"{name}_{kind}_from_genomic.fna.gz"
    archive = _fetch(directory + filename, Path(cache) / filename)
    records, header, sequence = [], None, []
    with gzip.open(archive, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(sequence)))
                header, sequence = line[1:].strip(), []
            else:
                sequence.append(line.strip())
    if header:
        records.append((header, "".join(sequence)))
    return records


def _select(records, gene=None, product=None):
    chosen = []
    for header, sequence in records:
        if gene and f"[gene={gene}]" not in header:
            continue
        if product and f"[product={product}]" not in header:
            continue
        chosen.append((header.split()[0], sequence))
    return chosen


def _reads(run, cache):
    with urllib.request.urlopen(ENA_FILEREPORT.format(run=run), timeout=300) as response:
        report = list(csv.DictReader(
            response.read().decode().splitlines(), delimiter="\t"
        ))[0]
    paths = []
    for url in report["fastq_ftp"].split(";"):
        paths.append(_fetch("https://" + url, Path(cache) / Path(url).name))
    if len(paths) != 2:
        raise RuntimeError(f"expected a read pair for {run}, got {len(paths)} file(s)")
    return paths


def _blast_identity(candidate, truth, work):
    """Best identity and aligned length of a candidate against one truth copy."""
    work = Path(work)
    query, subject = work / "candidate.fasta", work / "truth.fasta"
    query.write_text(f">candidate\n{candidate}\n")
    subject.write_text(f">truth\n{truth}\n")
    completed = subprocess.run(
        ["blastn", "-query", str(query), "-subject", str(subject), "-outfmt",
         "6 pident length", "-max_hsps", "1"],
        capture_output=True, text=True, check=False,
    )
    best = (0.0, 0)
    for line in completed.stdout.splitlines():
        pident, length = line.split("\t")[:2]
        if int(length) > best[1]:
            best = (float(pident), int(length))
    return best


def _write_bait(records, path, locus, source):
    """Write the first annotated record of a locus as a one-allele bait file."""
    if not records:
        raise RuntimeError(f"no {locus} record in the annotation of {source}")
    name, sequence = records[0]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 70):
            handle.write(sequence[start : start + 70] + "\n")
    return str(path)


def _reconstruct(repo_root, sample_dir, strain, locus, bait, draft, r1, r2, threads,
                 tag=""):
    """Run one reconstruction; return (report record or None, candidate, exit code)."""
    name = f"{locus}{tag}"
    out = Path(sample_dir) / name
    samplesheet = Path(sample_dir) / f"samples_{name}.tsv"
    samplesheet.write_text(
        "sample_id\tassembly_path\tr1_path\tr2_path\n"
        f"{strain}\t{draft}\t{r1}\t{r2}\n"
    )
    completed = subprocess.run(
        [sys.executable, "-m", "locus_recon.cli",
         "--samplesheet", str(samplesheet), "--main-output-dir", str(out),
         "--bait", bait, "--locus", locus,
         "--threads", str(threads), "--memory-per-sample", "32", "--no-progress"],
        cwd=repo_root, capture_output=True, text=True, check=False,
    )
    reports = list(out.glob(f"locus_recon_report_{locus}.tsv"))
    if not reports:
        return None, "", completed.returncode
    candidates = list(out.glob(f"{strain}/*_reconstructed.fasta"))
    sequence = read_fasta_sequences(str(candidates[0]))[0][1] if candidates else ""
    return read_tsv(reports[0])[0], sequence, completed.returncode


def _score(sequence, truth_records, work):
    """Compare a candidate to the annotated truth copies of the same genome.

    Returns the strict verdict used for calibration, the best identity and
    aligned length, and a relation string that separates the two ways a
    candidate can fail to be identical: disagreeing bases, or agreeing bases
    over a span that begins or ends somewhere else than the annotation does.
    """
    exact, best_identity, best_length, best_delta = "no", 0.0, 0, 0
    for _, truth_sequence in truth_records:
        if sequence and sequence in (
            truth_sequence, reverse_complement(truth_sequence)
        ):
            exact = "yes"
        identity, length = _blast_identity(sequence, truth_sequence, work)
        if identity * length > best_identity * best_length:
            best_identity, best_length = identity, length
            best_delta = len(sequence) - len(truth_sequence)
    if exact == "yes":
        relation = "identical to the annotated allele"
    elif best_identity >= 99.995 and best_length:
        relation = (
            f"identical over the {best_length} bp overlap; candidate "
            f"{best_delta:+d} bp against the annotated boundary"
        )
    elif best_length:
        relation = f"{best_identity:.2f}% identity over {best_length} bp"
    else:
        relation = "no alignment to the annotated allele"
    return exact, best_identity, best_length, relation


def _evidence(sequence, truth_records, identity, length, record):
    return (
        f"truth {len(truth_records[0][1])} bp x {len(truth_records)} copy/ies; "
        f"candidate {len(sequence)} bp at {identity:.2f}% over {length} bp; "
        f"catalogue identity {record.get('best_identity_pct', '')}%"
    )


def _draft(sample_dir, r1, r2, threads, name="spades", isolate=True):
    """Assemble one draft with SPAdes, reusing an existing one.

    ``--isolate`` is documented for high-coverage isolate data, so the
    subsampled arms are assembled in SPAdes' default mode instead.  That is a
    difference between the arms and it is stated in the deposit: each arm uses
    the settings one would actually use at that depth.
    """
    draft = Path(sample_dir) / name / "contigs.fasta"
    if not draft.is_file():
        run_checked(
            ["spades.py"] + (["--isolate"] if isolate else []) + [
                "-1", str(r1), "-2", str(r2),
                "-o", str(Path(sample_dir) / name),
                "-t", str(threads), "-m", "32",
            ], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return str(draft)


def stage_b(repo_root, outdir, manifest_path, workdir, threads):
    """Repeat the closed-genome audit on public data and tabulate the tiers."""
    manifest = read_tsv(manifest_path)
    workdir = Path(workdir)
    cache = workdir / "downloads"
    workdir.mkdir(parents=True, exist_ok=True)

    baits = {
        "gyrB": _write_bait(
            _select(_annotated_records(BAIT_ASSEMBLY, cache, "cds"), gene="gyrB"),
            Path(outdir) / "bait_gyrB_26695.fasta", "gyrB", BAIT_ASSEMBLY),
        "23S": _write_bait(
            _select(_annotated_records(BAIT_ASSEMBLY, cache, "rna"),
                    product="23S ribosomal RNA"),
            Path(outdir) / "bait_23S_26695.fasta", "23S", BAIT_ASSEMBLY),
    }

    rows = []
    for entry in manifest:
        strain, run = entry["strain"], entry["illumina_run"]
        sample_dir = workdir / strain
        sample_dir.mkdir(exist_ok=True)
        print(f"  [{strain}] reads {run}", flush=True)
        r1, r2 = _reads(run, cache)
        print(f"  [{strain}] draft", flush=True)
        draft = _draft(sample_dir, r1, r2, threads)

        truth = {
            "gyrB": _select(
                _annotated_records(entry["assembly_accession"], cache, "cds"),
                gene="gyrB"),
            "23S": _select(
                _annotated_records(entry["assembly_accession"], cache, "rna"),
                product="23S ribosomal RNA"),
        }

        for locus in ("gyrB", "23S"):
            # The class follows the annotation, not the result: a locus the
            # closed genome carries twice cannot be attributed to a single copy
            # from short reads, so withholding is correct behaviour.
            case_class = "supported" if len(truth[locus]) == 1 else "multi-copy"
            print(f"  [{strain}] reconstructing {locus}", flush=True)
            record, sequence, exit_code = _reconstruct(
                repo_root, sample_dir, strain, locus, baits[locus], draft,
                r1, r2, threads,
            )
            if record is None:
                rows.append([
                    "closed-genome audit", f"{strain} {locus}", case_class,
                    f"NO_REPORT(exit={exit_code})", "", "", "no",
                    "annotated closed-genome truth", "not assessed", "",
                ])
                continue
            exact, identity, length, relation = _score(
                sequence, truth[locus], sample_dir)
            rows.append([
                "closed-genome audit", f"{strain} {locus}", case_class,
                record.get("workflow_status", ""),
                record.get("result_disposition", ""),
                record.get("sequence_confidence") or record.get("qc_confidence", ""),
                exact,
                _evidence(sequence, truth[locus], identity, length, record),
                relation, record.get("qc_flags", ""),
            ])

    write_tsv(
        Path(outdir) / "tier_calibration_closed_genomes.tsv",
        ["stage", "case", "case_class", "workflow_status", "disposition",
         "sequence_confidence", "matches_truth", "evidence",
         "relation_to_annotation", "qc_flags"],
        rows,
    )
    return summarise(rows, stage="closed genomes"), rows


def select_divergent_loci(manifest, cache, outdir, wanted, min_bp=1200):
    """Choose single-copy loci that stress the tier through catalogue distance.

    The shortlist is prespecified: single-copy *H. pylori* genes with no known
    paralogue family, spanning housekeeping conservation to the mosaic
    virulence genes.  Selection among them is measured, not assumed: a locus is
    eligible only if the bait assembly and all five closed genomes annotate it
    exactly once, and the chosen loci are the eligible ones furthest from the
    bait allele.  Identity to the bait is what the tier has to survive, so it
    is the ranking key and it is deposited for every candidate.
    """
    bait_records = _annotated_records(BAIT_ASSEMBLY, cache, "cds")
    genomes = {
        entry["strain"]: _annotated_records(entry["assembly_accession"], cache, "cds")
        for entry in manifest
    }
    rows, eligible = [], []
    for gene, description in CANDIDATE_LOCI.items():
        bait = _select(bait_records, gene=gene)
        copies, identities = {}, []
        for strain, records in genomes.items():
            found = _select(records, gene=gene)
            copies[strain] = len(found)
            if len(bait) == 1 and len(found) == 1:
                identity, length = _blast_identity(
                    bait[0][1], found[0][1], Path(cache)
                )
                identities.append(identity)
        single_everywhere = (
            len(bait) == 1
            and all(count == 1 for count in copies.values())
            and len(bait[0][1]) >= min_bp
        )
        mean_identity = (
            round(sum(identities) / len(identities), 2) if identities else None
        )
        rows.append([
            gene, description, len(bait[0][1]) if bait else 0,
            len(bait), ";".join(f"{s}={c}" for s, c in sorted(copies.items())),
            "yes" if single_everywhere else "no",
            mean_identity if mean_identity is not None else "",
            round(min(identities), 2) if identities else "",
            round(max(identities), 2) if identities else "",
        ])
        if single_everywhere and mean_identity is not None:
            eligible.append((mean_identity, gene))
    write_tsv(
        Path(outdir) / "locus_selection.tsv",
        ["gene", "rationale", "bait_length_bp", "bait_copies", "copies_per_genome",
         "eligible", "mean_identity_to_bait_pct", "min_identity_pct",
         "max_identity_pct"],
        rows,
    )
    eligible.sort()
    return [gene for _, gene in eligible[:wanted]]


def stage_c_divergence(repo_root, outdir, manifest_path, workdir, threads, wanted):
    """Single-copy loci at greater catalogue distance, on the same drafts."""
    manifest = read_tsv(manifest_path)
    workdir = Path(workdir)
    cache = workdir / "downloads"
    loci = select_divergent_loci(manifest, cache, outdir, wanted)
    selected = {
        row["gene"]: float(row["mean_identity_to_bait_pct"])
        for row in read_tsv(Path(outdir) / "locus_selection.tsv")
        if row["mean_identity_to_bait_pct"]
    }
    print(f"  selected by measured distance to the bait: {', '.join(loci)}", flush=True)
    bait_records = _annotated_records(BAIT_ASSEMBLY, cache, "cds")
    baits = {
        locus: _write_bait(
            _select(bait_records, gene=locus),
            Path(outdir) / f"bait_{locus}_26695.fasta", locus, BAIT_ASSEMBLY)
        for locus in loci
    }

    rows = []
    for entry in manifest:
        strain = entry["strain"]
        sample_dir = workdir / strain
        r1, r2 = _reads(entry["illumina_run"], cache)
        draft = _draft(sample_dir, r1, r2, threads)
        annotation = _annotated_records(entry["assembly_accession"], cache, "cds")
        for locus in loci:
            truth = _select(annotation, gene=locus)
            case_class = (
                "divergent" if selected[locus] < DIVERGENT_IDENTITY_PCT
                else "supported"
            )
            print(f"  [{strain}] reconstructing {locus}", flush=True)
            record, sequence, exit_code = _reconstruct(
                repo_root, sample_dir, strain, locus, baits[locus], draft,
                r1, r2, threads,
            )
            if record is None:
                rows.append([
                    "divergent single-copy loci", f"{strain} {locus}", case_class,
                    f"NO_REPORT(exit={exit_code})", "", "", "no",
                    "annotated closed-genome truth", "not assessed", "",
                ])
                continue
            exact, identity, length, relation = _score(
                sequence, truth, sample_dir)
            rows.append([
                "divergent single-copy loci", f"{strain} {locus}", case_class,
                record.get("workflow_status", ""),
                record.get("result_disposition", ""),
                record.get("sequence_confidence") or record.get("qc_confidence", ""),
                exact, _evidence(sequence, truth, identity, length, record),
                relation, record.get("qc_flags", ""),
            ])
    write_tsv(
        Path(outdir) / "tier_calibration_divergent_loci.tsv",
        ["stage", "case", "case_class", "workflow_status", "disposition",
         "sequence_confidence", "matches_truth", "evidence",
         "relation_to_annotation", "qc_flags"],
        rows,
    )
    return summarise(rows, stage="divergent single-copy loci"), rows


def _read_stats(path):
    """Read count and read length of one FASTQ.gz, without loading it."""
    lines = int(subprocess.run(
        f"zcat {path} | wc -l", shell=True, capture_output=True, text=True, check=True
    ).stdout.split()[0])
    first = subprocess.run(
        f"zcat {path} | head -2 | tail -1", shell=True, capture_output=True,
        text=True, check=True,
    ).stdout.strip()
    return lines // 4, len(first)


def _subsample(source, destination, keep_every):
    """Keep every k-th read, deterministically and in step across the pair.

    The record test is ``(NR-1) % (4k) < 4``, not ``NR % (4k) < 4``: FASTQ
    records start at line 1, so the unshifted form straddles record
    boundaries and emits a file that is four lines per record but not four
    lines of the same record.  The result is accepted by gzip and rejected by
    every assembler, so the output is checked here rather than downstream.
    """
    destination = Path(destination)
    if destination.is_file() and destination.stat().st_size:
        return str(destination)
    subprocess.run(
        f"zcat {source} | awk '(NR-1)%({4 * keep_every})<4' | gzip -1 > "
        f"{destination}",
        shell=True, check=True,
    )
    lines = int(subprocess.run(
        f"zcat {destination} | wc -l", shell=True, capture_output=True,
        text=True, check=True).stdout.split()[0])
    first = subprocess.run(
        f"zcat {destination} | head -1", shell=True, capture_output=True,
        text=True, check=True).stdout
    if lines % 4 or not first.startswith("@"):
        raise RuntimeError(
            f"subsampled FASTQ is malformed: {lines} lines, first line "
            f"{first[:20]!r}"
        )
    return str(destination)


def stage_d_depth(repo_root, outdir, manifest_path, workdir, threads, depths,
                  locus="gyrB", isolate=False, filename=None, stage=None):
    """The same accepted locus at reduced depth, assembled and reconstructed.

    The depth series answers the other half of the acceptance question: the
    top tier was reached on these genomes at full depth, so at what depth does
    it stop being reached, and does the tool withhold rather than accept a
    degraded reconstruction.  Reads are subsampled by keeping every k-th pair,
    which is deterministic and preserves the pairing, and the draft is
    reassembled at each depth so the loss of depth is felt by the assembly as
    well as by the reconstruction.
    """
    manifest = read_tsv(manifest_path)
    workdir = Path(workdir)
    cache = workdir / "downloads"
    stage_label = stage or "depth series"
    bait = _write_bait(
        _select(_annotated_records(BAIT_ASSEMBLY, cache, "cds"), gene=locus),
        Path(outdir) / f"bait_{locus}_26695.fasta", locus, BAIT_ASSEMBLY)

    rows = []
    for entry in manifest:
        strain = entry["strain"]
        sample_dir = workdir / strain
        r1, r2 = _reads(entry["illumina_run"], cache)
        full_draft = _draft(sample_dir, r1, r2, threads)
        genome_bp = sum(
            len(sequence) for _, sequence in read_fasta_sequences(full_draft)
        )
        pairs, read_bp = _read_stats(r1)
        observed = pairs * 2 * read_bp / genome_bp
        truth = _select(
            _annotated_records(entry["assembly_accession"], cache, "cds"), gene=locus)
        print(f"  [{strain}] {observed:.0f}x observed over {genome_bp/1e6:.2f} Mb",
              flush=True)
        for target in depths:
            keep_every = max(2, round(observed / target))
            achieved = observed / keep_every
            suffix = "_isolate" if isolate else ""
            tag = f"_depth{target}x{suffix}"
            sub1 = _subsample(r1, sample_dir / f"sub{target}x_R1.fastq.gz", keep_every)
            sub2 = _subsample(r2, sample_dir / f"sub{target}x_R2.fastq.gz", keep_every)
            print(f"  [{strain}] {target}x arm: every {keep_every}th pair "
                  f"({achieved:.1f}x), assembling", flush=True)
            try:
                draft = _draft(sample_dir, sub1, sub2, threads,
                               name=f"spades{tag}", isolate=isolate)
            except RuntimeError as error:
                rows.append([
                    stage_label, f"{strain} {locus} {target}x{suffix}", "low-depth",
                    "ASSEMBLY_FAILED", "", "", "no", str(error)[:160],
                    "not assessed", "",
                ])
                continue
            record, sequence, exit_code = _reconstruct(
                repo_root, sample_dir, strain, locus, bait, draft, sub1, sub2,
                threads, tag=tag)
            if record is None:
                rows.append([
                    stage_label, f"{strain} {locus} {target}x{suffix}", "low-depth",
                    f"NO_REPORT(exit={exit_code})", "", "", "no",
                    f"target {target}x, achieved {achieved:.1f}x",
                    "not assessed", "",
                ])
                continue
            exact, identity, length, relation = _score(
                sequence, truth, sample_dir)
            rows.append([
                stage_label, f"{strain} {locus} {target}x{suffix}", "low-depth",
                record.get("workflow_status", ""),
                record.get("result_disposition", ""),
                record.get("sequence_confidence") or record.get("qc_confidence", ""),
                exact,
                (f"target {target}x, achieved {achieved:.1f}x; "
                 + _evidence(sequence, truth, identity, length, record)),
                relation, record.get("qc_flags", ""),
            ])
    write_tsv(
        Path(outdir) / (filename or "tier_calibration_depth_series.tsv"),
        ["stage", "case", "case_class", "workflow_status", "disposition",
         "sequence_confidence", "matches_truth", "evidence",
         "relation_to_annotation", "qc_flags"],
        rows,
    )
    return summarise(rows, stage="depth series"), rows


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate the confidence tiers against known truth.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--outdir", default="validation/tier-calibration")
    parser.add_argument("--manifest", help="Closed-genome manifest for stage B.")
    parser.add_argument("--workdir", default="tier_calibration_work")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--divergent-loci", type=int, default=0, metavar="N",
        help="Stage C: add the N eligible single-copy loci furthest from the bait.",
    )
    parser.add_argument(
        "--depth-isolate", action="store_true",
        help=("repeat the depth series with SPAdes --isolate at every depth. "
              "Reported as a control on the assembler mode, not as additional "
              "calibration cases: the primary series keeps the default mode."),
    )
    parser.add_argument(
        "--depth-series", default="", metavar="X,Y",
        help="Stage D: target depths for the subsampled arm, e.g. 15,8.",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = {"stage_a_benchmarks": stage_a(REPO_ROOT, outdir)}
    print("[A] deposited benchmarks: "
          f"{summary['stage_a_benchmarks']['overall']['cases']} cases, "
          f"{summary['stage_a_benchmarks']['overall']['false_accepts']} false accepts")
    for name, block in summary["stage_a_benchmarks"]["by_case_class"].items():
        print(f"    {name:10s} {block['reached_top_tier']}/{block['cases']} at "
              f"{TOP_TIER}, precision {block['top_tier_precision']}")

    inputs = [
        REPO_ROOT / "validation" / "reference_results" / "validation_results.tsv",
        REPO_ROOT / "validation" / "completeness-benchmark"
        / "completeness_benchmark_results.tsv",
    ]
    if args.manifest:
        block, rows = stage_b(
            REPO_ROOT, outdir, args.manifest, args.workdir, args.threads
        )
        summary["stage_b_closed_genomes"] = block
        print(f"[B] closed genomes: {block['overall']['cases']} cases, "
              f"{block['overall']['reached_top_tier']} at {TOP_TIER}, "
              f"{block['overall']['false_accepts']} false accepts")
        inputs.append(Path(args.manifest))

    if args.manifest and args.divergent_loci:
        block, divergent_rows = stage_c_divergence(
            REPO_ROOT, outdir, args.manifest, args.workdir, args.threads,
            args.divergent_loci,
        )
        summary["stage_c_divergent_loci"] = block
        print(f"[C] divergent loci: {block['overall']['cases']} cases, "
              f"{block['overall']['reached_top_tier']} at {TOP_TIER}, "
              f"{block['overall']['false_accepts']} false accepts")

    if args.manifest and args.depth_series:
        depths = [int(value) for value in args.depth_series.split(",")]
        block, depth_rows = stage_d_depth(
            REPO_ROOT, outdir, args.manifest, args.workdir, args.threads, depths,
        )
        summary["stage_d_depth_series"] = block
        print(f"[D] depth series: {block['overall']['cases']} cases, "
              f"{block['overall']['reached_top_tier']} at {TOP_TIER}, "
              f"{block['overall']['false_accepts']} false accepts")

        if args.depth_isolate:
            # Same subsampled libraries, assembled with --isolate at every
            # depth. Kept out of the calibration aggregate: it answers whether
            # the tier response follows depth or the assembler mode, and
            # counting it as further cases would double the same libraries.
            control, control_rows = stage_d_depth(
                REPO_ROOT, outdir, args.manifest, args.workdir, args.threads,
                depths, isolate=True,
                filename="tier_calibration_depth_series_isolate.tsv",
                stage="depth series, --isolate control",
            )
            summary["stage_d_isolate_control"] = control
            agreement = sum(
                1 for primary, repeated in zip(depth_rows, control_rows)
                if primary[5] == repeated[5]
            )
            summary["stage_d_isolate_control"]["tier_agreement_with_primary"] = (
                f"{agreement}/{len(control_rows)}"
            )
            print(f"[D-control] --isolate at every depth: "
                  f"{control['overall']['cases']} cases, "
                  f"{control['overall']['reached_top_tier']} at {TOP_TIER}, "
                  f"tier identical to the primary series in "
                  f"{agreement}/{len(control_rows)}")

    write_json(outdir / "TIER_CALIBRATION.json", summary)
    write_json(
        outdir / "PROVENANCE.json",
        build_provenance(
            workflow="tier_calibration",
            repo_root=REPO_ROOT,
            inputs=inputs,
            parameters={
                "top_tier": TOP_TIER,
                "top_disposition": TOP_DISPOSITION,
                "bait_assembly": BAIT_ASSEMBLY,
                "mock_case_classes": MOCK_CLASSES,
                "threads": args.threads,
                "stage_b_run": bool(args.manifest),
            },
            tools={
                "python": [sys.executable, "--version"],
                "blastn": ["blastn", "-version"],
                "spades": ["spades.py", "--version"],
            },
        ),
    )
    write_tsv(
        outdir / "INPUT_CHECKSUMS.tsv",
        ["path", "size_bytes", "sha256"],
        [[str(p), Path(p).stat().st_size, sha256_file(p)] for p in inputs],
    )
    print(f"\nWrote {outdir}/TIER_CALIBRATION.json")


if __name__ == "__main__":
    main()
