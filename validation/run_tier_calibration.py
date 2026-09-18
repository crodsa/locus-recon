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
TOP_TIER = "HIGH"
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


def stage_b(repo_root, outdir, manifest_path, workdir, threads):
    """Repeat the closed-genome audit on public data and tabulate the tiers."""
    manifest = read_tsv(manifest_path)
    workdir = Path(workdir)
    cache = workdir / "downloads"
    workdir.mkdir(parents=True, exist_ok=True)

    bait_cds = _annotated_records(BAIT_ASSEMBLY, cache, "cds")
    bait_rna = _annotated_records(BAIT_ASSEMBLY, cache, "rna")
    baits = {}
    for locus, records in (
        ("gyrB", _select(bait_cds, gene="gyrB")),
        ("23S", _select(bait_rna, product="23S ribosomal RNA")),
    ):
        if not records:
            raise RuntimeError(f"no {locus} record in the bait assembly {BAIT_ASSEMBLY}")
        path = Path(outdir) / f"bait_{locus}_26695.fasta"
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as handle:
            name, sequence = records[0]
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")
        baits[locus] = str(path)

    rows = []
    for entry in manifest:
        strain, run = entry["strain"], entry["illumina_run"]
        sample_dir = workdir / strain
        sample_dir.mkdir(exist_ok=True)
        print(f"  [{strain}] reads {run}", flush=True)
        r1, r2 = _reads(run, cache)

        draft = sample_dir / "spades" / "contigs.fasta"
        if not draft.is_file():
            print(f"  [{strain}] assembling", flush=True)
            run_checked([
                "spades.py", "--isolate", "-1", r1, "-2", r2,
                "-o", str(sample_dir / "spades"), "-t", str(threads), "-m", "32",
            ], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

        truth = {
            "gyrB": _select(
                _annotated_records(entry["assembly_accession"], cache, "cds"),
                gene="gyrB",
            ),
            "23S": _select(
                _annotated_records(entry["assembly_accession"], cache, "rna"),
                product="23S ribosomal RNA",
            ),
        }

        for locus in ("gyrB", "23S"):
            out = sample_dir / locus
            samplesheet = sample_dir / f"samples_{locus}.tsv"
            samplesheet.write_text(
                "sample_id\tassembly_path\tr1_path\tr2_path\n"
                f"{strain}\t{draft}\t{r1}\t{r2}\n"
            )
            print(f"  [{strain}] reconstructing {locus}", flush=True)
            completed = subprocess.run(
                [sys.executable, "-m", "locus_recon.cli",
                 "--samplesheet", str(samplesheet), "--main-output-dir", str(out),
                 "--bait", baits[locus], "--locus", locus,
                 "--threads", str(threads), "--memory-per-sample", "32",
                 "--no-progress"],
                cwd=repo_root, capture_output=True, text=True, check=False,
            )
            reports = list(out.glob(f"locus_recon_report_{locus}.tsv"))
            if not reports:
                rows.append([
                    "closed-genome audit", f"{strain} {locus}",
                    "supported" if len(truth[locus]) == 1 else "multi-copy",
                    f"NO_REPORT(exit={completed.returncode})", "", "", "no",
                    "annotated closed-genome truth",
                ])
                continue
            record = read_tsv(reports[0])[0]
            candidates = list(out.glob(f"{strain}/*_reconstructed.fasta"))
            sequence = (
                read_fasta_sequences(str(candidates[0]))[0][1] if candidates else ""
            )
            exact, best_identity, best_length = "no", 0.0, 0
            for _, truth_sequence in truth[locus]:
                if sequence and sequence in (truth_sequence, reverse_complement(truth_sequence)):
                    exact = "yes"
                identity, length = _blast_identity(sequence, truth_sequence, sample_dir)
                if identity * length > best_identity * best_length:
                    best_identity, best_length = identity, length
            rows.append([
                "closed-genome audit", f"{strain} {locus}",
                # The class follows the annotation, not the result: a locus the
                # closed genome carries twice cannot be attributed to a single
                # copy from short reads, so withholding is correct behaviour.
                "supported" if len(truth[locus]) == 1 else "multi-copy",
                record.get("workflow_status", ""),
                record.get("result_disposition", ""),
                record.get("sequence_confidence") or record.get("qc_confidence", ""),
                exact,
                (f"truth {len(truth[locus][0][1])} bp x {len(truth[locus])} copy/ies; "
                 f"candidate {len(sequence)} bp at {best_identity:.2f}% over "
                 f"{best_length} bp; catalogue identity "
                 f"{record.get('best_identity_pct', '')}%"),
            ])

    write_tsv(
        Path(outdir) / "tier_calibration_closed_genomes.tsv",
        ["stage", "case", "case_class", "workflow_status", "disposition",
         "sequence_confidence", "matches_truth", "evidence"],
        rows,
    )
    return summarise(rows, stage="closed genomes"), rows


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate the confidence tiers against known truth.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--outdir", default="validation/tier-calibration")
    parser.add_argument("--manifest", help="Closed-genome manifest for stage B.")
    parser.add_argument("--workdir", default="tier_calibration_work")
    parser.add_argument("--threads", type=int, default=8)
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
