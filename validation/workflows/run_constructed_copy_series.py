#!/usr/bin/env python3
"""Generate and analyse the frozen constructed 1/2/3-copy series."""

from __future__ import annotations

import argparse
import csv
import gzip
import random
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from locus_recon.depth_ratio import estimate_locus_copy_number
from validation.common import build_provenance, compact_summary, write_json


COPIES = (1, 2, 3)
DEPTHS = (20, 50, 100)
GC_LEVELS = (0.30, 0.50, 0.70)
BASE_SEED = 20260817
READ_LENGTH = 150
INSERT = 350
LEFT = 50_000
LOCUS_LENGTH = 4_000
RIGHT = 50_000


def _sequence(length: int, gc: float, rng: random.Random) -> str:
    weights = (("A", (1 - gc) / 2), ("C", gc / 2), ("G", gc / 2), ("T", (1 - gc) / 2))
    bases, cumulative = [], []
    running = 0.0
    for base, weight in weights:
        bases.append(base)
        running += weight
        cumulative.append(running)
    out = []
    for _ in range(length):
        value = rng.random()
        out.append(next(base for base, boundary in zip(bases, cumulative) if value <= boundary))
    return "".join(out)


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")


def _write_reads(path1: Path, path2: Path, source: str, depth: int, rng: random.Random) -> None:
    pairs = round(depth * len(source) / (2 * READ_LENGTH))
    quality = "I" * READ_LENGTH
    with path1.open("wb") as raw1, path2.open("wb") as raw2:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw1, mtime=0) as gz1, \
                gzip.GzipFile(filename="", mode="wb", fileobj=raw2, mtime=0) as gz2:
            for index in range(pairs):
                start = rng.randrange(0, len(source) - INSERT + 1)
                fragment = source[start:start + INSERT]
                read1 = fragment[:READ_LENGTH]
                read2 = _reverse_complement(fragment[-READ_LENGTH:])
                name = f"@pair_{index:07d}_start_{start}"
                gz1.write(f"{name}/1\n{read1}\n+\n{quality}\n".encode())
                gz2.write(f"{name}/2\n{read2}\n+\n{quality}\n".encode())


def _map(reference: Path, r1: Path, r2: Path, bam: Path, aligner: str, threads: int) -> None:
    if bam.exists() and Path(str(bam) + ".bai").exists():
        return
    index = Path(str(reference) + (".bwt.2bit.64" if aligner == "bwa-mem2" else ".bwt"))
    if not index.exists():
        subprocess.run([aligner, "index", str(reference)], check=True, capture_output=True)
    with bam.with_suffix(".mapping.log").open("w") as log:
        mapping = subprocess.Popen(
            [aligner, "mem", "-t", str(threads), str(reference), str(r1), str(r2)],
            stdout=subprocess.PIPE, stderr=log,
        )
        assert mapping.stdout is not None
        sorting = subprocess.run(
            ["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"],
            stdin=mapping.stdout, stdout=log, stderr=log,
        )
        mapping.stdout.close()
        mapping_code = mapping.wait()
    if mapping_code != 0 or sorting.returncode != 0:
        raise RuntimeError(f"mapping failed for {bam.parent.name}")
    subprocess.run(["samtools", "index", str(bam)], check=True)


_THREE_DECIMALS = (
    "ratio", "ratio_ci_low", "ratio_ci_high", "dosage_estimate", "absolute_error",
    "ambiguity_index",
)


def _write_results_table(path: Path, results: list[dict]) -> None:
    """Write the per-case table with fixed decimal places."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        for row in results:
            formatted = dict(row)
            for key in _THREE_DECIMALS:
                formatted[key] = f"{row[key]:.3f}"
            formatted["locus_gc"] = f"{row['locus_gc']:.2f}"
            writer.writerow(formatted)


def main() -> int:
    repo_root = REPO_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument(
        "--deposit-dir", type=Path,
        help="also write results.tsv, evaluation.json and run_summary.json here "
             "(validation/constructed-copy-series)",
    )
    args = parser.parse_args()
    if min(args.threads, args.n_boot) < 1:
        parser.error("threads and bootstrap count must be positive")
    for tool in ("samtools",):
        if not shutil.which(tool):
            raise FileNotFoundError(f"{tool} is required on PATH")
    aligner = "bwa-mem2" if shutil.which("bwa-mem2") else "bwa"
    if not shutil.which(aligner):
        raise FileNotFoundError("bwa-mem2 or bwa is required on PATH")

    args.work_dir.mkdir(parents=True, exist_ok=True)
    results = []
    inputs: list[Path] = []
    for gc_index, gc in enumerate(GC_LEVELS):
        genome_rng = random.Random(BASE_SEED + gc_index * 1000)
        left = _sequence(LEFT, 0.50, genome_rng)
        locus = _sequence(LOCUS_LENGTH, gc, genome_rng)
        right = _sequence(RIGHT, 0.50, genome_rng)
        reference_sequence = left + locus + right
        for copies in COPIES:
            source = left + locus * copies + right
            for depth in DEPTHS:
                case = f"copies_{copies}__depth_{depth}__gc_{gc:.2f}"
                case_dir = args.work_dir / "cases" / case
                case_dir.mkdir(parents=True, exist_ok=True)
                reference = case_dir / "collapsed_reference.fasta"
                source_path = case_dir / "source_truth.fasta"
                r1, r2 = case_dir / "reads_1.fastq.gz", case_dir / "reads_2.fastq.gz"
                _write_fasta(reference, "constructed_chromosome", reference_sequence)
                _write_fasta(source_path, f"constructed_truth_{copies}_copies", source)
                read_seed = BASE_SEED + gc_index * 10_000 + copies * 100 + depth
                _write_reads(r1, r2, source, depth, random.Random(read_seed))
                bam = case_dir / "reads_vs_collapsed_reference.sorted.bam"
                _map(reference, r1, r2, bam, aligner, args.threads)
                start, end = LEFT + 1, LEFT + LOCUS_LENGTH
                region = f"constructed_chromosome:{start}-{end}"
                result = estimate_locus_copy_number(
                    bam,
                    [region],
                    gene_spans=[(1, LOCUS_LENGTH)],
                    gene_length=LOCUS_LENGTH,
                    max_backbone_bootstrap_blocks=512,
                    n_boot=args.n_boot,
                    seed=BASE_SEED,
                    samtools="samtools",
                )
                results.append({
                    "case": case,
                    "truth_copies": copies,
                    "nominal_depth": depth,
                    "locus_gc": gc,
                    "seed": read_seed,
                    "ratio": result.ratio,
                    "ratio_ci_low": result.ratio_ci_low,
                    "ratio_ci_high": result.ratio_ci_high,
                    "depth_call": result.depth_call,
                    "dosage_estimate": result.dosage_estimate,
                    "absolute_error": round(abs(result.dosage_estimate - copies), 3),
                    "ambiguity_index": result.ambiguity_index,
                })
                inputs.extend([reference, source_path, r1, r2])
                print(
                    f"{case}\testimate={result.dosage_estimate:.3f}\t{result.depth_call}",
                    flush=True,
                )

    results_dir = args.work_dir / "results"
    results_dir.mkdir(exist_ok=True)
    table = results_dir / "constructed_copy_series.tsv"
    _write_results_table(table, results)
    write_json(results_dir / "constructed_copy_series.json", results)
    one_copy_false = sum(
        row["depth_call"] == "MULTICOPY_DEPTH" for row in results if row["truth_copies"] == 1
    )
    multi_missed = sum(
        row["depth_call"] != "MULTICOPY_DEPTH" for row in results if row["truth_copies"] > 1
    )
    errors = sorted(row["absolute_error"] for row in results)
    median_error = round(errors[len(errors) // 2], 3)
    evaluation = {
        "schema_version": "1.0",
        "cases": len(results),
        "one_copy_false_multicopy": one_copy_false,
        "two_or_three_copy_missed": multi_missed,
        "median_absolute_error": median_error,
        "criterion_zero_one_copy_false_multicopy": one_copy_false == 0,
        "criterion_all_two_three_copy_detected": multi_missed == 0,
        "criterion_median_absolute_error_le_0_25": median_error <= 0.25,
    }
    evaluation["all_prespecified_criteria_passed"] = all(
        evaluation[key] for key in evaluation if key.startswith("criterion_")
    )
    write_json(results_dir / "evaluation.json", evaluation)
    provenance = build_provenance(
        workflow="constructed-1-2-3-copy-series",
        repo_root=repo_root,
        inputs=[*inputs, table],
        parameters={
            "copies": list(COPIES),
            "depths": list(DEPTHS),
            "gc_levels": list(GC_LEVELS),
            "base_seed": BASE_SEED,
            "read_length": READ_LENGTH,
            "insert": INSERT,
            "backbone_length": LEFT + RIGHT,
            "locus_length": LOCUS_LENGTH,
            "threads": args.threads,
            "n_boot": args.n_boot,
            "max_backbone_bootstrap_blocks": 512,
        },
        tools={
            "python": [sys.executable, "--version"],
            "aligner": [aligner, "version"],
            "samtools": ["samtools", "--version"],
        },
    )
    provenance["complete"] = True
    provenance["criteria_passed"] = evaluation["all_prespecified_criteria_passed"]
    write_json(results_dir / "workflow_provenance.json", provenance)
    if args.deposit_dir:
        args.deposit_dir.mkdir(parents=True, exist_ok=True)
        _write_results_table(args.deposit_dir / "results.tsv", results)
        write_json(args.deposit_dir / "evaluation.json", evaluation)
        write_json(args.deposit_dir / "run_summary.json", compact_summary(
            provenance,
            criteria_passed=evaluation["all_prespecified_criteria_passed"],
            complete=True,
            note=("Every input is regenerated deterministically by the workflow; "
                  "workflow_provenance.json in its --work-dir holds the checksum "
                  "of each one."),
        ))
    return 0 if evaluation["all_prespecified_criteria_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
