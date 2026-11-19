#!/usr/bin/env python3
"""Sensitivity of the depth-ratio interval to the block size of the bootstrap.

The interval of the depth ratio comes from a block bootstrap with blocks of
500 bp.  This script repeats the estimate with blocks of 100, 250, 500 and
1,000 bp for the seven *aphA1* runs and for the 23S rRNA gene of the five
*Helicobacter pylori* isolates, through the library function
``locus_recon.depth_ratio.estimate_locus_copy_number`` and with the alignments,
regions, seed and number of resamples of the two validation workflows
(``workflows/run_aphA1_validation.py`` and ``workflows/run_low_copy_23S.py``).
Nothing in the package is changed: the block size is the ``block`` argument of
that function.

For each case and block size the table gives the ratio, its interval and the
width of the interval, the dosage and its interval, the number of blocks the
locus is cut into, and whether the interval contains the reference value.  The
reference is the copy number by quantitative PCR for MRSN57 and MRSN58 (McGann
et al. 2014), one copy for the four isolates with a tobramycin MIC of
0.5 mg/L, and two copies for the 23S rRNA gene (closed genomes).  The run
without a published value has no reference.

Usage:
    python validation/run_block_size_sensitivity.py \\
        --aphA1-results work/aphA1/results --aphA1-reference work/aphA1/MRSN56_GCA.fasta \\
        --low-copy-work work/low-copy-23S \\
        --output validation/block-size-sensitivity/block_size_sensitivity.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if (REPO_ROOT / "locus_recon").is_dir() and str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from locus_recon.depth_ratio import estimate_locus_copy_number  # noqa: E402

BLOCKS = (100, 250, 500, 1000)
DEFAULT_BLOCK = 500
SEED = 20260817
N_BOOT = 1000
MAX_BACKBONE_BLOCKS = 512
APHA1_REGION = "CP080452.1:250422-251237"
APHA1_LENGTH = 816
APHA1_EXCLUSIONS = ["CP080453.1", "CP080454.1", "CP080455.1", "CP080456.1"]
HEADER = [
    "dataset", "case", "run", "locus", "block_bp", "default_block",
    "locus_bp", "bootstrap_locus_blocks", "bootstrap_backbone_blocks_total",
    "bootstrap_backbone_blocks_used", "n_boot", "seed",
    "ratio", "ratio_ci_low", "ratio_ci_high", "ratio_ci_width",
    "depth_call", "dosage_estimate", "dosage_ci_low", "dosage_ci_high",
    "dosage_ci_width", "dosage_status",
    "reference_copies", "reference_low", "reference_high", "reference_source",
    "ratio_ci_contains_reference", "ratio_ci_overlaps_reference_range",
    "dosage_ci_contains_reference", "dosage_ci_overlaps_reference_range",
    "ratio_equals_workflow_record", "ratio_ci_equals_workflow_record",
]


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def qpcr(text):
    """Centre and range of a quantitative PCR value written as ``10 (+/-2)``."""
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*\(\+/-\s*(\d+(?:\.\d+)?)\)\s*", text)
    if match:
        centre, error = (float(value) for value in match.groups())
        return centre, centre - error, centre + error, "quantitative PCR, mean +/- SD (McGann et al. 2014)"
    try:
        value = float(text)
    except ValueError:
        return None
    return value, value, value, "one copy expected (tobramycin MIC 0.5 mg/L; McGann et al. 2014)"


def aphA1_cases(results_dir, reference, manifest):
    recorded = {row["run"]: row for row in read_tsv(Path(results_dir) / "aphA1_copy_number.tsv")}
    full = {entry["run"]: entry["result"]
            for entry in json.loads((Path(results_dir) / "aphA1_depth_ratio.json").read_text())}
    cases = []
    for row in read_tsv(manifest):
        run = row["illumina_run"]
        cases.append({
            "dataset": "aphA1", "case": row["sample"], "run": run, "locus": "aphA1",
            "locus_bp": APHA1_LENGTH,
            "bam": str(Path(results_dir) / f"{run}.reads_vs_MRSN56.sorted.bam"),
            "arguments": {
                "locus_regions": [APHA1_REGION], "gene_spans": [(1, APHA1_LENGTH)],
                "gene_length": APHA1_LENGTH, "reference_fasta": str(reference),
                "exclude_regions": APHA1_EXCLUSIONS,
            },
            "reference": qpcr(row["truth_qpcr_copies"]),
            "recorded": {"ratio": full[run].get("ratio"),
                         "ratio_ci_low": recorded[run]["ratio_ci_low"],
                         "ratio_ci_high": recorded[run]["ratio_ci_high"]},
        })
    return cases


def low_copy_cases(work_dir, locus="23S"):
    work_dir = Path(work_dir)
    cases = []
    for row in read_tsv(work_dir / "results" / "low_copy_results.tsv"):
        if row["locus"] != locus:
            continue
        regions = row["locus_regions"].split(";")
        spans = [tuple(int(value) for value in span.split("-"))
                 for span in row["region_gene_spans"].split(";")]
        gene_length = None
        for name, sequence in read_fasta(
                next((work_dir / "fixed_baits").glob("*_23S_gyrB.fasta"))):
            if name.split("|")[0] == locus:
                gene_length = len(sequence)
        truth = float(row["truth_copies"])
        cases.append({
            "dataset": "low-copy-23S", "case": row["strain"], "run": row["illumina_run"],
            "locus": locus, "locus_bp": gene_length,
            "bam": str(work_dir / "samples" / row["strain"]
                       / "reads_vs_illumina_draft.sorted.bam"),
            "arguments": {"locus_regions": regions, "gene_spans": spans,
                          "gene_length": gene_length},
            "reference": (truth, truth, truth, "copies annotated in the closed genome"),
            "recorded": {"ratio": row["ratio"], "ratio_ci_low": row["ratio_ci_low"],
                         "ratio_ci_high": row["ratio_ci_high"]},
        })
    return cases


def read_fasta(path):
    records, name, parts = [], None, []
    with open(path) as handle:
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


def number(value):
    return None if value in (None, "") else float(value)


def contains(low, high, value):
    if None in (low, high, value):
        return ""
    return "yes" if low <= value <= high else "no"


def overlaps(low, high, other_low, other_high):
    if None in (low, high, other_low, other_high):
        return ""
    return "yes" if low <= other_high and other_low <= high else "no"


def same(first, second):
    if first in (None, "") or second in (None, ""):
        return ""
    return "yes" if float(first) == float(second) else "no"


def evaluate(task):
    case, block = task
    arguments = dict(case["arguments"])
    regions = arguments.pop("locus_regions")
    result = estimate_locus_copy_number(
        case["bam"], regions, block=block,
        max_backbone_bootstrap_blocks=MAX_BACKBONE_BLOCKS,
        n_boot=N_BOOT, seed=SEED, samtools="samtools", **arguments,
    )
    row = result.as_row()
    ratio_low, ratio_high = number(row.get("ratio_ci_low")), number(row.get("ratio_ci_high"))
    dosage_low, dosage_high = number(row.get("dosage_ci_low")), number(row.get("dosage_ci_high"))
    reference = case["reference"] or (None, None, None, "no published value")
    recorded = case["recorded"]
    at_default = block == DEFAULT_BLOCK
    return {
        "dataset": case["dataset"], "case": case["case"], "run": case["run"],
        "locus": case["locus"], "block_bp": block,
        "default_block": "yes" if at_default else "no",
        "locus_bp": case["locus_bp"],
        "bootstrap_locus_blocks": row.get("bootstrap_locus_blocks", ""),
        "bootstrap_backbone_blocks_total": row.get("bootstrap_backbone_blocks_total", ""),
        "bootstrap_backbone_blocks_used": row.get("bootstrap_backbone_blocks_used", ""),
        "n_boot": N_BOOT, "seed": SEED,
        "ratio": row.get("ratio", ""), "ratio_ci_low": row.get("ratio_ci_low", ""),
        "ratio_ci_high": row.get("ratio_ci_high", ""),
        "ratio_ci_width": (
            "" if None in (ratio_low, ratio_high) else f"{ratio_high - ratio_low:.3f}"),
        "depth_call": row.get("depth_call", ""),
        "dosage_estimate": row.get("dosage_estimate", ""),
        "dosage_ci_low": row.get("dosage_ci_low", ""),
        "dosage_ci_high": row.get("dosage_ci_high", ""),
        "dosage_ci_width": (
            "" if None in (dosage_low, dosage_high) else f"{dosage_high - dosage_low:.3f}"),
        "dosage_status": row.get("dosage_status", ""),
        "reference_copies": "" if reference[0] is None else f"{reference[0]:g}",
        "reference_low": "" if reference[1] is None else f"{reference[1]:g}",
        "reference_high": "" if reference[2] is None else f"{reference[2]:g}",
        "reference_source": reference[3],
        "ratio_ci_contains_reference": contains(ratio_low, ratio_high, reference[0]),
        "ratio_ci_overlaps_reference_range": overlaps(
            ratio_low, ratio_high, reference[1], reference[2]),
        "dosage_ci_contains_reference": contains(dosage_low, dosage_high, reference[0]),
        "dosage_ci_overlaps_reference_range": overlaps(
            dosage_low, dosage_high, reference[1], reference[2]),
        "ratio_equals_workflow_record": (
            same(row.get("ratio"), recorded["ratio"]) if at_default else ""),
        "ratio_ci_equals_workflow_record": (
            ("yes" if same(row.get("ratio_ci_low"), recorded["ratio_ci_low"]) == "yes"
             and same(row.get("ratio_ci_high"), recorded["ratio_ci_high"]) == "yes" else "no")
            if at_default else ""),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--aphA1-results", type=Path,
                        help="output directory of run_aphA1_validation.py")
    parser.add_argument("--aphA1-reference", type=Path, help="reference FASTA of that run")
    parser.add_argument("--aphA1-manifest", type=Path,
                        default=REPO_ROOT / "validation" / "manifests" / "aphA1_copy_number.tsv")
    parser.add_argument("--low-copy-work", type=Path,
                        help="work directory of run_low_copy_23S.py")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1, help="cases evaluated at a time")
    args = parser.parse_args()

    cases = []
    if args.aphA1_results:
        cases += aphA1_cases(args.aphA1_results, args.aphA1_reference, args.aphA1_manifest)
    if args.low_copy_work:
        cases += low_copy_cases(args.low_copy_work)
    if not cases:
        parser.error("give --aphA1-results, --low-copy-work or both")
    tasks = [(case, block) for case in cases for block in BLOCKS]
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        rows = list(pool.map(evaluate, tasks))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} rows ({len(cases)} cases x {len(BLOCKS)} block sizes) "
          f"written to {args.output}")


if __name__ == "__main__":
    main()
