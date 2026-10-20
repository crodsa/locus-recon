#!/usr/bin/env python3
"""Compare Locus-Recon mock-run outputs with the deterministic truth set."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, Tuple


def read_fasta(path: Path) -> Dict[str, str]:
    records = {}
    header = None
    sequence = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence).upper()
                header = line[1:].split()[0]
                sequence = []
            elif line:
                sequence.append(line)
    if header is not None:
        records[header] = "".join(sequence).upper()
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def edit_distance(left: str, right: str) -> int:
    if len(left) < len(right):
        left, right = right, left
    previous = list(range(len(right) + 1))
    for left_index, left_base in enumerate(left, 1):
        current = [left_index]
        for right_index, right_base in enumerate(right, 1):
            current.append(min(
                current[-1] + 1,
                previous[right_index] + 1,
                previous[right_index - 1] + (left_base != right_base),
            ))
        previous = current
    return previous[-1]


def sequence_concordance(predicted: str, truths: Iterable[str]) -> Tuple[float, bool]:
    best_identity = 0.0
    exact = False
    for truth in truths:
        for oriented in (predicted, reverse_complement(predicted)):
            denominator = max(len(oriented), len(truth), 1)
            identity = 100.0 * (1.0 - edit_distance(oriented, truth) / denominator)
            best_identity = max(best_identity, identity)
            exact = exact or oriented == truth
    return best_identity, exact


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "yes", "1"}


def load_tsv(path: Path) -> list:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def evaluate(dataset_dir: Path, results_dir: Path, output_dir: Path) -> dict:
    expected_rows = load_tsv(dataset_dir / "truth" / "expected_results.tsv")
    result_rows = load_tsv(results_dir / "locus_recon_report_mockLocus.tsv")
    results_by_sample = {row["sample_id"]: row for row in result_rows}
    truth_alleles = read_fasta(dataset_dir / "truth" / "truth_alleles.fasta")

    evaluated_rows = []
    for expected in expected_rows:
        sample_id = expected["sample_id"]
        observed = results_by_sample.get(sample_id, {})
        status = observed.get("status", "MISSING")
        predicted_sequence = ""
        allele_path_text = observed.get("allele_file", "")
        if status == "SUCCESS" and allele_path_text:
            allele_path = Path(allele_path_text)
            if allele_path.exists():
                predicted_records = read_fasta(allele_path)
                predicted_sequence = next(iter(predicted_records.values()), "")

        allowed_ids = [
            allele_id for allele_id in (
                expected["truth_allele"], expected["secondary_allele"]
            ) if allele_id
        ]
        allowed_sequences = [truth_alleles[allele_id] for allele_id in allowed_ids]
        if predicted_sequence and allowed_sequences:
            identity, exact = sequence_concordance(predicted_sequence, allowed_sequences)
        else:
            identity, exact = 0.0, False

        mixture_detected = as_bool(observed.get("mixture_detected", ""))
        flags = observed.get("qc_flags", "")
        ambiguity_flagged = "AMBIGUOUS_BEST_HIT" in flags
        evaluated_rows.append({
            "sample_id": sample_id,
            "category": expected["category"],
            "nominal_depth": float(expected["nominal_depth"]),
            "mixture_fraction": float(expected["mixture_fraction"]),
            "expected_outcome": expected["expected_outcome"],
            "status": status,
            "qc_confidence": observed.get("qc_confidence", ""),
            "truth_identity_pct": identity,
            "exact_truth_match": exact,
            "mixture_detected": mixture_detected,
            "mixed_site_count": int(observed.get("mixed_site_count") or 0),
            "median_mixed_fraction": float(observed.get("median_mixed_fraction") or 0),
            "strand_biased_sites": int(observed.get("strand_biased_sites") or 0),
            "ambiguity_flagged": ambiguity_flagged,
            "remap_mean_depth": float(observed.get("remap_mean_depth") or 0),
            "remap_breadth_pct": float(observed.get("remap_breadth_pct") or 0),
            "message": observed.get("message", "result row missing"),
        })

    baseline = [row for row in evaluated_rows if row["category"] == "baseline"]
    fragmented = [row for row in evaluated_rows if row["category"] == "fragmented"]
    pure = [
        row for row in evaluated_rows
        if row["category"] in {"baseline", "fragmented", "depth_series"}
        and row["status"] == "SUCCESS"
    ]
    mixtures = [row for row in evaluated_rows if row["category"] == "mixture"]
    mixtures_10plus = [row for row in mixtures if row["mixture_fraction"] >= 0.10]
    paralog = next(row for row in evaluated_rows if row["category"] == "paralog")
    negative = next(row for row in evaluated_rows if row["category"] == "negative")
    depth_series = sorted(
        (row for row in evaluated_rows if row["category"] == "depth_series"),
        key=lambda row: row["nominal_depth"],
    )

    criteria = {
        "baseline_exact": bool(baseline) and all(
            row["status"] == "SUCCESS" and row["exact_truth_match"] for row in baseline
        ),
        "fragmented_exact": bool(fragmented) and all(
            row["status"] == "SUCCESS" and row["exact_truth_match"] for row in fragmented
        ),
        "no_false_mixture_in_successful_pure_controls": all(
            not row["mixture_detected"] for row in pure
        ),
        "mixtures_10pct_and_above_detected": bool(mixtures_10plus) and all(
            row["mixture_detected"] for row in mixtures_10plus
        ),
        "paralog_ambiguity_flagged": paralog["ambiguity_flagged"],
        "negative_not_reconstructed": negative["status"] != "SUCCESS",
    }
    summary = {
        "schema_version": "1.0",
        "criteria": criteria,
        "all_core_criteria_passed": all(criteria.values()),
        "mixture_detection": {
            "detected": sum(row["mixture_detected"] for row in mixtures),
            "total": len(mixtures),
            "by_fraction": {
                f"{row['mixture_fraction']:.2f}": row["mixture_detected"]
                for row in mixtures
            },
        },
        "depth_series": [
            {
                "nominal_depth": row["nominal_depth"],
                "status": row["status"],
                "exact_truth_match": row["exact_truth_match"],
                "observed_remap_depth": row["remap_mean_depth"],
                "breadth_pct": row["remap_breadth_pct"],
            }
            for row in depth_series
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "validation_results.tsv").open("w", newline="") as handle:
        columns = list(evaluated_rows[0])
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in evaluated_rows:
            writer.writerow(row)
    with (output_dir / "validation_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    validation_dir = Path(__file__).resolve().parent
    parser.add_argument(
        "--dataset-dir", type=Path, default=validation_dir / "mock_dataset"
    )
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir or args.results_dir / "validation"
    summary = evaluate(
        args.dataset_dir.resolve(), args.results_dir.resolve(), output_dir.resolve()
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    raise SystemExit(0 if summary["all_core_criteria_passed"] else 1)


if __name__ == "__main__":
    main()
