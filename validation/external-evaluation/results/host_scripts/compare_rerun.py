#!/usr/bin/env python3
"""Compare the tables of a second, independent run of one isolate with those of the main run.

    python compare_rerun.py MAIN_OUTDIR RERUN_OUTDIR STRAIN OUTPUT.tsv
"""
import csv
import sys
from pathlib import Path

COLUMNS = {
    "external_units.tsv": (("locus", "tool"), (
        "sequences_returned", "returned_bp", "returned_sha256", "comparison", "decision",
        "decision_secondary", "disposition", "tier", "qc_flags", "ariba_flag",
        "kma_template_coverage", "lowercase_bases", "uncalled_bases")),
    "external_copy_number.tsv": (("locus",), (
        "locus_regions", "gene_spans", "ratio", "ratio_ci_low", "ratio_ci_high", "depth_call",
        "dosage_estimate", "dosage_status", "graph_context_count", "graph_context_lower_bound",
        "graph_context_status", "copy_number_call", "copy_number_kind", "copy_number_status")),
    "block_size_sensitivity.tsv": (("locus", "block_bp"), (
        "ratio", "ratio_ci_low", "ratio_ci_high", "dosage_ci_low", "dosage_ci_high")),
}


def read(path, strain):
    with open(path, newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t") if row["strain"] == strain]


def main():
    first, second, strain, output = sys.argv[1:5]
    rows, different = [], 0
    for table, (keys, columns) in COLUMNS.items():
        main_rows = {tuple(row[key] for key in keys): row for row in read(Path(first) / table, strain)}
        rerun_rows = {tuple(row[key] for key in keys): row for row in read(Path(second) / table, strain)}
        for key in sorted(set(main_rows) | set(rerun_rows)):
            for column in columns:
                one = main_rows.get(key, {}).get(column, "missing")
                two = rerun_rows.get(key, {}).get(column, "missing")
                different += one != two
                rows.append({"strain": strain, "table": table, "row": " ".join(key),
                             "column": column, "main_run": one, "second_run": two,
                             "same": "yes" if one == two else "no"})
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} values compared for {strain}; {different} differ")


if __name__ == "__main__":
    main()
