#!/usr/bin/env python3
"""Compact record of one run of the simulated benchmark.

``validation/reference_results/`` holds three files of a run of
``run_mock_validation.py``: the batch report, the per-sample evaluation and
its summary.  This script copies them from a results directory.  In the batch
report every path under the results directory is written relative to it,
because an absolute path describes the machine and not the result.
"""

import argparse
import csv
import os
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--locus", default="mockLocus")
    args = parser.parse_args()

    results = args.results_dir.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    report = results / f"locus_recon_report_{args.locus}.tsv"
    with report.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    prefix = str(results) + os.sep
    changed = 0
    for row in rows[1:]:
        for index, value in enumerate(row):
            if value.startswith(prefix):
                row[index] = os.path.relpath(value, results)
                changed += 1
    with (args.out / report.name).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)
    for name in ("validation_results.tsv", "validation_summary.json"):
        shutil.copyfile(results / "validation" / name, args.out / name)
    print(f"{len(rows) - 1} report rows, {changed} path(s) made relative; "
          f"written to {args.out}")


if __name__ == "__main__":
    main()
