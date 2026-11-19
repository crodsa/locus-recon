#!/usr/bin/env python3
"""Apply the copy-number consensus of the installed release to stored results.

    python validation/relabel_consensus.py --out validation/consensus-relabel-1.0.2

Release 1.0.2 changed one label of ``reconcile_copy_number``: a graph context
count of two or more that the bootstrap interval of the depth ratio excludes
is returned with the status ``REVIEW`` and as a lower bound (1.0.1:
``SUPPORTED``).  The external evaluation and the targeted tests were run with
release 1.0.1 and their records are kept as they were written.  This script
reads, for every unit of those records, the stored depth result and the stored
graph result, calls ``reconcile_copy_number`` of the installed package on
them, and writes the stored consensus beside the recomputed one.

``reconcile_copy_number`` is a pure function of six fields of the depth
result (``reliable``, ``depth_call``, ``dosage_estimate``, ``dosage_status``,
``ratio_ci_low``, ``ratio_ci_high``) and three of the graph result
(``status``, ``context_count``, ``context_lower_bound``).  With release 1.0.1
installed the script must reproduce every stored consensus; ``--expect-identical``
makes it fail otherwise.  That run is the check that the stored fields are
sufficient (``check_with_1.0.1.txt`` in the output folder).

The change was made after the results of the evaluation were known.  The
table shows what it does to those results; it is not a test on new data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from locus_recon import VERSION  # noqa: E402
from locus_recon.copy_number import reconcile_copy_number  # noqa: E402
from validation.common import source_sha256  # noqa: E402

VALIDATION = REPO_ROOT / "validation"
FIELDS = ("copy_number_call", "copy_number_kind", "copy_number_method",
          "copy_number_status", "copy_number_lower_bound",
          "copy_number_upper_bound", "consensus_flags")


def _number(text):
    text = (text or "").strip()
    if text in ("", "nan", "NaN", "None"):
        return None
    value = float(text)
    return int(value) if value == int(value) else value


def _truth(text):
    return _number(text)


def _bool(text):
    return str(text).strip().lower() in ("true", "1", "yes")


def _depth(row, prefix=""):
    get = lambda name: row.get(prefix + name, "")  # noqa: E731
    estimate = get("dosage_estimate").strip()
    return SimpleNamespace(
        reliable=_bool(get("reliable")),
        depth_call=get("depth_call").strip(),
        dosage_estimate=None if estimate in ("", "nan", "None") else float(estimate),
        dosage_status=get("dosage_status").strip(),
        ratio_ci_low=float(get("ratio_ci_low")),
        ratio_ci_high=float(get("ratio_ci_high")),
    )


def _graph(row):
    status = (row.get("graph_context_status") or "").strip()
    if status in ("", "NOT_EVALUATED", "nan"):
        return None
    return SimpleNamespace(
        status=status,
        context_count=_number(row.get("graph_context_count")),
        context_lower_bound=_number(row.get("graph_context_lower_bound")),
    )


def _text(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return ";".join(value)
    if isinstance(value, float):
        return repr(round(value, 6)).rstrip("0").rstrip(".") if value != int(value) else str(int(value))
    return str(value)


def _same(stored, new):
    """Stored text against a recomputed value; numbers are compared as numbers."""
    stored = (stored or "").strip()
    if stored in ("nan", "NaN", "None"):
        stored = ""
    new = _text(new)
    if stored == new:
        return True
    try:
        return abs(float(stored) - float(new)) < 5e-4
    except ValueError:
        return False


def _against_truth(result, truth):
    """Scoring of the external evaluation (protocol 2.6), applied to a result."""
    if truth is None:
        return "truth not defined"
    if result.copy_number_kind == "INTEGER_CONTEXT_COUNT":
        return "correct" if result.copy_number_call == truth else "wrong"
    if result.copy_number_kind == "MEAN_DEPTH_DOSAGE":
        return "correct" if round(result.copy_number_call) == truth else "wrong"
    if result.copy_number_kind == "LOWER_BOUND":
        return "correct" if result.copy_number_lower_bound <= truth else "wrong"
    return "no value"


def _units():
    """Yield (source, unit, set, truth copies, depth, graph, stored consensus)."""
    path = VALIDATION / "external-evaluation/results/external_copy_number.tsv"
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row["ratio_ci_low"].strip():
                continue
            unit_set = "protocol set" if row["unit_set"] == "protocol" else "supplementary set"
            yield (f"external evaluation, {unit_set}", f"{row['strain']} {row['locus']}",
                   f"{row['unit_set']}; {row['stratum']}", _truth(row["truth_copies"]),
                   _depth(row), _graph(row), {f: row[f] for f in FIELDS},
                   row.get("reads_agree_with_genome", ""))
    for name, label in (("test4_graph_contexts.tsv", "targeted test 4"),
                        ("test5_overlapping_spans.tsv", "targeted test 5")):
        path = VALIDATION / "targeted-tests" / name
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if not row.get("dr_ratio_ci_low", "").strip():
                    continue
                truth = _truth(row.get("true_copies_per_genome") or row.get("truth_copies"))
                depth = _depth(row, "dr_")
                yield (label + ", with the graph", row["run_id"], row["scenario"], truth,
                       depth, _graph(row), {f: row["cn_" + f] for f in FIELDS}, "")
                yield (label + ", without the graph", row["run_id"], row["scenario"], truth,
                       depth, None, {f: row["cn_nograph_" + f] for f in FIELDS}, "")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expect-identical", action="store_true",
                        help="fail when a recomputed field differs from the stored one")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rows, differing = [], 0
    for source, unit, group, truth, depth, graph, stored, agree in _units():
        result = reconcile_copy_number(depth, graph)
        new = {f: getattr(result, f) for f in FIELDS}
        changed = [f for f in FIELDS if not _same(stored[f], new[f])]
        differing += bool(changed)
        stored_result = SimpleNamespace(
            copy_number_kind=stored["copy_number_kind"],
            copy_number_call=_number(stored["copy_number_call"]),
            copy_number_lower_bound=_number(stored["copy_number_lower_bound"]),
        )
        row = {"source": source, "unit": unit, "group": group,
               "truth_copies": _text(truth), "reads_agree_with_genome": agree,
               "ratio_ci_low": depth.ratio_ci_low, "ratio_ci_high": depth.ratio_ci_high,
               "depth_call": depth.depth_call,
               "graph_context_status": graph.status if graph else "NOT_EVALUATED",
               "graph_context_count": _text(graph.context_count) if graph else ""}
        for f in FIELDS:
            row["stored_" + f] = "" if stored[f].strip() in ("nan", "None") else stored[f].strip()
        for f in FIELDS:
            row["recomputed_" + f] = _text(new[f])
        row["fields_that_differ"] = ";".join(changed)
        row["stored_against_truth"] = _against_truth(stored_result, truth)
        row["recomputed_against_truth"] = _against_truth(result, truth)
        bound = result.copy_number_lower_bound
        row["recomputed_lower_bound_holds"] = (
            "" if bound is None or truth is None else ("yes" if bound <= truth else "no"))
        rows.append(row)

    table = args.out / "consensus_relabel.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary = {}
    for row in rows:
        for which in ("stored", "recomputed"):
            key = (row["source"].split(",")[0],
                   row["source"], which, row[which + "_copy_number_status"],
                   row[which + "_against_truth"])
            summary[key] = summary.get(key, 0) + 1
    with (args.out / "consensus_relabel_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["analysis", "source", "consensus", "copy_number_status",
                         "against_truth", "units"])
        for key in sorted(summary):
            writer.writerow([*key, summary[key]])

    record = {
        "software": {"name": "locus-recon", "version": VERSION,
                     "source_sha256": source_sha256(REPO_ROOT)},
        "units": len(rows),
        "units_with_a_field_that_differs": differing,
        "inputs": {
            str(p.relative_to(REPO_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (VALIDATION / "external-evaluation/results/external_copy_number.tsv",
                      VALIDATION / "targeted-tests/test4_graph_contexts.tsv",
                      VALIDATION / "targeted-tests/test5_overlapping_spans.tsv")},
    }
    name = "PROVENANCE.json" if not args.expect_identical else "check_identical.json"
    (args.out / name).write_text(json.dumps(record, indent=2) + "\n")
    print(f"locus-recon {VERSION}: {len(rows)} units, {differing} with a field that differs "
          f"from the stored consensus")
    if args.expect_identical and differing:
        raise SystemExit("stored consensus not reproduced")


if __name__ == "__main__":
    main()
