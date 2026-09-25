#!/usr/bin/env python3
"""Gate for numbers the release documentation still waits on.

A slot token such as ``[[CALIB-05]]`` or ``[[CALIB-06 top-tier precision]]``
stands where a number from the tier-calibration stages B to D belongs. It
cannot be mistaken for a value, and this script finds every one.

    python release/slot_gate.py            # exit 1 while any slot remains
    python release/slot_gate.py --list     # one row per distinct slot
    python release/slot_gate.py --values   # the values that close the slots

``--values`` reads ``validation/tier-calibration/`` after the four-stage run
and prints, slot by slot, the quantity to write in. It prints; it never edits
the documentation, because a sentence may need rewording once its number is
known. See ``release/SLOTS.md`` for what closes each slot.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SLOT_REGEX = re.compile(r"\[\[([A-Z][A-Z0-9]*)-(\d+)([^\]]*)\]\]")
SCANNED = ("README.md", "CHANGELOG.md", "CITATION.cff", "docs/**/*.md",
           "examples/**/*.md", "validation/**/*.md")
CALIBRATION = REPO_ROOT / "validation" / "tier-calibration"


def scanned_files(root: Path = REPO_ROOT) -> list[Path]:
    files: set[Path] = set()
    for pattern in SCANNED:
        files.update(path for path in root.glob(pattern) if path.is_file())
    return sorted(files)


def slot_scan(paths) -> list[dict]:
    """Every slot occurrence, scanning whole files so a wrapped token is found."""
    rows = []
    for path in paths:
        text = Path(path).read_text(encoding="utf-8")
        for match in SLOT_REGEX.finditer(text):
            rows.append({
                "file": str(Path(path).relative_to(REPO_ROOT)),
                "line": text.count("\n", 0, match.start()) + 1,
                "slot": f"{match.group(1)}-{match.group(2)}",
                "note": " ".join(match.group(3).split()),
            })
    return rows


def slot_summary(rows: list[dict]) -> list[dict]:
    """One entry per distinct slot, with its occurrences and first note."""
    summary: dict[str, dict] = {}
    for row in rows:
        entry = summary.setdefault(
            row["slot"], {"slot": row["slot"], "note": "", "where": []})
        entry["where"].append(f"{row['file']}:{row['line']}")
        entry["note"] = entry["note"] or row["note"]
    return [summary[key] for key in sorted(summary)]


def _read_tsv(directory: Path, name: str) -> list[dict]:
    path = directory / name
    if not path.is_file():
        raise SystemExit(f"missing {path}; run all four stages")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _catalogue_identity(row: dict) -> float | None:
    match = re.search(r"catalogue identity ([0-9.]+)%", row.get("evidence", ""))
    return float(match.group(1)) if match else None


def _locus(row: dict) -> str:
    return row["case"].split()[1]


def _flags(text: str) -> list[str]:
    """Split a qc_flags cell on the separators outside parentheses."""
    flags, depth, current = [], 0, []
    for character in text:
        depth += {"(": 1, ")": -1}.get(character, 0)
        if character == ";" and depth == 0:
            flags.append("".join(current).strip())
            current = []
        else:
            current.append(character)
    flags.append("".join(current).strip())
    return [flag for flag in flags if flag]


def _span_flags(row: dict) -> str:
    return "; ".join(flag for flag in _flags(row.get("qc_flags", ""))
                     if flag.startswith(("ALLELE_SPAN_CLIPPED", "LOCUS_SPLIT")))


def _truth(rows: list[dict]) -> str:
    """Exact count, then boundary shifts and base-level identities of the rest."""
    exact = sum(row["matches_truth"] == "yes" for row in rows)
    shifts, identities = [], []
    for row in rows:
        relation = row.get("relation_to_annotation", "")
        shift = re.search(r"candidate ([+-]\d+) bp against the annotated boundary", relation)
        if row["matches_truth"] == "yes":
            continue
        if shift:
            shifts.append(int(shift.group(1)))
        elif relation:
            identities.append(relation)
    parts = [f"{exact} of {len(rows)} exact"]
    if shifts:
        parts.append(f"{len(shifts)} identical over the overlap with the boundary "
                     f"{min(shifts):+d} to {max(shifts):+d} bp from the annotation")
    parts += identities
    return "; ".join(parts)


def _tiers(rows: list[dict]) -> str:
    counts = Counter(row["sequence_confidence"] or "none" for row in rows)
    return ", ".join(f"`{tier}` x {n}" for tier, n in sorted(counts.items()))


def _interval(entry: dict) -> str:
    interval = entry.get("top_tier_precision_ci95")
    if not interval:
        return f"{entry['top_tier_precision']:.2f}"
    return f"{entry['top_tier_precision']:.2f} ({interval[0]:.2f}-{interval[1]:.2f})"


def values(directory: Path = CALIBRATION) -> int:
    summary_path = directory / "TIER_CALIBRATION.json"
    summary = json.loads(summary_path.read_text())
    if "combined" not in summary:
        print("Stages still missing: " + ", ".join(summary.get("pending_stages", [])))
        return 1
    combined = summary["combined"]
    overall = combined["overall"]
    supported = combined["by_case_class"]["supported"]
    stage_b = _read_tsv(directory, "tier_calibration_closed_genomes.tsv")
    stage_c = _read_tsv(directory, "tier_calibration_divergent_loci.tsv")
    stage_d = _read_tsv(directory, "tier_calibration_depth_series.tsv")
    gyrb = [row for row in stage_b if _locus(row) == "gyrB"]
    rrna = [row for row in stage_b if _locus(row) == "23S"]
    out = []

    exact_b = sum(row["matches_truth"] == "yes" for row in stage_b)
    out.append(("CALIB-01", f"{exact_b} (of {len(stage_b)})"))
    exact_high = [row for row in gyrb
                  if row["sequence_confidence"] == "HIGH" and row["matches_truth"] == "yes"]
    out.append(("CALIB-02", f"{len(exact_high)} (of {len(gyrb)})"))
    identities = [value for value in map(_catalogue_identity, gyrb) if value is not None]
    out.append(("CALIB-03", f"{min(identities):.1f}-{max(identities):.1f}%"
                if identities else "no catalogue identity recorded"))
    not_exact = [row for row in stage_b if row["matches_truth"] != "yes"]
    out.append(("CALIB-04", "none" if not not_exact else "\n".join(
        f"  {row['case']}: `{row['sequence_confidence']}` / `{row['disposition']}`; "
        f"{row['relation_to_annotation']}; {row['evidence']}"
        + (f"; {_span_flags(row)}" if _span_flags(row) else "")
        for row in not_exact)))
    out.append(("CALIB-05", f"{supported['reached_top_tier']} of {supported['cases']}"))
    out.append(("CALIB-06", _interval(overall) if overall["top_tier_precision"]
                is not None else "no case reached the top tier"))
    out.append(("CALIB-07", f"{overall['false_accepts']} (in {overall['cases']} cases)"))
    withheld = overall["cases"] - overall["reached_top_tier"]
    out.append(("CALIB-08", f"{overall['withheld_but_exact']} of {withheld}"))
    base_level = [row for row in stage_b + stage_c + stage_d
                  if "% identity over" in row.get("relation_to_annotation", "")]
    out.append(("CALIB-09", "none" if not base_level else "\n".join(
        f"  {row['case']} ({row['stage']}): `{row['sequence_confidence']}`; "
        f"{row['relation_to_annotation']}; flags: {row['qc_flags']}"
        for row in base_level)))

    by_depth = defaultdict(list)
    for row in stage_d:
        by_depth[row["case"].split()[-1]].append(row)
    table = ["| depth | tiers | exact |", "|---|---|---|",
             f"| full | {_tiers(gyrb)} | "
             f"{sum(r['matches_truth'] == 'yes' for r in gyrb)} of {len(gyrb)} |"]
    for depth in sorted(by_depth, key=lambda label: -float(label.rstrip("x"))):
        rows = by_depth[depth]
        table.append(f"| ~{depth} | {_tiers(rows)} | "
                     f"{sum(r['matches_truth'] == 'yes' for r in rows)} of {len(rows)} |")
    out.append(("CALIB-10", "\n".join(table)))

    selection = {row["gene"]: row for row in _read_tsv(directory, "locus_selection.tsv")}
    by_locus = defaultdict(list)
    for row in stage_c:
        by_locus[_locus(row)].append(row)
    table = ["| locus | mean identity to the bait | bait length | reached `HIGH` | truth |",
             "|---|---|---|---|---|"]
    for locus, rows in sorted(by_locus.items()):
        chosen = selection.get(locus, {})
        truth = _truth(rows)
        table.append(
            f"| *{locus}* | {chosen.get('mean_identity_to_bait_pct', '?')}% | "
            f"{chosen.get('bait_length_bp', '?')} bp | "
            f"{sum(r['sequence_confidence'] == 'HIGH' for r in rows)} / {len(rows)} | "
            f"{truth} |")
    out.append(("CALIB-11", "\n".join(table)))
    out.append(("CALIB-12", f"{_tiers(rrna)}; consensus exact in "
                f"{sum(r['matches_truth'] == 'yes' for r in rrna)} of {len(rrna)}"))

    tables = (directory / "tier_calibration_tables.md").read_text()
    section = tables.split("## All four stages", 1)[-1].split("\n## ", 1)[0].strip()
    out.append(("CALIB-13", section))
    control = summary.get("stage_d_isolate_control", {})
    out.append(("CALIB-14", control.get("tier_agreement_with_primary",
                                        "control not run (--depth-isolate)")))

    for slot, value in out:
        separator = "\n" if "\n" in value else " "
        print(f"{slot}:{separator}{value}\n")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--list", action="store_true", help="inventory the slots")
    mode.add_argument("--values", action="store_true",
                      help="print the values that close the calibration slots")
    parser.add_argument("--calibration-dir", type=Path, default=CALIBRATION,
                        help="tier-calibration output directory (default: %(default)s)")
    args = parser.parse_args()
    if args.values:
        return values(args.calibration_dir)

    rows = slot_scan(scanned_files())
    if args.list:
        for entry in slot_summary(rows):
            print(f"{entry['slot']}\t{len(entry['where'])} occurrence(s)\t{entry['note']}")
            for where in entry["where"]:
                print(f"\t{where}")
        print(f"{len(slot_summary(rows))} distinct slot(s), {len(rows)} occurrence(s)")
        return 0
    if rows:
        first = "; ".join(f"{row['file']}:{row['line']} [[{row['slot']}]]"
                          for row in rows[:10])
        print(f"{len(rows)} unresolved slot occurrence(s) across "
              f"{len(slot_summary(rows))} slot(s); not releasable. First: {first}",
              file=sys.stderr)
        return 1
    print(f"No slot tokens in {len(scanned_files())} documentation files.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
