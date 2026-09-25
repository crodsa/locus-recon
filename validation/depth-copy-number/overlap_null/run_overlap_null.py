#!/usr/bin/env python3
"""Geometry-matched null for the additive dosage estimand.

The random-window control shows what a single-copy query returns when the
assembler represents it as one interval. It does not answer the question the
tcdB case raises: what does the same estimand return for a single-copy query
that the discovery step fragments into many mutually overlapping intervals,
which is the geometry that produces a dosage above one?

This script enriches for that geometry. It self-aligns the assembly once,
builds a per-base multiplicity track, and draws pseudo-queries centred on
repeat-rich positions, keeping those whose own discovery run returns at least
two intervals. Each retained query is single-copy in the sense that matters:
it is a stretch of the draft assembly, not a locus with independent evidence
for two copies. Any dosage above one that it returns is estimand behaviour,
not biology.

The retained queries are then summarised by a least-squares fit of dosage on
the overlapped fraction of the query, and the fit is evaluated, with a 95%
prediction interval, at the overlapped fraction and dosage of the case under
review, read from its evaluation record.

Usage:
    python validation/depth-copy-number/overlap_null/run_overlap_null.py \\
        --assembly LIBA6656_ST154.assembly.fasta \\
        --bam ERR467623_vs_LIBA6656_draft.sorted.bam \\
        --locus-contigs .11940_5_90.1,.11940_5_90.29,.11940_5_90.56,.11940_5_90.59,\\
.11940_5_90.63,.11940_5_90.76,.11940_5_90.77,.11940_5_90.8 \\
        --case-evaluation validation/liba6656-depth/evaluation.json \\
        --outdir validation/depth-copy-number/overlap_null
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from locus_recon.depth_ratio import _query_segments, query_copy_profile  # noqa: E402

MIN_HIT_BP = 300
MIN_HIT_IDENTITY = 90.0
WINDOW_STEP = 250
EDGE_BP = 100


def read_fasta(path):
    records, name = {}, None
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name is not None:
                records[name].append(line.strip())
    return {key: "".join(value) for key, value in records.items()}


def median(values):
    ordered = sorted(values)
    n = len(ordered)
    if not n:
        return 0.0
    middle = n // 2
    return ordered[middle] if n % 2 else (ordered[middle - 1] + ordered[middle]) / 2


def _betacf(a, b, x):
    """Continued fraction for the regularised incomplete beta function."""
    tiny, eps = 1e-300, 3e-16
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c if abs(1.0 + aa / c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c if abs(1.0 + aa / c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _incomplete_beta(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_quantile(probability, df):
    """Student t quantile by bisection on the exact CDF."""
    def cdf(t):
        tail = 0.5 * _incomplete_beta(df / 2.0, 0.5, df / (df + t * t))
        return 1.0 - tail if t >= 0 else tail

    low, high = -60.0, 60.0
    for _ in range(200):
        middle = (low + high) / 2
        if cdf(middle) < probability:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def fit_null(rows, case_overlap, case_dosage):
    """Least-squares fit of dosage on overlap, evaluated at the case."""
    x = [row["overlap_fraction"] for row in rows]
    y = [row["dosage"] for row in rows]
    n = len(x)
    if n < 3:
        raise RuntimeError("fewer than three multi-region pseudo-queries; no fit")
    mean_x, mean_y = sum(x) / n, sum(y) / n
    sxx = sum((value - mean_x) ** 2 for value in x)
    sxy = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    residuals = [b - (intercept + slope * a) for a, b in zip(x, y)]
    sse = sum(value * value for value in residuals)
    sst = sum((value - mean_y) ** 2 for value in y)
    residual_sd = math.sqrt(sse / (n - 2))
    prediction = intercept + slope * case_overlap
    prediction_se = residual_sd * math.sqrt(
        1 + 1 / n + (case_overlap - mean_x) ** 2 / sxx
    )
    t = t_quantile(0.975, n - 2)
    return {
        "fit_intercept": round(intercept, 4),
        "fit_slope": round(slope, 4),
        "fit_r2": round(1 - sse / sst, 4),
        "residual_sd": round(residual_sd, 4),
        "null_prediction_at_tcdB_overlap": round(prediction, 3),
        "null_prediction_interval_95": [
            round(prediction - t * prediction_se, 3),
            round(prediction + t * prediction_se, 3),
        ],
        "z_of_observed": round((case_dosage - prediction) / prediction_se, 3),
        "overlap_at_which_null_reaches_1.5": round((1.5 - intercept) / slope, 3),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--assembly", type=Path, required=True)
    parser.add_argument("--bam", type=Path, required=True,
                        help="reads mapped to --assembly, coordinate-sorted and indexed")
    parser.add_argument("--locus-contigs", required=True,
                        help="comma-separated contigs carrying the locus; excluded")
    parser.add_argument("--case-evaluation", type=Path, required=True,
                        help="evaluation.json of the case (dosage and overlap fraction)")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--query-length", type=int, default=7104)
    parser.add_argument("--seed", type=int, default=20260903)
    parser.add_argument("--target", type=int, default=400,
                        help="stop after this many multi-region pseudo-queries")
    parser.add_argument("--max-tried", type=int, default=3000)
    parser.add_argument("--min-regions", type=int, default=2)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    locus = {name for name in args.locus_contigs.split(",") if name}
    case = json.loads(args.case_evaluation.read_text())
    case_dosage = float(case["dosage_estimate"])
    case_overlap = float(case["gene_overlap_fraction"])
    qlen = args.query_length

    assembly = read_fasta(args.assembly)
    depth_output = subprocess.run(
        ["samtools", "depth", "-a", "-q", "20", "-Q", "20", str(args.bam)],
        capture_output=True, text=True, check=True,
    )
    depth = {}
    for line in depth_output.stdout.splitlines():
        contig, _position, value = line.split("\t")
        depth.setdefault(contig, []).append(int(value))
    backbone = [
        value for contig, values in depth.items()
        if contig not in locus and len(values) >= 5000
        for value in values[EDGE_BP:-EDGE_BP]
    ]
    backbone_median = median(backbone)

    rows, tried = [], 0
    with tempfile.TemporaryDirectory(prefix="overlap-null-") as temporary:
        database = f"{temporary}/assembly"
        subprocess.run(
            ["makeblastdb", "-in", str(args.assembly), "-dbtype", "nucl",
             "-out", database], check=True, capture_output=True,
        )
        # One self-alignment locates repeat-rich positions.
        self_hits = subprocess.run(
            ["blastn", "-query", str(args.assembly), "-db", database,
             "-evalue", "1e-10", "-num_threads", str(args.threads),
             "-outfmt", "6 qseqid sseqid pident length qstart qend"],
            capture_output=True, text=True, check=True,
        )
        multiplicity = {contig: [0] * len(seq) for contig, seq in assembly.items()}
        for line in self_hits.stdout.splitlines():
            fields = line.split("\t")
            if int(fields[3]) < MIN_HIT_BP or float(fields[2]) < MIN_HIT_IDENTITY:
                continue
            track = multiplicity.get(fields[0])
            if track is None:
                continue
            for index in range(int(fields[4]) - 1, min(int(fields[5]), len(track))):
                track[index] += 1

        windows = []
        for contig, sequence in assembly.items():
            if contig in locus or len(sequence) < qlen + 4 * EDGE_BP:
                continue
            track = multiplicity[contig]
            for start in range(EDGE_BP, len(sequence) - qlen - EDGE_BP, WINDOW_STEP):
                if sequence[start:start + qlen].upper().count("N") == 0:
                    rich = sum(1 for value in track[start:start + qlen] if value >= 2) / qlen
                    windows.append((contig, start, rich))
        random.Random(args.seed).shuffle(windows)
        print(f"repeat-rich candidate windows: {len(windows)}", flush=True)

        query = Path(temporary) / "query.fasta"
        for contig, start, _rich in windows:
            if len(rows) >= args.target or tried >= args.max_tried:
                break
            tried += 1
            query.write_text(f">query\n{assembly[contig][start:start + qlen]}\n")
            discovery = subprocess.run(
                ["blastn", "-query", str(query), "-db", database, "-evalue", "1e-10",
                 "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send"],
                capture_output=True, text=True, check=True,
            )
            hits = []
            for line in discovery.stdout.splitlines():
                fields = line.split("\t")
                if int(fields[3]) < MIN_HIT_BP or float(fields[2]) < MIN_HIT_IDENTITY:
                    continue
                sstart, send = int(fields[6]), int(fields[7])
                hits.append({"sseqid": fields[1], "sstart": min(sstart, send),
                             "send": max(sstart, send),
                             "qstart": int(fields[4]), "qend": int(fields[5])})
            merged, by_contig = [], {}
            for hit in hits:
                by_contig.setdefault(hit["sseqid"], []).append(hit)
            for group in by_contig.values():
                group.sort(key=lambda item: item["sstart"])
                current = dict(group[0])
                for hit in group[1:]:
                    if hit["sstart"] <= current["send"] + 20:
                        current["send"] = max(current["send"], hit["send"])
                        current["qstart"] = min(current["qstart"], hit["qstart"])
                        current["qend"] = max(current["qend"], hit["qend"])
                    else:
                        merged.append(current)
                        current = dict(hit)
                merged.append(current)
            if len(merged) < args.min_regions:
                continue
            spans, ratios = [], []
            for hit in merged:
                values = depth.get(hit["sseqid"], [])[hit["sstart"] - 1:hit["send"]]
                if not values:
                    continue
                spans.append((min(hit["qstart"], hit["qend"]),
                              max(hit["qstart"], hit["qend"])))
                ratios.append(median(values) / backbone_median)
            if len(spans) < args.min_regions:
                continue
            dosage, _lower, _upper, covered, max_mult = query_copy_profile(
                _query_segments(spans, qlen), ratios, qlen)
            total = sum(end - begin + 1 for begin, end in spans)
            rows.append({
                "contig": contig, "start": start, "n_regions": len(spans),
                "coverage_fraction": covered / qlen,
                "overlap_fraction": (total - covered) / qlen,
                "max_multiplicity": max_mult,
                "dosage": dosage,
            })

    args.outdir.mkdir(parents=True, exist_ok=True)
    # The fit uses full precision; the table is rounded for reading.
    digits = {"coverage_fraction": 4, "overlap_fraction": 4,
              "max_multiplicity": 3, "dosage": 3}
    with (args.outdir / "null_dose_overlap_control.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {key: round(value, digits[key]) if key in digits else value
             for key, value in row.items()}
            for row in rows
        )
    overlaps = [round(row["overlap_fraction"], 4) for row in rows]
    summary = {
        "n_windows_tried": tried,
        "n_multi_region": len(rows),
        "overlap_fraction_range": [min(overlaps), max(overlaps)],
        **fit_null(rows, case_overlap, case_dosage),
        "tcdB_overlap_fraction": round(case_overlap, 4),
        "tcdB_dosage": case_dosage,
        "note": (
            "Single-copy pseudo-queries drawn from the same draft assembly and "
            "run through the same discovery and dosage steps. Dosage rises with "
            "the fraction of the query resolved by more than one interval. The "
            "tcdB overlap fraction lies beyond the observed control range, so "
            "the null at that geometry is extrapolated."
        ),
    }
    (args.outdir / "null_dose_overlap_control.json").write_text(
        json.dumps(summary, indent=1) + "\n"
    )
    print(f"written: {len(rows)} multi-region pseudo-queries from {tried} windows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
