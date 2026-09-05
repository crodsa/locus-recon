#!/usr/bin/env python3
# WITHDRAWN CRITERION (v1.1). The statistics this script computes — the
# dispersion index H and the window P90 — were removed from
# `locus_recon.depth_ratio` in v1.1 because their thresholds were calibrated
# against single-copy backbone only and their sensitivity to real partial
# amplifications was never measured. The script is kept unchanged as the
# record of how the published thresholds were derived, and is referenced by
# `../PRESPECIFIED_CRITERIA.md`. It no longer describes shipped behaviour and
# is not executed by any test.
"""How finely the length-weighted estimator resolves adjacent low copy numbers.

Post hoc. Predictions (b), (e) and (f) were written down in advance and are
evaluated by `estimator_predictions.py`; this script was written after those
results were known and answers a different question, so nothing here is
presented as a test that could have falsified the estimator.

The question is what the estimator can and cannot distinguish near two copies.
Geometries with a known number of copies are cut out of single-copy backbone
sequence and passed to the released function:

    1.0 copies   one gene length, four fragments
    2.0 copies   two gene lengths, four fragments each
    2.5 copies   two gene lengths plus one half gene, four fragments each
                 and two more over the half

Each geometry is then read as a draw from the sampling distribution of the
estimator at that truth, and the spacing between the three distributions is the
resolution the estimator actually achieves on this genome. An observed estimate
can be placed inside them, which turns a point estimate into an empirical
two-sided p value for each candidate truth.

The dispersion measured here is that of the backbone, not of the locus, so the
intervals describe depth variation of single-copy sequence in this library and
carry no information about locus-specific mapping behaviour.

    python copy_number_resolution.py --bam full_map.sorted.bam \
        --backbone-depth 62.0 --observed 2.258 --exclude contig_A

Depth filters are the module's: base quality >= 20, mapping quality >= 20.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

from locus_recon.depth_ratio import _median, length_weighted_copies

sys.path.insert(0, str(Path(__file__).resolve().parent))
from estimator_predictions import GENE, _quantile, backbone_depth, tile

WINDOW = 200            # window size of the module's dispersion profile


def labelled_stretches(depth: dict[str, list[int]]) -> list[tuple[str, int, list[int]]]:
    """Non-overlapping single-copy stretches one gene long, with their origin."""
    out: list[tuple[str, int, list[int]]] = []
    for contig, buf in depth.items():
        for s in range(0, len(buf) - GENE + 1, GENE):
            out.append((contig, s + 1, buf[s:s + GENE]))
    return out


def window_profile(block: list[int], backbone: float) -> dict:
    """The module's own two window statistics over one stretch."""
    w = [_median(block[i:i + WINDOW]) / backbone
         for i in range(0, len(block) - WINDOW + 1, WINDOW)]
    p25, p50, p75 = (_quantile(w, p) for p in (0.25, 0.50, 0.75))
    return {"windows": len(w), "H": (p75 - p25) / p50 if p50 else 0.0,
            "p90": _quantile(w, 0.90)}


def geometries(blocks: list[tuple[str, int, list[int]]], backbone: float,
               copies: float) -> list[tuple[float, int]]:
    """Estimates over disjoint constructions of `copies` copies of the gene."""
    whole, half = int(copies), copies - int(copies)
    per_construction = whole + (1 if half else 0)
    full_spans = tile(GENE, 4)
    half_spans = tile(GENE // 2, 2) if half else []
    out: list[tuple[float, int]] = []
    for i in range(0, len(blocks) - per_construction + 1, per_construction):
        spans: list[tuple[int, int]] = []
        ratios: list[float] = []
        for k in range(whole):
            blk = blocks[i + k][2]
            spans.extend(full_spans)
            ratios.extend(_median(blk[a - 1:b]) / backbone for a, b in full_spans)
        if half:
            blk = blocks[i + whole][2]
            spans.extend(half_spans)
            ratios.extend(_median(blk[a - 1:b]) / backbone for a, b in half_spans)
        out.append((length_weighted_copies(spans, ratios, GENE), i))
    return out


def describe(est: list[float]) -> dict:
    return {"n": len(est), "median": statistics.median(est),
            "p2.5": _quantile(est, 0.025), "p25": _quantile(est, 0.25),
            "p75": _quantile(est, 0.75), "p97.5": _quantile(est, 0.975),
            "min": min(est), "max": max(est)}


def two_sided_p(est: list[float], observed: float) -> float:
    n = len(est)
    below = sum(1 for x in est if x <= observed) / n
    above = sum(1 for x in est if x >= observed) / n
    return min(1.0, 2 * min(below, above))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--backbone-depth", type=float, required=True)
    ap.add_argument("--observed", type=float,
                    help="an estimate to place inside the distributions")
    ap.add_argument("--exclude", action="append", default=[])
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--json")
    args = ap.parse_args()

    depth = backbone_depth(args.bam, set(args.exclude), args.samtools)
    blocks = labelled_stretches(depth)
    print(f"single-copy stretches of {GENE} bp: {len(blocks)}", file=sys.stderr)

    res: dict = {"truths": {}}
    draws: dict[float, list[float]] = {}
    drawn: dict[float, list[tuple[float, int]]] = {}
    for copies in (1.0, 2.0, 2.5):
        pairs = geometries(blocks, args.backbone_depth, copies)
        drawn[copies] = pairs
        est = [e for e, _ in pairs]
        draws[copies] = est
        d = describe(est)
        if args.observed is not None:
            d["two_sided_p_for_observed"] = two_sided_p(est, args.observed)
        res["truths"][f"{copies:.1f}"] = d
        line = (f"truth {copies:.1f}  n = {d['n']:>3}  median {d['median']:.3f}"
                f"  central 95% {d['p2.5']:.2f}-{d['p97.5']:.2f}"
                f"  range {d['min']:.2f}-{d['max']:.2f}")
        if args.observed is not None:
            line += f"  p = {d['two_sided_p_for_observed']:.3f}"
        print(line)

    res["one_versus_two_gap"] = {
        "max_at_one_copy": max(draws[1.0]), "min_at_two_copies": min(draws[2.0]),
        "disjoint": max(draws[1.0]) < min(draws[2.0])}
    cut = 2.25
    res["two_versus_two_and_a_half"] = {
        "midpoint": cut,
        "at_two_copies_above_midpoint":
            sum(1 for x in draws[2.0] if x >= cut) / len(draws[2.0]),
        "at_two_and_a_half_below_midpoint":
            sum(1 for x in draws[2.5] if x < cut) / len(draws[2.5])}
    if args.observed is not None:
        res["observed"] = args.observed

    worst, idx = max(drawn[1.0])
    contig, start, block = blocks[idx]
    res["largest_at_one_copy"] = {
        "estimate": worst, "contig": contig, "start": start,
        "fragment_ratios": [_median(block[a - 1:b]) / args.backbone_depth
                            for a, b in tile(GENE, 4)],
        "window_profile": window_profile(block, args.backbone_depth)}
    w = res["largest_at_one_copy"]
    print(f"largest at one copy: {worst:.3f} on {contig} at {start}, "
          f"fragment ratios "
          f"{', '.join(f'{r:.2f}' for r in w['fragment_ratios'])}; "
          f"window H {w['window_profile']['H']:.3f} and P90 "
          f"{w['window_profile']['p90']:.3f} over "
          f"{w['window_profile']['windows']} windows")

    g = res["one_versus_two_gap"]
    o = res["two_versus_two_and_a_half"]
    print(f"one against two: {g['max_at_one_copy']:.3f} at most against "
          f"{g['min_at_two_copies']:.3f} at least, "
          f"{'disjoint' if g['disjoint'] else 'overlapping'}")
    print(f"two against 2.5 at a midpoint of {cut}: "
          f"{o['at_two_copies_above_midpoint']:.1%} of two-copy geometries land "
          f"above it and {o['at_two_and_a_half_below_midpoint']:.1%} of "
          f"2.5-copy geometries land below it")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
