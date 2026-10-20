#!/usr/bin/env python3
"""Falsifiable predictions (a), (b), (e) and (f) for the length-weighted estimator.

The four predictions were written down before any of them was evaluated, in
`PRESPECIFIED_CRITERIA.md`, and are reproduced here through the released
function itself, so that the number in the manuscript and the number a reader
obtains come from the same code path:

    from locus_recon.depth_ratio import length_weighted_copies

Prediction (a) is read off the deposited *aphA1* output. The remaining three
geometries are cut out of single-copy backbone sequence, which fixes the truth
without needing an organism that carries a known duplication.

    (a)  one fragment spanning the whole gene              truth = the ratio
    (b)  one fragment spanning a 4,000 bp segment          truth = 1 copy
    (e)  one 7,104 bp stretch cut into n = 2..9 fragments  truth = 1 copy
    (f)  two 7,104 bp stretches, each cut into quarters    truth = 2 copies

Prediction (a) requires a deviation from the ratio of at most 0.02 on the five
prespecified single-fragment cases, (b) a median inside [0.90, 1.10] with
P99 <= 1.65, (e) the same median bound at every n together with |slope| < 0.02
copies per fragment, and (f) a median inside [1.80, 2.20]. Nothing here is
fitted and no bound is read off the output.

Prediction (f) states two copies, so each copy is passed as its own tiling of
the gene: the estimator is a statement about one gene length, and a span running
past `gene_length` is clipped by design.

    python estimator_predictions.py --bam full_map.sorted.bam \
        --backbone-depth 62.0 --exclude contig_A --exclude contig_B \
        --external-validation ../external_validation.json

Depth filters are the module's: base quality >= 20, mapping quality >= 20.
"""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys

from locus_recon.depth_ratio import _median, length_weighted_copies

PRESPECIFIED_SINGLE_FRAGMENT = ("MRSN3361", "MRSN3363", "MRSN56",
                                "MRSN57", "MRSN58")
SEGMENT = 4000          # calibration scale of prediction (b)
GENE = 7104            # one copy of tcdB in LIBA-6656
MIN_CONTIG = 5000


def _quantile(v: list[float], p: float) -> float:
    if not v:
        return 0.0
    s = sorted(v)
    i = (len(s) - 1) * p
    lo = int(i)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (i - lo)


def backbone_depth(bam: str, exclude: set[str], samtools: str) -> dict[str, list[int]]:
    """Per-base depth of every single-copy backbone contig, contig by contig."""
    lens: dict[str, int] = {}
    idx = subprocess.run([samtools, "idxstats", bam],
                         capture_output=True, text=True, check=True).stdout
    for line in idx.strip().split("\n"):
        f = line.split("\t")
        if f[0] != "*":
            lens[f[0]] = int(f[1])
    keep = [c for c, L in lens.items() if c not in exclude and L >= MIN_CONTIG]
    if not keep:
        sys.exit("no backbone contig survived the exclusions")

    bed = "\n".join(f"{c}\t0\t{lens[c]}" for c in keep) + "\n"
    out = subprocess.run([samtools, "depth", "-a", "-q", "20", "-Q", "20",
                          "-b", "/dev/stdin", bam],
                         input=bed, capture_output=True, text=True,
                         check=True).stdout

    depth: dict[str, list[int]] = {}
    cur: str | None = None
    buf: list[int] = []
    for line in out.split("\n"):
        if not line:
            continue
        c, _, d = line.split("\t")
        if c != cur:
            if cur is not None:
                depth[cur] = buf
            cur, buf = c, []
        buf.append(int(d))
    if cur is not None:
        depth[cur] = buf
    return depth


def tile(length: int, n: int) -> list[tuple[int, int]]:
    """`n` consecutive fragments covering 1..length, remainder on the last one."""
    spans: list[tuple[int, int]] = []
    step = length // n
    pos = 1
    for k in range(n):
        end = length if k == n - 1 else pos + step - 1
        spans.append((pos, end))
        pos = end + 1
    return spans


def prediction_a(path: str) -> dict:
    """One fragment spanning the whole gene, truth the ratio itself."""
    with open(path) as fh:
        runs = json.load(fh)
    rows = []
    for rec in runs:
        if rec["strain"] not in PRESPECIFIED_SINGLE_FRAGMENT:
            continue
        g = rec["gene"]
        spans = [tuple(int(x) for x in s.split("-"))
                 for s in g["region_gene_spans"].split(";")]
        ratios = [float(x) for x in g["region_ratios"].split(";")]
        est = length_weighted_copies(spans, ratios, g["gene_length"])
        rows.append([rec["strain"], g["ratio"], est, abs(est - g["ratio"])])
    return {"cases": rows, "max_deviation": max(r[3] for r in rows),
            "pass": len(rows) == len(PRESPECIFIED_SINGLE_FRAGMENT)
                    and all(r[3] <= 0.02 for r in rows)}


def prediction_b(depth: dict[str, list[int]], backbone: float) -> dict:
    """One fragment, one 4 kb segment, truth one copy."""
    est: list[float] = []
    for buf in depth.values():
        for s in range(0, len(buf) - SEGMENT + 1, SEGMENT):
            r = _median(buf[s:s + SEGMENT]) / backbone
            if r > 0:
                est.append(length_weighted_copies([(1, SEGMENT)], [r], SEGMENT))
    med, p99 = _quantile(est, 0.50), _quantile(est, 0.99)
    return {"segments": len(est), "median": med, "p99": p99, "max": max(est),
            "pass": 0.90 <= med <= 1.10 and p99 <= 1.65}


def stretches(depth: dict[str, list[int]]) -> list[list[int]]:
    """Non-overlapping single-copy stretches one gene long."""
    out: list[list[int]] = []
    for buf in depth.values():
        for s in range(0, len(buf) - GENE + 1, GENE):
            out.append(buf[s:s + GENE])
    return out


def prediction_e(blocks: list[list[int]], backbone: float) -> dict:
    """One gene cut into n fragments, truth one copy, for n = 2..9."""
    rows: list[list[float]] = []
    for n in range(2, 10):
        spans = tile(GENE, n)
        est = [length_weighted_copies(
            spans, [_median(blk[a - 1:b]) / backbone for a, b in spans], GENE)
            for blk in blocks]
        rows.append([n, statistics.median(est)])
    xs = [r[0] for r in rows]
    ys = [r[1] for r in rows]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    slope = (sum((x - mx) * (y - my) for x, y in rows)
             / sum((x - mx) ** 2 for x in xs))
    return {"stretches": len(blocks), "median_by_fragment_count": rows,
            "slope_copies_per_fragment": slope,
            "pass": all(0.90 <= m <= 1.10 for _, m in rows) and abs(slope) < 0.02}


def prediction_f(blocks: list[list[int]], backbone: float) -> dict:
    """Two genes, each cut into quarters, truth two copies."""
    spans = tile(GENE, 4)
    est: list[float] = []
    for i in range(0, len(blocks) - 1, 2):
        geometry: list[tuple[int, int]] = []
        ratios: list[float] = []
        for blk in (blocks[i], blocks[i + 1]):
            geometry.extend(spans)
            ratios.extend(_median(blk[a - 1:b]) / backbone for a, b in spans)
        est.append(length_weighted_copies(geometry, ratios, GENE))
    med = statistics.median(est)
    return {"geometries": len(est), "median": med,
            "p25": _quantile(est, 0.25), "p75": _quantile(est, 0.75),
            "min": min(est), "max": max(est),
            "pass": 1.80 <= med <= 2.20}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True, help="reads mapped to the draft assembly")
    ap.add_argument("--backbone-depth", type=float, required=True,
                    help="median depth of the single-copy backbone")
    ap.add_argument("--exclude", action="append", default=[],
                    help="contig to keep out of the backbone; repeatable")
    ap.add_argument("--external-validation",
                    help="deposited aphA1 module output, for prediction (a)")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--json", help="write the results to this file")
    args = ap.parse_args()

    depth = backbone_depth(args.bam, set(args.exclude), args.samtools)
    print(f"backbone contigs {len(depth)}, "
          f"{sum(len(v) for v in depth.values()):,} bp", file=sys.stderr)

    blocks = stretches(depth)
    res = {"b": prediction_b(depth, args.backbone_depth),
           "e": prediction_e(blocks, args.backbone_depth),
           "f": prediction_f(blocks, args.backbone_depth)}
    keys = "bef"
    if args.external_validation:
        res["a"] = prediction_a(args.external_validation)
        keys = "abef"
        a = res["a"]
        print(f"(a) {len(a['cases'])} single-fragment cases")
        for strain, ratio, est, dev in a["cases"]:
            print(f"    {strain:<10} ratio {ratio:>7.3f}  estimator {est:>7.3f}"
                  f"  deviation {dev:.3f}")
        print(f"    largest deviation {a['max_deviation']:.3f}"
              f"   {'PASS' if a['pass'] else 'FALSIFIED'}")
    res["all_pass"] = all(res[k]["pass"] for k in keys)

    b, e, f = res["b"], res["e"], res["f"]
    print(f"(b) {b['segments']} segments of {SEGMENT} bp, one fragment each\n"
          f"    median {b['median']:.4f}  P99 {b['p99']:.4f}  max {b['max']:.4f}"
          f"   {'PASS' if b['pass'] else 'FALSIFIED'}")
    print(f"(e) {e['stretches']} stretches of {GENE} bp")
    for n, m in e["median_by_fragment_count"]:
        print(f"    n = {n}  median {m:.4f}")
    print(f"    slope {e['slope_copies_per_fragment']:+.5f} copies per fragment"
          f"   {'PASS' if e['pass'] else 'FALSIFIED'}")
    print(f"(f) {f['geometries']} two-copy geometries\n"
          f"    median {f['median']:.4f}  IQR {f['p25']:.2f}-{f['p75']:.2f}"
          f"   {'PASS' if f['pass'] else 'FALSIFIED'}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=2)
    return 0 if res["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
