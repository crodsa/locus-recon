#!/usr/bin/env python3
# WITHDRAWN CRITERION (v1.1). The statistics this script computes — the
# dispersion index H and the window P90 — were removed from
# `locus_recon.depth_ratio` in v1.1 because their thresholds were calibrated
# against single-copy backbone only and their sensitivity to real partial
# amplifications was never measured. The script is kept unchanged as the
# record of how the published thresholds were derived, and is referenced by
# `../PRESPECIFIED_CRITERIA.md`. It no longer describes shipped behaviour and
# is not executed by any test.
"""Null distributions of the two windowed statistics over single-copy sequence.

Both thresholds used by `locus_recon.depth_ratio` come from this script:
the dispersion index H = (P75 - P25) / P50 and the window P90, computed over
non-overlapping 4,000 bp segments (20 windows of 200 bp) of a backbone that
contains no target locus.

The segments must be single-copy. Contigs carrying the locus under study, short
contigs and any replicon with its own copy number are excluded, mirroring the
backbone definition used by the module itself.

    python null_distribution.py --bam reads_vs_draft.sorted.bam \
        --backbone-depth 62.0 --exclude contig_A --exclude contig_B

Depth filters are the module's: base quality >= 20, mapping quality >= 20.
"""
from __future__ import annotations

import argparse
import subprocess
import sys

WINDOW = 200
N_WINDOWS = 20          # 20 x 200 bp = 4 kb, the calibration scale
MIN_CONTIG = 5000


def _median(v: list[int]) -> float:
    if not v:
        return 0.0
    s = sorted(v)
    n = len(s)
    return float(s[n // 2]) if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


def _quantile(v: list[float], p: float) -> float:
    if not v:
        return 0.0
    s = sorted(v)
    i = (len(s) - 1) * p
    lo = int(i)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (i - lo)


def segments(bam: str, exclude: set[str], samtools: str) -> list[list[float]]:
    """Window-ratio vectors, one per non-overlapping 4 kb single-copy segment."""
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
    depth = subprocess.run([samtools, "depth", "-a", "-q", "20", "-Q", "20",
                            "-b", "/dev/stdin", bam],
                           input=bed, capture_output=True, text=True,
                           check=True).stdout

    out: list[list[float]] = []
    cur: str | None = None
    buf: list[int] = []

    def flush(b: list[int], backbone: float) -> None:
        step = WINDOW * N_WINDOWS
        for s in range(0, len(b) - step + 1, step):
            wr = [_median(b[s + i * WINDOW:s + (i + 1) * WINDOW]) / backbone
                  for i in range(N_WINDOWS)]
            wr = [x for x in wr if x > 0]
            if len(wr) == N_WINDOWS:
                out.append(wr)

    rows: list[tuple[str, int]] = []
    for line in depth.split("\n"):
        if not line:
            continue
        c, _, d = line.split("\t")
        rows.append((c, int(d)))

    backbone = _median([d for _, d in rows])
    if backbone <= 0:
        sys.exit("backbone median depth is zero")

    for c, d in rows:
        if c != cur:
            if cur is not None:
                flush(buf, backbone)
            cur, buf = c, []
        buf.append(d)
    if cur is not None:
        flush(buf, backbone)
    print(f"backbone median depth: {backbone:.1f}", file=sys.stderr)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True, help="reads mapped to the draft assembly")
    ap.add_argument("--exclude", action="append", default=[],
                    help="contig to keep out of the backbone; repeatable")
    ap.add_argument("--samtools", default="samtools")
    args = ap.parse_args()

    segs = segments(args.bam, set(args.exclude), args.samtools)
    if not segs:
        sys.exit("no 4 kb segment was fully covered")

    h = [(_quantile(s, 0.75) - _quantile(s, 0.25)) / _quantile(s, 0.50)
         for s in segs if _quantile(s, 0.50) > 0]
    p90 = [_quantile(s, 0.90) for s in segs]

    for name, v in (("H", h), ("window P90", p90)):
        print(f"{name}: {len(v)} single-copy 4 kb segments")
        for p in (0.50, 0.75, 0.90, 0.95, 0.99):
            print(f"  P{int(p * 100):<2} = {_quantile(v, p):.3f}")
        print(f"  max = {max(v):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
