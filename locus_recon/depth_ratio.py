"""
depth_ratio.py -- normalised locus depth ratio for Locus-Recon.

Why this exists
---------------
Locus-Recon adjudicates copy structure from *allelic* evidence: bidirectionally
supported mixed sites, frame anomalies, graph branches. That evidence is silent
when the copies are identical. A tandem amplification of an antimicrobial
resistance gene produces no mixed site at all, yet it is emphatically not a
single-copy locus, and the phenotype depends on how many copies are present.

The only signal left in short reads is depth. When a repeat collapses in the
draft assembly, every copy's reads pile onto the single collapsed
representative, so

    copies ~= depth(locus) / depth(single-copy backbone)

This module computes that ratio with the same quality filters as the rest of
the evidence table, with a block bootstrap confidence interval and --
importantly -- with an explicit ambiguity diagnostic that tells you when the
ratio cannot be trusted.

Design notes that matter
------------------------
1. **No new alignment is needed.** The standard workflow already maps the
   paired reads back to the draft assembly during read recruitment. That BAM is
   exactly what the backbone median needs. Pass it in; do not realign.

2. **The denominator must use the same filters as the numerator.** A backbone
   computed without the Q20/MQ20 filters inflates itself relative to a filtered
   locus and biases the ratio downward.

3. **MAPQ filtering cuts both ways.** If the locus collapsed to ONE copy in the
   draft, reads from all physical copies map there unambiguously and MAPQ is
   high: the ratio is meaningful. If the assembler kept two or more copies,
   reads distribute ambiguously, MAPQ collapses to 0, and the filtered depth
   *understates* the truth. `ambiguity_index` (filtered/unfiltered depth at the
   locus) detects this. Below `ambiguity_warn` the ratio is reported but marked
   unreliable rather than silently wrong.

4. **Fragmented loci require explicit coverage accounting.** The ratio is per
   collapsed representative. When discovery returns several non-overlapping
   assembly intervals that jointly cover the whole query, `dosage_estimate`
   weights each region's ratio by the query bases that region covers,

       copies = sum_r (aligned_gene_bp_r * ratio_r) / gene_length

   which reduces exactly to `ratio` when a single fragment spans the gene. Where
   query spans OVERLAP the sum is not a copy count: a single-copy query resolved
   into mutually overlapping intervals returns a value above one from the
   geometry alone. The value is still reported, with `dosage_status` set to
   `REVIEW_OVERLAPPING_SPANS` and `geometry_expected_dosage` giving what the same
   spans return at unit ratios. An empirical control on 87 single-copy
   pseudo-queries from a real assembly shows the returned dosage rising with the
   overlapped fraction (slope 1.54, R^2 0.64), so the geometric figure is a lower
   bound on the inflation. Incomplete query coverage sets
   `ESTIMATED_PARTIAL_COVERAGE`. The raw aligned, union-covered and overlapping
   lengths remain available for audit. See `validation/depth-copy-number/` and
   `validation/depth-copy-number/overlap_null/`.

5. **Depth is autocorrelated.** Bootstrapping individual bases would give an
   absurdly tight interval. Resampling is done over contiguous blocks.

6. **Non-integer ratios are informative, not noise.** An amplification carried
   by a subpopulation gives a fractional mean copy number. Report it as mean
   copies per cell; do not round it into a lie.

What this module does not answer
--------------------------------
The ratio measures **collapse**: the pile-up of reads onto a representative the
assembler failed to separate. It is not, on its own, a total copy number. Where the
assembler resolved the copies onto different contigs there is no pile-up, and a
ratio near 1 is the correct answer to the question the ratio asks while saying
nothing about how many copies exist. `dosage_estimate` combines query coverage
with collapsed depth only for calibration-eligible geometries; otherwise its
status and flags state why no total dosage is reported.

It also says nothing about **which** copy carries which allele. Assigning a
short read to one physical copy of a repeat requires the read to overlap a
position at which that copy differs from every other copy; whether such
positions exist is a property of the locus, not of the estimator, and is
computable from a closed reference before any sample is sequenced.
"""
from __future__ import annotations

import random
import re
import statistics
import subprocess
from dataclasses import dataclass, asdict, field
from pathlib import Path

try:  # Optional acceleration; the algorithm and seed are unchanged without it.
    import numpy as _np
except ImportError:  # pragma: no cover - exercised in the minimal runtime
    _np = None

__all__ = [
    "DepthRatioResult",
    "estimate_locus_copy_number",
    "depth_array",
    "block_bootstrap_ci",
    "query_copy_profile",
]

# Excluded SAM flags: unmapped, secondary, QC fail, optical/PCR duplicate,
# supplementary. Matches what the remap evidence table already discards.
DEFAULT_EXCL_FLAGS = "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY"

# Lower bound the aggregate interval must clear before single copy is rejected.
#
# The bound is asserted rather than fitted to these data, and it is asserted
# deliberately rather than by inertia: applying this rule to
# 300 single-copy 7,104 bp segments drawn from the backbones of two real
# alignments (150 from ERR467623 vs the LIBA-6656 draft, 150 from SRR15734215
# vs MRSN56) rejected single copy 0 times. The 95th percentile of the interval
# lower bound over those segments was 1.145 and 1.094 respectively; the largest
# value observed was 1.365. The margin between that maximum and 1.5 is the
# constant's whole justification, and it is not decorative -- see below.
#
# The margin is needed because the interval is narrower than the segment-to-
# segment spread of true single-copy depth. Over the same 300 segments the
# nominal 95% interval covered the true value of 1.0 in 35% and 41% of cases,
# and excluded 1.0 from above in 30% and 24%. The block bootstrap resamples
# 500 bp blocks within one locus, so it measures how precisely that locus's own
# median ratio is pinned down; it does not measure how much 7 kb regions differ
# from one another for reasons -- coverage waves, GC, mappability -- that are
# fixed under the frozen alignment and therefore invisible to resampling. A
# rule of `ci_low > 1.0` would inherit that mismatch and reject single copy on
# roughly a quarter of single-copy loci. 1.5 absorbs it.
#
# Consequences for anyone editing this number: lowering it is not a
# sensitivity improvement, it is an unquantified rise in the false rejection
# rate, and it must not be done on the argument that the interval is a 95%
# interval. Calibration per sample is available (call_threshold_calibration)
# and can only raise the effective threshold, never lower it.
CALL_THRESHOLD_DEFAULT = 1.5
CALL_THRESHOLD_NULL_SEGMENTS = 300
CALL_THRESHOLD_NULL_REJECTIONS = 0


@dataclass
class DepthRatioResult:
    """One locus, one sample."""

    locus_regions: list[str]
    locus_median_depth: float
    locus_median_depth_unfiltered: float
    backbone_median_depth: float
    backbone_mad: float
    ratio: float
    ratio_ci_low: float
    ratio_ci_high: float
    assembly_copies: int
    # Copies, weighted by the gene bases each assembly fragment resolves. Over
    # 512 single-copy stretches cut into n = 2..9 pieces the estimator returns
    # 0.984-0.989, and over 256 two-copy geometries a median of 1.985
    # (P25/P75 1.81/2.17). See `validation/depth-copy-number/`.
    # None when the locus spans several regions and no gene coordinates were
    # supplied: there is then no correct way to compute it, and None is
    # preferable to a wrong number.
    copies_estimate: float | None
    ambiguity_index: float
    gc_locus: float | None
    single_copy_rejected: bool
    reliable: bool
    notes: list[str] = field(default_factory=list)
    # `assembly_copies` and `copies_estimate` are retained as compatibility
    # aliases. Their explicit replacements name the estimands rather than
    # implying that an assembly interval is a biological copy.
    assembly_region_count: int = 0
    dosage_estimate: float | None = None
    dosage_status: str = "NOT_ESTIMATED"
    dosage_flags: list[str] = field(default_factory=list)
    depth_call: str = "DEPTH_NOT_EVALUATED"
    bootstrap_block_size: int = 0
    bootstrap_locus_blocks: int = 0
    bootstrap_backbone_blocks_total: int = 0
    bootstrap_backbone_blocks_used: int = 0
    bootstrap_resamples: int = 0
    bootstrap_seed: int = 0
    # --- provenance of the aggregate call threshold ---------------------------- #
    # The default is asserted, not fitted (see CALL_THRESHOLD_DEFAULT). These
    # fields record which value a run used and where it came from, so a reader
    # is not left to infer that it was derived from these data.
    call_threshold_used: float = CALL_THRESHOLD_DEFAULT
    call_threshold_source: str = "DEFAULT_ASSERTED"
    # --- length-weighted copy estimate --------------------------------------- #
    gene_length: int = 0
    gene_covered_bp: int = 0
    gene_total_aligned_bp: int = 0
    gene_union_covered_bp: int = 0
    gene_overlap_bp: int = 0
    gene_overlap_fraction: float = 0.0
    region_gene_spans: list[tuple[int, int]] = field(default_factory=list)
    region_ratios: list[float] = field(default_factory=list)
    # --- per-query-base dosage --------------------------------------------- #
    # `dosage_estimate` is the mean multiplicity over COVERED query bases.
    # The bracket is deterministic given the geometry: lower assigns zero to
    # uncovered bases, upper assigns the largest observed multiplicity. The CI
    # is stochastic, from resampling every region and the backbone in blocks.
    dosage_lower: float | None = None
    dosage_upper: float | None = None
    dosage_ci_low: float | None = None
    dosage_ci_high: float | None = None
    query_coverage_fraction: float = 0.0
    query_multiplicity_max: float = 0.0
    # The dosage the SAME span geometry returns when every region is single
    # copy: the mean number of spans covering each covered query base. It is
    # exact and parameter-free, and it is a LOWER bound on geometric inflation,
    # because sequence represented several times in an assembly also attracts
    # multi-mapped reads that raise the per-region ratios themselves. See
    # `validation/depth-copy-number/overlap_null/`.
    geometry_expected_dosage: float | None = None

    def as_row(self) -> dict:
        d = asdict(self)
        d["locus_regions"] = ";".join(self.locus_regions)
        d["notes"] = "; ".join(self.notes)
        d["dosage_flags"] = ";".join(self.dosage_flags)
        d["region_gene_spans"] = ";".join(f"{a}-{b}" for a, b in self.region_gene_spans)
        d["region_ratios"] = ";".join(f"{x:.3f}" for x in self.region_ratios)
        return d

    def summary_line(self) -> str:
        rel = "" if self.reliable else "  [UNRELIABLE: see notes]"
        n_regions = self.assembly_region_count or self.assembly_copies
        frag = "region" if n_regions == 1 else "regions"
        head = (f"depth ratio {self.ratio:.2f} "
                f"(95% CI {self.ratio_ci_low:.2f}-{self.ratio_ci_high:.2f}), "
                f"{n_regions} assembly {frag}")
        if self.dosage_estimate is not None:
            head += f", dosage {self.dosage_estimate:.2f} ({self.dosage_status})"
        return f"{head} -> {self.depth_call}{rel}"


# --------------------------------------------------------------------------- #
# depth extraction
# --------------------------------------------------------------------------- #
def depth_array(
    bam: str | Path,
    region: str | None = None,
    *,
    min_bq: int = 20,
    min_mq: int = 20,
    excl_flags: str = DEFAULT_EXCL_FLAGS,
    samtools: str = "samtools",
) -> list[int]:
    """Per-base depth as a list. `region` is contig or contig:start-end."""
    cmd = [samtools, "depth", "-a", "-q", str(min_bq), "-Q", str(min_mq),
           "--excl-flags", excl_flags]
    if region:
        cmd += ["-r", region]
    cmd.append(str(bam))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools depth failed:\n{proc.stderr[:1000]}")
    out = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            out.append(int(parts[2]))
    return out


def _parse_region(region: str) -> tuple[str, int | None, int | None]:
    m = re.match(r"^([^:]+)(?::(\d+)-(\d+))?$", region)
    if not m:
        raise ValueError(f"malformed region: {region}")
    contig, start, end = m.group(1), m.group(2), m.group(3)
    return contig, int(start) if start else None, int(end) if end else None


def _contig_lengths(bam: str | Path, samtools: str = "samtools") -> dict[str, int]:
    proc = subprocess.run([samtools, "idxstats", str(bam)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools idxstats failed:\n{proc.stderr[:1000]}")
    lengths = {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 2 and f[0] != "*":
            lengths[f[0]] = int(f[1])
    return lengths


def _region_length(region: str) -> int:
    """Return the inclusive length of a fully specified region."""
    _contig, start, end = _parse_region(region)
    if start is None or end is None:
        return 0
    return abs(end - start) + 1


# --------------------------------------------------------------------------- #
# statistics
# --------------------------------------------------------------------------- #
def _median(xs) -> float:
    return float(statistics.median(xs)) if xs else 0.0


def _mad(xs) -> float:
    if not xs:
        return 0.0
    med = statistics.median(xs)
    return float(statistics.median([abs(x - med) for x in xs]))


def block_bootstrap_ci(
    locus: list[int],
    backbone: list[int],
    *,
    block: int = 500,
    n_boot: int = 1000,
    seed: int = 20260817,
    alpha: float = 0.05,
    max_backbone_blocks: int = 512,
) -> tuple[float, float]:
    """Percentile CI for median(locus)/median(backbone), resampling blocks.

    Depth at adjacent bases is highly correlated, so per-base resampling would
    produce a spuriously narrow interval. Blocks preserve that correlation.
    """
    if not locus or not backbone:
        return (0.0, 0.0)
    rng = random.Random(seed)

    if block <= 0 or n_boot <= 0 or max_backbone_blocks <= 0:
        raise ValueError("block, n_boot, and max_backbone_blocks must be positive")

    def blocks(xs: list[int]) -> list[list[int]]:
        return [xs[i:i + block] for i in range(0, len(xs), block)] or [xs]

    lb = blocks(locus)
    bb = _evenly_spaced_blocks(blocks(backbone), max_backbone_blocks)

    def resampled_median(source: list[list[int]]) -> float:
        indices = [rng.randrange(len(source)) for _ in range(len(source))]
        if _np is not None:
            arrays = [_np.asarray(source[index], dtype=_np.int32) for index in indices]
            return float(_np.median(_np.concatenate(arrays)))
        values = [value for index in indices for value in source[index]]
        return float(statistics.median(values))

    ratios = []
    for _ in range(n_boot):
        lm = resampled_median(lb)
        bm = resampled_median(bb)
        if bm > 0:
            ratios.append(lm / bm)
    if not ratios:
        return (0.0, 0.0)
    ratios.sort()
    lo = ratios[max(0, int(len(ratios) * alpha / 2) - 1)]
    hi = ratios[min(len(ratios) - 1, int(len(ratios) * (1 - alpha / 2)))]
    return (float(lo), float(hi))


def _evenly_spaced_blocks(
    blocks: list[list[int]], max_blocks: int,
) -> list[list[int]]:
    """Deterministically cap bootstrap cost while retaining genome-wide span.

    The point denominator still uses every eligible base. Only uncertainty
    resampling is capped. First and last blocks are retained and intervening
    blocks are selected at even index intervals, avoiding result-dependent or
    random selection of backbone regions.
    """
    if max_blocks <= 0:
        raise ValueError("max_blocks must be positive")
    if len(blocks) <= max_blocks:
        return blocks
    if max_blocks == 1:
        return [blocks[len(blocks) // 2]]
    indices = [round(i * (len(blocks) - 1) / (max_blocks - 1)) for i in range(max_blocks)]
    return [blocks[index] for index in indices]


# --------------------------------------------------------------------------- #
# GC handling
# --------------------------------------------------------------------------- #
def _quantile(xs, p: float) -> float:
    """Linear-interpolation quantile on an unsorted iterable."""
    v = sorted(xs)
    if not v:
        return 0.0
    if len(v) == 1:
        return float(v[0])
    i = (len(v) - 1) * p
    lo = int(i)
    hi = min(lo + 1, len(v) - 1)
    return float(v[lo] + (v[hi] - v[lo]) * (i - lo))






def _read_fasta(path: str | Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def _gc(seq: str) -> float:
    seq = seq.upper()
    n = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / n if n else 0.0


def _span_coverage(
    gene_spans: list[tuple[int, int]], gene_length: int,
) -> tuple[int, int, int]:
    """Total aligned, union-covered and overlapping query bases."""
    clipped: list[tuple[int, int]] = []
    for a, b in gene_spans:
        lo, hi = max(1, min(a, b)), min(gene_length, max(a, b))
        if hi >= lo:
            clipped.append((lo, hi))
    total = sum(hi - lo + 1 for lo, hi in clipped)
    if not clipped:
        return 0, 0, 0
    merged: list[list[int]] = []
    for lo, hi in sorted(clipped):
        if not merged or lo > merged[-1][1] + 1:
            merged.append([lo, hi])
        else:
            merged[-1][1] = max(merged[-1][1], hi)
    union = sum(hi - lo + 1 for lo, hi in merged)
    return total, union, total - union


def _query_segments(
    gene_spans: list[tuple[int, int]], gene_length: int,
) -> list[tuple[int, tuple[int, ...]]]:
    """Cut the query into maximal stretches of constant fragment membership.

    Returns `(n_bases, region_indices)` per stretch. Membership is a property
    of the geometry alone, so this is computed once and reused across bootstrap
    resamples, where only the ratios change.
    """
    if gene_length <= 0:
        raise ValueError("gene_length must be positive")
    clipped: list[tuple[int, int] | None] = []
    cuts = {1, gene_length + 1}
    for a, b in gene_spans:
        lo, hi = max(1, min(a, b)), min(gene_length, max(a, b))
        if hi >= lo:
            clipped.append((lo, hi))
            cuts.add(lo)
            cuts.add(hi + 1)
        else:
            clipped.append(None)
    points = sorted(c for c in cuts if 1 <= c <= gene_length + 1)
    segments = []
    for start, stop in zip(points, points[1:]):
        members = tuple(
            i for i, span in enumerate(clipped)
            if span is not None and span[0] <= start and stop - 1 <= span[1]
        )
        segments.append((stop - start, members))
    return segments


def _unit_ratio_dosage(spans, gene_length: int) -> float:
    """Dosage the same span geometry returns when every region is single copy.

    Mean number of spans covering each covered query base. Exact, parameter
    free, and a lower bound on geometric inflation: sequence represented
    several times also attracts multi-mapped reads, which raises the per-region
    ratios on top of this.
    """
    cover = [0] * gene_length
    for start, end in spans:
        lo = max(1, int(start))
        hi = min(gene_length, int(end))
        for i in range(lo - 1, hi):
            cover[i] += 1
    covered = [c for c in cover if c > 0]
    return (sum(covered) / len(covered)) if covered else 0.0


def query_copy_profile(
    segments: list[tuple[int, tuple[int, ...]]],
    region_ratios: list[float],
    gene_length: int,
) -> tuple[float | None, float, float, int, float]:
    """Total locus dosage as the mean copy multiplicity across query bases.

    The estimand is a per-query-base quantity. Under the collapse model each
    assembly fragment's ratio is the number of genome copies collapsed into it,
    so a query base resolved into k fragments carries the SUM of their ratios:

        R(j) = SUM over fragments covering query base j of ratio_r
        dosage = mean of R(j) over covered query bases

    Overlap between query spans is therefore not double counting. It is the
    signature of copies the assembler resolved into separate fragments, and
    summing is what recovers their total. Where the spans are disjoint and
    complete this reduces exactly to `length_weighted_copies`, and a single
    full-length fragment reduces it to `ratio`, so neither of those
    quantities is displaced.

    Incomplete coverage does not withhold the quantity, it widens it. The
    bracket returned assigns the uncovered bases zero copies (lower) and the
    largest multiplicity observed anywhere on the query (upper). When coverage
    is complete the three numbers coincide.

    Returns `(dosage, lower, upper, covered_bp, max_multiplicity)`.
    """
    if gene_length <= 0:
        raise ValueError("gene_length must be positive")
    total = 0.0
    covered = 0
    max_multiplicity = 0.0
    for n_bases, members in segments:
        multiplicity = sum(region_ratios[i] for i in members)
        if members:
            covered += n_bases
            max_multiplicity = max(max_multiplicity, multiplicity)
        total += n_bases * multiplicity
    lower = total / gene_length
    upper = (total + (gene_length - covered) * max_multiplicity) / gene_length
    dosage = (total / covered) if covered else None
    return dosage, lower, upper, covered, max_multiplicity


def block_bootstrap_dosage_ci(
    per_region_depth: list[list[int]],
    backbone: list[int],
    segments: list[tuple[int, tuple[int, ...]]],
    gene_length: int,
    *,
    block: int = 500,
    n_boot: int = 1000,
    seed: int = 20260817,
    alpha: float = 0.05,
    max_backbone_blocks: int = 512,
) -> tuple[float, float]:
    """Percentile CI for the dosage, resampling every region and the backbone.

    Rescaling the ratio interval would carry only the backbone's uncertainty.
    Each fragment contributes its own depth variability to the dosage, so each
    is resampled in blocks alongside the backbone within a single iteration.
    """
    if not backbone or not per_region_depth or not any(per_region_depth):
        return (0.0, 0.0)
    if block <= 0 or n_boot <= 0 or max_backbone_blocks <= 0:
        raise ValueError("block, n_boot, and max_backbone_blocks must be positive")
    rng = random.Random(seed)

    def blocks(xs: list[int]) -> list[list[int]]:
        return [xs[i:i + block] for i in range(0, len(xs), block)] or [xs]

    region_blocks = [blocks(d) if d else [[0]] for d in per_region_depth]
    bb = _evenly_spaced_blocks(blocks(backbone), max_backbone_blocks)

    def resampled_median(source: list[list[int]]) -> float:
        indices = [rng.randrange(len(source)) for _ in range(len(source))]
        if _np is not None:
            arrays = [_np.asarray(source[i], dtype=_np.int32) for i in indices]
            return float(_np.median(_np.concatenate(arrays)))
        return float(statistics.median(
            [v for i in indices for v in source[i]]))

    draws = []
    for _ in range(n_boot):
        bm = resampled_median(bb)
        if bm <= 0:
            continue
        ratios = [resampled_median(rb) / bm for rb in region_blocks]
        dosage, _lo, _hi, _cov, _mx = query_copy_profile(
            segments, ratios, gene_length)
        if dosage is not None:
            draws.append(dosage)
    if not draws:
        return (0.0, 0.0)
    draws.sort()
    lo = draws[max(0, int(len(draws) * alpha / 2) - 1)]
    hi = draws[min(len(draws) - 1, int(len(draws) * (1 - alpha / 2)))]
    return (float(lo), float(hi))


# --------------------------------------------------------------------------- #
# main entry point
# --------------------------------------------------------------------------- #
def length_weighted_copies(
    gene_spans: list[tuple[int, int]],
    region_ratios: list[float],
    gene_length: int,
) -> float:
    """Copies of the locus, weighting each assembly fragment by how much of the
    gene it actually spans.

        copies = SUM_r ( aligned_gene_bp_r * ratio_r ) / gene_length

    `gene_spans` are QUERY coordinates from the discovery BLAST -- how much of
    the bait each assembly fragment covers -- not the padded subject intervals.

    The two halves of the evidence are combined only when the query-coordinate
    spans are complete and non-overlapping: length captures what the assembler
    resolved, while depth captures what it collapsed. Neither works alone.

    A single-fragment locus reduces this to `ratio` exactly, so the estimator is
    a strict generalisation and cannot change any single-region result.

    Retained for the disjoint, complete geometry it was written for, where it
    equals the mean per-query-base multiplicity. `query_copy_profile` is the
    general form and is what the estimator now calls: it handles overlap by
    summing multiplicities, which is what the collapse model implies, and
    handles incomplete coverage by bracketing rather than withholding.
    """
    if gene_length <= 0:
        raise ValueError("gene_length must be positive")
    total = 0.0
    for (a, b), r in zip(gene_spans, region_ratios):
        lo, hi = max(1, min(a, b)), min(gene_length, max(a, b))
        if hi >= lo:
            total += (hi - lo + 1) * r
    return total / gene_length


def estimate_locus_copy_number(
    bam: str | Path,
    locus_regions: list[str],
    *,
    gene_spans: list[tuple[int, int]] | None = None,
    gene_length: int | None = None,
    reference_fasta: str | Path | None = None,
    exclude_regions: list[str] | None = None,
    min_bq: int = 20,
    min_mq: int = 20,
    min_contig_len: int = 5000,
    contig_end_trim: int = 100,
    block: int = 500,
    max_backbone_bootstrap_blocks: int = 512,
    n_boot: int = 1000,
    seed: int = 20260817,
    call_threshold: float = CALL_THRESHOLD_DEFAULT,
    call_threshold_calibration: float | None = None,
    ambiguity_warn: float = 0.70,
    samtools: str = "samtools",
) -> DepthRatioResult:
    """Estimate copies per collapsed representative for one locus.

    Parameters
    ----------
    bam
        Alignment of the sample's reads against the DRAFT ASSEMBLY -- the BAM
        the recruitment step already produces. Not the candidate remap.
    locus_regions
        One or more `contig:start-end` intervals for the locus, i.e. the
        discovery intervals that passed bait alignment.
    gene_spans, gene_length
        QUERY coordinates from the discovery BLAST: for each entry of
        `locus_regions`, in the same order, the `(start, end)` stretch of the
        bait that fragment aligns to, and the full length of one copy of the
        gene. Required to turn the ratio into a copy number when the assembler
        split the locus across several fragments; without them `copies_estimate`
        is left at None rather than guessed. Ignored (and inferred) when the
        locus is a single region, where the estimate equals `ratio` by
        construction. See `validation/depth-copy-number/`.
    exclude_regions
        Regions to keep out of the backbone: rRNA operons, IS elements, known
        repeat families. Excluding them matters because their own collapse
        would drag the backbone median upward and mask a real amplification.
    """
    bam = Path(bam)
    notes: list[str] = []
    excl = exclude_regions or []
    contig_lengths = _contig_lengths(bam, samtools=samtools)

    # ---- numerator -------------------------------------------------------- #
    locus_depth: list[int] = []
    locus_depth_raw: list[int] = []
    per_region_depth: list[list[int]] = []       # copies_estimate needs it per region
    for reg in locus_regions:
        d = depth_array(bam, reg, min_bq=min_bq, min_mq=min_mq, samtools=samtools)
        per_region_depth.append(d)
        locus_depth += d
        locus_depth_raw += depth_array(bam, reg, min_bq=0, min_mq=0,
                                       samtools=samtools)
    d_locus = _median(locus_depth)
    d_locus_raw = _median(locus_depth_raw)

    ambiguity = (d_locus / d_locus_raw) if d_locus_raw > 0 else 0.0
    reliable = True
    if ambiguity < ambiguity_warn:
        reliable = False
        notes.append(
            f"only {ambiguity:.0%} of raw depth at the locus survives the MAPQ filter; "
            "reads are mapping ambiguously, which happens when the assembler kept "
            "more than one copy. The ratio understates the true copy number."
        )

    # ---- denominator ------------------------------------------------------ #
    parsed_excl: list[tuple[str, int, int]] = []
    for reg in list(locus_regions) + excl:
        c, s, e = _parse_region(reg)
        parsed_excl.append((c, s or 1, e or 10**12))

    ref = _read_fasta(reference_fasta) if reference_fasta else {}
    gc_locus = None
    if ref:
        chunks = []
        for reg in locus_regions:
            c, s, e = _parse_region(reg)
            if c in ref and s and e:
                chunks.append(ref[c][s - 1:e])
        if chunks:
            gc_locus = _gc("".join(chunks))

    # The denominator is the whole backbone: every contig above the length
    # floor, ends trimmed, locus and excluded intervals removed. A GC-matched
    # denominator is deliberately not offered: it has not been validated against
    # a held-out panel with matched backgrounds, and on the panel where a
    # matched background was available it changed no verdict. `gc_locus` is
    # reported as a description of the locus rather than as a correction
    # applied to it.
    backbone = []
    for contig, length in contig_lengths.items():
        if length < min_contig_len:
            continue
        lo, hi = contig_end_trim + 1, length - contig_end_trim
        if hi <= lo:
            continue
        vals = depth_array(bam, f"{contig}:{lo}-{hi}",
                           min_bq=min_bq, min_mq=min_mq, samtools=samtools)
        # drop excluded intervals falling on this contig
        spans = [(s, e) for c, s, e in parsed_excl if c == contig]
        if spans:
            keep = []
            for i, v in enumerate(vals, start=lo):
                if not any(s <= i <= e for s, e in spans):
                    keep.append(v)
            vals = keep
        backbone += vals

    d_backbone = _median(backbone)
    mad = _mad(backbone)
    if d_backbone <= 0:
        notes.append("backbone median depth is zero; ratio undefined")
        return DepthRatioResult(
            locus_regions=list(locus_regions), locus_median_depth=d_locus,
            locus_median_depth_unfiltered=d_locus_raw, backbone_median_depth=0.0,
            backbone_mad=0.0, ratio=0.0, ratio_ci_low=0.0, ratio_ci_high=0.0,
            assembly_copies=len(locus_regions), copies_estimate=None,
            ambiguity_index=ambiguity, gc_locus=gc_locus,
            single_copy_rejected=False, reliable=False, notes=notes,
            assembly_region_count=len(locus_regions), dosage_estimate=None,
            dosage_status="NOT_ESTIMATED", depth_call="DEPTH_INDETERMINATE",
            bootstrap_block_size=block,
            bootstrap_locus_blocks=(len(locus_depth) + block - 1) // block,
            bootstrap_backbone_blocks_total=(len(backbone) + block - 1) // block,
            bootstrap_backbone_blocks_used=min(
                (len(backbone) + block - 1) // block, max_backbone_bootstrap_blocks
            ),
            bootstrap_resamples=n_boot, bootstrap_seed=seed)

    ratio = d_locus / d_backbone
    ci_low, ci_high = block_bootstrap_ci(locus_depth, backbone, block=block,
                                         n_boot=n_boot, seed=seed,
                                         max_backbone_blocks=max_backbone_bootstrap_blocks)

    if mad > 0 and mad / d_backbone > 0.35:
        notes.append(
            f"backbone depth is highly dispersed (MAD/median = {mad/d_backbone:.2f}); "
            "uneven coverage weakens the ratio")

    n_asm = len(locus_regions)

    # ---- dosage, weighted by resolved query length ------------------------- #
    # `assembly_region_count` is a geometry fact. `dosage_estimate` combines
    # that geometry with depth: the mean copy multiplicity per query base.
    # Incomplete coverage widens the bracket rather than withholding the
    # quantity; overlap raises the multiplicity rather than double-counting.
    region_ratios = [_median(d) / d_backbone for d in per_region_depth]
    spans: list[tuple[int, int]] = list(gene_spans or [])
    g_len = gene_length or 0
    dosage: float | None = None
    dosage_lower: float | None = None
    dosage_upper: float | None = None
    dosage_ci_low: float | None = None
    dosage_ci_high: float | None = None
    query_coverage_fraction = 0.0
    query_multiplicity_max = 0.0
    dosage_flags: list[str] = []
    geometry_expected_dosage: float | None = None
    total_aligned_bp = union_covered_bp = overlap_bp = 0

    if spans and len(spans) != n_asm:
        notes.append(
            f"gene_spans has {len(spans)} entries for {n_asm} locus regions; "
            "they must correspond one to one. Copy estimate not computed."
        )
        spans = []
        dosage_flags.append("GENE_SPAN_COUNT_MISMATCH")
    if not spans and n_asm == 1:
        # A single fragment spans the whole thing by definition; this is the
        # branch that makes the estimator reduce to `ratio` exactly.
        c, s, e = _parse_region(locus_regions[0])
        if s and e:
            spans, g_len = [(1, e - s + 1)], g_len or (e - s + 1)
    if spans and g_len > 0:
        total_aligned_bp, union_covered_bp, overlap_bp = _span_coverage(spans, g_len)
        segments = _query_segments(spans, g_len)
        raw, lower, upper, covered_bp, mult_max = query_copy_profile(
            segments, region_ratios, g_len)
        query_coverage_fraction = round(covered_bp / g_len, 4)
        query_multiplicity_max = round(mult_max, 3)
        if raw is not None:
            dosage = round(raw, 3)
            dosage_lower = round(lower, 3)
            dosage_upper = round(upper, 3)
            ci = block_bootstrap_dosage_ci(
                per_region_depth, backbone, segments, g_len,
                block=block, n_boot=n_boot, seed=seed,
                max_backbone_blocks=max_backbone_bootstrap_blocks)
            dosage_ci_low, dosage_ci_high = round(ci[0], 3), round(ci[1], 3)
        if union_covered_bp < g_len:
            dosage_flags.append("INCOMPLETE_QUERY_COVERAGE")
            notes.append(
                f"the discovery intervals cover {union_covered_bp}/{g_len} unique "
                f"query bases ({query_coverage_fraction:.1%}); dosage is the mean "
                f"multiplicity over covered bases, and the bracket "
                f"[{dosage_lower}, {dosage_upper}] assigns the {g_len - union_covered_bp} "
                "uncovered bases zero copies and the largest multiplicity observed "
                "on the query respectively."
            )
        if overlap_bp > 0:
            dosage_flags.append("OVERLAPPING_QUERY_SPANS")
            geometry_expected_dosage = round(
                _unit_ratio_dosage(gene_spans, g_len), 3)
            notes.append(
                f"query-coordinate spans overlap by {overlap_bp} bp "
                f"({overlap_bp/g_len:.1%} of the gene); under the collapse model a "
                "query base resolved into several fragments carries the sum of "
                f"their ratios, so overlap raises the per-base multiplicity (max "
                f"{query_multiplicity_max}) rather than counting bases twice. "
                f"The same geometry with every region at single copy would return "
                f"{geometry_expected_dosage}, so the dosage is NOT a copy count "
                "here; that figure is itself a lower bound on the inflation, "
                "because multiply represented sequence also attracts multi-mapped "
                "reads that raise the per-region ratios."
            )
        zero_regions = [i for i, r in enumerate(region_ratios) if r == 0.0]
        if zero_regions:
            dosage_flags.append("ZERO_RATIO_REGIONS")
            notes.append(
                f"{len(zero_regions)} of {len(region_ratios)} locus regions have zero "
                "filtered depth, so the query bases they resolve contribute nothing "
                "to the dosage. The estimate is biased downward by that amount."
            )
    elif n_asm > 1:
        dosage_flags.append("MISSING_QUERY_GEOMETRY")
        notes.append(
            f"the locus is split across {n_asm} assembly fragments, so the global "
            "median cannot be turned into a copy number without the gene "
            "coordinates each fragment covers. Pass gene_spans and gene_length "
            "(the discovery BLAST query coordinates) to get a copy estimate. "
            "Multiplying the ratio by the fragment count is NOT a valid "
            "substitute -- it returns ~n on single-copy loci split into n pieces."
        )

    if not reliable:
        dosage_flags.append("LOW_MAPPING_UNIQUENESS")
    # Overlap IS a review state. A single-copy query whose bases are resolved
    # by several overlapping intervals returns a dosage above one from the
    # geometry alone, and an empirical control on this genome shows the
    # returned value rising with the overlapped fraction on sequence known to
    # be single copy (validation/depth-copy-number/overlap_null/). The estimand
    # does not separate that inflation from real copy number, so the status
    # says so rather than presenting the number as a count.
    if dosage is None:
        dosage_status = "NOT_ESTIMATED"
    elif not reliable:
        dosage_status = "LOWER_BOUND"
    elif overlap_bp > 0:
        dosage_status = "REVIEW_OVERLAPPING_SPANS"
    elif union_covered_bp < g_len:
        dosage_status = "ESTIMATED_PARTIAL_COVERAGE"
    else:
        dosage_status = "ESTIMATED"

    # Low ambiguity annotates the ratio; it does not withhold the verdict. The
    # bias runs one way only: filtered depth can never exceed raw depth, and the
    # single-copy backbone is high-MAPQ by construction, so it is not depressed.
    # Requiring `reliable` before rejecting single copy would therefore use that
    # downward bias to suppress precisely the conclusion the bias protects -- a
    # locus whose repeated flanks collapse under the MAPQ filter would be called
    # single-copy consistent while its own lower confidence bound sits far above
    # the threshold. Low ambiguity marks the ratio as a lower bound instead.
    # The threshold may be raised by a per-sample null but never lowered: a
    # calibration that came out below the shipped constant would mean this
    # sample's backbone is quieter than the two used to justify 1.5, and
    # trading that quiet for easier rejection is exactly the move the constant
    # exists to prevent.
    threshold_used = call_threshold
    threshold_source = ("DEFAULT_ASSERTED" if call_threshold == CALL_THRESHOLD_DEFAULT
                        else "USER_SUPPLIED")
    if call_threshold_calibration is not None:
        if call_threshold_calibration > call_threshold:
            threshold_used = float(call_threshold_calibration)
            threshold_source = "CALIBRATED_RAISED"
            notes.append(
                f"call threshold raised from {call_threshold:.2f} to "
                f"{threshold_used:.2f} by this sample's single-copy null")
        else:
            threshold_source = "CALIBRATED_FLOOR_HELD"
            notes.append(
                f"this sample's single-copy null gave "
                f"{call_threshold_calibration:.2f}, below the configured "
                f"{call_threshold:.2f}; the configured value was kept")

    rejected = bool(ci_low > threshold_used)

    # An interval that excludes 1 is not by itself evidence against single
    # copy. On single-copy segments of two real backbones this interval
    # excluded 1 in 59% and 65% of cases (above and below combined), because it
    # quantifies within-locus sampling only. Saying so in the record keeps a
    # reader from reading "CI 1.111-1.349" as "more than one copy" when the
    # verdict says otherwise.
    if ratio > 0 and not (ci_low <= 1.0 <= ci_high):
        side = "above" if ci_low > 1.0 else "below"
        notes.append(
            f"the interval [{ci_low:.3f}, {ci_high:.3f}] excludes 1 from {side}, "
            "which does not on its own indicate a departure from single copy: it "
            "covers within-locus sampling noise only, and single-copy loci "
            "commonly produce intervals that exclude 1. The verdict, not the "
            "interval, carries the copy-number statement.")

    if not [k for k in range(0, 64) if ci_low <= k <= ci_high]:
        notes.append(
            "the interval contains no integer; it is narrower than the spacing "
            "between adjacent copy numbers and cannot be read as an interval "
            "for an integer copy count.")

    if rejected and not reliable:
        notes.append(
            f"the ratio ({ratio:.2f}) is a LOWER BOUND: only {ambiguity:.0%} of raw "
            "depth at the locus survives the MAPQ filter, and that filter can only "
            "remove depth, never add it. The extra copies are established; how many "
            "there are is not."
        )

    depth_call = "MULTICOPY_DEPTH" if rejected else "SINGLE_COPY_COMPATIBLE"

    if 0 < ratio < 0.6:
        notes.append("ratio well below 1: possible partial deletion, "
                     "subpopulation loss, or a locus present in only some cells")

    return DepthRatioResult(
        locus_regions=list(locus_regions),
        locus_median_depth=round(d_locus, 2),
        locus_median_depth_unfiltered=round(d_locus_raw, 2),
        backbone_median_depth=round(d_backbone, 2),
        backbone_mad=round(mad, 2),
        ratio=round(ratio, 3),
        ratio_ci_low=round(ci_low, 3),
        ratio_ci_high=round(ci_high, 3),
        assembly_copies=n_asm,
        copies_estimate=dosage,
        ambiguity_index=round(ambiguity, 3),
        gc_locus=round(gc_locus, 4) if gc_locus is not None else None,
        single_copy_rejected=rejected,
        reliable=reliable,
        notes=notes,
        call_threshold_used=round(threshold_used, 3),
        call_threshold_source=threshold_source,
        assembly_region_count=n_asm,
        dosage_estimate=dosage,
        dosage_lower=dosage_lower,
        dosage_upper=dosage_upper,
        dosage_ci_low=dosage_ci_low,
        dosage_ci_high=dosage_ci_high,
        query_coverage_fraction=query_coverage_fraction,
        query_multiplicity_max=query_multiplicity_max,
        dosage_status=dosage_status,
        geometry_expected_dosage=geometry_expected_dosage,
        dosage_flags=dosage_flags,
        depth_call=depth_call,
        bootstrap_block_size=block,
        bootstrap_locus_blocks=(len(locus_depth) + block - 1) // block,
        bootstrap_backbone_blocks_total=(len(backbone) + block - 1) // block,
        bootstrap_backbone_blocks_used=min(
            (len(backbone) + block - 1) // block, max_backbone_bootstrap_blocks
        ),
        bootstrap_resamples=n_boot,
        bootstrap_seed=seed,
        gene_length=g_len,
        gene_covered_bp=union_covered_bp,
        gene_total_aligned_bp=total_aligned_bp,
        gene_union_covered_bp=union_covered_bp,
        gene_overlap_bp=overlap_bp,
        gene_overlap_fraction=round(overlap_bp / g_len, 6) if g_len else 0.0,
        region_gene_spans=spans,
        region_ratios=[round(r, 3) for r in region_ratios],
    )


# --------------------------------------------------------------------------- #
# CLI (also usable standalone, before integration)
# --------------------------------------------------------------------------- #
def _cli() -> int:
    import argparse
    import csv
    import sys

    ap = argparse.ArgumentParser(
        description="Normalised locus depth ratio (copies per collapsed representative)")
    ap.add_argument("--bam", required=True, help="reads mapped to the DRAFT assembly")
    ap.add_argument("--locus", required=True, action="append",
                    help="contig:start-end; repeat for multiple discovery intervals")
    ap.add_argument("--reference", help="draft assembly FASTA (reports locus GC)")
    ap.add_argument("--exclude", action="append", default=[],
                    help="region to keep out of the backbone; repeatable")
    ap.add_argument("--gene-span", action="append", default=[], metavar="START-END",
                    help="query coordinates of the bait covered by each --locus "
                         "interval, in the same order and repeated once per "
                         "interval; required with --gene-length to obtain a copy "
                         "estimate for a locus split across several fragments")
    ap.add_argument("--gene-length", type=int,
                    help="length of one copy of the gene (bp); pair with --gene-span")
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--min-mq", type=int, default=20)
    ap.add_argument("--call-threshold", type=float, default=1.5)
    ap.add_argument("--seed", type=int, default=20260817)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--max-backbone-bootstrap-blocks", type=int, default=512,
                    help="deterministic cap on 500 bp backbone blocks resampled "
                         "for the confidence interval; the point denominator "
                         "still uses every eligible base")
    ap.add_argument("--tsv", help="write one-row TSV here")
    args = ap.parse_args()

    spans: list[tuple[int, int]] | None = None
    if args.gene_span:
        if len(args.gene_span) != len(args.locus):
            ap.error("--gene-span must be given once per --locus, in the same order")
        spans = []
        for sp in args.gene_span:
            try:
                a, b = sp.replace(",", "").split("-")
                spans.append((int(a), int(b)))
            except ValueError:
                ap.error(f"--gene-span expects START-END, got {sp!r}")
    if spans and not args.gene_length:
        ap.error("--gene-span requires --gene-length")

    res = estimate_locus_copy_number(
        args.bam, args.locus,
        gene_spans=spans, gene_length=args.gene_length,
        reference_fasta=args.reference,
        exclude_regions=args.exclude, min_bq=args.min_bq, min_mq=args.min_mq,
        call_threshold=args.call_threshold,
        seed=args.seed, n_boot=args.n_boot,
        max_backbone_bootstrap_blocks=args.max_backbone_bootstrap_blocks)

    print(res.summary_line())
    for n in res.notes:
        print(f"  note: {n}")

    if args.tsv:
        row = res.as_row()
        with open(args.tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
            w.writeheader()
            w.writerow(row)
        print(f"wrote {args.tsv}", file=sys.stderr)
    return 0


def main() -> int:
    """Console-script entry point (`locus-recon-depth-ratio`)."""
    return _cli()


if __name__ == "__main__":
    raise SystemExit(main())
