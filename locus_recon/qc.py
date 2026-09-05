"""
Allele quality control: bait database profiling, metric computation,
confidence classification, and QC report formatting.
"""

import math
from collections import Counter
from typing import List, Optional, Tuple

from .utils import (
    DESCRIPTIVE_FLAG_PREFIXES,
    PER_BASE_THRESHOLDS,
    SPAN_CLIP_THRESHOLDS,
    QC_THRESHOLDS,
    REMAP_THRESHOLDS,
    STOP_CODONS,
    UNCERTAINTY_THRESHOLDS,
    log,
)
from .io import read_fasta_sequences

VALID_DNA = set("ACGTRYSWKMBDHVN")


# ==============================================================================
# Percentile helper
# ==============================================================================


def _percentile(sorted_vals: List[float], p: float) -> float:
    """Compute the p-th percentile (0–100) from a pre-sorted list.

    Args:
        sorted_vals: Sorted list of numeric values.
        p:           Percentile to compute (0–100).

    Returns:
        Interpolated percentile value.
    """
    if not sorted_vals:
        return 0.0
    k = (len(sorted_vals) - 1) * (p / 100.0)
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return float(sorted_vals[int(k)])
    return sorted_vals[int(f)] * (c - k) + sorted_vals[int(c)] * (k - f)


# ==============================================================================
# Bait database profiling
# ==============================================================================


def profile_bait_database(bait_path: str) -> dict:
    """Analyse all known alleles in a bait FASTA to build an expected profile.

    Computes length distribution statistics (median, IQR, mode, stdev) and
    mean GC content. These values drive QC classification thresholds.

    Args:
        bait_path: Path to the multi-FASTA bait file (all known alleles).

    Returns:
        Dict with keys: n_alleles, lengths, len_min, len_max, len_median,
        len_mean, len_stdev, len_q1, len_q3, len_iqr, len_mode,
        gc_mean, gc_stdev.
    """
    log.info("Profiling bait database (expected allele characteristics)...")
    sequences = read_fasta_sequences(bait_path)
    if not sequences:
        raise ValueError(f"Bait FASTA is empty or unreadable: {bait_path}")

    seen_headers = set()
    for header, sequence in sequences:
        if header in seen_headers:
            raise ValueError(f"Bait FASTA contains duplicate sequence ID {header!r}.")
        if not sequence:
            raise ValueError(f"Bait FASTA record {header!r} has an empty sequence.")
        invalid = sorted(set(sequence) - VALID_DNA)
        if invalid:
            raise ValueError(
                f"Bait FASTA record {header!r} contains invalid DNA symbols: "
                + ", ".join(invalid)
            )
        seen_headers.add(header)

    lengths: List[int] = sorted(len(seq) for _, seq in sequences)
    gc_values: List[float] = []
    for _, seq in sequences:
        if seq:
            gc_values.append(100.0 * (seq.count("G") + seq.count("C")) / len(seq))

    n = len(lengths)
    mean_len = sum(lengths) / n
    var_len  = sum((x - mean_len) ** 2 for x in lengths) / n if n > 1 else 0.0
    stdev_len = math.sqrt(var_len)

    q1  = _percentile(lengths, 25)
    q3  = _percentile(lengths, 75)
    iqr = q3 - q1

    mode_len = Counter(lengths).most_common(1)[0][0]
    mode_mod3 = Counter(length % 3 for length in lengths).most_common(1)[0][0]

    frame_stop_counts = {
        frame: sum(_count_internal_stops(sequence, frame) for _, sequence in sequences)
        for frame in range(3)
    }
    minimum_stops = min(frame_stop_counts.values())
    expected_frames = [
        frame for frame, stop_count in frame_stop_counts.items()
        if stop_count == minimum_stops
    ]
    expected_frame = expected_frames[0]

    gc_mean  = sum(gc_values) / len(gc_values) if gc_values else 0.0
    gc_var   = (
        sum((g - gc_mean) ** 2 for g in gc_values) / len(gc_values)
        if len(gc_values) > 1 else 0.0
    )
    gc_stdev = math.sqrt(gc_var)

    profile = {
        "n_alleles":   n,
        "lengths":     lengths,
        "len_min":     lengths[0],
        "len_max":     lengths[-1],
        "len_median":  _percentile(lengths, 50),
        "len_mean":    mean_len,
        "len_stdev":   stdev_len,
        "len_q1":      q1,
        "len_q3":      q3,
        "len_iqr":     iqr,
        "len_mode":    mode_len,
        "len_mod3_mode": mode_mod3,
        "expected_coding_frame": expected_frame,
        "expected_coding_frames": expected_frames,
        "frame_stop_counts": frame_stop_counts,
        "gc_mean":     gc_mean,
        "gc_stdev":    gc_stdev,
    }

    log.info(f"  Bait alleles : {n}")
    log.info(f"  Length range : {lengths[0]}-{lengths[-1]} bp")
    log.info(f"  Length median: {profile['len_median']:.0f} bp "
             f"(Q1={q1:.0f}, Q3={q3:.0f}, IQR={iqr:.0f})")
    log.info(f"  Length mode  : {mode_len} bp")
    log.info(f"  Modal len %% 3: {mode_mod3}")
    log.info(
        "  Plausible coding frame(s): "
        + ", ".join(f"+{frame}" for frame in expected_frames)
    )
    log.info(f"  Length stdev : {stdev_len:.1f} bp")
    log.info(f"  GC content   : {gc_mean:.1f}% +/- {gc_stdev:.1f}%")
    if iqr == 0:
        log.info("  Note: IQR=0 (all alleles same length). "
                 "Will use stdev as fallback tolerance.")

    return profile


# ==============================================================================
# QC assessment
# ==============================================================================


def classify_catalogue_status(
    validation_identity: float,
    exact_known_allele: Optional[bool] = None,
) -> dict:
    """Describe the reconstruction's relationship to the supplied catalogue.

    This is deliberately independent of sequence confidence.  Sequence
    confidence answers whether the reported bases are supported by the reads;
    catalogue status answers whether the resulting sequence is already known.
    A divergent sequence may be an accurate reconstruction of a genuinely novel
    allele, so catalogue status never constrains the confidence tier.

    Args:
        validation_identity: % identity to the nearest catalogue allele.
        exact_known_allele:  True when the reconstruction is byte-equivalent to
                             a catalogue record (100% identity, 100% query
                             coverage, equal aligned/query/subject lengths).
                             When None, exactness is inferred from identity
                             alone.

    Returns:
        Dict with 'catalogue_status', 'nearest_allele_identity' and
        'catalogue_review_recommended'.
    """
    th = QC_THRESHOLDS
    if exact_known_allele is None:
        exact_known_allele = validation_identity >= 100.0

    if exact_known_allele:
        status = "EXACT_MATCH"
    elif validation_identity < th["low"]["min_identity"]:
        status = "HIGHLY_DIVERGENT_FROM_REFERENCE"
    elif validation_identity < th["high"]["min_identity"]:
        status = "DIVERGENT_FROM_REFERENCE"
    else:
        status = "NONEXACT_MATCH"

    return {
        "catalogue_status": status,
        "nearest_allele_identity": validation_identity,
        "catalogue_review_recommended": status != "EXACT_MATCH",
    }


def assess_allele_quality(
    allele_seq: str,
    bait_profile: dict,
    validation_identity: float,
    validation_qcov: float,
    locus: str,
    remap_metrics: Optional[dict] = None,
    ambiguity_metrics: Optional[dict] = None,
    expect_cds: bool = True,
    span_metrics: Optional[dict] = None,
    exact_known_allele: Optional[bool] = None,
) -> dict:
    """Compute quality metrics for a reconstructed allele vs the bait profile.

    Args:
        allele_seq:           Reconstructed allele sequence (uppercase).
        bait_profile:         Dict from profile_bait_database().
        validation_identity:  % identity from allele-vs-bait BLAST.
        validation_qcov:      % query coverage from allele-vs-bait BLAST.
        locus:                Locus name (informational).
        remap_metrics:        Dict from remap_reads_to_final_allele(), or None.
        ambiguity_metrics:    Dict from summarize_hit_ambiguity(), or None.
        expect_cds:           If False, skip CDS-specific checks (--noncoding_locus).
        span_metrics:         Dict from describe_full_query_span(), optionally
                              carrying a 'continuation' entry from
                              find_span_continuation(), or None.

    Returns:
        QC metrics dict including 'confidence' (HIGH/MEDIUM/LOW/SUSPECT) and
        'flags' list.
    """
    allele_len = len(allele_seq)

    median   = bait_profile["len_median"]
    iqr      = bait_profile["len_iqr"]
    stdev    = bait_profile["len_stdev"]
    mode     = bait_profile["len_mode"]

    length_delta     = allele_len - int(round(median))
    length_delta_abs = abs(length_delta)

    if stdev > 0.5:
        length_zscore = (allele_len - bait_profile["len_mean"]) / stdev
    else:
        length_zscore = 0.0 if length_delta_abs <= 1 else 999.0

    eff_tolerance    = iqr if iqr > 0 else max(stdev, 3.0)
    length_within_1x = length_delta_abs <= eff_tolerance * 1.0
    length_within_2x = length_delta_abs <= eff_tolerance * 2.0
    length_within_3x = length_delta_abs <= eff_tolerance * 3.0

    n_count    = allele_seq.count("N")
    n_fraction = n_count / allele_len if allele_len > 0 else 0.0
    gc_count   = allele_seq.count("G") + allele_seq.count("C")
    gc_pct     = 100.0 * gc_count / allele_len if allele_len > 0 else 0.0
    gc_deviation = abs(gc_pct - bait_profile["gc_mean"])

    if expect_cds:
        expected_frames = bait_profile.get(
            "expected_coding_frames",
            [bait_profile.get("expected_coding_frame", 0)],
        )
        expected_frame = expected_frames[0]
        expected_mod3 = bait_profile.get("len_mod3_mode", mode % 3)
        internal_stops = min(
            _count_internal_stops(allele_seq, frame) for frame in expected_frames
        )
        has_start_codon = allele_seq[:3] in ("ATG", "GTG", "TTG")
        has_stop_codon  = allele_seq[-3:] in STOP_CODONS if len(allele_seq) >= 3 else False
        length_mod3     = allele_len % 3
        frame_disrupted = length_mod3 != expected_mod3
    else:
        expected_frame = 0
        expected_frames = [0]
        expected_mod3 = 0
        internal_stops  = 0
        has_start_codon = False
        has_stop_codon  = False
        length_mod3     = 0
        frame_disrupted = False

    remap_metrics     = remap_metrics or {}
    ambiguity_metrics = ambiguity_metrics or {}

    confidence, flags = _classify_confidence(
        length_within_1x=length_within_1x,
        length_within_2x=length_within_2x,
        length_within_3x=length_within_3x,
        length_delta_abs=length_delta_abs,
        identity=validation_identity,
        qcov=validation_qcov,
        n_fraction=n_fraction,
        gc_deviation=gc_deviation,
        internal_stops=internal_stops,
        length_mod3=length_mod3,
        expected_mod3=expected_mod3,
        frame_disrupted=frame_disrupted,
        remap_metrics=remap_metrics,
        ambiguity_metrics=ambiguity_metrics,
        expect_cds=expect_cds,
        span_metrics=span_metrics,
    )

    catalogue = classify_catalogue_status(
        validation_identity=validation_identity,
        exact_known_allele=exact_known_allele,
    )

    return {
        "allele_length":     allele_len,
        "expected_median":   median,
        "expected_mode":     mode,
        "expected_range":    f"{bait_profile['len_min']}-{bait_profile['len_max']}",
        "expected_iqr_range": f"{median - eff_tolerance:.0f}-{median + eff_tolerance:.0f}",
        "length_delta":      length_delta,
        "length_delta_abs":  length_delta_abs,
        "length_zscore":     length_zscore,
        "length_within_iqr": length_within_1x,
        "effective_tolerance": eff_tolerance,
        "validation_identity": validation_identity,
        "validation_qcov":   validation_qcov,
        "n_count":           n_count,
        "n_fraction":        n_fraction,
        "gc_pct":            gc_pct,
        "gc_expected":       bait_profile["gc_mean"],
        "gc_deviation":      gc_deviation,
        "internal_stops":    internal_stops,
        "has_start_codon":   has_start_codon,
        "has_stop_codon":    has_stop_codon,
        "length_mod3":       length_mod3,
        "expected_length_mod3": expected_mod3,
        "expected_coding_frame": expected_frame,
        "evaluated_coding_frames": expected_frames,
        "frame_disrupted": frame_disrupted,
        "expect_cds":        expect_cds,
        "remap_mean_depth":       remap_metrics.get("mean_depth", 0.0),
        "remap_min_depth":        remap_metrics.get("min_depth", 0),
        "remap_max_depth":        remap_metrics.get("max_depth", 0),
        "remap_breadth_pct":      remap_metrics.get("breadth_pct", 0.0),
        "remap_pct_bases_lt5":    remap_metrics.get("pct_bases_lt5", 100.0),
        "remap_pct_bases_lt10":   remap_metrics.get("pct_bases_lt10", 100.0),
        "remap_mapped_reads":     remap_metrics.get("mapped_reads", 0),
        "remap_mapped_pct":       remap_metrics.get("mapped_pct", 0.0),
        "per_base_available":      remap_metrics.get("per_base_available", False),
        "uncertain_base_count":    remap_metrics.get("uncertain_base_count", 0),
        "uncertain_base_fraction": remap_metrics.get("uncertain_base_fraction", 0.0),
        "interior_uncertain_base_count": remap_metrics.get(
            "interior_uncertain_base_count", 0),
        "interior_uncertain_base_fraction": remap_metrics.get(
            "interior_uncertain_base_fraction", 0.0),
        "tier_uncertain_base_fraction": remap_metrics.get(
            "tier_uncertain_base_fraction", 0.0),
        "interior_metric_used":    remap_metrics.get("interior_metric_used", True),
        "interior_length":         remap_metrics.get("interior_length", 0),
        "internal_zero_depth_positions": remap_metrics.get(
            "internal_zero_depth_positions", 0),
        "terminal_margin_bp":      remap_metrics.get("terminal_margin_bp", 0),
        "low_coverage_position_count": remap_metrics.get(
            "low_coverage_position_count", 0),
        "mean_base_quality":       remap_metrics.get("mean_base_quality", 0.0),
        "mean_mapping_quality":    remap_metrics.get("mean_mapping_quality", 0.0),
        "strand_balance_pct":      remap_metrics.get("strand_balance_pct", 0.0),
        "candidate_mixed_sites":   remap_metrics.get("candidate_mixed_sites", 0),
        "mixed_site_count":        remap_metrics.get("mixed_site_count", 0),
        "strand_biased_sites":     remap_metrics.get("strand_biased_sites", 0),
        "reference_discordant_sites": remap_metrics.get("reference_discordant_sites", 0),
        "max_alt_fraction":        remap_metrics.get("max_alt_fraction", 0.0),
        "median_mixed_fraction":   remap_metrics.get("median_mixed_fraction", 0.0),
        "mixture_detected":        remap_metrics.get("mixture_detected", False),
        "second_best_bitscore":   ambiguity_metrics.get("second_best_bitscore", 0.0),
        "second_best_identity":   ambiguity_metrics.get("second_best_identity", 0.0),
        "bitscore_margin":        ambiguity_metrics.get("bitscore_margin", 999_999.0),
        "hit_is_ambiguous":       ambiguity_metrics.get("ambiguous", False),
        # 'confidence' is the sequence confidence: it is decided by read
        # support, sequence integrity and span completeness only.
        # 'sequence_confidence' is the explicit alias; 'confidence' is retained
        # for backwards compatibility with existing parsers.
        "confidence":  confidence,
        "sequence_confidence": confidence,
        "flags":       flags,
        **catalogue,
    }


def _count_internal_stops(seq: str, frame: int = 0) -> int:
    """Count internal stop codons in one defined reading frame.

    Args:
        seq: DNA sequence (uppercase).

    Returns:
        Number of internal stops in the requested frame.  The final complete
        codon is excluded because it may be the locus' expected terminal stop.
    """
    if frame not in (0, 1, 2):
        raise ValueError(f"Reading frame must be 0, 1, or 2; found {frame}.")
    codons = [seq[i : i + 3] for i in range(frame, len(seq) - 2, 3)]
    return sum(1 for codon in codons[:-1] if codon in STOP_CODONS)


def _classify_confidence(
    length_within_1x: bool,
    length_within_2x: bool,
    length_within_3x: bool,
    length_delta_abs: float,
    identity: float,
    qcov: float,
    n_fraction: float,
    gc_deviation: float,
    internal_stops: int,
    length_mod3: int,
    expected_mod3: int,
    frame_disrupted: bool,
    remap_metrics: dict,
    ambiguity_metrics: dict,
    expect_cds: bool,
    span_metrics: Optional[dict] = None,
) -> Tuple[str, List[str]]:
    """Assign a sequence confidence tier and generate explanatory flag strings.

    The tier reflects read support, sequence integrity and span completeness.
    Catalogue relationship is annotated here but excluded from the tier via
    DESCRIPTIVE_FLAG_PREFIXES; see classify_catalogue_status().

    Args:
        See assess_allele_quality() for parameter descriptions.

    Returns:
        (confidence_str, flags_list) where confidence is one of
        HIGH / MEDIUM / LOW / SUSPECT.
    """
    flags: List[str] = []
    th       = QC_THRESHOLDS
    remap_th = REMAP_THRESHOLDS

    if not length_within_1x:
        if length_within_2x:
            flags.append(f"LENGTH_MARGINAL (delta={length_delta_abs:.0f} bp, outside 1xIQR)")
        elif length_within_3x:
            flags.append(f"LENGTH_DEVIANT (delta={length_delta_abs:.0f} bp, outside 2xIQR)")
        else:
            flags.append(f"LENGTH_ANOMALOUS (delta={length_delta_abs:.0f} bp, outside 3xIQR)")

    # Catalogue relationship is reported, not scored.  Divergence from the
    # nearest curated allele is a statement about the reference catalogue, not
    # about whether the reported bases are supported by the reads, so these
    # flags are descriptive and never constrain the tier.  The identity
    # thresholds are unchanged from earlier releases; they are reused here only
    # to annotate how far the reconstruction sits from the nearest allele.
    if identity < th["low"]["min_identity"]:
        flags.append(
            f"CATALOGUE_HIGHLY_DIVERGENT ({identity:.1f}% to nearest allele)"
        )
    elif identity < th["medium"]["min_identity"]:
        flags.append(f"CATALOGUE_DIVERGENT ({identity:.1f}% to nearest allele)")
    elif identity < th["high"]["min_identity"]:
        flags.append(
            f"CATALOGUE_DIVERGENCE_BELOW_HIGH ({identity:.1f}% to nearest allele)"
        )

    if qcov < th["low"]["min_qcov"]:
        flags.append(f"LOW_COVERAGE ({qcov:.1f}%)")
    elif qcov < th["medium"]["min_qcov"]:
        flags.append(f"REDUCED_COVERAGE ({qcov:.1f}%)")
    elif qcov < th["high"]["min_qcov"]:
        flags.append(f"COVERAGE_BELOW_HIGH ({qcov:.1f}%)")

    if n_fraction > th["low"]["max_n_fraction"]:
        flags.append(f"HIGH_N_CONTENT ({n_fraction*100:.2f}%)")
    elif n_fraction > th["medium"]["max_n_fraction"]:
        flags.append(f"ELEVATED_N_CONTENT ({n_fraction*100:.2f}%)")
    elif n_fraction > th["high"]["max_n_fraction"]:
        flags.append(f"SOME_N_CONTENT ({n_fraction*100:.2f}%)")

    if gc_deviation > th["low"]["max_gc_deviation"]:
        flags.append(f"GC_ANOMALY (delta={gc_deviation:.1f}%)")
    elif gc_deviation > th["medium"]["max_gc_deviation"]:
        flags.append(f"GC_SHIFT (delta={gc_deviation:.1f}%)")
    elif gc_deviation > th["high"]["max_gc_deviation"]:
        flags.append(f"GC_BELOW_HIGH (delta={gc_deviation:.1f}%)")

    if expect_cds:
        if internal_stops > 0:
            flags.append(f"INTERNAL_STOPS (min {internal_stops} in best frame)")
        if frame_disrupted:
            flags.append(
                f"FRAME_LENGTH_SHIFT (length % 3 = {length_mod3}; "
                f"bait mode = {expected_mod3})"
            )

    mean_depth    = remap_metrics.get("mean_depth", 0.0)
    breadth_pct   = remap_metrics.get("breadth_pct", 0.0)
    pct_bases_lt5 = remap_metrics.get("pct_bases_lt5", 100.0)

    if mean_depth < remap_th["low"]["min_mean_depth"]:
        flags.append(f"LOW_REMAP_DEPTH (mean={mean_depth:.1f}x)")
    elif mean_depth < remap_th["medium"]["min_mean_depth"]:
        flags.append(f"REDUCED_REMAP_DEPTH (mean={mean_depth:.1f}x)")
    elif mean_depth < remap_th["high"]["min_mean_depth"]:
        flags.append(f"REMAP_DEPTH_BELOW_HIGH (mean={mean_depth:.1f}x)")

    if breadth_pct < remap_th["low"]["min_breadth"]:
        flags.append(f"LOW_REMAP_BREADTH ({breadth_pct:.1f}%)")
    elif breadth_pct < remap_th["medium"]["min_breadth"]:
        flags.append(f"REDUCED_REMAP_BREADTH ({breadth_pct:.1f}%)")
    elif breadth_pct < remap_th["high"]["min_breadth"]:
        flags.append(f"REMAP_BREADTH_BELOW_HIGH ({breadth_pct:.1f}%)")

    # Coverage patchiness is reported but is descriptive only. On a short
    # extracted-and-remapped locus, pct_bases_lt5 is largely determined by read
    # length, mean depth and locus length rather than by evidence quality, so it
    # does not constrain the tier. Whether the affected bases are actually
    # resolved is decided by the per-base certainty model below.
    if pct_bases_lt5 > remap_th["low"]["max_pct_bases_lt5"]:
        flags.append(f"PATCHY_REMAP_SUPPORT ({pct_bases_lt5:.1f}% <5x, descriptive)")
    elif pct_bases_lt5 > remap_th["medium"]["max_pct_bases_lt5"]:
        flags.append(f"MODEST_REMAP_PATCHINESS ({pct_bases_lt5:.1f}% <5x, descriptive)")
    elif pct_bases_lt5 > remap_th["high"]["max_pct_bases_lt5"]:
        flags.append(
            f"REMAP_PATCHINESS_BELOW_HIGH ({pct_bases_lt5:.1f}% <5x, descriptive)"
        )

    # Per-base certainty: the tier-affecting replacement for patchiness.
    if remap_metrics.get("per_base_available", False):
        unc_th = UNCERTAINTY_THRESHOLDS
        unc_fraction = remap_metrics.get("tier_uncertain_base_fraction", 0.0)
        interior_used = remap_metrics.get("interior_metric_used", True)
        scope = "interior" if interior_used else "full-length"
        if unc_fraction > unc_th["low"]:
            flags.append(
                f"WIDESPREAD_UNCERTAIN_BASES ({unc_fraction*100:.2f}% {scope})"
            )
        elif unc_fraction > unc_th["medium"]:
            flags.append(
                f"ELEVATED_UNCERTAIN_BASES ({unc_fraction*100:.2f}% {scope})"
            )
        elif unc_fraction > unc_th["high"]:
            flags.append(
                f"UNCERTAIN_BASES_BELOW_HIGH ({unc_fraction*100:.2f}% {scope})"
            )

        internal_gap = remap_metrics.get("internal_zero_depth_positions", 0)
        if internal_gap > 0:
            flags.append(
                f"INTERNAL_ZERO_DEPTH_GAP ({internal_gap} interior position(s) "
                f"with no read support)"
            )

    if ambiguity_metrics.get("ambiguous", False):
        margin = ambiguity_metrics.get("bitscore_margin", 0.0)
        flags.append(f"AMBIGUOUS_BEST_HIT (bitscore margin={margin:.1f})")

    if remap_metrics.get("per_base_available", False):
        per_base_th = PER_BASE_THRESHOLDS
        mean_base_quality = remap_metrics.get("mean_base_quality", 0.0)
        mean_mapping_quality = remap_metrics.get("mean_mapping_quality", 0.0)
        strand_balance = remap_metrics.get("strand_balance_pct", 0.0)
        mixed_sites = remap_metrics.get("mixed_site_count", 0)
        candidate_sites = remap_metrics.get("candidate_mixed_sites", 0)
        strand_biased_sites = remap_metrics.get("strand_biased_sites", 0)
        discordant_sites = remap_metrics.get("reference_discordant_sites", 0)

        if mean_base_quality < per_base_th["low"]["min_mean_base_quality"]:
            flags.append(f"LOW_BASE_QUALITY (mean={mean_base_quality:.1f})")
        elif mean_base_quality < per_base_th["medium"]["min_mean_base_quality"]:
            flags.append(f"REDUCED_BASE_QUALITY (mean={mean_base_quality:.1f})")
        elif mean_base_quality < per_base_th["high"]["min_mean_base_quality"]:
            flags.append(f"BASE_QUALITY_BELOW_HIGH (mean={mean_base_quality:.1f})")

        if mean_mapping_quality < per_base_th["low"]["min_mean_mapping_quality"]:
            flags.append(f"LOW_MAPPING_QUALITY (mean={mean_mapping_quality:.1f})")
        elif mean_mapping_quality < per_base_th["medium"]["min_mean_mapping_quality"]:
            flags.append(f"REDUCED_MAPPING_QUALITY (mean={mean_mapping_quality:.1f})")
        elif mean_mapping_quality < per_base_th["high"]["min_mean_mapping_quality"]:
            flags.append(f"MAPPING_QUALITY_BELOW_HIGH (mean={mean_mapping_quality:.1f})")

        if strand_balance < per_base_th["low"]["min_strand_balance"]:
            flags.append(f"SEVERE_STRAND_IMBALANCE (balance={strand_balance:.1f}%)")
        elif strand_balance < per_base_th["medium"]["min_strand_balance"]:
            flags.append(f"STRAND_IMBALANCE (balance={strand_balance:.1f}%)")
        elif strand_balance < per_base_th["high"]["min_strand_balance"]:
            flags.append(f"STRAND_BALANCE_BELOW_HIGH ({strand_balance:.1f}%)")

        if remap_metrics.get("mixture_detected", False):
            fraction = remap_metrics.get("median_mixed_fraction", 0.0)
            flags.append(
                f"ALLELE_MIXTURE ({mixed_sites} bidirectional sites, "
                f"median alt fraction={fraction:.3f})"
            )
        if candidate_sites and strand_biased_sites:
            flags.append(
                f"STRAND_BIASED_ALT_SUPPORT ({strand_biased_sites}/"
                f"{candidate_sites} candidate sites)"
            )
        if discordant_sites:
            flags.append(f"REFERENCE_DISCORDANCE ({discordant_sites} near-fixed sites)")

    # A contig boundary inside the locus truncates the allele.  The per-base
    # model scores whether each reported base is supported and is silent on
    # bases that were never reported, so completeness is asserted here or not
    # at all.
    span_metrics = span_metrics or {}
    clipped_bp = span_metrics.get("clipped_bp", 0)
    if clipped_bp >= SPAN_CLIP_THRESHOLDS["min_flag_bp"]:
        ends = []
        if span_metrics.get("clipped_start_bp"):
            ends.append(f"{span_metrics['clipped_start_bp']} bp at contig start")
        if span_metrics.get("clipped_end_bp"):
            ends.append(f"{span_metrics['clipped_end_bp']} bp at contig end")
        flags.append(
            f"ALLELE_SPAN_CLIPPED_AT_CONTIG_END ({', '.join(ends)}; "
            f"allele is truncated by {clipped_bp} bp)"
        )
        continuation = span_metrics.get("continuation")
        if continuation:
            flags.append(
                f"LOCUS_SPLIT_ACROSS_CONTIGS "
                f"({continuation['recovered_bp']} of the {clipped_bp} lost bp align to "
                f"{continuation['sseqid']}, sharing {continuation['overlap_bp']} bp "
                f"with the reported contig)"
            )

    # Descriptive flags are reported for transparency but are excluded from
    # every tier decision below.
    tier_flags = [f for f in flags if not f.startswith(DESCRIPTIVE_FLAG_PREFIXES)]

    if not tier_flags:
        # HIGH may now carry descriptive coverage notes.
        return "HIGH", flags

    has_length_flag = any(f.startswith("LENGTH_") for f in flags)

    # SUSPECT is reserved for evidence that the sequence may misrepresent the
    # underlying biology -- not for merely thin coverage. Patchy remap support
    # is therefore absent from this set, while INTERNAL_ZERO_DEPTH_GAP is in it,
    # because an unsupported interior position can indicate a chimeric join.
    has_severe_flag = any(
        f.startswith((
            "LOW_COVERAGE", "HIGH_N_CONTENT", "GC_ANOMALY",
            "INTERNAL_STOPS", "LOW_REMAP_DEPTH", "LOW_REMAP_BREADTH",
            "AMBIGUOUS_BEST_HIT", "FRAME_LENGTH_SHIFT",
            "LOW_BASE_QUALITY", "LOW_MAPPING_QUALITY", "SEVERE_STRAND_IMBALANCE",
            "ALLELE_MIXTURE", "STRAND_BIASED_ALT_SUPPORT", "REFERENCE_DISCORDANCE",
            "INTERNAL_ZERO_DEPTH_GAP",
            "LOCUS_SPLIT_ACROSS_CONTIGS",
        ))
        for f in tier_flags
    )

    # A length difference from the bait profile, with no severe read-evidence
    # flag, caps the tier at MEDIUM.  The cap does not depend on catalogue
    # identity: gating it on identity would keep catalogue distance as a tier
    # determinant, only in the upward direction.  The identity threshold is
    # retained solely to label the case as a candidate novel allele.
    if has_length_flag and not has_severe_flag:
        if identity >= th["novel_override_identity"]:
            flags.append(
                f"POSSIBLE_NOVEL_ALLELE (length differs but identity >= "
                f"{th['novel_override_identity']:.0f}% -- may be genuine indel variant)"
            )
        return "MEDIUM", flags

    if has_severe_flag:
        return "SUSPECT", flags

    n_flags = len([f for f in tier_flags if not f.startswith("POSSIBLE_")])
    if n_flags <= 1 and length_within_2x:
        return "MEDIUM", flags
    if n_flags <= 2 and length_within_3x:
        return "LOW", flags
    return "SUSPECT", flags


# ==============================================================================
# QC report formatting
# ==============================================================================


def format_qc_report(
    sample_id: str,
    locus: str,
    qc: dict,
    bait_profile: dict,
) -> str:
    """Format a human-readable QC scorecard for one reconstructed allele.

    Args:
        sample_id:    Sample identifier.
        locus:        Locus name.
        qc:           Dict from assess_allele_quality().
        bait_profile: Dict from profile_bait_database().

    Returns:
        Multi-line string suitable for printing to console or writing to file.
    """
    conf = qc["confidence"]
    tag = {
        "HIGH": "[PASS]",
        "MEDIUM": "[REVIEW]",
        "LOW": "[WARN]",
        "SUSPECT": "[HOLD]",
    }.get(conf, "[????]")

    W   = 64
    bar = "+" + "-" * (W - 2) + "+"

    lines = [
        bar,
        f"|  ALLELE QC REPORT: {sample_id} / {locus}",
        bar,
        f"|  Sequence confidence: {tag} {conf}",
        f"|  Catalogue status   : {qc.get('catalogue_status', 'NOT_ASSESSED')}",
        bar,
        "|  LENGTH",
        f"|    Reconstructed : {qc['allele_length']:>7} bp",
        f"|    Expected (med): {qc['expected_median']:>7.0f} bp",
        f"|    Expected (mod): {qc['expected_mode']:>7} bp",
        f"|    Known range   : {qc['expected_range']:>15}",
        f"|    IQR band      : {qc['expected_iqr_range']:>15}",
    ]
    delta_sign = "+" if qc["length_delta"] > 0 else ""
    lines += [
        f"|    Delta         : {delta_sign}{qc['length_delta']:>7} bp  "
        f"(Z-score={qc['length_zscore']:+.2f})",
        f"|    Within IQR    : {'YES' if qc['length_within_iqr'] else 'NO <<'}",
        bar,
        "|  CATALOGUE RELATIONSHIP (does not set sequence confidence)",
        f"|    Identity      : {qc['validation_identity']:>7.1f}%",
        f"|    Query coverage: {qc['validation_qcov']:>7.1f}%",
        f"|    2nd-hit ident.: {qc['second_best_identity']:>7.1f}%",
        f"|    Bitscore gap  : {qc['bitscore_margin']:>7.1f}",
        bar,
        "|  REMAP SUPPORT TO FINAL ALLELE",
        f"|    Mapped reads  : {qc['remap_mapped_reads']:>7}",
        f"|    Mean depth    : {qc['remap_mean_depth']:>7.1f}x",
        f"|    Min depth     : {qc['remap_min_depth']:>7}x",
        f"|    Max depth     : {qc['remap_max_depth']:>7}x",
        f"|    Breadth       : {qc['remap_breadth_pct']:>7.1f}%",
        f"|    Bases <5x     : {qc['remap_pct_bases_lt5']:>7.1f}%",
        f"|    Bases <10x    : {qc['remap_pct_bases_lt10']:>7.1f}%",
        f"|    (coverage patchiness is descriptive only)",
        bar,
        "|  PER-BASE CERTAINTY  (tier-affecting)",
        f"|    Low-coverage   : {qc['low_coverage_position_count']:>7} positions (<5x)",
        f"|    Uncertain      : {qc['uncertain_base_count']:>7} positions ({qc['uncertain_base_fraction']*100:.2f}%)",
        f"|    Terminal margin: {qc['terminal_margin_bp']:>7} bp (mean read length)",
        f"|    Interior       : {qc['interior_length']:>7} bp"
        + ("" if qc["interior_metric_used"] else "  [too short; full-length used]"),
        f"|    Interior uncert: {qc['interior_uncertain_base_count']:>7} positions ({qc['interior_uncertain_base_fraction']*100:.2f}%)",
        f"|    Interior gaps  : {qc['internal_zero_depth_positions']:>7} positions (0x)",
        bar,
        "|  QUALITY-FILTERED PER-BASE EVIDENCE",
        f"|    Mean base Q   : {qc['mean_base_quality']:>7.1f}",
        f"|    Mean mapping Q: {qc['mean_mapping_quality']:>7.1f}",
        f"|    Strand balance: {qc['strand_balance_pct']:>6.1f}%  (100%=balanced)",
        f"|    Candidate sites: {qc['candidate_mixed_sites']:>5}",
        f"|    Mixed sites   : {qc['mixed_site_count']:>7}  (both strands)",
        f"|    Strand-biased : {qc['strand_biased_sites']:>7}",
        f"|    Max alt fract.: {qc['max_alt_fraction']:>7.3f}",
        f"|    Median mixture: {qc['median_mixed_fraction']:>7.3f}",
        f"|    Mixture call  : {'YES <<' if qc['mixture_detected'] else 'NO':>7}",
        bar,
        "|  SEQUENCE COMPOSITION",
        f"|    GC content    : {qc['gc_pct']:>7.1f}%  "
        f"(expected: {qc['gc_expected']:.1f}% +/- {bait_profile['gc_stdev']:.1f}%)",
        f"|    GC deviation  : {qc['gc_deviation']:>7.1f} pp",
        f"|    N bases       : {qc['n_count']:>7}   ({qc['n_fraction']*100:.3f}%)",
    ]

    if qc["expect_cds"]:
        frame_tag = (
            "(matches bait mode)" if not qc["frame_disrupted"]
            else "(differs from bait mode) <<"
        )
        lines += [
            bar,
            "|  CDS INTEGRITY",
            "|    Bait-supported frame(s): "
            + ", ".join(f"+{frame}" for frame in qc["evaluated_coding_frames"]),
            f"|    Length % 3    : {qc['length_mod3']:>7}   {frame_tag}",
            f"|    Bait mode % 3 : {qc['expected_length_mod3']:>7}",
            f"|    Internal stops: {qc['internal_stops']:>7}   (in inferred frame)",
            "|    Note: MLST definitions may be internal gene fragments;",
            "|          start/stop codons are therefore not required.",
        ]
    else:
        lines += [bar, "|  CDS INTEGRITY", "|    Skipped (--noncoding_locus enabled)"]

    if qc["flags"]:
        lines += [bar, "|  FLAGS"] + [f"|    >> {flag}" for flag in qc["flags"]]

    lines += [bar, "|  INTERPRETATION"]
    if conf == "HIGH":
        lines += [
            "|    Similarity, remap support, and structure all look consistent.",
            "|    -> Strong computational support; curate before submission.",
        ]
    elif conf == "MEDIUM":
        lines += [
            "|    Some deviations were detected, but support is still plausible.",
            "|    -> Review flagged metrics before submission.",
        ]
    elif conf == "LOW":
        lines += [
            "|    Multiple deviations from the expected profile were detected.",
            "|    -> Manual inspection of coverage and local assembly is required.",
        ]
    else:
        lines += [
            "|    Evidence suggests ambiguity, weak remap support, or misassembly.",
            "|    -> Do NOT submit without further investigation.",
        ]
    lines.append(bar)
    return "\n".join(lines)


def write_qc_report_file(report_text: str, output_path: str) -> None:
    """Write the QC report card to a text file.

    Args:
        report_text: Formatted QC report string.
        output_path: Destination file path.
    """
    with open(output_path, "w") as f:
        f.write(report_text)
        f.write("\n")
    log.debug(f"  QC report written to: {output_path}")
