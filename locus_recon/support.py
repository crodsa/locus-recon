"""Quality-filtered per-base read support and mixture detection."""

from __future__ import annotations

import math
import subprocess
from collections import defaultdict
from typing import Dict, Iterable, Tuple

from .utils import BASE_CERTAINTY_THRESHOLDS, PipelineStepError, log

CANONICAL_BASES = {"A", "C", "G", "T"}


def assess_base_certainty(
    depth: int,
    dominant_forward: int,
    dominant_reverse: int,
    alt_fraction: float,
    mean_base_quality: float,
    mean_mapping_quality: float,
    credible_mixed_site: bool,
    min_base_quality: float,
    min_mapping_quality: float,
    thresholds: dict = None,
) -> dict:
    """Decide whether the base reported at one position is adequately supported.

    This is deliberately *not* a depth cutoff. Low filtered depth is compatible
    with a certain base when the available observations are concordant, of
    adequate quality and seen from both strands; conversely a high-depth
    position with conflicting evidence is uncertain. The distinction between
    ``is_low_coverage`` and ``is_uncertain`` is the point of the function.

    Args:
        depth:                Quality-filtered depth at the position.
        dominant_forward:     Forward observations supporting the reported base.
        dominant_reverse:     Reverse observations supporting the reported base.
        alt_fraction:         Fraction of filtered observations supporting the
                              most frequent alternative base.
        mean_base_quality:    Mean base quality of retained observations.
        mean_mapping_quality: Mean mapping quality of retained observations.
        credible_mixed_site:  Whether the site meets the credible mixed-site
                              definition already used for mixture detection.
        min_base_quality:     Configured base-quality intake filter.
        min_mapping_quality:  Configured mapping-quality intake filter.
        thresholds:           Override for BASE_CERTAINTY_THRESHOLDS.

    Returns:
        dict with keys ``uncertain`` (bool), ``reason_codes`` (list of str) and
        ``low_coverage`` (bool).
    """
    th = thresholds or BASE_CERTAINTY_THRESHOLDS
    reasons = []

    low_coverage = depth < th["low_coverage_depth"]

    # A. No usable evidence.
    if depth <= 0:
        return {
            "uncertain": True,
            "reason_codes": ["ZERO_DEPTH"],
            "low_coverage": True,
        }
    if depth < th["min_filtered_depth"]:
        reasons.append("INSUFFICIENT_FILTERED_DEPTH")

    # C. Directionally weak evidence: strand restriction only matters while the
    #    depth is low enough that it materially weakens the call.
    if (depth < th["bidirectional_required_below_depth"]
            and (dominant_forward == 0 or dominant_reverse == 0)):
        reasons.append("NO_BIDIRECTIONAL_SUPPORT")

    # B. Conflicting evidence.
    if alt_fraction >= th["max_alt_fraction"]:
        reasons.append("EXCESS_ALTERNATIVE_SUPPORT")
    if credible_mixed_site:
        reasons.append("CREDIBLE_MIXED_SITE")

    # D. Low-quality evidence. Retained observations already passed the intake
    #    filters, so this fires only when the site mean falls below them.
    if mean_base_quality < min_base_quality:
        reasons.append("LOW_BASE_QUALITY")
    if mean_mapping_quality < min_mapping_quality:
        reasons.append("LOW_MAPPING_QUALITY")

    return {
        "uncertain": bool(reasons),
        "reason_codes": reasons,
        "low_coverage": low_coverage,
    }


def estimate_mean_read_length(
    samtools_path: str, sorted_bam: str, max_reads: int = 20000
) -> int:
    """Return the mean aligned read length in a BAM, or 0 if undeterminable.

    The value defines the terminal margin over which remap depth is
    structurally limited (a read cannot start before the first position of the
    reconstructed allele), so it is a property of the library rather than a
    tunable parameter.
    """
    try:
        completed = subprocess.run(
            [samtools_path, "view", sorted_bam],
            capture_output=True, text=True, check=True,
        )
    except (subprocess.CalledProcessError, OSError):
        return 0
    total = count = 0
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 10 or fields[9] == "*":
            continue
        total += len(fields[9])
        count += 1
        if count >= max_reads:
            break
    return int(round(total / count)) if count else 0


def _empty_strands() -> dict:
    return {"forward": 0, "reverse": 0}


def parse_pileup_bases(read_bases: str, reference_base: str) -> dict:
    """Parse samtools-mpileup base syntax into allele and strand counts.

    Matches (``.`` and ``,``), mismatches, deletion placeholders, read-boundary
    markers, and inline insertion/deletion descriptions are supported. Inline
    indels are represented as ``+SEQUENCE`` or ``-SEQUENCE`` alleles anchored
    at the current reference position.
    """
    reference_base = reference_base.upper()
    counts: Dict[str, dict] = defaultdict(_empty_strands)
    observations = 0
    index = 0
    last_strand = "forward"

    while index < len(read_bases):
        symbol = read_bases[index]
        if symbol == "^":
            index += 2  # Skip '^' and the following encoded mapping quality.
            continue
        if symbol == "$":
            index += 1
            continue
        if symbol in "+-":
            sign = symbol
            index += 1
            length_start = index
            while index < len(read_bases) and read_bases[index].isdigit():
                index += 1
            if length_start == index:
                continue
            event_length = int(read_bases[length_start:index])
            sequence = read_bases[index:index + event_length].upper()
            index += event_length
            event = sign + sequence
            counts[event][last_strand] += 1
            continue

        allele = None
        strand = "forward"
        if symbol == ".":
            allele = reference_base
        elif symbol == ",":
            allele = reference_base
            strand = "reverse"
        elif symbol in "ACGTN":
            allele = symbol
        elif symbol in "acgtn":
            allele = symbol.upper()
            strand = "reverse"
        elif symbol == "*":
            allele = "*"
        elif symbol == "#":
            allele = "*"
            strand = "reverse"
        elif symbol in "<>":
            # Reference skips consume a quality character but are not useful
            # nucleotide evidence for bacterial loci.
            observations += 1
            index += 1
            continue

        if allele is not None:
            counts[allele][strand] += 1
            observations += 1
            last_strand = strand
        index += 1

    return {"counts": dict(counts), "observations": observations}


def _mean_phred(encoded: str) -> float:
    if not encoded or encoded == "*":
        return 0.0
    return sum(max(0, ord(char) - 33) for char in encoded) / len(encoded)


def _log_hypergeom_probability(a: int, b: int, c: int, d: int) -> float:
    row_1 = a + b
    row_2 = c + d
    col_1 = a + c
    total = row_1 + row_2
    return (
        math.lgamma(row_1 + 1)
        + math.lgamma(row_2 + 1)
        + math.lgamma(col_1 + 1)
        + math.lgamma(total - col_1 + 1)
        - math.lgamma(total + 1)
        - math.lgamma(a + 1)
        - math.lgamma(row_1 - a + 1)
        - math.lgamma(col_1 - a + 1)
        - math.lgamma(row_2 - col_1 + a + 1)
    )


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Return a two-sided Fisher exact p-value for a 2x2 count table."""
    if min(a, b, c, d) < 0:
        raise ValueError("Fisher exact counts cannot be negative.")
    total = a + b + c + d
    if total == 0:
        return 1.0

    row_1 = a + b
    row_2 = c + d
    col_1 = a + c
    minimum = max(0, col_1 - row_2)
    maximum = min(row_1, col_1)
    observed_log_p = _log_hypergeom_probability(a, b, c, d)
    probability = 0.0
    for candidate_a in range(minimum, maximum + 1):
        candidate_b = row_1 - candidate_a
        candidate_c = col_1 - candidate_a
        candidate_d = row_2 - candidate_c
        candidate_log_p = _log_hypergeom_probability(
            candidate_a, candidate_b, candidate_c, candidate_d
        )
        if candidate_log_p <= observed_log_p + 1e-12:
            probability += math.exp(candidate_log_p)
    return min(1.0, probability)


def analyze_pileup_lines(
    lines: Iterable[str],
    allele_length: int,
    mixture_min_fraction: float,
    mixture_min_sites: int,
    mixture_min_alt_depth: int,
    strand_bias_pvalue: float = 0.01,
    min_base_quality: float = 20.0,
    min_mapping_quality: float = 20.0,
    terminal_margin: int = 0,
) -> Tuple[dict, list]:
    """Summarize quality-filtered mpileup records and return per-site rows."""
    rows = []
    total_depth = total_forward = total_reverse = 0
    base_quality_sum = mapping_quality_sum = 0.0
    quality_observations = mapping_observations = 0
    covered = bases_lt5 = bases_lt10 = 0
    candidate_mixed_sites = credible_mixed_sites = 0
    strand_biased_sites = reference_discordant_sites = 0
    max_alt_fraction = 0.0
    credible_fractions = []
    low_coverage_positions = 0
    uncertain_positions = 0
    uncertainty_reason_counts = defaultdict(int)

    for raw_line in lines:
        fields = raw_line.rstrip("\n").split("\t")
        if len(fields) < 4:
            continue
        contig, position_text, reference_base, depth_text = fields[:4]
        position = int(position_text)
        reported_depth = int(depth_text)
        read_bases = fields[4] if len(fields) > 4 else ""
        base_qualities = fields[5] if len(fields) > 5 else ""
        mapping_qualities = fields[6] if len(fields) > 6 else ""
        # With ``-aa``, samtools represents uncovered positions as depth 0
        # followed by ``*`` placeholders. They are not deletion observations.
        parsed = (
            parse_pileup_bases(read_bases, reference_base)
            if reported_depth > 0
            else {"counts": {}, "observations": 0}
        )
        counts = parsed["counts"]
        depth = parsed["observations"] if reported_depth else 0

        # mpileup's quality strings have one character per read observation.
        mean_base_quality = _mean_phred(base_qualities)
        mean_mapping_quality = _mean_phred(mapping_qualities)
        if base_qualities and base_qualities != "*":
            base_quality_sum += sum(max(0, ord(char) - 33) for char in base_qualities)
            quality_observations += len(base_qualities)
        if mapping_qualities and mapping_qualities != "*":
            mapping_quality_sum += sum(
                max(0, ord(char) - 33) for char in mapping_qualities
            )
            mapping_observations += len(mapping_qualities)

        ref = reference_base.upper()
        ref_counts = counts.get(ref, _empty_strands())
        nucleotide_counts = [
            value for key, value in counts.items()
            if not key.startswith(("+", "-"))
        ]
        forward = sum(value["forward"] for value in nucleotide_counts)
        reverse = sum(value["reverse"] for value in nucleotide_counts)
        total_forward += forward
        total_reverse += reverse
        total_depth += depth
        if depth > 0:
            covered += 1
        if depth < 5:
            bases_lt5 += 1
        if depth < 10:
            bases_lt10 += 1

        alternatives = {
            allele: strand_counts
            for allele, strand_counts in counts.items()
            if allele != ref and allele != "N"
        }
        if alternatives:
            alt_base, alt_counts = max(
                alternatives.items(),
                key=lambda item: item[1]["forward"] + item[1]["reverse"],
            )
        else:
            alt_base, alt_counts = "", _empty_strands()

        alt_depth = alt_counts["forward"] + alt_counts["reverse"]
        alt_fraction = (alt_depth / depth) if depth else 0.0
        max_alt_fraction = max(max_alt_fraction, alt_fraction)
        strand_balance = (
            200.0 * min(forward, reverse) / (forward + reverse)
            if (forward + reverse)
            else 0.0
        )
        strand_pvalue = fisher_exact_two_sided(
            ref_counts["forward"], ref_counts["reverse"],
            alt_counts["forward"], alt_counts["reverse"],
        ) if alt_depth else 1.0

        candidate = (
            depth > 0
            and alt_depth >= mixture_min_alt_depth
            and mixture_min_fraction <= alt_fraction <= 1.0 - mixture_min_fraction
        )
        strand_biased = candidate and (
            alt_counts["forward"] == 0
            or alt_counts["reverse"] == 0
            or strand_pvalue < strand_bias_pvalue
        )
        credible = candidate and not strand_biased
        reference_discordant = depth > 0 and alt_fraction > 1.0 - mixture_min_fraction

        candidate_mixed_sites += int(candidate)
        credible_mixed_sites += int(credible)
        strand_biased_sites += int(strand_biased)
        reference_discordant_sites += int(reference_discordant)
        if credible:
            credible_fractions.append(alt_fraction)

        certainty = assess_base_certainty(
            depth=depth,
            dominant_forward=ref_counts["forward"],
            dominant_reverse=ref_counts["reverse"],
            alt_fraction=alt_fraction,
            mean_base_quality=mean_base_quality,
            mean_mapping_quality=mean_mapping_quality,
            credible_mixed_site=credible,
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
        )
        low_coverage_positions += int(certainty["low_coverage"])
        uncertain_positions += int(certainty["uncertain"])
        for code in certainty["reason_codes"]:
            uncertainty_reason_counts[code] += 1

        rows.append({
            "contig": contig,
            "position": position,
            "reference_base": ref,
            "depth": depth,
            "forward_depth": forward,
            "reverse_depth": reverse,
            "strand_balance_pct": strand_balance,
            "ref_forward": ref_counts["forward"],
            "ref_reverse": ref_counts["reverse"],
            "alt_base": alt_base,
            "alt_forward": alt_counts["forward"],
            "alt_reverse": alt_counts["reverse"],
            "alt_fraction": alt_fraction,
            "mean_base_quality": mean_base_quality,
            "mean_mapping_quality": mean_mapping_quality,
            "strand_bias_pvalue": strand_pvalue,
            "candidate_mixed_site": candidate,
            "credible_mixed_site": credible,
            "strand_biased_site": strand_biased,
            "reference_discordant_site": reference_discordant,
            "is_low_coverage": certainty["low_coverage"],
            "is_uncertain": certainty["uncertain"],
            "uncertainty_reason": ";".join(certainty["reason_codes"]),
        })

    observed_positions = len(rows)
    missing_positions = max(0, allele_length - observed_positions)
    bases_lt5 += missing_positions
    bases_lt10 += missing_positions
    # Positions absent from the pileup have no evidence at all, so they are
    # both low coverage and uncertain.
    low_coverage_positions += missing_positions
    uncertain_positions += missing_positions
    if missing_positions:
        uncertainty_reason_counts["ZERO_DEPTH"] += missing_positions
    denominator = allele_length or 1

    # Terminal margin: within one read length of either end, remap depth is
    # structurally limited because no read can begin outside the reconstructed
    # allele. Interior positions carry the evidence that is diagnostic of
    # reconstruction quality; both the interior and full-length fractions are
    # reported.
    margin = max(0, int(terminal_margin))
    interior_rows = [
        row for row in rows
        if min(row["position"], allele_length + 1 - row["position"]) > margin
    ]
    interior_length = len(interior_rows)
    interior_uncertain = sum(1 for row in interior_rows if row["is_uncertain"])
    internal_zero_depth = sum(1 for row in interior_rows if row["depth"] <= 0)
    interior_fraction = (
        interior_uncertain / interior_length if interior_length else 0.0
    )
    full_fraction = uncertain_positions / denominator
    # For loci short relative to the read length the interior is too small to
    # carry the decision; fall back to the full-length fraction so that the
    # metric can never become vacuous.
    interior_usable = interior_length >= max(50, int(0.2 * allele_length))
    tier_fraction = interior_fraction if interior_usable else full_fraction
    overall_strand_balance = (
        200.0 * min(total_forward, total_reverse) / (total_forward + total_reverse)
        if total_forward + total_reverse else 0.0
    )
    sorted_fractions = sorted(credible_fractions)
    if sorted_fractions:
        middle = len(sorted_fractions) // 2
        if len(sorted_fractions) % 2:
            median_fraction = sorted_fractions[middle]
        else:
            median_fraction = (
                sorted_fractions[middle - 1] + sorted_fractions[middle]
            ) / 2.0
    else:
        median_fraction = 0.0

    summary = {
        "per_base_available": True,
        "allele_len": allele_length,
        "mean_depth": total_depth / denominator,
        "min_depth": min((row["depth"] for row in rows), default=0) if not missing_positions else 0,
        "max_depth": max((row["depth"] for row in rows), default=0),
        "breadth_pct": 100.0 * covered / denominator,
        "bases_lt5": bases_lt5,
        "bases_lt10": bases_lt10,
        "pct_bases_lt5": 100.0 * bases_lt5 / denominator,
        "pct_bases_lt10": 100.0 * bases_lt10 / denominator,
        # Coverage and certainty are reported separately and never conflated.
        "low_coverage_position_count": low_coverage_positions,
        "low_coverage_position_fraction": low_coverage_positions / denominator,
        "uncertain_base_count": uncertain_positions,
        "uncertain_base_fraction": uncertain_positions / denominator,
        "terminal_margin_bp": margin,
        "interior_length": interior_length,
        "interior_uncertain_base_count": interior_uncertain,
        "interior_uncertain_base_fraction": interior_fraction,
        "interior_metric_used": interior_usable,
        "tier_uncertain_base_fraction": tier_fraction,
        "internal_zero_depth_positions": internal_zero_depth,
        "uncertainty_reasons": dict(uncertainty_reason_counts),
        "mean_base_quality": (
            base_quality_sum / quality_observations if quality_observations else 0.0
        ),
        "mean_mapping_quality": (
            mapping_quality_sum / mapping_observations
            if mapping_observations else 0.0
        ),
        "forward_observations": total_forward,
        "reverse_observations": total_reverse,
        "strand_balance_pct": overall_strand_balance,
        "candidate_mixed_sites": candidate_mixed_sites,
        "mixed_site_count": credible_mixed_sites,
        "strand_biased_sites": strand_biased_sites,
        "reference_discordant_sites": reference_discordant_sites,
        "max_alt_fraction": max_alt_fraction,
        "median_mixed_fraction": median_fraction,
        "mixture_detected": credible_mixed_sites >= mixture_min_sites,
        "mixture_min_fraction": mixture_min_fraction,
        "mixture_min_sites": mixture_min_sites,
        "mixture_min_alt_depth": mixture_min_alt_depth,
    }
    return summary, rows


def write_support_table(rows: list, output_path: str) -> None:
    """Write one transparent evidence row per reconstructed-allele position."""
    columns = [
        "contig", "position", "reference_base", "depth",
        "forward_depth", "reverse_depth", "strand_balance_pct",
        "ref_forward", "ref_reverse", "alt_base", "alt_forward", "alt_reverse",
        "alt_fraction", "mean_base_quality", "mean_mapping_quality",
        "strand_bias_pvalue", "candidate_mixed_site", "credible_mixed_site",
        "strand_biased_site", "reference_discordant_site",
        "is_low_coverage", "is_uncertain", "uncertainty_reason",
    ]
    with open(output_path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            values = []
            for column in columns:
                value = row[column]
                if isinstance(value, bool):
                    values.append(str(value).lower())
                elif isinstance(value, float):
                    values.append(f"{value:.6g}")
                else:
                    values.append(str(value))
            handle.write("\t".join(values) + "\n")


def collect_per_base_support(
    samtools_path: str,
    allele_fasta: str,
    sorted_bam: str,
    output_path: str,
    allele_length: int,
    min_base_quality: int,
    min_mapping_quality: int,
    mixture_min_fraction: float,
    mixture_min_sites: int,
    mixture_min_alt_depth: int,
) -> dict:
    """Run quality-filtered mpileup, write the evidence table, and summarize it."""
    command = [
        samtools_path, "mpileup", "-aa", "-s", "-d", "1000000",
        "-Q", str(min_base_quality), "-q", str(min_mapping_quality),
        "-f", allele_fasta, sorted_bam,
    ]
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        detail = exc.stderr.strip()
        raise PipelineStepError(f"samtools mpileup failed: {detail}") from exc

    mean_read_length = estimate_mean_read_length(samtools_path, sorted_bam)

    summary, rows = analyze_pileup_lines(
        completed.stdout.splitlines(),
        allele_length=allele_length,
        mixture_min_fraction=mixture_min_fraction,
        mixture_min_sites=mixture_min_sites,
        mixture_min_alt_depth=mixture_min_alt_depth,
        min_base_quality=float(min_base_quality),
        min_mapping_quality=float(min_mapping_quality),
        terminal_margin=mean_read_length,
    )
    summary.update({
        "min_base_quality": min_base_quality,
        "min_mapping_quality": min_mapping_quality,
        "mean_read_length": mean_read_length,
        "support_file": output_path,
    })
    write_support_table(rows, output_path)
    log.debug(
        "  Per-base support: mean depth %.1fx, breadth %.1f%%, mixed sites %d",
        summary["mean_depth"], summary["breadth_pct"], summary["mixed_site_count"],
    )
    return summary
