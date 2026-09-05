import pytest

from locus_recon.support import (
    analyze_pileup_lines,
    assess_base_certainty,
    fisher_exact_two_sided,
    parse_pileup_bases,
)


def test_pileup_parser_handles_strands_markers_and_indels():
    parsed = parse_pileup_bases("^F.,Aa$+2tt,-1c", "C")
    counts = parsed["counts"]
    assert parsed["observations"] == 5
    assert counts["C"] == {"forward": 1, "reverse": 2}
    assert counts["A"] == {"forward": 1, "reverse": 1}
    assert counts["+TT"]["reverse"] == 1
    assert counts["-C"]["reverse"] == 1


def test_fisher_exact_matches_reference_example():
    assert fisher_exact_two_sided(1, 9, 11, 3) == pytest.approx(0.00275946, rel=1e-5)
    assert fisher_exact_two_sided(0, 0, 0, 0) == 1.0


def test_quality_filtered_support_detects_bidirectional_mixture():
    qualities = "I" * 10
    mappings = "J" * 10
    lines = [
        f"allele\t1\tA\t10\t.....,,,,,\t{qualities}\t{mappings}",
        f"allele\t2\tA\t10\t....,,GGgg\t{qualities}\t{mappings}",
        f"allele\t3\tC\t10\t....,,TTtt\t{qualities}\t{mappings}",
        f"allele\t4\tG\t10\t.....,,,,,\t{qualities}\t{mappings}",
    ]
    summary, rows = analyze_pileup_lines(
        lines,
        allele_length=4,
        mixture_min_fraction=0.1,
        mixture_min_sites=2,
        mixture_min_alt_depth=2,
    )
    assert summary["mixture_detected"]
    assert summary["mixed_site_count"] == 2
    assert summary["median_mixed_fraction"] == pytest.approx(0.4)
    assert summary["mean_base_quality"] == 40.0
    assert summary["mean_mapping_quality"] == 41.0
    assert rows[1]["alt_forward"] == 2
    assert rows[1]["alt_reverse"] == 2


def test_one_strand_alt_support_is_not_a_credible_mixture_site():
    qualities = "I" * 10
    lines = [f"allele\t1\tA\t10\t....,,,,GG\t{qualities}\t{qualities}"]
    summary, _ = analyze_pileup_lines(
        lines,
        allele_length=1,
        mixture_min_fraction=0.1,
        mixture_min_sites=1,
        mixture_min_alt_depth=2,
    )
    assert summary["candidate_mixed_sites"] == 1
    assert summary["mixed_site_count"] == 0
    assert summary["strand_biased_sites"] == 1
    assert not summary["mixture_detected"]


def test_zero_depth_mpileup_placeholder_is_not_a_strand_observation():
    summary, rows = analyze_pileup_lines(
        ["allele\t1\tA\t0\t*\t*\t*"],
        allele_length=1,
        mixture_min_fraction=0.1,
        mixture_min_sites=1,
        mixture_min_alt_depth=2,
    )
    assert rows[0]["depth"] == 0
    assert rows[0]["forward_depth"] == 0
    assert rows[0]["reverse_depth"] == 0
    assert rows[0]["alt_base"] == ""
    assert summary["forward_observations"] == 0


def test_significant_strand_bias_prevents_a_credible_mixture_call():
    # The alternative is present on both strands, but its strand ratio is
    # incompatible with the reference observations.
    bases = "." * 20 + "," * 20 + "G" * 20 + "g"
    qualities = "I" * len(bases)
    summary, rows = analyze_pileup_lines(
        [f"allele\t1\tA\t61\t{bases}\t{qualities}\t{qualities}"],
        allele_length=1,
        mixture_min_fraction=0.1,
        mixture_min_sites=1,
        mixture_min_alt_depth=2,
    )
    assert rows[0]["strand_bias_pvalue"] < 0.01
    assert rows[0]["strand_biased_site"]
    assert not rows[0]["credible_mixed_site"]
    assert not summary["mixture_detected"]


# ==============================================================================
# Per-base certainty model
#
# These cases fix the intended semantics of the coverage/certainty split: low
# filtered depth is a coverage statement, not an uncertainty verdict. Each case
# is asserted both on the certainty function directly and, where it can be
# expressed as a pileup record, end to end through analyze_pileup_lines().
# ==============================================================================


def _certainty(depth, fwd, rev, alt_fraction=0.0, credible=False,
               base_q=40.0, map_q=41.0):
    return assess_base_certainty(
        depth=depth,
        dominant_forward=fwd,
        dominant_reverse=rev,
        alt_fraction=alt_fraction,
        mean_base_quality=base_q,
        mean_mapping_quality=map_q,
        credible_mixed_site=credible,
        min_base_quality=20.0,
        min_mapping_quality=20.0,
    )


def _one_site(read_bases, reference_base="A", depth=None, allele_length=1):
    n = depth if depth is not None else len(read_bases)
    quals = "I" * n if n else "*"
    maps = "J" * n if n else "*"
    bases = read_bases if n else "*"
    return analyze_pileup_lines(
        [f"allele\t1\t{reference_base}\t{n}\t{bases}\t{quals}\t{maps}"],
        allele_length=allele_length,
        mixture_min_fraction=0.1,
        mixture_min_sites=1,
        mixture_min_alt_depth=2,
    )


def test_case1_low_depth_bidirectional_concordant_base_is_certain():
    """Depth 4 seen from both strands: low coverage, but the base is resolved."""
    verdict = _certainty(depth=4, fwd=2, rev=2)
    assert verdict["low_coverage"] is True
    assert verdict["uncertain"] is False
    assert verdict["reason_codes"] == []

    _, rows = _one_site("..,,")
    assert rows[0]["is_low_coverage"] is True
    assert rows[0]["is_uncertain"] is False
    assert rows[0]["uncertainty_reason"] == ""


def test_case2_low_depth_single_strand_base_is_uncertain():
    """Depth 4 from one strand only: strand-restricted evidence is not enough."""
    verdict = _certainty(depth=4, fwd=4, rev=0)
    assert verdict["low_coverage"] is True
    assert verdict["uncertain"] is True
    assert "NO_BIDIRECTIONAL_SUPPORT" in verdict["reason_codes"]

    _, rows = _one_site("....")
    assert rows[0]["is_uncertain"] is True
    assert "NO_BIDIRECTIONAL_SUPPORT" in rows[0]["uncertainty_reason"]


def test_case3_high_depth_conflicting_base_is_uncertain_but_not_low_coverage():
    """Depth 20 split 50/50: ample coverage, unresolved base."""
    verdict = _certainty(depth=20, fwd=5, rev=5, alt_fraction=0.5, credible=True)
    assert verdict["low_coverage"] is False
    assert verdict["uncertain"] is True
    assert "EXCESS_ALTERNATIVE_SUPPORT" in verdict["reason_codes"]
    assert "CREDIBLE_MIXED_SITE" in verdict["reason_codes"]

    _, rows = _one_site(".....,,,,,GGGGGggggg")
    assert rows[0]["is_low_coverage"] is False
    assert rows[0]["is_uncertain"] is True


def test_case4_high_depth_concordant_base_is_certain():
    verdict = _certainty(depth=20, fwd=10, rev=10)
    assert verdict["low_coverage"] is False
    assert verdict["uncertain"] is False

    _, rows = _one_site("..........,,,,,,,,,,")
    assert rows[0]["is_low_coverage"] is False
    assert rows[0]["is_uncertain"] is False


def test_case5_zero_depth_is_uncertain_and_low_coverage():
    verdict = _certainty(depth=0, fwd=0, rev=0)
    assert verdict["low_coverage"] is True
    assert verdict["uncertain"] is True
    assert verdict["reason_codes"] == ["ZERO_DEPTH"]

    _, rows = _one_site("*", depth=0)
    assert rows[0]["is_uncertain"] is True
    assert rows[0]["uncertainty_reason"] == "ZERO_DEPTH"


def test_case6_low_depth_minor_alternative_is_uncertain():
    """Depth 4 with one alternative observation.

    At this depth a 25% alternative fraction cannot be distinguished from a
    genuine minor allele, so the position is reported uncertain even though the
    site does not qualify as a credible mixed site.
    """
    verdict = _certainty(depth=4, fwd=2, rev=1, alt_fraction=0.25, credible=False)
    assert verdict["uncertain"] is True
    assert verdict["reason_codes"] == ["EXCESS_ALTERNATIVE_SUPPORT"]

    _, rows = _one_site("..,G")
    assert rows[0]["credible_mixed_site"] is False
    assert rows[0]["is_uncertain"] is True
    assert "EXCESS_ALTERNATIVE_SUPPORT" in rows[0]["uncertainty_reason"]


def test_single_filtered_observation_cannot_corroborate_itself():
    verdict = _certainty(depth=1, fwd=1, rev=0)
    assert verdict["uncertain"] is True
    assert "INSUFFICIENT_FILTERED_DEPTH" in verdict["reason_codes"]


def test_terminal_margin_separates_interior_from_edge_uncertainty():
    """Edge positions are excluded from the interior metric, interior are not."""
    quals, maps = "I" * 4, "J" * 4
    deep_q, deep_m = "I" * 20, "J" * 20
    lines = []
    for pos in range(1, 21):
        if pos in (1, 20):
            # single-strand, low depth: uncertain, but terminal
            lines.append(f"allele\t{pos}\tA\t4\t....\t{quals}\t{maps}")
        elif pos == 10:
            lines.append(f"allele\t{pos}\tA\t4\t....\t{quals}\t{maps}")
        else:
            lines.append(
                f"allele\t{pos}\tA\t20\t..........,,,,,,,,,,\t{deep_q}\t{deep_m}"
            )
    summary, _ = analyze_pileup_lines(
        lines, allele_length=20, mixture_min_fraction=0.1,
        mixture_min_sites=2, mixture_min_alt_depth=2, terminal_margin=3,
    )
    assert summary["uncertain_base_count"] == 3
    assert summary["terminal_margin_bp"] == 3
    assert summary["interior_length"] == 14
    assert summary["interior_uncertain_base_count"] == 1


def test_interior_metric_falls_back_when_locus_is_short():
    """A locus shorter than twice the read length has no usable interior."""
    quals, maps = "I" * 4, "J" * 4
    lines = [f"allele\t{p}\tA\t4\t....\t{quals}\t{maps}" for p in range(1, 11)]
    summary, _ = analyze_pileup_lines(
        lines, allele_length=10, mixture_min_fraction=0.1,
        mixture_min_sites=2, mixture_min_alt_depth=2, terminal_margin=150,
    )
    assert summary["interior_length"] == 0
    assert summary["interior_metric_used"] is False
    assert summary["tier_uncertain_base_fraction"] == pytest.approx(1.0)
