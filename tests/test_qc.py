import pytest

from locus_recon.io import classify_result_disposition
from locus_recon.qc import (
    assess_allele_quality,
    classify_catalogue_status,
    profile_bait_database,
)
from locus_recon.utils import DESCRIPTIVE_FLAG_PREFIXES


def _strong_remap():
    return {
        "mean_depth": 30.0,
        "min_depth": 20,
        "max_depth": 45,
        "breadth_pct": 100.0,
        "pct_bases_lt5": 0.0,
        "pct_bases_lt10": 0.0,
        "mapped_reads": 200,
        "mapped_pct": 99.0,
    }


def test_profile_and_high_confidence_qc(tmp_path):
    bait = tmp_path / "bait.fa"
    bait.write_text(
        ">allele_1\nATGAAACCCGGG\n"
        ">allele_2\nATGAAACCCGGG\n"
        ">allele_3\nATGAAACCCGGG\n"
    )
    profile = profile_bait_database(str(bait))
    qc = assess_allele_quality(
        "ATGAAACCCGGG", profile, 100.0, 100.0, "abc",
        remap_metrics=_strong_remap(),
        ambiguity_metrics={"ambiguous": False, "bitscore_margin": 100.0},
    )
    assert profile["len_mod3_mode"] == 0
    assert qc["confidence"] == "HIGH"
    assert qc["length_delta"] == 0


def test_modal_nonzero_length_modulo_is_not_called_frameshift(tmp_path):
    bait = tmp_path / "bait.fa"
    bait.write_text(
        ">a1\nATGAAACCCGGGA\n>a2\nATGAAACCCGGGT\n>a3\nATGAAACCCGGGC\n"
    )
    profile = profile_bait_database(str(bait))
    qc = assess_allele_quality(
        "ATGAAACCCGGGA", profile, 100.0, 100.0, "fragment",
        remap_metrics=_strong_remap(),
        ambiguity_metrics={"ambiguous": False, "bitscore_margin": 100.0},
    )
    assert profile["len_mod3_mode"] == 1
    assert not qc["frame_disrupted"]
    assert not any(flag.startswith("FRAME_LENGTH_SHIFT") for flag in qc["flags"])


def _identical_bait(tmp_path):
    """Three identical catalogue records, so length and GC are unremarkable."""
    bait = tmp_path / "bait.fa"
    bait.write_text(
        ">allele_1\nATGAAACCCGGG\n"
        ">allele_2\nATGAAACCCGGG\n"
        ">allele_3\nATGAAACCCGGG\n"
    )
    return profile_bait_database(str(bait))


def _assess_catalogue(profile, identity, **kwargs):
    remap = _strong_remap()
    remap.update(kwargs.pop("remap", {}))
    ambiguity = {"ambiguous": False, "bitscore_margin": 100.0}
    ambiguity.update(kwargs.pop("ambiguity", {}))
    return assess_allele_quality(
        "ATGAAACCCGGG", profile, identity, 100.0, "abc",
        remap_metrics=remap,
        ambiguity_metrics=ambiguity,
        **kwargs,
    )


def test_case_a_exact_catalogue_match_is_high(tmp_path):
    """Case A: exact catalogue match with strong read support -> HIGH/PASS."""
    qc = _assess_catalogue(_identical_bait(tmp_path), 100.0, exact_known_allele=True)
    assert qc["confidence"] == "HIGH"
    assert qc["sequence_confidence"] == "HIGH"
    assert qc["catalogue_status"] == "EXACT_MATCH"
    assert not qc["catalogue_review_recommended"]
    assert classify_result_disposition("SUCCESS", qc["confidence"]) == "PASS"


@pytest.mark.parametrize(
    "identity,expected_status,expected_flag",
    [
        (98.0, "NONEXACT_MATCH", None),
        (96.0, "DIVERGENT_FROM_REFERENCE", "CATALOGUE_DIVERGENCE_BELOW_HIGH"),
        (90.0, "DIVERGENT_FROM_REFERENCE", "CATALOGUE_DIVERGENT"),
        (80.0, "HIGHLY_DIVERGENT_FROM_REFERENCE", "CATALOGUE_HIGHLY_DIVERGENT"),
    ],
)
def test_case_b_catalogue_divergence_alone_does_not_downgrade(
    tmp_path, identity, expected_status, expected_flag
):
    """Case B: identical read evidence, divergent nearest allele -> still HIGH.

    Divergence from the nearest curated allele is reported as catalogue status
    and must not reduce sequence confidence on its own, however far the nearest
    allele is: at 80% identity the tier is still set by the read evidence.
    """
    profile = _identical_bait(tmp_path)
    qc = _assess_catalogue(profile, identity, exact_known_allele=False)
    assert qc["confidence"] == "HIGH"
    assert qc["catalogue_status"] == expected_status
    assert qc["catalogue_review_recommended"]
    assert qc["nearest_allele_identity"] == identity
    assert classify_result_disposition("SUCCESS", qc["confidence"]) == "PASS"
    if expected_flag:
        assert any(f.startswith(expected_flag) for f in qc["flags"])
    assert not any(f.startswith("IDENTITY_") for f in qc["flags"])
    assert not any(f.startswith("REDUCED_IDENTITY") for f in qc["flags"])


def test_case_b_catalogue_flags_are_never_tier_flags(tmp_path):
    """Catalogue annotations must be excluded from the tier-deciding flag set."""
    qc = _assess_catalogue(_identical_bait(tmp_path), 80.0, exact_known_allele=False)
    catalogue_flags = [f for f in qc["flags"] if f.startswith("CATALOGUE_")]
    assert catalogue_flags
    tier_flags = [
        f for f in qc["flags"] if not f.startswith(DESCRIPTIVE_FLAG_PREFIXES)
    ]
    assert not [f for f in tier_flags if f.startswith("CATALOGUE_")]


def test_case_c_divergence_plus_ambiguity_is_suspect_for_ambiguity(tmp_path):
    """Case C: divergent AND ambiguous -> SUSPECT, attributed to ambiguity."""
    profile = _identical_bait(tmp_path)
    qc = _assess_catalogue(
        profile, 90.0, exact_known_allele=False,
        ambiguity={"ambiguous": True, "bitscore_margin": 1.0},
    )
    assert qc["confidence"] == "SUSPECT"
    assert any(f.startswith("AMBIGUOUS_BEST_HIT") for f in qc["flags"])
    assert classify_result_disposition("SUCCESS", qc["confidence"]) == "HOLD"


def test_case_d_divergence_plus_uncertain_bases_downgrades_on_read_evidence(
    tmp_path,
):
    """Case D: divergent AND unresolved bases -> downgraded on read evidence."""
    profile = _identical_bait(tmp_path)
    qc = _assess_catalogue(
        profile, 90.0, exact_known_allele=False,
        remap={
            "per_base_available": True,
            "tier_uncertain_base_fraction": 0.20,
            "interior_uncertain_base_fraction": 0.20,
            "interior_metric_used": True,
        },
    )
    assert qc["confidence"] != "HIGH"
    assert any(f.startswith("WIDESPREAD_UNCERTAIN_BASES") for f in qc["flags"])
    assert qc["catalogue_status"] == "DIVERGENT_FROM_REFERENCE"


def test_length_deviation_caps_at_medium_independently_of_identity(tmp_path):
    """The length cap must not depend on catalogue identity in either direction.

    Pre-correction, identity >= 97% rescued a length-flagged reconstruction to
    MEDIUM while the same evidence at lower identity fell through to LOW or
    SUSPECT.  Both now resolve to MEDIUM.
    """
    bait = tmp_path / "bait.fa"
    bait.write_text(
        ">a1\nATGAAACCCGGG\n>a2\nATGAAACCCGGG\n>a3\nATGAAACCCGGG\n"
    )
    profile = profile_bait_database(str(bait))
    tiers = set()
    for identity in (99.0, 90.0):
        qc = assess_allele_quality(
            "ATGAAACCCGGGATGAAACCCGGGATGAAACCCGGG", profile, identity, 100.0,
            "abc",
            remap_metrics=_strong_remap(),
            ambiguity_metrics={"ambiguous": False, "bitscore_margin": 100.0},
            exact_known_allele=False,
        )
        assert any(f.startswith("LENGTH_") for f in qc["flags"])
        tiers.add(qc["confidence"])
    assert tiers == {"MEDIUM"}


def test_profile_rejects_invalid_bait_symbols(tmp_path):
    bait = tmp_path / "bait.fa"
    bait.write_text(">bad\nACGTZ\n")
    with pytest.raises(ValueError, match="invalid DNA symbols"):
        profile_bait_database(str(bait))


def test_bidirectional_allele_mixture_is_suspect(tmp_path):
    bait = tmp_path / "bait.fa"
    bait.write_text(">a1\nATGAAACCCGGG\n>a2\nATGAAACCCGGG\n")
    profile = profile_bait_database(str(bait))
    remap = {
        **_strong_remap(),
        "per_base_available": True,
        "mean_base_quality": 38.0,
        "mean_mapping_quality": 55.0,
        "strand_balance_pct": 95.0,
        "candidate_mixed_sites": 4,
        "mixed_site_count": 4,
        "strand_biased_sites": 0,
        "reference_discordant_sites": 0,
        "max_alt_fraction": 0.22,
        "median_mixed_fraction": 0.20,
        "mixture_detected": True,
    }
    qc = assess_allele_quality(
        "ATGAAACCCGGG", profile, 100.0, 100.0, "abc",
        remap_metrics=remap,
        ambiguity_metrics={"ambiguous": False, "bitscore_margin": 100.0},
    )
    assert qc["confidence"] == "SUSPECT"
    assert any(flag.startswith("ALLELE_MIXTURE") for flag in qc["flags"])


# ==============================================================================
# Confidence tier and coverage
#
# fraction_below_5x is not a veto on HIGH; the tier-affecting evidence statistic
# is the interior uncertain-base fraction.
# ==============================================================================


def _bait(tmp_path, sequence="ATGAAACCCGGG"):
    bait = tmp_path / "bait.fa"
    bait.write_text(f">a1\n{sequence}\n>a2\n{sequence}\n>a3\n{sequence}\n")
    return profile_bait_database(str(bait))


def _remap_with_certainty(pct_bases_lt5=0.0, interior_fraction=0.0,
                          internal_gap=0, mean_depth=30.0, breadth=100.0):
    return {
        "mean_depth": mean_depth,
        "min_depth": 20,
        "max_depth": 45,
        "breadth_pct": breadth,
        "pct_bases_lt5": pct_bases_lt5,
        "pct_bases_lt10": 0.0,
        "mapped_reads": 200,
        "mapped_pct": 99.0,
        "per_base_available": True,
        "mean_base_quality": 38.0,
        "mean_mapping_quality": 42.0,
        "strand_balance_pct": 98.0,
        "tier_uncertain_base_fraction": interior_fraction,
        "interior_uncertain_base_fraction": interior_fraction,
        "interior_metric_used": True,
        "internal_zero_depth_positions": internal_gap,
    }


def _assess(tmp_path, remap, sequence="ATGAAACCCGGG"):
    return assess_allele_quality(
        sequence, _bait(tmp_path, sequence), 100.0, 100.0, "abc",
        remap_metrics=remap,
        ambiguity_metrics={"ambiguous": False, "bitscore_margin": 100.0},
    )


def test_patchy_coverage_alone_no_longer_vetoes_high(tmp_path):
    """17% of bases below 5x, but every interior base resolved -> HIGH."""
    qc = _assess(tmp_path, _remap_with_certainty(pct_bases_lt5=17.2))
    assert qc["confidence"] == "HIGH"
    # The observation is still reported, marked descriptive.
    assert any(f.startswith("PATCHY_REMAP_SUPPORT") for f in qc["flags"])
    assert any("descriptive" in f for f in qc["flags"])


def test_descriptive_patchiness_is_no_longer_a_severe_flag(tmp_path):
    """PATCHY_REMAP_SUPPORT must not by itself produce SUSPECT."""
    qc = _assess(tmp_path, _remap_with_certainty(pct_bases_lt5=40.0))
    assert qc["confidence"] != "SUSPECT"


def test_interior_uncertainty_is_graded_and_tier_affecting(tmp_path):
    below_high = _assess(tmp_path, _remap_with_certainty(interior_fraction=0.03))
    elevated = _assess(tmp_path, _remap_with_certainty(interior_fraction=0.08))
    widespread = _assess(tmp_path, _remap_with_certainty(interior_fraction=0.30))

    assert any(f.startswith("UNCERTAIN_BASES_BELOW_HIGH") for f in below_high["flags"])
    assert any(f.startswith("ELEVATED_UNCERTAIN_BASES") for f in elevated["flags"])
    assert any(f.startswith("WIDESPREAD_UNCERTAIN_BASES") for f in widespread["flags"])
    assert below_high["confidence"] != "HIGH"


def test_uncertainty_within_budget_still_reaches_high(tmp_path):
    qc = _assess(tmp_path, _remap_with_certainty(interior_fraction=0.005))
    assert qc["confidence"] == "HIGH"
    assert not any("UNCERTAIN_BASES" in f for f in qc["flags"])


def test_internal_zero_depth_gap_is_severe(tmp_path):
    """An unsupported interior position can indicate a chimeric join."""
    qc = _assess(tmp_path, _remap_with_certainty(internal_gap=4))
    assert any(f.startswith("INTERNAL_ZERO_DEPTH_GAP") for f in qc["flags"])
    assert qc["confidence"] == "SUSPECT"


def test_low_remap_depth_remains_tier_affecting(tmp_path):
    """Recalibration must not weaken the genuine depth advisories."""
    qc = _assess(tmp_path, _remap_with_certainty(mean_depth=4.0))
    assert any(f.startswith("LOW_REMAP_DEPTH") for f in qc["flags"])
    assert qc["confidence"] == "SUSPECT"


def _span_bait(tmp_path):
    bait = tmp_path / "bait.fa"
    bait.write_text(
        ">allele_1\nATGAAACCCGGG\n"
        ">allele_2\nATGAAACCCGGG\n"
        ">allele_3\nATGAAACCCGGG\n"
    )
    return profile_bait_database(str(bait))


def _assess_span(tmp_path, span_metrics):
    return assess_allele_quality(
        "ATGAAACCCGGG", _span_bait(tmp_path), 100.0, 100.0, "abc",
        remap_metrics=_strong_remap(),
        ambiguity_metrics={"ambiguous": False, "bitscore_margin": 100.0},
        span_metrics=span_metrics,
    )


def test_unclipped_span_leaves_high_confidence_intact(tmp_path):
    qc = _assess_span(tmp_path, {"clipped_bp": 0, "clipped_start_bp": 0,
                            "clipped_end_bp": 0, "continuation": None})
    assert qc["confidence"] == "HIGH"
    assert not any("CLIPPED" in f or "SPLIT" in f for f in qc["flags"])


def test_omitted_span_metrics_behave_as_unclipped(tmp_path):
    # Backward compatibility: callers predating the span layer must be unaffected.
    assert _assess_span(tmp_path, None)["confidence"] == "HIGH"


def test_sub_threshold_clip_is_treated_as_alignment_jitter(tmp_path):
    qc = _assess_span(tmp_path, {"clipped_bp": 3, "clipped_start_bp": 0,
                            "clipped_end_bp": 3, "continuation": None})
    assert qc["confidence"] == "HIGH"
    assert not any("CLIPPED" in f for f in qc["flags"])


def test_contig_end_clip_blocks_high_confidence(tmp_path):
    # An otherwise flawless allele that is missing sequence is not HIGH: every
    # per-base metric describes only the bases that were reported.
    qc = _assess_span(tmp_path, {"clipped_bp": 317, "clipped_start_bp": 0,
                            "clipped_end_bp": 317, "continuation": None})
    assert qc["confidence"] != "HIGH"
    flag = next(f for f in qc["flags"] if f.startswith("ALLELE_SPAN_CLIPPED_AT_CONTIG_END"))
    assert "317 bp at contig end" in flag
    assert not any(f.startswith("LOCUS_SPLIT_ACROSS_CONTIGS") for f in qc["flags"])


def test_split_locus_is_escalated_to_suspect(tmp_path):
    # The lost sequence is demonstrably present in the local assembly, so the
    # reported allele is an arbitrary fragment of a locus the data resolves.
    qc = _assess_span(tmp_path, {
        "clipped_bp": 317, "clipped_start_bp": 0, "clipped_end_bp": 317,
        "continuation": {"sseqid": "NODE_2_length_784_cov_947.591231",
                         "recovered_bp": 317, "overlap_bp": 77},
    })
    assert qc["confidence"] == "SUSPECT"
    flag = next(f for f in qc["flags"] if f.startswith("LOCUS_SPLIT_ACROSS_CONTIGS"))
    assert "NODE_2_length_784_cov_947.591231" in flag
    assert "77 bp" in flag


def test_both_contig_ends_clipped_are_summed_and_named(tmp_path):
    qc = _assess_span(tmp_path, {"clipped_bp": 60, "clipped_start_bp": 25,
                            "clipped_end_bp": 35, "continuation": None})
    flag = next(f for f in qc["flags"] if f.startswith("ALLELE_SPAN_CLIPPED_AT_CONTIG_END"))
    assert "25 bp at contig start" in flag and "35 bp at contig end" in flag
    assert "truncated by 60 bp" in flag
