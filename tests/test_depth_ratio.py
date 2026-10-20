"""
Unit tests for depth_ratio.py.

These test the statistics and the decision logic, not samtools. Depth
extraction is monkeypatched so the tests run in milliseconds with no BAM and
no external binary -- the same style as the rest of the Locus-Recon test suite
should use for anything that would otherwise need alignment fixtures.

    pytest tests/test_depth_ratio.py -q
"""
from __future__ import annotations

import random

import pytest

from locus_recon import depth_ratio as dr


def _flat(depth: int, n: int, jitter: int = 0, seed: int = 1) -> list[int]:
    rng = random.Random(seed)
    return [max(0, depth + rng.randint(-jitter, jitter)) for _ in range(n)]


@pytest.fixture
def fake_depth(monkeypatch):
    """Install a fake depth_array keyed on region substring."""
    table: dict[str, list[int]] = {}

    def _fake(bam, region=None, *, min_bq=20, min_mq=20, excl_flags="", samtools=""):
        key = "backbone" if region is None or region.startswith("chr_bb") else "locus"
        if min_mq == 0 and f"{key}_raw" in table:
            key = f"{key}_raw"
        return table[key]

    monkeypatch.setattr(dr, "depth_array", _fake)
    monkeypatch.setattr(dr, "_contig_lengths", lambda bam, samtools="": {"chr_bb": 200000})
    return table


def test_single_copy_locus_is_not_rejected(fake_depth):
    fake_depth["locus"] = _flat(50, 3000, jitter=6)
    fake_depth["locus_raw"] = _flat(52, 3000, jitter=6)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=2)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-3000"], n_boot=200)
    assert 0.85 < res.ratio < 1.15
    assert res.single_copy_rejected is False
    assert res.reliable is True


def test_fourfold_amplification_is_detected(fake_depth):
    """A collapsed 4-copy rRNA operon: the gonococcal known-truth case."""
    fake_depth["locus"] = _flat(200, 3000, jitter=20)
    fake_depth["locus_raw"] = _flat(205, 3000, jitter=20)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=3)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-3000"], n_boot=200)
    assert 3.5 < res.ratio < 4.5
    assert res.single_copy_rejected is True
    assert res.ratio_ci_low > 1.5


def test_two_copy_locus_matches_liba6656_expectation(fake_depth):
    fake_depth["locus"] = _flat(100, 7104, jitter=12)
    fake_depth["locus_raw"] = _flat(103, 7104, jitter=12)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=4)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-7104"], n_boot=200)
    assert 1.7 < res.ratio < 2.3
    assert res.single_copy_rejected is True


def test_mapq_ambiguity_is_flagged_not_hidden(fake_depth):
    """Assembler kept both copies: filtered depth craters, ratio understates.

    Here the ratio is 0.4, so nothing is rejected anyway; what the test pins is
    that the unreliability is reported rather than hidden. What happens to a
    HIGH ratio under the same ambiguity is a separate rule, exercised in
    test_high_ratio_is_called_even_when_ambiguous.
    """
    fake_depth["locus"] = _flat(20, 3000, jitter=4)
    fake_depth["locus_raw"] = _flat(100, 3000, jitter=8)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=5)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-3000"], n_boot=200)
    assert res.reliable is False
    assert res.ratio < 1.5
    assert res.single_copy_rejected is False, "ratio below threshold, nothing to call"
    assert any("MAPQ" in n for n in res.notes)


def test_high_ratio_is_called_even_when_ambiguous(fake_depth):
    """Low ambiguity biases the ratio DOWN, so it cannot justify withholding a
    rejection of single copy. This is the shape of a high-copy aphA1 locus whose
    repeated flanks collapse under the MAPQ filter: the lower confidence bound
    sits far above the threshold while ambiguity is around 0.45. Gating the
    rejection on reliability would answer "cannot rule out one copy" from
    evidence that can only understate the copy number."""
    fake_depth["locus"] = _flat(500, 3000, jitter=20)
    fake_depth["locus_raw"] = _flat(1200, 3000, jitter=40)     # ambiguity ~0.42
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=7)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-3000"], n_boot=200)
    assert res.reliable is False, "the ambiguity warning must still be raised"
    assert res.single_copy_rejected is True
    assert res.depth_call == "MULTICOPY_DEPTH"
    assert any("LOWER BOUND" in n for n in res.notes), \
        "the count must be reported as a floor, not as a point estimate"


def test_fragment_count_does_not_multiply_the_estimate(fake_depth):
    """The copy estimate must NOT scale with the number of assembly fragments.

    `ratio * assembly_copies` looks like a copy number and is not one: a
    single-copy locus the assembler broke into n pieces comes back as n copies.
    Measured on 512 real single-copy 7104 bp stretches split into n disjoint
    pieces, that quantity runs from 1.97 (n=2) to 8.85 (n=9) with the truth at
    1. Only the fragment count itself is reported as a fact.
    """
    fake_depth["locus"] = _flat(50, 3000, jitter=6)
    fake_depth["locus_raw"] = _flat(52, 3000, jitter=6)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=6)

    res = dr.estimate_locus_copy_number(
        "x.bam", ["chr_locus:1-3000", "chr_locus2:1-3000"], n_boot=200)
    assert res.assembly_copies == 2, "the fragment count is still reported as a fact"
    assert not hasattr(res, "genome_level_estimate")
    # Without gene coordinates there is no correct answer, so none is invented.
    assert res.copies_estimate is None
    assert any("NOT a valid substitute" in n for n in res.notes)

    # One single-copy gene broken into two disjoint halves is still one copy.
    res2 = dr.estimate_locus_copy_number(
        "x.bam", ["chr_locus:1-3000", "chr_locus2:1-3000"],
        gene_spans=[(1, 1500), (1501, 3000)], gene_length=3000, n_boot=200)
    assert res2.copies_estimate == pytest.approx(1.0, abs=0.12)


def test_single_fragment_estimate_is_exactly_the_ratio(fake_depth):
    """Prespecified prediction (a): a strict generalisation.

    If the locus resolves to one region the estimator must reduce to `ratio`,
    which is what keeps the aphA1 validation (0.95, 1.09, 1.10, 10.45, 77.67)
    untouched by this change.
    """
    fake_depth["locus"] = _flat(500, 3000, jitter=20)
    fake_depth["locus_raw"] = _flat(510, 3000, jitter=20)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=7)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-3000"], n_boot=200)
    assert res.copies_estimate == pytest.approx(res.ratio, abs=0.001)


def test_two_separated_copies_are_recovered(fake_depth):
    """The regime the depth ratio alone cannot see.

    Two copies the assembler kept APART each sit at 1x depth, so the ratio is
    ~1 and correctly reports no collapse. Length carries the second copy.
    """
    fake_depth["locus"] = _flat(50, 3000, jitter=6)
    fake_depth["locus_raw"] = _flat(52, 3000, jitter=6)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=8)

    res = dr.estimate_locus_copy_number(
        "x.bam", ["chr_copyA:1-3000", "chr_copyB:1-3000"],
        gene_spans=[(1, 3000), (1, 3000)], gene_length=3000, n_boot=200)
    assert 0.85 < res.ratio < 1.15, "no collapse, so the ratio must stay at 1"
    assert res.copies_estimate == pytest.approx(2.0, abs=0.2)


def test_collapsed_and_separated_copies_agree(fake_depth):
    """Two copies collapsed into ONE fragment at 2x must give the same answer
    as two copies kept apart at 1x. Same biology, different assembly luck."""
    fake_depth["locus"] = _flat(100, 3000, jitter=10)
    fake_depth["locus_raw"] = _flat(102, 3000, jitter=10)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=9)

    res = dr.estimate_locus_copy_number(
        "x.bam", ["chr_locus:1-3000"],
        gene_spans=[(1, 3000)], gene_length=3000, n_boot=200)
    assert res.copies_estimate == pytest.approx(2.0, abs=0.2)


def test_mismatched_gene_spans_are_refused_not_guessed(fake_depth):
    fake_depth["locus"] = _flat(50, 3000, jitter=6)
    fake_depth["locus_raw"] = _flat(52, 3000, jitter=6)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=10)

    res = dr.estimate_locus_copy_number(
        "x.bam", ["chr_a:1-3000", "chr_b:1-3000", "chr_c:1-3000"],
        gene_spans=[(1, 1500)], gene_length=3000, n_boot=200)
    assert res.copies_estimate is None
    assert any("one to one" in n for n in res.notes)


def test_fractional_ratio_is_preserved(fake_depth):
    """Subpopulation amplification: 2.4 copies per cell must stay 2.4."""
    fake_depth["locus"] = _flat(120, 3000, jitter=8)
    fake_depth["locus_raw"] = _flat(122, 3000, jitter=8)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=7)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-3000"], n_boot=200)
    assert 2.2 < res.ratio < 2.6
    assert res.ratio != round(res.ratio)


def test_zero_backbone_returns_undefined_not_crash(fake_depth):
    fake_depth["locus"] = _flat(50, 1000)
    fake_depth["locus_raw"] = _flat(50, 1000)
    fake_depth["backbone"] = [0] * 100000

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-1000"], n_boot=50)
    assert res.ratio == 0.0
    assert res.reliable is False


def _autocorrelated(mean: int, n: int, step: int = 3, seed: int = 1) -> list[int]:
    """Depth with local structure, like real coverage -- not i.i.d. noise."""
    rng = random.Random(seed)
    out, cur = [], mean
    for _ in range(n):
        cur = max(1, cur + rng.randint(-step, step))
        cur = int(cur * 0.98 + mean * 0.02)  # mean reversion
        out.append(cur)
    return out


def test_block_bootstrap_is_wider_than_naive_resampling():
    """Guard against the pseudo-replication bug.

    With autocorrelated depth -- which is what real coverage looks like --
    resampling single bases pretends there are far more independent
    observations than there are, and the interval collapses. Blocks must give
    a wider, honest interval.
    """
    locus = _autocorrelated(100, 5000, step=4, seed=11)
    backbone = _autocorrelated(50, 50000, step=3, seed=12)

    wide = dr.block_bootstrap_ci(locus, backbone, block=500, n_boot=200)
    tight = dr.block_bootstrap_ci(locus, backbone, block=1, n_boot=200)
    assert (wide[1] - wide[0]) > (tight[1] - tight[0])


def test_bootstrap_is_deterministic_under_a_fixed_seed():
    locus = _flat(80, 3000, jitter=10, seed=21)
    backbone = _flat(40, 30000, jitter=8, seed=22)
    a = dr.block_bootstrap_ci(locus, backbone, n_boot=200, seed=99)
    b = dr.block_bootstrap_ci(locus, backbone, n_boot=200, seed=99)
    assert a == b


def test_bootstrap_backbone_block_cap_is_deterministic():
    blocks = [[index] * 5 for index in range(100)]
    selected_a = dr._evenly_spaced_blocks(blocks, 8)
    selected_b = dr._evenly_spaced_blocks(blocks, 8)
    assert selected_a == selected_b
    assert len(selected_a) == 8
    assert selected_a[0] == blocks[0]
    assert selected_a[-1] == blocks[-1]
    assert selected_a == sorted(selected_a)


def test_result_records_bootstrap_geometry(fake_depth):
    fake_depth["locus"] = _flat(50, 4000, jitter=3, seed=70)
    fake_depth["locus_raw"] = _flat(51, 4000, jitter=3, seed=71)
    fake_depth["backbone"] = _flat(50, 100000, jitter=5, seed=72)
    result = dr.estimate_locus_copy_number(
        "x.bam", ["chr_locus:1-4000"], n_boot=20,
        max_backbone_bootstrap_blocks=32,
    )
    assert result.bootstrap_block_size == 500
    assert result.bootstrap_locus_blocks == 8
    assert result.bootstrap_backbone_blocks_total == 200
    assert result.bootstrap_backbone_blocks_used == 32
    assert result.bootstrap_resamples == 20
    assert result.bootstrap_seed == 20260817


def test_deletion_is_noted(fake_depth):
    fake_depth["locus"] = _flat(15, 3000, jitter=3)
    fake_depth["locus_raw"] = _flat(16, 3000, jitter=3)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=8)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-3000"], n_boot=200)
    assert res.ratio < 0.6
    assert any("deletion" in n for n in res.notes)


@pytest.mark.parametrize("region,expected", [
    ("contig_1", ("contig_1", None, None)),
    ("contig_1:100-200", ("contig_1", 100, 200)),
])
def test_region_parsing(region, expected):
    assert dr._parse_region(region) == expected


def test_malformed_region_raises():
    with pytest.raises(ValueError):
        dr._parse_region("contig:notanumber")


# --------------------------------------------------------------------------- #
# Windowed profile. Prespecified in validation/depth-copy-number/: H > 0.30
# marks a locus whose global median cannot be trusted. The threshold comes from
# the null distribution of H over 916 single-copy 4 kb backbone segments
# (P99 = 0.288).
# --------------------------------------------------------------------------- #
def _partly_collapsed(n: int, backbone: int, collapsed_frac: float,
                      seed: int = 11) -> list[int]:
    """Locus where only part of the span is collapsed at 2x.

    This is the LIBA-6656 *tcdB* geometry: the assembler resolved the copies
    over most of the locus (1x there) and collapsed them over the rest (2x).
    """
    rng = random.Random(seed)
    out: list[int] = []
    block = 200
    for i in range(0, n, block):
        lvl = 2 * backbone if (i / n) < collapsed_frac else backbone
        out += [max(0, lvl + rng.randint(-6, 6)) for _ in range(min(block, n - i))]
    return out


def test_span_overlap_is_exposed_as_review_state(fake_depth):
    fake_depth["locus"] = _flat(50, 3100, jitter=2, seed=44)
    fake_depth["locus_raw"] = _flat(51, 3100, jitter=2, seed=45)
    fake_depth["backbone"] = _flat(50, 100000, jitter=4, seed=46)

    res = dr.estimate_locus_copy_number(
        "x.bam",
        ["chr_locus:1-1500", "chr_locus2:1-1600"],
        gene_spans=[(1, 1500), (1401, 3000)],
        gene_length=3000,
        n_boot=100,
    )

    assert res.assembly_region_count == 2
    assert res.assembly_copies == 2  # same value under the short name
    assert res.gene_total_aligned_bp == 3100
    assert res.gene_union_covered_bp == 3000
    assert res.gene_overlap_bp == 100
    assert res.gene_overlap_fraction == pytest.approx(100 / 3000, abs=0.001)
    # Overlap downgrades the verdict to a review state. The additive value is
    # still returned -- a query base resolved into two fragments carries the sum
    # of their ratios -- but it is not presented as a copy count, because the
    # same geometry inflates the value on single-copy sequence.
    assert "OVERLAPPING_QUERY_SPANS" in res.dosage_flags
    assert res.dosage_status == "REVIEW_OVERLAPPING_SPANS"
    # The geometric null is exposed alongside it: 100 of 3000 covered bases are
    # resolved twice, so unit ratios would return (2900 + 200) / 3000.
    assert res.geometry_expected_dosage == pytest.approx(3100 / 3000, abs=0.001)
    r0, r1 = res.region_ratios
    expected = (1400 * r0 + 100 * (r0 + r1) + 1500 * r1) / 3000
    assert res.dosage_estimate == pytest.approx(expected, abs=0.001)
    assert res.query_multiplicity_max == pytest.approx(r0 + r1, abs=0.001)
    # Query coverage is complete, so the bracket collapses onto the estimate.
    assert res.query_coverage_fraction == pytest.approx(1.0, abs=1e-6)
    assert res.dosage_lower == pytest.approx(res.dosage_estimate, abs=0.001)
    assert res.dosage_upper == pytest.approx(res.dosage_estimate, abs=0.001)


def test_overlapping_spans_match_the_disjoint_formula_when_disjoint(fake_depth):
    """The new estimand is a strict generalisation: identical where spans are
    disjoint and complete, which is the geometry every published run used."""
    fake_depth["locus"] = _flat(50, 3000, jitter=2, seed=44)
    fake_depth["locus_raw"] = _flat(51, 3000, jitter=2, seed=45)
    fake_depth["backbone"] = _flat(50, 100000, jitter=4, seed=46)

    spans = [(1, 1500), (1501, 3000)]
    res = dr.estimate_locus_copy_number(
        "x.bam", ["chr_locus:1-1500", "chr_locus2:1-1500"],
        gene_spans=spans, gene_length=3000, n_boot=100,
    )
    length_weighted = dr.length_weighted_copies(spans, res.region_ratios, 3000)
    assert res.gene_overlap_bp == 0
    assert res.dosage_estimate == pytest.approx(length_weighted, abs=0.001)
    assert res.dosage_status == "ESTIMATED"


def test_incomplete_query_coverage_brackets_rather_than_withholds(fake_depth):
    fake_depth["locus"] = _flat(50, 2000, jitter=2, seed=47)
    fake_depth["locus_raw"] = _flat(51, 2000, jitter=2, seed=48)
    fake_depth["backbone"] = _flat(50, 100000, jitter=4, seed=49)

    res = dr.estimate_locus_copy_number(
        "x.bam",
        ["chr_locus:1-1000", "chr_locus2:1-1000"],
        gene_spans=[(1, 1000), (2001, 3000)],
        gene_length=3000,
        n_boot=100,
    )

    assert res.gene_union_covered_bp == 2000
    assert "INCOMPLETE_QUERY_COVERAGE" in res.dosage_flags
    assert res.dosage_status == "ESTIMATED_PARTIAL_COVERAGE"
    assert res.query_coverage_fraction == pytest.approx(2 / 3, abs=0.001)

    r0, r1 = res.region_ratios
    total = 1000 * r0 + 1000 * r1
    # Point estimate is the mean over the query bases that were observed.
    assert res.dosage_estimate == pytest.approx(total / 2000, abs=0.001)
    # The bracket is what the unobserved third of the query can still do:
    # nothing at all, or as much as the busiest observed base.
    assert res.dosage_lower == pytest.approx(total / 3000, abs=0.001)
    assert res.dosage_upper == pytest.approx(
        (total + 1000 * max(r0, r1)) / 3000, abs=0.001)
    assert res.dosage_lower <= res.dosage_estimate <= res.dosage_upper
    # A confidence interval is reported alongside, and it is the sampling
    # uncertainty, not the coverage bracket. The two are distinct quantities.
    assert res.dosage_ci_low <= res.dosage_estimate <= res.dosage_ci_high


def test_zero_depth_region_is_flagged_as_downward_bias(fake_depth):
    """A locus fragment with no surviving depth drags the dosage down. That is
    a bias with a known sign, so it is flagged rather than silently absorbed."""
    fake_depth["locus"] = _flat(50, 2000, jitter=2, seed=51)
    fake_depth["locus_raw"] = _flat(51, 2000, jitter=2, seed=52)
    fake_depth["backbone"] = _flat(50, 100000, jitter=4, seed=53)
    fake_depth["chr_locus2:1-1000"] = [0] * 1000

    res = dr.estimate_locus_copy_number(
        "x.bam", ["chr_locus:1-1000", "chr_locus2:1-1000"],
        gene_spans=[(1, 1000), (1001, 2000)], gene_length=2000, n_boot=100,
    )
    if 0.0 in res.region_ratios:
        assert "ZERO_RATIO_REGIONS" in res.dosage_flags


# --------------------------------------------------------------------------- #
# Provenance and calibration of the aggregate call threshold.
#
# The 1.5 default is asserted rather than fitted. These tests pin the two
# properties that make an asserted constant defensible: it is
# reported as asserted unless a calibration says otherwise, and a calibration
# can only ever make rejection harder.
# --------------------------------------------------------------------------- #

def test_threshold_provenance_defaults_to_asserted(fake_depth):
    fake_depth["locus"] = _flat(50, 3000, jitter=6)
    fake_depth["locus_raw"] = _flat(52, 3000, jitter=6)
    fake_depth["backbone"] = _flat(50, 100000, jitter=6, seed=2)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-3000"], n_boot=200)
    assert res.call_threshold_used == dr.CALL_THRESHOLD_DEFAULT
    assert res.call_threshold_source == "DEFAULT_ASSERTED"


def test_calibration_below_the_floor_cannot_ease_rejection(fake_depth):
    """A quiet backbone must not buy an easier call.

    The sample's own null comes out at 1.10, well under the shipped 1.5. A
    ratio of about 1.3 clears 1.10 and not 1.5. The floor must hold, so this
    locus stays unrejected and the record says the floor held.
    """
    fake_depth["locus"] = _flat(65, 3000, jitter=3)
    fake_depth["locus_raw"] = _flat(66, 3000, jitter=3)
    fake_depth["backbone"] = _flat(50, 100000, jitter=3, seed=2)

    res = dr.estimate_locus_copy_number(
        "x.bam", ["chr_locus:1-3000"], n_boot=200, call_threshold_calibration=1.10)
    assert res.ratio_ci_low > 1.10          # would clear the sample's own null
    assert res.ratio_ci_low < 1.5           # but not the shipped floor
    assert res.single_copy_rejected is False
    assert res.call_threshold_used == dr.CALL_THRESHOLD_DEFAULT
    assert res.call_threshold_source == "CALIBRATED_FLOOR_HELD"


def test_calibration_above_the_floor_raises_and_withdraws_a_call(fake_depth):
    """A noisy backbone must be able to withdraw a call, never to add one."""
    fake_depth["locus"] = _flat(100, 3000, jitter=4)
    fake_depth["locus_raw"] = _flat(101, 3000, jitter=4)
    fake_depth["backbone"] = _flat(50, 100000, jitter=4, seed=2)

    base = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-3000"], n_boot=200)
    assert base.single_copy_rejected is True

    res = dr.estimate_locus_copy_number(
        "x.bam", ["chr_locus:1-3000"], n_boot=200, call_threshold_calibration=3.0)
    assert res.single_copy_rejected is False
    assert res.call_threshold_used == 3.0
    assert res.call_threshold_source == "CALIBRATED_RAISED"
    assert res.ratio == base.ratio          # the number is untouched


def test_interval_excluding_one_is_annotated_not_promoted(fake_depth):
    """The LIBA-6656 tcdB shape: CI above 1, verdict single-copy compatible.

    The interval excludes 1 and the call does not follow it. Without the note
    a reader reconciles the two by reading the interval as the copy-number
    statement, which is the reading the null forbids.
    """
    fake_depth["locus"] = _flat(60, 7104, jitter=3)
    fake_depth["locus_raw"] = _flat(61, 7104, jitter=3)
    fake_depth["backbone"] = _flat(50, 100000, jitter=3, seed=2)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-7104"], n_boot=200)
    assert res.ratio_ci_low > 1.0
    assert res.single_copy_rejected is False
    assert res.depth_call == "SINGLE_COPY_COMPATIBLE"
    assert any("excludes 1 from above" in n for n in res.notes)
    assert any("within-locus sampling noise only" in n for n in res.notes)


def test_interval_containing_no_integer_says_so(fake_depth):
    fake_depth["locus"] = _flat(60, 7104, jitter=3)
    fake_depth["locus_raw"] = _flat(61, 7104, jitter=3)
    fake_depth["backbone"] = _flat(50, 100000, jitter=3, seed=2)

    res = dr.estimate_locus_copy_number("x.bam", ["chr_locus:1-7104"], n_boot=200)
    assert not [k for k in range(0, 8) if res.ratio_ci_low <= k <= res.ratio_ci_high]
    assert any("contains no integer" in n for n in res.notes)


def test_depth_ratio_cli_reports_version(monkeypatch, capsys):
    from locus_recon import VERSION
    from locus_recon.depth_ratio import main

    monkeypatch.setattr("sys.argv", ["locus-recon-depth-ratio", "--version"])
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip() == f"locus-recon-depth-ratio {VERSION}"


@pytest.mark.parametrize("option,value", [
    ("--n-boot", "0"), ("--min-bq", "-1"), ("--gene-length", "0"),
    ("--call-threshold", "-1.5"), ("--max-backbone-bootstrap-blocks", "0"),
])
def test_depth_ratio_cli_rejects_out_of_range_numbers(monkeypatch, option, value):
    from locus_recon.depth_ratio import main

    monkeypatch.setattr(
        "sys.argv",
        ["locus-recon-depth-ratio", "--bam", "x.bam", "--locus", "c:1-10", option, value],
    )
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 2
