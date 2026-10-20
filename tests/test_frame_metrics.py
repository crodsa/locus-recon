"""Reading-frame metrics, bait validation, and length-profile reporting."""

import random

import pytest

from locus_recon.io import reverse_complement
from locus_recon.qc import (
    MIN_LENGTH_PROFILE_N,
    assess_allele_quality,
    count_internal_stops_six_frames,
    format_qc_report,
    profile_bait_database,
    _count_internal_stops,
)

REMAP = {
    "mean_depth": 90.0, "breadth_pct": 100.0,
    "pct_bases_lt5": 0.0, "mapped_reads": 5000,
}


STOP_CODON_TRIPLETS = ("TAA", "TAG", "TGA")


def _orf(codons: int = 120, seed: int = 0) -> str:
    """A clean reading frame whose reverse complement has none.

    Built from a fixed seed so the fixture is deterministic: a tandem repeat
    would read cleanly on both strands and could not distinguish a six-frame
    scan from a forward-only one.
    """
    alphabet = [
        first + second + third
        for first in "ACGT" for second in "ACGT" for third in "ACGT"
        if first + second + third not in STOP_CODON_TRIPLETS
    ]
    body = "".join(random.Random(seed).choices(alphabet, k=codons))
    return "ATG" + body + "TAA"


def _write_bait(tmp_path, records, name="bait.fasta"):
    path = tmp_path / name
    path.write_text("".join(f">{header}\n{sequence}\n" for header, sequence in records))
    return str(path)


def test_six_frames_find_a_reverse_strand_reading_frame():
    antisense = reverse_complement(_orf())
    forward_only = min(_count_internal_stops(antisense, frame) for frame in range(3))
    stops, label, counts = count_internal_stops_six_frames(antisense)
    assert forward_only > 0          # the 1.0.0 forward-only view
    assert (stops, label) == (0, "-1")
    assert len(counts) == 6


def test_profile_reports_an_antisense_bait_catalogue(tmp_path):
    orf = _orf()
    bait = _write_bait(
        tmp_path,
        [(f"allele_{i}", reverse_complement(orf[: len(orf) - 3 * i])) for i in range(6)],
    )
    profile = profile_bait_database(bait)
    assert profile["bait_appears_antisense"] is True
    assert profile["expected_coding_frame_labels"] == ["-1"]
    assert set(profile["six_frame_stop_counts"]) == {"+1", "+2", "+3", "-1", "-2", "-3"}


def test_antisense_bait_does_not_raise_internal_stops(tmp_path):
    """The false positive this test pins: a clean ORF read on the wrong strand."""
    orf = _orf()
    bait = _write_bait(
        tmp_path,
        [(f"allele_{i}", reverse_complement(orf[: len(orf) - 3 * i])) for i in range(6)],
    )
    profile = profile_bait_database(bait)
    qc = assess_allele_quality(
        allele_seq=reverse_complement(orf), bait_profile=profile,
        validation_identity=99.0, validation_qcov=100.0, locus="demo",
        remap_metrics=REMAP, expect_cds=True,
    )
    assert qc["internal_stops"] == 0
    assert qc["coding_frame_used"] == "-1"
    assert not any(flag.startswith("INTERNAL_STOPS") for flag in qc["flags"])


def test_long_identifier_is_rejected_before_makeblastdb(tmp_path):
    bait = _write_bait(tmp_path, [("A" * 77, _orf()), ("short", _orf())])
    with pytest.raises(ValueError, match="parse_seqids"):
        profile_bait_database(bait)


def test_gapped_alignment_is_named_as_such(tmp_path):
    bait = _write_bait(tmp_path, [("aligned_1", "ACGT---ACGT"), ("aligned_2", "ACGTACGTAC")])
    with pytest.raises(ValueError, match="gapped alignment"):
        profile_bait_database(bait)


def test_paralogous_fragment_is_reported(tmp_path):
    orf = _orf()
    records = [(f"allele_{i}", orf[: len(orf) - 3 * i]) for i in range(5)]
    records.append(("paralog_fragment", "GTCAAGGTCAAGCCTTGCAGGTTGCA" * 4))
    profile = profile_bait_database(_write_bait(tmp_path, records))
    assert profile["unrelated_records"] == ["paralog_fragment"]


def test_length_flags_carry_the_profile_size(tmp_path):
    orf = _orf()
    small = profile_bait_database(_write_bait(
        tmp_path, [(f"allele_{i}", orf[: len(orf) - 30 * i]) for i in range(3)],
        name="small.fasta",
    ))
    assert small["n_alleles"] < MIN_LENGTH_PROFILE_N
    assert small["length_profile_underpowered"] is True
    qc = assess_allele_quality(
        allele_seq=orf[:200], bait_profile=small,
        validation_identity=99.0, validation_qcov=90.0, locus="demo",
        remap_metrics=REMAP, expect_cds=True,
    )
    length_flags = [flag for flag in qc["flags"] if flag.startswith("LENGTH_")]
    assert length_flags and f"profile n={small['n_alleles']}" in length_flags[0]
    assert any(
        flag.startswith("CATALOGUE_LENGTH_PROFILE_UNDERPOWERED")
        for flag in qc["flags"]
    )
    assert qc["length_profile_n"] == small["n_alleles"]


def test_frame_length_shift_is_attributed_to_a_contig_clip(tmp_path):
    orf = _orf()
    profile = profile_bait_database(_write_bait(
        tmp_path, [(f"allele_{i}", orf[: len(orf) - 3 * i]) for i in range(6)],
    ))
    truncated = orf[:1000] + "A"          # length % 3 != bait mode
    qc = assess_allele_quality(
        allele_seq=truncated, bait_profile=profile,
        validation_identity=99.0, validation_qcov=60.0, locus="demo",
        remap_metrics=REMAP, expect_cds=True,
        span_metrics={"clipped_bp": 1200, "clipped_start_bp": 0,
                      "clipped_end_bp": 1200, "unclaimed_query": 1200},
    )
    assert any(flag.startswith("FRAME_LENGTH_SHIFT_TRUNCATED") for flag in qc["flags"])
    assert not any(flag.startswith("FRAME_LENGTH_SHIFT (") for flag in qc["flags"])


def test_placeholder_frameshift_is_named_as_the_source_of_stops(tmp_path):
    """A placeholder whose length is not a multiple of three manufactures stops."""
    orf = _orf(codons=200)
    profile = profile_bait_database(_write_bait(
        tmp_path, [(f"allele_{i}", orf[: len(orf) - 3 * i]) for i in range(6)],
        name="frameshift_bait.fasta",
    ))
    midpoint = 3 * (len(orf) // 6)
    scaffolded = orf[:midpoint] + "N" * 100 + orf[midpoint:]
    qc = assess_allele_quality(
        allele_seq=scaffolded, bait_profile=profile,
        validation_identity=99.0, validation_qcov=100.0, locus="demo",
        remap_metrics=REMAP, expect_cds=True,
    )
    assert qc["internal_stops"] > 0                      # as delivered
    assert qc["internal_stops_placeholder_closed"] == 0  # without the placeholder
    assert any(flag.startswith("INTERNAL_STOPS") for flag in qc["flags"])
    explanation = [
        flag for flag in qc["flags"]
        if flag.startswith("PLACEHOLDER_FRAMESHIFT_EXPLAINS_STOPS")
    ]
    assert explanation and "100 bp of interior placeholder" in explanation[0]


def test_genuine_stops_are_not_attributed_to_a_placeholder(tmp_path):
    """Stops that survive excision are reported without the explanation."""
    orf = _orf(codons=200)
    profile = profile_bait_database(_write_bait(
        tmp_path, [(f"allele_{i}", orf[: len(orf) - 3 * i]) for i in range(6)],
        name="genuine_bait.fasta",
    ))
    midpoint = 3 * (len(orf) // 6)
    damaged = orf[:midpoint] + "N" * 99 + "TAA" + orf[midpoint + 3:]
    qc = assess_allele_quality(
        allele_seq=damaged, bait_profile=profile,
        validation_identity=99.0, validation_qcov=100.0, locus="demo",
        remap_metrics=REMAP, expect_cds=True,
    )
    assert qc["internal_stops_placeholder_closed"] >= 1
    assert not any(
        flag.startswith("PLACEHOLDER_FRAMESHIFT_EXPLAINS_STOPS")
        for flag in qc["flags"]
    )


def test_qc_report_names_the_six_frame_labels_used_in_the_table(tmp_path):
    """The readable scorecard and the batch table must name the same frame."""
    antisense = reverse_complement(_orf())
    bait = _write_bait(tmp_path, [(f"a{i}", antisense) for i in range(3)])
    profile = profile_bait_database(bait)
    qc = assess_allele_quality(
        antisense, profile, 100.0, 100.0, "locus",
        remap_metrics=REMAP,
        ambiguity_metrics={"ambiguous": False, "bitscore_margin": 100.0},
    )
    report = format_qc_report("sample", "locus", qc, profile)

    assert qc["coding_frame_used"].startswith("-")
    assert f"Frame used    : {qc['coding_frame_used']:>7}" in report
    assert "Bait frame(s) : " + ", ".join(profile["expected_coding_frame_labels"]) in report
    cds = report[report.index("CDS INTEGRITY"):]
    assert not any(label in cds for label in ("+0", "-0"))
