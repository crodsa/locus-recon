import csv
import gzip

import pytest

from locus_recon.io import (
    classify_result_disposition,
    fastq_has_records,
    merge_intervals,
    parse_samplesheet,
    read_fasta_region,
    reverse_complement,
    write_allele_catalogs,
    write_batch_report,
)


def test_samplesheet_accepts_comments_header_and_relative_paths(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "# Locus-Recon input\n"
        "sample_id\tassembly_path\tr1_path\tr2_path\n"
        "isolate-01\tdata/a.fa\treads/a_R1.fastq.gz\treads/a_R2.fastq.gz\n"
    )
    samples = parse_samplesheet(str(sheet))
    assert samples[0][0] == "isolate-01"
    assert samples[0][1] == str(tmp_path / "data" / "a.fa")


@pytest.mark.parametrize("sample_id", ["../escape", "bad/id", "white space", ""])
def test_samplesheet_rejects_unsafe_sample_ids(tmp_path, sample_id):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"{sample_id}\ta.fa\tr1.fq\tr2.fq\n")
    with pytest.raises(ValueError, match="unsafe sample ID"):
        parse_samplesheet(str(sheet))


def test_samplesheet_rejects_duplicate_ids(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("s1\ta.fa\tr1.fq\tr2.fq\ns1\tb.fa\tr1.fq\tr2.fq\n")
    with pytest.raises(ValueError, match="duplicate"):
        parse_samplesheet(str(sheet))


def test_fastq_record_detection_handles_empty_gzip(tmp_path):
    empty = tmp_path / "empty.fastq.gz"
    populated = tmp_path / "reads.fastq.gz"
    with gzip.open(empty, "wt"):
        pass
    with gzip.open(populated, "wt") as handle:
        handle.write("@read1\nACGT\n+\nIIII\n")
    assert not fastq_has_records(str(empty))
    assert fastq_has_records(str(populated))


def test_fasta_region_and_reverse_complement(tmp_path):
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(">contig_1 description\nAACCGGTT\n")
    assert read_fasta_region(str(fasta), "contig_1", 2, 6) == "CCGG"
    assert reverse_complement("ACGTRYMKBDHVN") == "NBDHVMKRYACGT"


def test_merge_intervals_is_contig_aware():
    merged = merge_intervals([
        ("c1", 10, 20), ("c1", 18, 30), ("c2", 1, 3), ("c1", 40, 50),
    ])
    assert merged == [("c1", 10, 30), ("c1", 40, 50), ("c2", 1, 3)]


def test_batch_report_contains_traceability_fields(tmp_path):
    report = tmp_path / "report.tsv"
    write_batch_report([
        {
            "sample_id": "s1", "status": "SUCCESS", "allele_length": 500,
            "nearest_allele": "aroE_17", "exact_known_allele": True,
            "best_identity": 100.0, "validation_coverage": 100.0,
            "qc_confidence": "HIGH",
        }
    ], str(report), "aroE", run_date="2026-01-01 00:00:00")
    with report.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["nearest_allele"] == "aroE_17"
    assert row["exact_known_allele"] == "true"
    assert row["validation_coverage_pct"] == "100.0"


@pytest.mark.parametrize(
    ("workflow_status", "confidence", "expected"),
    [
        ("SUCCESS", "HIGH", "PASS"),
        ("SUCCESS", "MEDIUM", "REVIEW"),
        ("SUCCESS", "LOW", "REVIEW"),
        ("SUCCESS", "SUSPECT", "HOLD"),
        ("SKIP", "", "HOLD"),
        ("FAIL", "", "HOLD"),
    ],
)
def test_result_disposition_is_separate_from_workflow_completion(
    workflow_status, confidence, expected
):
    assert classify_result_disposition(workflow_status, confidence) == expected


def test_batch_report_writes_workflow_status_and_disposition(tmp_path):
    report = tmp_path / "report.tsv"
    write_batch_report([
        {
            "sample_id": "s1",
            "status": "SUCCESS",
            "allele_length": 500,
            "best_identity": 100.0,
            "validation_coverage": 100.0,
            "qc_confidence": "SUSPECT",
        }
    ], str(report), "aroE", run_date="2026-01-01 00:00:00")
    with report.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["status"] == "SUCCESS"  # documented compatibility alias
    assert row["workflow_status"] == "SUCCESS"
    assert row["result_disposition"] == "HOLD"


def _write_candidate(path, sample_id):
    path.write_text(f">{sample_id}\nACGT\n")
    return str(path)


def _fasta_ids(path):
    return [
        line[1:].strip()
        for line in path.read_text().splitlines()
        if line.startswith(">")
    ]


def test_catalogues_separate_accepted_review_hold_and_all_candidates(tmp_path):
    results = []
    for confidence in ("HIGH", "MEDIUM", "LOW", "SUSPECT"):
        sample_id = confidence.lower()
        results.append({
            "sample_id": sample_id,
            "status": "SUCCESS",
            "qc_confidence": confidence,
            "allele_file": _write_candidate(tmp_path / f"{sample_id}.fa", sample_id),
        })
    results.append({
        "sample_id": "skipped",
        "status": "SKIP",
        "qc_confidence": "",
        "allele_file": _write_candidate(tmp_path / "skipped.fa", "skipped"),
    })

    catalogues = write_allele_catalogs(results, str(tmp_path), "aroE")

    assert _fasta_ids(tmp_path / "aroE_accepted_alleles.fasta") == ["high"]
    assert _fasta_ids(tmp_path / "aroE_review_required_candidates.fasta") == [
        "medium", "low",
    ]
    assert _fasta_ids(tmp_path / "aroE_hold_candidates.fasta") == ["suspect"]
    assert _fasta_ids(tmp_path / "aroE_reconstructed_candidates.fasta") == [
        "high", "medium", "low", "suspect",
    ]
    assert _fasta_ids(tmp_path / "aroE_reconstructed_alleles.fasta") == ["high"]
    assert catalogues["accepted"]["count"] == 1
    assert catalogues["review"]["count"] == 2
    assert catalogues["hold"]["count"] == 1
    assert catalogues["all_candidates"]["count"] == 4
