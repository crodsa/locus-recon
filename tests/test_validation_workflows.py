from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from locus_recon import VERSION
from validation.common import (
    build_provenance,
    require_files,
    run_checked,
    sha256_file,
    write_json,
)
from validation.workflows.run_low_copy_23S import (
    LEGACY_RESULT_FIELDS,
    _audit_truth,
    _require_legacy_depth_compatibility,
)


def test_sha256_and_json_are_deterministic(tmp_path: Path):
    payload = tmp_path / "payload.txt"
    payload.write_text("locus-recon\n")
    expected = hashlib.sha256(b"locus-recon\n").hexdigest()
    assert sha256_file(payload) == expected

    destination = tmp_path / "record.json"
    write_json(destination, {"b": 2, "a": 1})
    assert destination.read_text() == '{\n  "a": 1,\n  "b": 2\n}\n'


def test_require_files_names_every_missing_input(tmp_path: Path):
    present = tmp_path / "present"
    present.write_text("x")
    with pytest.raises(FileNotFoundError) as exc:
        require_files([present, tmp_path / "missing_a", tmp_path / "missing_b"])
    assert "missing_a" in str(exc.value)
    assert "missing_b" in str(exc.value)


def test_run_checked_fails_loudly_and_records_command():
    with pytest.raises(RuntimeError) as exc:
        run_checked([sys.executable, "-c", "import sys; sys.exit(7)"])
    assert "exit 7" in str(exc.value)
    assert sys.executable in str(exc.value)


def test_provenance_hashes_inputs_and_records_parameters(tmp_path: Path):
    input_file = tmp_path / "input.txt"
    input_file.write_text("evidence")
    record = build_provenance(
        workflow="unit-test",
        repo_root=Path(__file__).resolve().parents[1],
        inputs=[input_file],
        parameters={"seed": 17},
        tools={"python": [sys.executable, "--version"]},
    )
    assert record["schema_version"] == "1.0"
    assert record["workflow"] == "unit-test"
    assert record["parameters"] == {"seed": 17}
    assert record["inputs"][0]["sha256"] == sha256_file(input_file)
    # Inputs outside the repository are named, not located: the checksum
    # identifies them, and an absolute path only describes the machine.
    assert record["inputs"][0]["path"] == "input.txt"
    assert record["tools"]["python"]["command"] == [Path(sys.executable).name, "--version"]
    assert record["tools"]["python"]["version"]
    software = record["software"]
    assert software["name"] == "locus-recon"
    assert software["version"] == VERSION
    assert len(software["source_sha256"]) == 64
    # A provenance record is always produced. Inside a checkout the commit is
    # the full SHA; outside one it is a recorded value or an explicit sentinel.
    # Asserting a 40-char SHA unconditionally made this test a statement about
    # the machine it ran on rather than about build_provenance.
    commit = record["git_commit"]
    assert isinstance(commit, str) and commit
    assert len(commit) == 40 or commit.startswith("UNAVAILABLE") or "recorded" in commit


def test_provenance_survives_outside_a_git_checkout(tmp_path: Path):
    """Source archives ship validation/ but no .git.

    Querying git outside a checkout must therefore not raise: a workflow that
    failed here would do so after writing its results and before writing
    provenance.json, leaving a reproduced run undocumented.
    """
    input_file = tmp_path / "input.txt"
    input_file.write_text("evidence")
    record = build_provenance(
        workflow="unit-test",
        repo_root=tmp_path,          # a directory, not a checkout
        inputs=[input_file],
        parameters={"seed": 17},
        tools={"python": [sys.executable, "--version"]},
    )
    assert record["git_commit"].startswith("UNAVAILABLE")

    (tmp_path / "COMMIT").write_text("67f1c81\n")
    recorded = build_provenance(
        workflow="unit-test",
        repo_root=tmp_path,
        inputs=[input_file],
        parameters={"seed": 17},
        tools={"python": [sys.executable, "--version"]},
    )
    assert recorded["git_commit"].startswith("67f1c81")


def test_provenance_json_round_trip(tmp_path: Path):
    path = tmp_path / "p.json"
    write_json(path, {"schema_version": "1.0", "complete": True})
    assert json.loads(path.read_text())["complete"] is True


def test_committed_constructed_copy_series_matches_prespecified_criteria():
    root = Path(__file__).resolve().parents[1]
    with (root / "validation/constructed-copy-series/results.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    evaluation = json.loads(
        (root / "validation/constructed-copy-series/evaluation.json").read_text()
    )

    assert len(rows) == 27
    assert all(
        row["depth_call"] == "SINGLE_COPY_COMPATIBLE"
        for row in rows if row["truth_copies"] == "1"
    )
    assert all(
        row["depth_call"] == "MULTICOPY_DEPTH"
        for row in rows if row["truth_copies"] in {"2", "3"}
    )
    assert evaluation["all_prespecified_criteria_passed"] is True


def test_committed_liba_depth_result_reports_the_additive_dosage():
    """The committed artifact must agree with what the released estimator returns.

    A committed artifact that disagrees with the estimator would let the test
    pin the file without exercising the code, so these assertions are the
    values regenerated from the BAM. The aggregate ratio and the additive
    dosage are asserted together because they answer different questions and
    disagreeing is their expected behaviour on a fragmented locus. Overlap
    holds the dosage in a review state: an empirical control on single-copy
    sequence with the same geometry reproduces the value, so it is reported
    with its geometric null and not as a count.
    """
    root = Path(__file__).resolve().parents[1]
    evaluation = json.loads(
        (root / "validation/liba6656-depth/evaluation.json").read_text()
    )

    assert evaluation["schema_version"] == "1.0"
    assert evaluation["depth_call"] == "SINGLE_COPY_COMPATIBLE"
    assert evaluation["dosage_status"] == "REVIEW_OVERLAPPING_SPANS"
    assert evaluation["dosage_estimate"] == 2.287
    # 83.3% of the query is resolved by more than one interval, so the additive
    # value is not a copy count and the record carries the geometric null.
    assert evaluation["geometry_expected_dosage"] == 1.843
    assert (evaluation["dosage_ci_low"], evaluation["dosage_ci_high"]) == (2.069, 2.702)
    assert evaluation["query_coverage_fraction"] == 0.9876
    assert evaluation["dosage_flags"] == [
        "INCOMPLETE_QUERY_COVERAGE",
        "OVERLAPPING_QUERY_SPANS",
        "ZERO_RATIO_REGIONS",
    ]
    assert "profile_call" not in evaluation


def test_committed_aphA1_panel_retains_full_class_and_qpcr_concordance():
    root = Path(__file__).resolve().parents[1]
    with (root / "validation/aphA1-copy-number/results.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    evaluation = json.loads(
        (root / "validation/aphA1-copy-number/evaluation.json").read_text()
    )

    assert len(rows) == 7
    assert sum(row["depth_call"] == "SINGLE_COPY_COMPATIBLE" for row in rows) == 4
    assert sum(row["depth_call"] == "MULTICOPY_DEPTH" for row in rows) == 3
    assert evaluation["class_concordance"] is True
    assert evaluation["quantitative_qpcr_intervals_contain_estimate"] is True
    assert hashlib.sha256(
        (root / "validation/aphA1-copy-number/results.tsv").read_bytes()
    ).hexdigest() == "0b50d0084b695f83414eb0a4c3d924e785751bf98d10700463e0cd550f6dbe4c"
    assert hashlib.sha256(
        (root / "validation/aphA1-copy-number/evaluation.json").read_bytes()
    ).hexdigest() == "c6f185f0726b688ea9069c0f049d5ba861ba11da9fb336b1b72401ded3ef5e5d"
    assert hashlib.sha256(
        (root / "validation/aphA1-copy-number/run_summary.json").read_bytes()
    ).hexdigest() == "63fe7d8a2723ff85468848a37c114916e127c40ace44244ffd70cce4a2da12fc"


def test_committed_low_copy_panel_meets_prespecified_two_copy_boundary():
    root = Path(__file__).resolve().parents[1]
    with (root / "validation/low-copy-23S/results.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    evaluation = json.loads(
        (root / "validation/low-copy-23S/evaluation.json").read_text()
    )

    assert len(rows) == 10
    rrna = [row for row in rows if row["locus"] == "23S"]
    controls = [row for row in rows if row["locus"] == "gyrB"]
    assert len(rrna) == len(controls) == 5
    assert sum(row["depth_call"] == "MULTICOPY_DEPTH" for row in rrna) >= 4
    assert not any(row["depth_call"] == "MULTICOPY_DEPTH" for row in controls)
    assert all(row["graph_context_count"] == "2" for row in rrna)
    assert all(row["copy_number_call"] == "2" for row in rrna)
    assert all(row["graph_context_count"] == "1" for row in controls)
    assert all(row["copy_number_call"] == "1" for row in controls)
    assert evaluation["criterion_legacy_depth_fields_unchanged"] is True
    assert evaluation["criterion_all_two_copy_consensus_calls_correct"] is True
    assert evaluation["criterion_all_single_copy_control_consensus_calls_correct"] is True
    assert evaluation["all_primary_criteria_passed"] is True


def test_low_copy_legacy_compatibility_audit_names_changed_field(tmp_path: Path):
    baseline = tmp_path / "baseline.tsv"
    candidate = tmp_path / "candidate.tsv"
    row = {field: "unchanged" for field in LEGACY_RESULT_FIELDS}
    row.update({"strain": "strain-1", "locus": "23S", "ratio": "2.0"})
    with baseline.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    row["ratio"] = "2.1"
    with candidate.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    with pytest.raises(RuntimeError, match="strain-1/23S:ratio"):
        _require_legacy_depth_compatibility(baseline, candidate)


def test_low_copy_truth_audit_ignores_circular_plasmid(tmp_path: Path):
    chromosome = SeqRecord(
        Seq("A" * 5000), id="CP_TEST_CHR.1", description="strain chromosome, complete genome"
    )
    chromosome.annotations.update({"topology": "circular", "molecule_type": "DNA"})
    chromosome.features = [
        SeqFeature(SimpleLocation(0, 5000), type="source"),
        SeqFeature(
            SimpleLocation(0, 2900), type="rRNA",
            qualifiers={"product": ["23S ribosomal RNA"]},
        ),
        SeqFeature(
            SimpleLocation(1000, 3900), type="rRNA",
            qualifiers={"product": ["23S ribosomal RNA"]},
        ),
        SeqFeature(
            SimpleLocation(2000, 4300), type="CDS", qualifiers={"gene": ["gyrB"]}
        ),
    ]
    plasmid = SeqRecord(
        Seq("C" * 1000), id="CP_TEST_PLS.1", description="plasmid pTest, complete sequence"
    )
    plasmid.annotations.update({"topology": "circular", "molecule_type": "DNA"})
    plasmid.features = [
        SeqFeature(
            SimpleLocation(0, 1000), type="source", qualifiers={"plasmid": ["pTest"]}
        )
    ]

    truth_dir = tmp_path / "truth"
    truth_dir.mkdir()
    gbff = truth_dir / "GCF_TEST_ASM_TEST_genomic.gbff.gz"
    import gzip
    from Bio import SeqIO

    with gzip.open(gbff, "wt") as handle:
        SeqIO.write([chromosome, plasmid], handle, "genbank")
    row = {
        "assembly_ftp": "https://example.invalid/GCF_TEST_ASM_TEST",
        "assembly_accession": "GCF_TEST.1",
        "strain": "test",
        "biosample": "SAMNTEST",
        "truth_23S_copies": "2",
        "truth_control_copies": "1",
    }

    audit, observed_gbff = _audit_truth(row, truth_dir)

    assert observed_gbff == gbff
    assert audit["chromosome_accession"] == "CP_TEST_CHR.1"
    assert audit["full_length_23S_copies"] == 2
    assert audit["other_deposited_replicons"] == [{
        "accession": "CP_TEST_PLS.1",
        "length": 1000,
        "topology": "circular",
        "classified_as_plasmid": True,
    }]
