import csv
import hashlib
import json
import sys
from pathlib import Path

import pytest

from locus_recon import VERSION
from validation.common import (
    build_provenance,
    require_files,
    run_checked,
    sha256_file,
    write_json,
)
from validation.evaluate_mock_results import evaluate
from validation.generate_mock_dataset import build_dataset


ROOT = Path(__file__).resolve().parents[1]
COMMITTED_DATASET = ROOT / "validation" / "mock_dataset"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_mock_dataset_is_deterministically_regenerated(tmp_path):
    regenerated = tmp_path / "mock_dataset"
    build_dataset(regenerated)
    assert sha256(regenerated / "dataset_manifest.json") == sha256(
        COMMITTED_DATASET / "dataset_manifest.json"
    )
    assert sha256(regenerated / "reads" / "mixture_20pct_B_R1.fastq.gz") == sha256(
        COMMITTED_DATASET / "reads" / "mixture_20pct_B_R1.fastq.gz"
    )


def test_complete_genome_is_the_exact_locus_truth():
    manifest = json.loads((COMMITTED_DATASET / "dataset_manifest.json").read_text())
    genome = "".join(
        line.strip()
        for line in (
            COMMITTED_DATASET / "truth" / "complete_genome_A.fasta"
        ).read_text().splitlines()
        if not line.startswith(">")
    )
    alleles = {}
    current = None
    for line in (
        COMMITTED_DATASET / "truth" / "truth_alleles.fasta"
    ).read_text().splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            alleles[current] = ""
        else:
            alleles[current] += line.strip()
    start = manifest["locus_start_0_based"]
    end = manifest["locus_end_0_based_exclusive"]
    assert genome[start:end] == alleles["mockLocus_A"]
    assert end - start == manifest["allele_length"]


def test_truth_evaluator_core_criteria(tmp_path):
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    report_path = results_dir / "locus_recon_report_mockLocus.tsv"
    truth_fasta = COMMITTED_DATASET / "truth" / "truth_alleles.fasta"
    with (COMMITTED_DATASET / "truth" / "expected_results.tsv").open() as handle:
        expected = list(csv.DictReader(handle, delimiter="\t"))

    columns = [
        "sample_id", "status", "allele_file", "qc_confidence", "qc_flags",
        "mixture_detected", "mixed_site_count", "median_mixed_fraction",
        "strand_biased_sites", "remap_mean_depth", "remap_breadth_pct", "message",
    ]
    with report_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in expected:
            negative = row["category"] == "negative"
            mixture = row["category"] == "mixture"
            paralog = row["category"] == "paralog"
            writer.writerow({
                "sample_id": row["sample_id"],
                "status": "SKIP" if negative else "SUCCESS",
                "allele_file": "" if negative else str(truth_fasta),
                "qc_confidence": "SUSPECT" if mixture or paralog else "HIGH",
                "qc_flags": "AMBIGUOUS_BEST_HIT" if paralog else "",
                "mixture_detected": str(mixture).lower(),
                "mixed_site_count": 8 if mixture else 0,
                "median_mixed_fraction": row["mixture_fraction"],
                "strand_biased_sites": 0,
                "remap_mean_depth": row["nominal_depth"],
                "remap_breadth_pct": 100,
                "message": "synthetic evaluator test",
            })

    summary = evaluate(COMMITTED_DATASET, results_dir, results_dir / "validation")
    assert summary["all_core_criteria_passed"]
    assert all(summary["criteria"].values())


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
    assert record["inputs"][0]["path"] == "input.txt"
    assert record["tools"]["python"]["command"] == [Path(sys.executable).name, "--version"]
    assert record["tools"]["python"]["version"]
    software = record["software"]
    assert software["name"] == "locus-recon"
    assert software["version"] == VERSION
    assert len(software["source_sha256"]) == 64
    commit = record["git_commit"]
    assert isinstance(commit, str) and commit
    assert len(commit) == 40 or commit.startswith("UNAVAILABLE") or "recorded" in commit


def test_provenance_survives_outside_a_git_checkout(tmp_path: Path):
    input_file = tmp_path / "input.txt"
    input_file.write_text("evidence")
    record = build_provenance(
        workflow="unit-test",
        repo_root=tmp_path,
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
