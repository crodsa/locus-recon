import csv
import hashlib
import json
from pathlib import Path

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
