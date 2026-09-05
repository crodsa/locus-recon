import json

import pytest

from locus_recon import VERSION
from locus_recon.copy_number_cli import main
from locus_recon.depth_ratio import DepthRatioResult


def _depth_result() -> DepthRatioResult:
    return DepthRatioResult(
        locus_regions=["contig:1-100"],
        locus_median_depth=240.0,
        locus_median_depth_unfiltered=245.0,
        backbone_median_depth=100.0,
        backbone_mad=5.0,
        ratio=2.4,
        ratio_ci_low=2.1,
        ratio_ci_high=2.7,
        assembly_copies=1,
        copies_estimate=2.4,
        ambiguity_index=0.98,
        gc_locus=None,
        single_copy_rejected=True,
        reliable=True,
        assembly_region_count=1,
        dosage_estimate=2.4,
        dosage_status="ESTIMATED",
        depth_call="MULTICOPY_DEPTH",
        dosage_flags=[],
    )


def test_depth_only_cli_writes_additive_json_and_tsv(tmp_path, monkeypatch, capsys):
    bam = tmp_path / "reads.bam"
    bam.write_bytes(b"bam-placeholder")
    output_json = tmp_path / "result.json"
    output_tsv = tmp_path / "result.tsv"
    observed = {}

    def fake_estimate(path, **parameters):
        observed["path"] = path
        observed["parameters"] = parameters
        return _depth_result()

    monkeypatch.setattr(
        "locus_recon.copy_number_cli.estimate_locus_copy_number", fake_estimate,
    )
    monkeypatch.setattr(
        "locus_recon.copy_number_cli._tool_version", lambda executable: "test-version",
    )

    assert main([
        "--bam", str(bam),
        "--locus", "contig:1-100",
        "--output-json", str(output_json),
        "--output-tsv", str(output_tsv),
        "--n-boot", "10",
    ]) == 0

    result = json.loads(output_json.read_text())
    assert result["software_version"] == VERSION
    assert result["ratio"] == 2.4
    assert result["dosage_estimate"] == 2.4
    assert result["graph_context_status"] == "NOT_EVALUATED"
    assert result["copy_number_call"] == 2.4
    assert result["copy_number_kind"] == "MEAN_DEPTH_DOSAGE"
    assert result["copy_number_method"] == "DEPTH_ONLY"
    assert result["provenance"]["input_sha256"]["bam"]
    assert result["provenance"]["parameters"]["depth"]["n_boot"] == 10
    assert output_tsv.read_text().splitlines()[0].startswith("schema_version\t")
    assert observed["parameters"]["locus_regions"] == ["contig:1-100"]
    assert "copy_number_call=2.4" in capsys.readouterr().out


def test_cli_requires_gfa_and_bait_together(tmp_path):
    bam = tmp_path / "reads.bam"
    gfa = tmp_path / "graph.gfa"
    bam.write_bytes(b"bam-placeholder")
    gfa.write_text("S\t1\tAAAA\n")

    with pytest.raises(SystemExit) as exc:
        main([
            "--bam", str(bam),
            "--locus", "contig:1-100",
            "--gfa", str(gfa),
            "--output-json", str(tmp_path / "result.json"),
        ])
    assert exc.value.code == 2


def test_cli_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert f"locus-recon-copy-number {VERSION}" in capsys.readouterr().out
