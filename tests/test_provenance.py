import json
from argparse import Namespace

from locus_recon.provenance import sha256_file, write_run_manifest


def test_manifest_records_checksums_and_resolved_inputs(tmp_path):
    bait = tmp_path / "bait.fa"
    sheet = tmp_path / "samples.tsv"
    assembly = tmp_path / "assembly.fa"
    r1 = tmp_path / "r1.fq"
    r2 = tmp_path / "r2.fq"
    bait.write_text(">a1\nACGT\n")
    sheet.write_text("s1\tassembly.fa\tr1.fq\tr2.fq\n")
    assembly.write_text(">c1\nACGT\n")
    r1.write_text("@r\nAC\n+\nII\n")
    r2.write_text("@r\nGT\n+\nII\n")
    manifest_path = tmp_path / "manifest.json"
    args = Namespace(bait=str(bait), samplesheet=str(sheet), locus="abc")

    write_run_manifest(
        str(manifest_path), args,
        [("s1", str(assembly), str(r1), str(r2))],
        tools={}, bait_profile={"n_alleles": 1, "lengths": [4]},
        run_date="2026-01-01 00:00:00",
    )
    manifest = json.loads(manifest_path.read_text())
    assert manifest["bait"]["sha256"] == sha256_file(str(bait))
    assert manifest["bait"]["profile"]["n_alleles"] == 1
    assert "lengths" not in manifest["bait"]["profile"]
    assert manifest["samples"][0]["sample_id"] == "s1"

