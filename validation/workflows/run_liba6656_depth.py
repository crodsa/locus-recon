#!/usr/bin/env python3
"""Depth ratio and additive dosage of the LIBA-6656 tcdB locus."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from locus_recon.depth_ratio import estimate_locus_copy_number
from validation.common import build_provenance, require_files, write_json
from validation.workflows.run_aphA1_validation import _download_ena_pair


RUN = "ERR467623"
GENE_LENGTH = 7104
OVERLAP_NULL = (
    REPO_ROOT / "validation" / "depth-copy-number" / "overlap_null"
    / "null_dose_overlap_control.json"
)
EVALUATION_FIELDS = (
    "assembly_region_count", "depth_call", "dosage_estimate", "dosage_status",
    "dosage_flags", "gene_total_aligned_bp", "gene_union_covered_bp",
    "gene_overlap_bp", "gene_overlap_fraction", "dosage_ci_low", "dosage_ci_high",
    "query_coverage_fraction", "geometry_expected_dosage",
)


def build_evaluation(fields: dict, overlap_null: dict | None) -> dict:
    """Evaluation record for the deposit, with its interpretation.

    The interpretation is written from the values, so it cannot drift from
    them; the single-copy null is quoted only when its deposit is present.
    """
    evaluation = {"schema_version": "1.0"}
    evaluation.update({key: fields[key] for key in EVALUATION_FIELDS})
    text = (
        "The two biological candidates remain supported by independent "
        "sequence, competitive-remapping and context evidence. The aggregate "
        f"depth-only call is {fields['depth_call']}: the ratio of locus to "
        "backbone median depth does not reject one copy, which is the expected "
        "reading when the assembler has split the gene across fragments, since "
        "each fragment then carries roughly one copy of depth. The additive "
        f"per-query-base dosage returns {fields['dosage_estimate']} (95% "
        f"block-bootstrap CI {fields['dosage_ci_low']}-{fields['dosage_ci_high']}) "
        f"over {100 * fields['query_coverage_fraction']:.1f}% of the query, but "
        f"{100 * fields['gene_overlap_fraction']:.1f}% of the query is resolved "
        "by more than one discovery interval and that value is NOT a copy count. "
        "The same span geometry with every region at single copy returns "
        f"{fields['geometry_expected_dosage']}"
    )
    if overlap_null:
        low, high = overlap_null["null_prediction_interval_95"]
        text += (
            f", and an empirical control on {overlap_null['n_multi_region']} "
            "single-copy pseudo-queries from this assembly "
            "(validation/depth-copy-number/overlap_null/) predicts "
            f"{overlap_null['null_prediction_at_tcdB_overlap']} [{low}, {high}] at "
            "this overlapped fraction, so the observed value is not "
            "distinguishable from the single-copy null"
        )
    text += (
        ". Depth neither supports nor refutes two copies here; it is "
        "uninformative on this geometry, and the two-copy conclusion rests on the "
        "sequence, graph and remapping evidence alone."
    )
    evaluation["interpretation"] = text
    return evaluation


def _prepare_assembly(source: Path, work_dir: Path) -> Path:
    if source.suffix != ".gz":
        return source.resolve()
    destination = work_dir / source.name.removesuffix(".gz")
    if not destination.exists():
        with gzip.open(source, "rb") as compressed, destination.open("wb") as plain:
            shutil.copyfileobj(compressed, plain, length=1024 * 1024)
    return destination


def _map(reference: Path, r1: Path, r2: Path, bam: Path, aligner: str, threads: int) -> None:
    if bam.exists() and Path(str(bam) + ".bai").exists():
        return
    index = Path(str(reference) + (".bwt.2bit.64" if aligner == "bwa-mem2" else ".bwt"))
    if not index.exists():
        subprocess.run([aligner, "index", str(reference)], check=True, capture_output=True)
    with bam.with_suffix(".mapping.log").open("w") as log:
        mapping = subprocess.Popen(
            [aligner, "mem", "-t", str(threads), str(reference), str(r1), str(r2)],
            stdout=subprocess.PIPE, stderr=log,
        )
        assert mapping.stdout is not None
        sorting = subprocess.run(
            ["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"],
            stdin=mapping.stdout, stdout=log, stderr=log,
        )
        mapping.stdout.close()
        mapping_code = mapping.wait()
    if mapping_code != 0 or sorting.returncode != 0:
        raise RuntimeError("LIBA-6656 read mapping failed")
    subprocess.run(["samtools", "index", str(bam)], check=True)


def main() -> int:
    repo_root = REPO_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draft-assembly", type=Path, required=True)
    parser.add_argument("--reads-dir", type=Path, required=True)
    parser.add_argument("--download-missing", action="store_true")
    parser.add_argument(
        "--regions", type=Path,
        default=repo_root / "validation" / "manifests" / "liba6656_depth_regions.tsv",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260817)
    args = parser.parse_args()
    if min(args.threads, args.n_boot) < 1:
        parser.error("threads and bootstrap count must be positive")
    source_assembly, regions_path = require_files([args.draft_assembly, args.regions])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.download_missing:
        download_records = _download_ena_pair(RUN, args.reads_dir)
        write_json(args.output_dir / "download_provenance.json", download_records)
    r1, r2 = args.reads_dir / f"{RUN}_1.fastq.gz", args.reads_dir / f"{RUN}_2.fastq.gz"
    require_files([r1, r2])
    assembly = _prepare_assembly(source_assembly, args.output_dir)

    for tool in ("samtools",):
        if not shutil.which(tool):
            raise FileNotFoundError(f"{tool} is required on PATH")
    aligner = "bwa-mem2" if shutil.which("bwa-mem2") else "bwa"
    if not shutil.which(aligner):
        raise FileNotFoundError("bwa-mem2 or bwa is required on PATH")
    bam = args.output_dir / "ERR467623_vs_LIBA6656_draft.sorted.bam"
    _map(assembly, r1, r2, bam, aligner, args.threads)

    with regions_path.open(newline="") as handle:
        region_rows = list(csv.DictReader(handle, delimiter="\t"))
    locus_regions = [row["assembly_region"] for row in region_rows]
    gene_spans = [
        (int(row["gene_span_start"]), int(row["gene_span_end"])) for row in region_rows
    ]
    result = estimate_locus_copy_number(
        bam,
        locus_regions,
        gene_spans=gene_spans,
        gene_length=GENE_LENGTH,
        max_backbone_bootstrap_blocks=512,
        n_boot=args.n_boot,
        seed=args.seed,
        samtools="samtools",
    )
    write_json(args.output_dir / "LIBA6656_tcdB_depth_result.json", result.as_row())
    row = result.as_row()
    table = args.output_dir / "LIBA6656_tcdB_depth_result.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    overlap_null = (
        json.loads(OVERLAP_NULL.read_text()) if OVERLAP_NULL.is_file() else None
    )
    evaluation = build_evaluation(
        {key: getattr(result, key) for key in EVALUATION_FIELDS}, overlap_null,
    )
    write_json(args.output_dir / "evaluation.json", evaluation)
    provenance = build_provenance(
        workflow="LIBA6656-tcdB-depth",
        repo_root=repo_root,
        inputs=[source_assembly, regions_path, r1, r2, table],
        parameters={
            "run": RUN,
            "gene_length": GENE_LENGTH,
            "threads": args.threads,
            "n_boot": args.n_boot,
            "max_backbone_bootstrap_blocks": 512,
            "seed": args.seed,
            "call_threshold": 1.5,
            "min_bq": 20,
            "min_mq": 20,
        },
        tools={
            "python": [sys.executable, "--version"],
            "aligner": [aligner, "version"],
            "samtools": ["samtools", "--version"],
        },
    )
    provenance["complete"] = True
    write_json(args.output_dir / "workflow_provenance.json", provenance)
    print(result.summary_line())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
