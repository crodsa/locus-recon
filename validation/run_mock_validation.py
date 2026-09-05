#!/usr/bin/env python3
"""Generate, run, and evaluate the complete Locus-Recon mock benchmark."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from common import build_provenance, require_files, write_json
from generate_mock_dataset import build_dataset


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    validation_dir = Path(__file__).resolve().parent
    parser.add_argument(
        "--dataset-dir", type=Path, default=validation_dir / "mock_dataset"
    )
    parser.add_argument(
        "--results-dir", type=Path, default=validation_dir / "mock_results"
    )
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--memory-gb", type=int, default=4)
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    dataset_dir = args.dataset_dir.resolve()
    results_dir = args.results_dir.resolve()
    if args.regenerate or not dataset_dir.exists():
        build_dataset(dataset_dir, force=args.regenerate and dataset_dir.exists())
    if results_dir.exists():
        if not args.force:
            parser.error(f"{results_dir} exists; use --force to replace it")
        shutil.rmtree(results_dir)

    root_dir = validation_dir.parent
    command = [
        sys.executable, "-m", "locus_recon.cli",
        "--samplesheet", str(dataset_dir / "samples.tsv"),
        "--main-output-dir", str(results_dir),
        "--bait", str(dataset_dir / "mock_locus_alleles.fasta"),
        "--locus", "mockLocus",
        "--threads", str(args.threads),
        "--parallel-samples", "1",
        "--memory-per-sample", str(args.memory_gb),
        "--min-bitscore-margin", "80",
        "--min-remap-depth", "3",
        "--min-remap-breadth", "90",
        "--min-base-quality", "20",
        "--min-mapping-quality", "20",
        "--mixture-min-fraction", "0.05",
        "--mixture-min-sites", "2",
        "--mixture-min-alt-depth", "2",
        "--no-progress",
    ]
    print("Running:", " ".join(command), flush=True)
    pipeline = subprocess.run(command, cwd=root_dir, check=False)
    if pipeline.returncode != 0:
        print(
            f"mock reconstruction failed with exit {pipeline.returncode}; "
            "evaluation was not run",
            file=sys.stderr,
        )
        raise SystemExit(pipeline.returncode)

    report = results_dir / "locus_recon_report_mockLocus.tsv"
    manifest = results_dir / "run_manifest.json"
    require_files([report, manifest])

    evaluation = subprocess.run([
        sys.executable,
        str(validation_dir / "evaluate_mock_results.py"),
        "--dataset-dir", str(dataset_dir),
        "--results-dir", str(results_dir),
    ], cwd=root_dir, check=False)
    if evaluation.returncode != 0:
        raise SystemExit(evaluation.returncode)

    summary = results_dir / "validation" / "validation_summary.json"
    require_files([summary])
    if not json.loads(summary.read_text()).get("all_core_criteria_passed"):
        raise SystemExit("mock evaluation completed without passing all core criteria")

    provenance = build_provenance(
        workflow="deterministic-mock-end-to-end",
        repo_root=root_dir,
        inputs=[
            dataset_dir / "dataset_manifest.json",
            dataset_dir / "samples.tsv",
            dataset_dir / "mock_locus_alleles.fasta",
            report,
            manifest,
            summary,
        ],
        parameters={
            "threads": args.threads,
            "memory_gb": args.memory_gb,
            "regenerate": args.regenerate,
        },
        tools={
            "python": [sys.executable, "--version"],
            "blastn": ["blastn", "-version"],
            "bwa_mem2": ["bwa-mem2", "version"],
            "samtools": ["samtools", "--version"],
            "spades": ["spades.py", "--version"],
        },
    )
    provenance["complete"] = True
    write_json(results_dir / "validation" / "workflow_provenance.json", provenance)


if __name__ == "__main__":
    main()
