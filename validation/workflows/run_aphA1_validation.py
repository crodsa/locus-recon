#!/usr/bin/env python3
"""Reproduce the seven-run aphA1 depth and copy-number validation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from locus_recon.depth_ratio import estimate_locus_copy_number
from validation.common import build_provenance, require_files, write_json


EXPECTED_RUNS = 7
DEFAULT_GENE_REGION = "CP080452.1:250422-251237"
DEFAULT_EXCLUSIONS = ["CP080453.1", "CP080454.1", "CP080455.1", "CP080456.1"]
ENA_REPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def _md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_ena_pair(run: str, reads_dir: Path) -> list[dict]:
    query = urllib.parse.urlencode({
        "accession": run,
        "result": "read_run",
        "fields": "fastq_ftp,fastq_md5",
        "format": "tsv",
    })
    with urllib.request.urlopen(f"{ENA_REPORT}?{query}", timeout=60) as response:
        lines = response.read().decode().strip().splitlines()
    if len(lines) != 2:
        raise RuntimeError(f"ENA returned {len(lines) - 1} records for {run}, expected one")
    report = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    urls = report["fastq_ftp"].split(";")
    md5s = report["fastq_md5"].split(";")
    if len(urls) != 2 or len(md5s) != 2:
        raise RuntimeError(f"{run} is not represented by exactly two ENA FASTQ files")
    records = []
    reads_dir.mkdir(parents=True, exist_ok=True)
    for mate, (remote, expected_md5) in enumerate(zip(urls, md5s), start=1):
        destination = reads_dir / f"{run}_{mate}.fastq.gz"
        if not destination.exists():
            temporary = destination.with_suffix(".fastq.gz.part")
            request = urllib.request.Request(
                f"https://{remote}", headers={"User-Agent": "locus-recon-validation/1.0"}
            )
            with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            temporary.replace(destination)
        observed_md5 = _md5(destination)
        if observed_md5 != expected_md5:
            raise RuntimeError(
                f"MD5 mismatch for {destination}: expected {expected_md5}, observed {observed_md5}"
            )
        records.append({
            "run": run,
            "mate": mate,
            "url": f"https://{remote}",
            "ena_md5": expected_md5,
            "local_md5": observed_md5,
        })
    return records


def _load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != EXPECTED_RUNS:
        raise ValueError(f"aphA1 manifest must contain {EXPECTED_RUNS} runs, found {len(rows)}")
    runs = [row["illumina_run"] for row in rows]
    if len(set(runs)) != EXPECTED_RUNS:
        raise ValueError("aphA1 manifest contains duplicate Illumina runs")
    return rows


def _aligner() -> str:
    for program in ("bwa-mem2", "bwa"):
        if shutil.which(program):
            return program
    raise FileNotFoundError("bwa-mem2 or bwa is required on PATH")


def _ensure_reference_index(reference: Path, aligner: str) -> None:
    mem2_index = Path(str(reference) + ".bwt.2bit.64")
    bwa_index = Path(str(reference) + ".bwt")
    expected = mem2_index if aligner == "bwa-mem2" else bwa_index
    if not expected.exists():
        subprocess.run([aligner, "index", str(reference)], check=True)


def _map_reads(
    *, aligner: str, samtools: str, reference: Path, r1: Path, r2: Path,
    bam: Path, threads: int,
) -> None:
    log = bam.with_suffix(".mapping.log")
    with log.open("w") as log_handle:
        align = subprocess.Popen(
            [aligner, "mem", "-t", str(threads), str(reference), str(r1), str(r2)],
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )
        assert align.stdout is not None
        sort = subprocess.run(
            [samtools, "sort", "-@", str(threads), "-o", str(bam), "-"],
            stdin=align.stdout,
            stdout=log_handle,
            stderr=log_handle,
        )
        align.stdout.close()
        align_code = align.wait()
    if align_code != 0 or sort.returncode != 0:
        bam.unlink(missing_ok=True)
        raise RuntimeError(
            f"mapping failed for {r1.name}: aligner={align_code}, samtools={sort.returncode}; "
            f"see {log}"
        )
    subprocess.run([samtools, "index", str(bam)], check=True)


def _qpcr_interval(text: str) -> tuple[float, float] | None:
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*\(\+/-\s*(\d+(?:\.\d+)?)\)\s*", text)
    if not match:
        return None
    centre, error = (float(value) for value in match.groups())
    return centre - error, centre + error


def main() -> int:
    repo_root = REPO_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reads-dir", type=Path, required=True)
    parser.add_argument(
        "--manifest", type=Path,
        default=repo_root / "validation" / "manifests" / "aphA1_copy_number.tsv",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--aphA1-region", default=DEFAULT_GENE_REGION)
    parser.add_argument("--exclude", action="append", default=None)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--force-remap", action="store_true")
    parser.add_argument(
        "--download-missing", action="store_true",
        help="download missing paired FASTQ files from ENA and verify their MD5 values",
    )
    args = parser.parse_args()

    reference, manifest = require_files([args.reference, args.manifest])
    rows = _load_manifest(manifest)
    reads: list[Path] = []
    download_records: list[dict] = []
    for row in rows:
        run = row["illumina_run"]
        if args.download_missing:
            download_records.extend(_download_ena_pair(run, args.reads_dir))
        reads.extend([args.reads_dir / f"{run}_1.fastq.gz", args.reads_dir / f"{run}_2.fastq.gz"])
    require_files(reads)
    if args.download_missing:
        write_json(args.output_dir / "download_provenance.json", download_records)

    if args.threads < 1 or args.n_boot < 1:
        parser.error("--threads and --n-boot must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    aligner, samtools = _aligner(), "samtools"
    if not shutil.which(samtools):
        raise FileNotFoundError("samtools is required on PATH")
    _ensure_reference_index(reference, aligner)
    exclusions = args.exclude if args.exclude is not None else DEFAULT_EXCLUSIONS

    records: list[dict] = []
    summary_rows: list[dict] = []
    for truth in rows:
        run = truth["illumina_run"]
        r1, r2 = args.reads_dir / f"{run}_1.fastq.gz", args.reads_dir / f"{run}_2.fastq.gz"
        bam = args.output_dir / f"{run}.reads_vs_MRSN56.sorted.bam"
        if args.force_remap:
            bam.unlink(missing_ok=True)
            Path(str(bam) + ".bai").unlink(missing_ok=True)
        if not bam.exists() or not Path(str(bam) + ".bai").exists():
            _map_reads(
                aligner=aligner, samtools=samtools, reference=reference,
                r1=r1, r2=r2, bam=bam, threads=args.threads,
            )

        result = estimate_locus_copy_number(
            bam,
            [args.aphA1_region],
            gene_spans=[(1, 816)],
            gene_length=816,
            reference_fasta=reference,
            exclude_regions=exclusions,
            max_backbone_bootstrap_blocks=512,
            n_boot=args.n_boot,
            seed=args.seed,
            samtools=samtools,
        )
        observed = result.dosage_estimate
        qpcr = _qpcr_interval(truth["truth_qpcr_copies"])
        within_qpcr = qpcr is not None and observed is not None and qpcr[0] <= observed <= qpcr[1]
        record = {
            "run": run,
            "strain": truth["sample"],
            "expected_ratio_class": truth["expected_ratio_class"],
            "truth_qpcr_copies": truth["truth_qpcr_copies"],
            "truth_wgs_coverage_fold": truth["truth_wgs_coverage_fold"],
            "tobramycin_mic_ug_ml": truth["tobramycin_mic_ug_ml"],
            "qpcr_interval": list(qpcr) if qpcr else None,
            "estimate_within_qpcr_interval": within_qpcr if qpcr else None,
            "result": result.as_row(),
        }
        records.append(record)
        summary_rows.append({
            "run": run,
            "strain": truth["sample"],
            "expected_ratio_class": truth["expected_ratio_class"],
            "truth_qpcr_copies": truth["truth_qpcr_copies"],
            "dosage_estimate": result.dosage_estimate,
            "ratio_ci_low": result.ratio_ci_low,
            "ratio_ci_high": result.ratio_ci_high,
            "depth_call": result.depth_call,
            "dosage_status": result.dosage_status,
            "ambiguity_index": result.ambiguity_index,
            "estimate_within_qpcr_interval": within_qpcr if qpcr else "",
        })
        print(
            f"{run}\t{truth['sample']}\tdosage={result.dosage_estimate}\t"
            f"{result.depth_call}", flush=True
        )

    write_json(args.output_dir / "aphA1_depth_ratio.json", records)
    table = args.output_dir / "aphA1_copy_number.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    classes_pass = all(
        (row["depth_call"] == "MULTICOPY_DEPTH") == (row["expected_ratio_class"] != "~1x")
        for row in summary_rows
    )
    qpcr_pass = all(
        row["estimate_within_qpcr_interval"] is True
        for row in summary_rows if row["truth_qpcr_copies"] in {"10 (+/-2)", "75 (+/-14)"}
    )
    evaluation = {
        "schema_version": "1.0",
        "runs": len(summary_rows),
        "class_concordance": classes_pass,
        "quantitative_qpcr_intervals_contain_estimate": qpcr_pass,
        "all_prespecified_criteria_passed": classes_pass and qpcr_pass,
    }
    write_json(args.output_dir / "evaluation.json", evaluation)
    provenance = build_provenance(
        workflow="aphA1-seven-run-copy-number",
        repo_root=repo_root,
        inputs=[reference, manifest, *reads, table, args.output_dir / "aphA1_depth_ratio.json"],
        parameters={
            "aphA1_region": args.aphA1_region,
            "exclude": exclusions,
            "threads": args.threads,
            "n_boot": args.n_boot,
            "max_backbone_bootstrap_blocks": 512,
            "seed": args.seed,
            "min_bq": 20,
            "min_mq": 20,
            "call_threshold": 1.5,
        },
        tools={
            "python": [sys.executable, "--version"],
            "aligner": [aligner, "version"],
            "samtools": [samtools, "--version"],
        },
    )
    provenance["complete"] = evaluation["all_prespecified_criteria_passed"]
    write_json(args.output_dir / "workflow_provenance.json", provenance)
    return 0 if evaluation["all_prespecified_criteria_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
