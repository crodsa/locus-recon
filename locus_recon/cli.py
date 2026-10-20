"""Command-line interface for reproducible Locus-Recon batch runs."""

import argparse
import os
import re
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from typing import List, Tuple

from . import VERSION, PROGRAM
from .utils import (
    DependencyError,
    PipelineStepError,
    check_dependencies,
    log,
    run_command,
    setup_logging,
)
from .io import (
    fastq_has_records,
    parse_samplesheet,
    write_allele_catalogs,
    write_batch_report,
)
from .qc import profile_bait_database
from .pipeline import process_sample, _empty_result
from .provenance import write_run_manifest


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("cannot be negative")
    return parsed


def _percentage(value: str) -> float:
    parsed = float(value)
    if not 0.0 <= parsed <= 100.0:
        raise argparse.ArgumentTypeError("must be between 0 and 100")
    return parsed


def _mixture_fraction(value: str) -> float:
    parsed = float(value)
    if not 0.0 < parsed < 0.5:
        raise argparse.ArgumentTypeError("must be greater than 0 and less than 0.5")
    return parsed


def _safe_name(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", value):
        raise argparse.ArgumentTypeError(
            "must start with a letter or number and contain only letters, "
            "numbers, dots, underscores, or hyphens"
        )
    return value


# ==============================================================================
# Parallel worker entry-point
# ==============================================================================


def _worker(args_tuple: tuple) -> dict:
    """Unpack arguments and call process_sample in a worker process.

    Using a top-level function (not a lambda or nested def) is required for
    ProcessPoolExecutor pickling.

    Args:
        args_tuple: (sample_id, assembly, r1, r2, global_args, tools, bait_profile)

    Returns:
        Per-sample result dict.
    """
    sample_id, assembly, r1, r2, global_args, tools, bait_profile = args_tuple
    return process_sample(
        sample_id, assembly, r1, r2, global_args, tools, bait_profile,
        parallel_worker=True,
    )


# ==============================================================================
# CLI entry point
# ==============================================================================


def main() -> None:
    """Parse arguments and run Locus-Recon in batch mode."""

    parser = argparse.ArgumentParser(
        prog="locus-recon",
        description=(
            f"{PROGRAM} v{VERSION} -- reconstruct a target locus from fragmented "
            "assemblies and the paired reads used to build them"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    # ── Required ─────────────────────────────────────────────────────────────
    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "--samplesheet", required=True, metavar="TSV",
        help=(
            "Four-column TSV: sample_id, assembly_path, r1_path, r2_path. "
            "A header is optional; relative paths are resolved from the TSV directory."
        ),
    )
    req.add_argument(
        "--main-output-dir", "--main_output_dir", dest="main_output_dir",
        required=True, metavar="DIR",
        help="Root output directory.",
    )
    req.add_argument(
        "-b", "--bait", required=True, metavar="FASTA",
        help="Multi-FASTA of curated known alleles for the target locus (e.g. from PubMLST).",
    )
    req.add_argument(
        "-l", "--locus", required=True, type=_safe_name, metavar="NAME",
        help="Locus name (e.g. 'aroE').",
    )

    # ── Runtime ───────────────────────────────────────────────────────────────
    run = parser.add_argument_group("runtime")
    run.add_argument(
        "-t", "--threads", type=_positive_int, default=8, metavar="N",
        help="Threads per sample (passed to BWA, SPAdes, samtools).",
    )
    run.add_argument(
        "-j", "--parallel-samples", "--parallel_samples",
        dest="parallel_samples", type=_positive_int, default=1, metavar="N",
        help=(
            "Number of samples to process in parallel. "
            "Each sample uses --threads threads, so total CPU = threads × parallel_samples. "
            "In parallel mode each sample writes its own log; the shared batch log "
            "is written by the main process only."
        ),
    )
    run.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print DEBUG messages to stderr.",
    )
    run.add_argument(
        "--no-progress", "--no_progress", dest="no_progress", action="store_true",
        help="Disable tqdm progress bars.",
    )
    run.add_argument(
        "--cleanup", action="store_true",
        help="Remove intermediate files after a successful sample.",
    )
    run.add_argument(
        "--memory-per-sample", type=_positive_int, default=16, metavar="GB",
        help="Maximum RAM passed to each SPAdes process.",
    )
    run.add_argument(
        "--dry-run", action="store_true",
        help="Validate inputs, dependencies, bait database, and manifest without processing reads.",
    )

    # ── BLAST / mapping thresholds ────────────────────────────────────────────
    filt = parser.add_argument_group("BLAST / mapping thresholds")
    filt.add_argument(
        "--min-discovery-identity", type=_percentage, default=80.0, metavar="PCT",
        help="Minimum identity for short bait HSPs used to find target regions.",
    )
    filt.add_argument(
        "--min-hsp-length", type=_positive_int, default=50, metavar="BP",
        help=(
            "Minimum bait HSP length during target discovery. Unlike query-coverage "
            "filtering, this permits recovery when a locus is split across contigs."
        ),
    )
    filt.add_argument(
        "--min-identity", "--min_identity", dest="min_identity",
        type=_percentage, default=85.0, metavar="PCT",
        help="Minimum identity to the nearest known allele during final validation.",
    )
    filt.add_argument(
        "--min-coverage", "--min_coverage", dest="min_coverage",
        type=_percentage, default=80.0, metavar="PCT",
        help="Minimum final reconstructed-allele query coverage during validation.",
    )
    filt.add_argument(
        "--region-padding", "--region_padding", dest="region_padding",
        type=int, default=300, metavar="BP",
        help="Padding (bp) added around BLAST target regions before read extraction.",
    )
    filt.add_argument(
        "--min-bitscore-margin", "--min_bitscore_margin",
        dest="min_bitscore_margin", type=float, default=10.0, metavar="SCORE",
        help="Minimum bitscore gap between best and second-best local BLAST hit.",
    )

    # ── Remap validation ──────────────────────────────────────────────────────
    remap = parser.add_argument_group("remap validation")
    remap.add_argument(
        "--min-remap-depth", "--min_remap_depth", dest="min_remap_depth",
        type=float, default=5.0, metavar="X",
        help="Minimum mean depth when remapping extracted reads to the final allele.",
    )
    remap.add_argument(
        "--min-remap-breadth", "--min_remap_breadth",
        dest="min_remap_breadth", type=_percentage, default=90.0, metavar="PCT",
        help="Minimum breadth of coverage (%%) when remapping to the final allele.",
    )

    evidence = parser.add_argument_group("per-base evidence and mixture detection")
    evidence.add_argument(
        "--min-base-quality", type=_nonnegative_int, default=20, metavar="PHRED",
        help="Minimum base quality included in depth and allele-fraction evidence.",
    )
    evidence.add_argument(
        "--min-mapping-quality", type=_nonnegative_int, default=20, metavar="PHRED",
        help="Minimum read mapping quality included in per-base evidence.",
    )
    evidence.add_argument(
        "--mixture-min-fraction", type=_mixture_fraction, default=0.05, metavar="FRACTION",
        help="Minimum non-reference fraction considered at a candidate mixed site.",
    )
    evidence.add_argument(
        "--mixture-min-sites", type=_positive_int, default=2, metavar="N",
        help="Minimum bidirectionally supported variant sites required to call a mixture.",
    )
    evidence.add_argument(
        "--mixture-min-alt-depth", type=_positive_int, default=2, metavar="N",
        help="Minimum quality-filtered alternative observations at a mixed site.",
    )

    # ── Locus options ─────────────────────────────────────────────────────────
    locus = parser.add_argument_group("locus options")
    locus.add_argument(
        "--noncoding-locus", "--noncoding_locus",
        dest="noncoding_locus", action="store_true",
        help="Disable reading-frame integrity checks for a non-coding target.",
    )

    args = parser.parse_args()
    args.use_progress = not args.no_progress

    if args.region_padding < 0:
        parser.error("--region-padding cannot be negative")
    if args.min_remap_depth < 0:
        parser.error("--min-remap-depth cannot be negative")
    if args.min_bitscore_margin < 0:
        parser.error("--min-bitscore-margin cannot be negative")

    args.samplesheet = os.path.abspath(os.path.expanduser(args.samplesheet))
    args.bait = os.path.abspath(os.path.expanduser(args.bait))
    args.main_output_dir = os.path.abspath(os.path.expanduser(args.main_output_dir))

    # ── Setup ─────────────────────────────────────────────────────────────────
    run_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    os.makedirs(args.main_output_dir, exist_ok=True)
    log_file_path = os.path.join(
        args.main_output_dir, f"locus_recon_batch_{args.locus}.log"
    )
    setup_logging(log_file_path, args.verbose)

    log.info("=" * 60)
    log.info(f"  {PROGRAM} v{VERSION} -- Batch run")
    log.info(f"  Locus : {args.locus}")
    log.info(f"  Bait  : {args.bait}")
    log.info(f"  Date  : {run_date}")
    log.info("=" * 60)

    if not os.path.isfile(args.bait):
        log.error(f"Bait file not found: {args.bait}")
        sys.exit(1)

    # ── Bait profile ──────────────────────────────────────────────────────────
    try:
        bait_profile = profile_bait_database(args.bait)
    except ValueError as exc:
        parser.error(str(exc))

    print(f"\n{'=' * 62}")
    print(f"  BAIT DATABASE PROFILE: {args.locus}")
    print(f"{'=' * 62}")
    print(f"  Alleles in DB  : {bait_profile['n_alleles']}")
    print(f"  Length range   : {bait_profile['len_min']}-{bait_profile['len_max']} bp")
    print(f"  Length median  : {bait_profile['len_median']:.0f} bp")
    print(f"  Length mode    : {bait_profile['len_mode']} bp")
    print(f"  IQR            : {bait_profile['len_q1']:.0f}-{bait_profile['len_q3']:.0f} "
          f"(range={bait_profile['len_iqr']:.0f} bp)")
    print(f"  Length stdev   : {bait_profile['len_stdev']:.1f} bp")
    print(f"  GC content     : {bait_profile['gc_mean']:.1f}% "
          f"+/- {bait_profile['gc_stdev']:.1f}%")
    eff_tol = (bait_profile["len_iqr"]
               if bait_profile["len_iqr"] > 0
               else max(bait_profile["len_stdev"], 3.0))
    print(f"  Eff. tolerance : +/- {eff_tol:.0f} bp (for HIGH confidence)")
    print(f"{'=' * 62}\n")

    # ── Dependency check ──────────────────────────────────────────────────────
    # ── Read samplesheet ──────────────────────────────────────────────────────
    try:
        samples = parse_samplesheet(args.samplesheet)
    except FileNotFoundError:
        log.error(f"Samplesheet not found: {args.samplesheet}")
        sys.exit(1)
    except ValueError as exc:
        parser.error(str(exc))

    if not samples:
        log.error("Samplesheet is empty.")
        sys.exit(1)

    log.info(f"Loaded {len(samples)} sample(s).")

    # Validate inputs before checking executables so malformed batch files fail
    # quickly even on a login node where the conda environment is not active.

    if args.parallel_samples > 1:
        log.info(
            f"Parallel mode: {args.parallel_samples} workers × "
            f"{args.threads} threads = "
            f"{args.parallel_samples * args.threads} logical CPUs max."
        )

    # ── Validate input files ──────────────────────────────────────────────────
    results_by_id = {}
    valid_samples: List[Tuple[str, str, str, str]] = []

    for sample_id, assembly, r1, r2 in samples:
        missing = [p for p in (assembly, r1, r2) if not os.path.isfile(p)]
        if missing:
            log.error(f"  Missing files for {sample_id}: {', '.join(missing)}")
            r = _empty_result(sample_id)
            r["message"] = f"Missing input files: {', '.join(missing)}"
            results_by_id[sample_id] = r
        elif os.path.getsize(assembly) == 0:
            r = _empty_result(sample_id)
            r["message"] = f"Assembly file is empty: {assembly}"
            log.error(f"  {sample_id}: {r['message']}")
            results_by_id[sample_id] = r
        elif not fastq_has_records(r1) or not fastq_has_records(r2):
            r = _empty_result(sample_id)
            r["message"] = "One or both FASTQ files are empty or unreadable."
            log.error(f"  {sample_id}: {r['message']}")
            results_by_id[sample_id] = r
        else:
            valid_samples.append((sample_id, assembly, r1, r2))

    try:
        tools, aligner_name = check_dependencies(args.use_progress)
    except DependencyError as exc:
        log.error(f"Dependency check failed: {exc}")
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    args.aligner_name = aligner_name

    bait_db_dir = os.path.join(args.main_output_dir, "bait_database")
    os.makedirs(bait_db_dir, exist_ok=True)
    args.bait_db = os.path.join(bait_db_dir, args.locus)
    bait_db_log = os.path.join(bait_db_dir, "makeblastdb.log")
    try:
        with open(bait_db_log, "w") as log_handle:
            run_command(
                [tools["makeblastdb"], "-in", args.bait, "-dbtype", "nucl",
                 "-parse_seqids", "-out", args.bait_db],
                log_handle,
                description="makeblastdb (shared bait database)",
            )
    except PipelineStepError as exc:
        log.error(f"Could not build the shared bait database: {exc}")
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    manifest_path = os.path.join(args.main_output_dir, "run_manifest.json")
    write_run_manifest(
        manifest_path=manifest_path,
        args=args,
        samples=samples,
        tools=tools,
        bait_profile=bait_profile,
        run_date=run_date,
    )

    if args.dry_run:
        if len(valid_samples) != len(samples):
            print("Dry run failed: one or more samples have missing inputs.", file=sys.stderr)
            sys.exit(1)
        print(
            f"Dry run passed: {len(samples)} sample(s), bait profile, and all "
            f"external dependencies validated. Manifest: {manifest_path}"
        )
        return

    # ── Run samples ───────────────────────────────────────────────────────────
    batch_start = time.monotonic()

    if args.parallel_samples <= 1 or len(valid_samples) == 1:
        # Sequential (default) — unchanged behaviour
        for i, (sample_id, assembly, r1, r2) in enumerate(valid_samples, 1):
            log.info(f"[{i}/{len(valid_samples)}] Starting: {sample_id}")
            r = process_sample(
                sample_id, assembly, r1, r2, args, tools, bait_profile,
                parallel_worker=False,
            )
            results_by_id[sample_id] = r
            log.info(f"[{i}/{len(valid_samples)}] Done: {sample_id} -> {r['status']}")
    else:
        # Parallel mode
        n_workers = min(args.parallel_samples, len(valid_samples))
        work_items = [
            (sid, asm, r1, r2, args, tools, bait_profile)
            for sid, asm, r1, r2 in valid_samples
        ]
        completed = 0
        total = len(work_items)

        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            future_to_id = {
                pool.submit(_worker, item): item[0]
                for item in work_items
            }
            for future in as_completed(future_to_id):
                sample_id = future_to_id[future]
                completed += 1
                try:
                    r = future.result()
                except Exception as exc:
                    log.error(
                        f"Worker for {sample_id} raised an exception: {exc}"
                    )
                    r = _empty_result(sample_id)
                    r["message"] = f"Worker exception: {exc}"
                results_by_id[sample_id] = r
                log.info(
                    f"[{completed}/{total}] {sample_id} -> "
                    f"{r['status']} ({r.get('elapsed_sec', 0):.1f}s)"
                )

    batch_elapsed = time.monotonic() - batch_start
    all_results = [results_by_id[sample_id] for sample_id, *_ in samples]

    # ── Batch report ──────────────────────────────────────────────────────────
    success = [r for r in all_results if r["status"] == "SUCCESS"]
    skipped = [r for r in all_results if r["status"] == "SKIP"]
    failed  = [r for r in all_results if r["status"] == "FAIL"]

    report_path = os.path.join(
        args.main_output_dir, f"locus_recon_report_{args.locus}.tsv"
    )
    write_batch_report(all_results, report_path, args.locus, run_date=run_date)

    # Disposition-aware candidate collections: accepted (PASS), review, hold,
    # and every technically successful candidate.
    catalogues = write_allele_catalogs(
        all_results, args.main_output_dir, args.locus
    )
    accepted = catalogues["accepted"]
    review = catalogues["review"]
    hold = catalogues["hold"]
    all_candidates = catalogues["all_candidates"]

    qc_counts = Counter(r.get("qc_confidence", "") for r in success)

    # ── Logging summary ───────────────────────────────────────────────────────
    log.info("=" * 60)
    log.info(f"  BATCH COMPLETE -- {args.locus}")
    log.info(f"  Total time : {timedelta(seconds=int(batch_elapsed))}")
    log.info(f"  Success    : {len(success)}")
    log.info(f"  Skipped    : {len(skipped)}")
    log.info(f"  Failed     : {len(failed)}")
    if qc_counts:
        log.info(
            f"  QC tiers   : "
            f"HIGH={qc_counts.get('HIGH', 0)}  "
            f"MEDIUM={qc_counts.get('MEDIUM', 0)}  "
            f"LOW={qc_counts.get('LOW', 0)}  "
            f"SUSPECT={qc_counts.get('SUSPECT', 0)}"
        )
    log.info(f"  Report     : {report_path}")
    log.info(
        "  Catalogues : accepted=%d review=%d hold=%d all=%d",
        accepted["count"], review["count"], hold["count"],
        all_candidates["count"],
    )
    log.info(f"  Log        : {log_file_path}")
    log.info("=" * 60)

    # ── Console summary ───────────────────────────────────────────────────────
    print(f"\n{'=' * 62}")
    print(f"  {PROGRAM} v{VERSION} -- BATCH SUMMARY")
    print(f"  Locus: {args.locus}")
    print(f"{'=' * 62}")
    print(f"  Total time : {timedelta(seconds=int(batch_elapsed))}")
    print(f"  Success    : {len(success)}")
    print(f"  Skipped    : {len(skipped)}")
    print(f"  Failed     : {len(failed)}")

    if qc_counts:
        print("\n  QC CONFIDENCE BREAKDOWN:")
        print(
            f"    [PASS]   HIGH    : {qc_counts.get('HIGH', 0):>3}  "
            "(strong computational support)"
        )
        print(f"    [REVIEW] MEDIUM  : {qc_counts.get('MEDIUM', 0):>3}  (review recommended)")
        print(f"    [WARN]   LOW     : {qc_counts.get('LOW', 0):>3}  (manual inspection needed)")
        print(
            f"    [HOLD]   SUSPECT : {qc_counts.get('SUSPECT', 0):>3}  "
            "(severe evidence issue -- do NOT submit)"
        )

    if all_results:
        print("\n  PER-SAMPLE RESULTS:")
        for r in all_results:
            if r["status"] == "SUCCESS":
                conf = r.get("qc_confidence", "?")
                tag  = {
                    "HIGH": "[PASS]  ", "MEDIUM": "[REVIEW]",
                    "LOW": "[WARN]  ", "SUSPECT": "[HOLD]  ",
                }.get(conf, "[????]  ")
                delta = r.get("qc_length_delta", 0)
                delta_str = f"{delta:+d}" if delta != 0 else "0"
                nearest = r.get("nearest_allele", "") or "none"
                print(
                    f"    {tag} {r['sample_id']:<20} "
                    f"{conf:<8} "
                    f"len={r['allele_length']} bp (delta={delta_str})  "
                    f"nearest={nearest} id={r['best_identity']:.1f}%"
                )
            elif r["status"] == "SKIP":
                print(f"    [SKIP]   {r['sample_id']:<20} {r['message']}")
            else:
                print(f"    [FAIL]   {r['sample_id']:<20} {r['message']}")

    print("\n  OUTPUT FILES:")
    print(f"    Report  : {report_path}")
    print(
        f"    Accepted candidates : {accepted['path']}  "
        f"({accepted['count']} PASS)"
    )
    print(
        f"    Review candidates   : {review['path']}  "
        f"({review['count']} REVIEW)"
    )
    print(
        f"    Hold candidates     : {hold['path']}  "
        f"({hold['count']} HOLD)"
    )
    print(
        f"    All candidates      : {all_candidates['path']}  "
        f"({all_candidates['count']} successful reconstructions)"
    )
    print(f"    Log     : {log_file_path}")
    print(f"    Manifest: {manifest_path}")
    print(f"{'=' * 62}")

    if not success and samples:
        sys.exit(1)


if __name__ == "__main__":
    main()
