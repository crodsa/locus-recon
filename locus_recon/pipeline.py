"""
Per-sample pipeline orchestration: process_sample() and cleanup helpers.
"""

import argparse
import os
import shutil
import time

from .utils import (
    add_sample_log_handler, log, PipelineStepError, run_command,
    configure_sample_logger,
    BLAST_FMT_STEP1, BLAST_FMT_STEP5,
)
from .io import (
    classify_result_disposition,
    get_fasta_lengths, merge_intervals, count_fasta_sequences,
    get_fasta_total_length, read_single_fasta_sequence, write_single_fasta,
    fastq_has_records, read_fasta_region, reverse_complement,
)
from .mapping import (
    run_bwa_to_sorted_bam, collect_mapping_stats,
    write_read_name_list_from_bam, select_alignments_by_read_names,
    remap_reads_to_final_allele,
)
from .blast import (
    collapse_overlapping_hits,
    describe_full_query_span,
    find_span_continuation,
    parse_discovery_blast_hit,
    rank_local_blast_hit,
    summarize_hit_ambiguity,
)
from .assembly_graph import assess_placeholder_junctions
from .qc import (
    assess_allele_quality, format_qc_report, write_qc_report_file,
)


def build_spades_command(
    spades: str,
    out_dir: str,
    threads: int,
    memory_gb: int,
    *,
    r1: str = None,
    r2: str = None,
    singletons: str = None,
    single_cell: bool = False,
) -> list:
    """Return the SPAdes command for one local assembly.

    SPAdes runs in its default mode, or with ``--sc`` for the documented
    retry, and never with ``--careful``: each reported base is checked by
    remapping the recruited reads instead of by SPAdes' mismatch corrector, so
    the command is the same whichever SPAdes version is installed.
    """
    command = [spades]
    if single_cell:
        command.append("--sc")
    if singletons:
        command += ["-s", singletons]
    if r1 and r2:
        command += ["-1", r1, "-2", r2]
    command += [
        "-o", out_dir, "-t", str(threads), "-m", str(memory_gb),
        "--phred-offset", "33",
    ]
    return command


# ---------------------------------------------------------------------------
# Empty result template (avoids repeating dict keys in multiple code paths)
# ---------------------------------------------------------------------------

def _empty_result(sample_id: str) -> dict:
    return {
        "sample_id":          sample_id,
        "status":             "FAIL",
        "workflow_status":    "FAIL",
        "result_disposition": "HOLD",
        "allele_file":        "",
        "allele_length":      0,
        "best_identity":      0.0,
        "validation_coverage": 0.0,
        "nearest_allele":     "",
        "exact_known_allele": False,
        "sequence_confidence": "",
        "catalogue_status":   "",
        "nearest_allele_identity": 0.0,
        "catalogue_review_recommended": False,
        "contigs_hit":        0,
        "target_regions":     0,
        "extracted_read_names": 0,
        "spades_mode":         "",
        "mapped_reads":       0,
        "elapsed_sec":        0,
        "message":            "",
        "qc_confidence":      "",
        "qc_length_delta":    0,
        "qc_length_zscore":   0.0,
        "qc_gc_pct":          0.0,
        "qc_gc_deviation":    0.0,
        "qc_n_count":         0,
        "qc_internal_stops":  0,
        "qc_coding_frame_used": "",
        "qc_internal_stops_placeholder_closed": 0,
        "qc_length_mod3":     0,
        "qc_length_profile_n": 0,
        "placeholder_runs":   0,
        "placeholder_bp":     0,
        "placeholder_junction_support": "NOT_ASSESSED",
        "qc_flags":           "",
        "span_clipped_bp":       0,
        "span_clipped_start_bp": 0,
        "span_clipped_end_bp":   0,
        "span_continuation_contig":    "",
        "span_continuation_bp":        0,
        "span_continuation_overlap_bp": 0,
        "remap_mean_depth":   0.0,
        "remap_breadth_pct":  0.0,
        "remap_pct_bases_lt5": 0.0,
        "uncertain_base_count":   0,
        "uncertain_base_fraction": 0.0,
        "interior_uncertain_base_count": 0,
        "interior_uncertain_base_fraction": 0.0,
        "internal_zero_depth_positions": 0,
        "terminal_margin_bp":     0,
        "mean_base_quality":  0.0,
        "mean_mapping_quality": 0.0,
        "strand_balance_pct": 0.0,
        "candidate_mixed_sites": 0,
        "mixed_site_count":   0,
        "strand_biased_sites": 0,
        "reference_discordant_sites": 0,
        "max_alt_fraction":   0.0,
        "median_mixed_fraction": 0.0,
        "mixture_detected":   False,
        "bitscore_margin":    0.0,
    }


# ==============================================================================
# Core per-sample pipeline
# ==============================================================================


def process_sample(
    sample_id: str,
    assembly: str,
    r1: str,
    r2: str,
    global_args: argparse.Namespace,
    tools: dict,
    bait_profile: dict,
    parallel_worker: bool = False,
) -> dict:
    """Execute the full Locus-Recon pipeline for a single isolate.

    Steps:
        1. BLAST bait vs assembly to identify target contigs.
        2. Map all reads to the assembly (BWA/BWA-mem2).
        3. Extract target-region reads + mate rescue.
        4. Local assembly with SPAdes.
        5. Extract allele from local assembly (bitscore-ranked best hit).
        6. Validate reconstructed allele by BLAST vs bait + remap support.
        7. QC scorecard (length, GC, N-content, CDS integrity, remap depth).

    Args:
        sample_id:       Unique sample identifier.
        assembly:        Path to draft genome FASTA.
        r1:              Path to forward Illumina reads.
        r2:              Path to reverse Illumina reads.
        global_args:     Parsed CLI arguments (argparse.Namespace).
        tools:           Tool path mapping from check_dependencies().
        bait_profile:    Dict from profile_bait_database().
        parallel_worker: If True, reconfigure the module logger to write only
                         to this sample's log file (avoids multi-process contention
                         on the shared batch log).

    Returns:
        Dict with status, QC metrics, and file paths.
    """
    t_start = time.monotonic()
    result  = _empty_result(sample_id)

    # ── Per-sample output directory ──────────────────────────────────────────
    sample_out_dir  = os.path.join(global_args.main_output_dir, sample_id)
    os.makedirs(sample_out_dir, exist_ok=True)
    sample_log_path = os.path.join(
        sample_out_dir, f"{sample_id}_{global_args.locus}_recon.log"
    )

    # In parallel workers, redirect Python logging to the sample log so that
    # concurrent writes to the shared batch log file are avoided.
    with open(sample_log_path, "w"):
        pass
    sample_log_handler = None
    if parallel_worker:
        configure_sample_logger(sample_log_path, mode="a")
    else:
        sample_log_handler = add_sample_log_handler(sample_log_path)

    log.info("=" * 60)
    log.info(f"  SAMPLE: {sample_id}")
    log.info("=" * 60)

    local_assembly = os.path.join(sample_out_dir, "assembly_input.fasta")
    threads_str  = str(global_args.threads)
    db_path      = os.path.join(sample_out_dir, "assembly_db")
    blast_out    = os.path.join(sample_out_dir, "blast_step1.tsv")
    target_bed   = os.path.join(sample_out_dir, "target_regions.bed")

    try:
        if os.path.exists(local_assembly):
            os.remove(local_assembly)
        shutil.copy2(assembly, local_assembly)
        contig_lengths = get_fasta_lengths(local_assembly)
        if not contig_lengths:
            raise ValueError("Assembly FASTA contains no sequence records.")

        with open(sample_log_path, "a") as slh:

            # ── STEP 1: BLAST bait vs assembly ──────────────────────────────
            log.info("  STEP 1/7: Identifying target contigs (BLAST)...")
            run_command(
                [tools["makeblastdb"], "-in", local_assembly,
                 "-dbtype", "nucl", "-out", db_path],
                slh, description="makeblastdb (assembly)",
            )
            run_command(
                [tools["blastn"], "-query", global_args.bait,
                 "-db", db_path, "-out", blast_out,
                 "-outfmt", BLAST_FMT_STEP1,
                 "-evalue", "1e-10", "-num_threads", threads_str],
                slh, description="blastn (bait vs assembly)",
            )

            target_regions = []
            target_contigs = set()
            hits_total = hits_passed = 0

            with open(blast_out) as bfh:
                for line in bfh:
                    parts = line.strip().split("\t")
                    if len(parts) < 11:
                        continue
                    hits_total += 1
                    hit = parse_discovery_blast_hit(parts)
                    if (hit["pident"] < global_args.min_discovery_identity
                            or hit["length"] < global_args.min_hsp_length):
                        continue
                    hits_passed += 1
                    sseqid = hit["sseqid"]
                    sstart, send = hit["sstart"], hit["send"]
                    target_contigs.add(sseqid)
                    raw_start = min(sstart, send) - 1
                    raw_end   = max(sstart, send)
                    pad       = global_args.region_padding
                    target_regions.append((
                        sseqid,
                        max(0, raw_start - pad),
                        min(contig_lengths.get(sseqid, raw_end), raw_end + pad),
                    ))

            if not target_contigs:
                msg = f"No contigs matched bait ({hits_total} raw, 0 passed)."
                log.warning(f"  SKIP {sample_id}: {msg}")
                result["status"]  = "SKIP"
                result["message"] = msg
                return result

            target_regions = merge_intervals(target_regions)
            with open(target_bed, "w") as fout:
                for contig, start, end in target_regions:
                    fout.write(f"{contig}\t{start}\t{end}\n")

            result["contigs_hit"] = len(target_contigs)
            result["target_regions"] = len(target_regions)
            log.info(
                f"    {len(target_contigs)} target contig(s), "
                f"{hits_passed} HSPs passed, "
                f"{len(target_regions)} merged interval(s)."
            )

            # ── STEP 2: Map reads to assembly ────────────────────────────────
            log.info("  STEP 2/7: Mapping reads to assembly...")
            run_command(
                [tools[global_args.aligner_name], "index", local_assembly],
                slh, description=f"{global_args.aligner_name} index",
            )
            sorted_bam = os.path.join(sample_out_dir, "full_map.sorted.bam")
            bwa_log     = os.path.join(sample_out_dir, "bwa_pipe.log")
            run_bwa_to_sorted_bam(
                tools, global_args.aligner_name, global_args.threads,
                local_assembly, r1, r2, global_args.locus,
                bwa_log, sorted_bam, slh, global_args.use_progress,
            )
            run_command(
                [tools["samtools"], "index", "-@", threads_str, sorted_bam],
                slh, description="samtools index",
            )
            flagstat_path = os.path.join(sample_out_dir, "mapping_flagstat.txt")
            mapping_stats = collect_mapping_stats(
                tools["samtools"], sorted_bam, flagstat_path, slh, global_args.threads,
            )
            result["mapped_reads"] = mapping_stats.get("mapped_reads", 0)

            # ── STEP 3: Extract target reads + mate rescue ───────────────────
            log.info("  STEP 3/7: Extracting reads from target regions...")
            target_bam = os.path.join(sample_out_dir, "target_region.bam")
            run_command(
                [tools["samtools"], "view", "-@", threads_str, "-b", "-h",
                 "-L", target_bed, "-o", target_bam, sorted_bam],
                slh, description="samtools view -L",
            )
            read_names_path = os.path.join(sample_out_dir, "target_read_names.txt")
            n_seed = write_read_name_list_from_bam(
                tools["samtools"], target_bam, read_names_path,
            )
            result["extracted_read_names"] = n_seed
            if n_seed == 0:
                msg = "No reads mapped to target regions."
                log.warning(f"  SKIP {sample_id}: {msg}")
                result["status"] = "SKIP"
                result["message"] = msg
                return result
            rescued_bam = os.path.join(sample_out_dir, "target_region_rescued.bam")
            select_alignments_by_read_names(
                tools, global_args.threads, sorted_bam,
                read_names_path, rescued_bam, slh,
            )
            namesorted_bam = os.path.join(sample_out_dir, "target_namesorted.bam")
            run_command(
                [tools["samtools"], "sort", "-n", "-@", threads_str,
                 "-o", namesorted_bam, rescued_bam],
                slh, description="samtools sort -n",
            )

            r1_gz = os.path.join(sample_out_dir, "target_R1.fq.gz")
            r2_gz = os.path.join(sample_out_dir, "target_R2.fq.gz")
            s_gz  = os.path.join(sample_out_dir, "target_S.fq.gz")
            run_command(
                [tools["samtools"], "fastq", "-@", threads_str,
                 "-1", r1_gz, "-2", r2_gz, "-s", s_gz, namesorted_bam],
                slh, description="samtools fastq",
            )

            has_pairs = fastq_has_records(r1_gz) and fastq_has_records(r2_gz)
            has_singletons = fastq_has_records(s_gz)
            if not has_pairs and not has_singletons:
                msg = "No reads were extracted after mate rescue."
                log.warning(f"  SKIP {sample_id}: {msg}")
                result["status"]  = "SKIP"
                result["message"] = msg
                return result
            log.info(f"    Extracted target read names: {n_seed}")

            # ── STEP 4: Local assembly ────────────────────────────────────────
            log.info("  STEP 4/7: Local assembly with SPAdes...")
            spades_dir = os.path.join(sample_out_dir, "spades_local")
            local_asm  = os.path.join(spades_dir, "scaffolds.fasta")
            if os.path.isdir(spades_dir):
                shutil.rmtree(spades_dir)
            spades_inputs = {
                "r1": r1_gz if has_pairs else None,
                "r2": r2_gz if has_pairs else None,
                "singletons": s_gz if has_singletons else None,
            }
            spades_cmd = build_spades_command(
                tools["spades.py"], spades_dir, global_args.threads,
                global_args.memory_per_sample, **spades_inputs,
            )
            try:
                run_command(spades_cmd, slh, description="SPAdes")
                result["spades_mode"] = "standard"
            except PipelineStepError:
                spades_log_path = os.path.join(spades_dir, "spades.log")
                spades_log_text = ""
                if os.path.isfile(spades_log_path):
                    with open(spades_log_path) as spades_log_handle:
                        spades_log_text = spades_log_handle.read()
                if "Invalid kmer coverage histogram" not in spades_log_text:
                    raise
                log.warning(
                    "    Standard SPAdes could not model the local coverage; "
                    "retrying with --sc for non-uniform/low-depth reads."
                )
                shutil.rmtree(spades_dir)
                retry_cmd = build_spades_command(
                    tools["spades.py"], spades_dir, global_args.threads,
                    global_args.memory_per_sample, single_cell=True,
                    **spades_inputs,
                )
                try:
                    run_command(retry_cmd, slh, description="SPAdes --sc retry")
                    result["spades_mode"] = "single-cell-retry"
                except PipelineStepError:
                    retry_log_path = os.path.join(spades_dir, "spades.log")
                    retry_log_text = ""
                    if os.path.isfile(retry_log_path):
                        with open(retry_log_path) as retry_log_handle:
                            retry_log_text = retry_log_handle.read()
                    if "Invalid kmer coverage histogram" in retry_log_text:
                        msg = (
                            "Local assembly could not model the low/non-uniform "
                            "coverage in standard or --sc mode."
                        )
                        log.warning(f"  SKIP {sample_id}: {msg}")
                        result["status"] = "SKIP"
                        result["message"] = msg
                        return result
                    raise

            if not os.path.exists(local_asm) or os.path.getsize(local_asm) == 0:
                msg = "SPAdes produced no scaffolds."
                log.warning(f"  SKIP {sample_id}: {msg}")
                result["status"]  = "SKIP"
                result["message"] = msg
                return result
            log.info(
                f"    SPAdes: {count_fasta_sequences(local_asm)} scaffold(s), "
                f"{get_fasta_total_length(local_asm):,} bp total."
            )

            # ── STEP 5: Extract allele from local assembly ────────────────────
            log.info("  STEP 5/7: Extracting allele from local assembly...")
            local_db       = os.path.join(sample_out_dir, "local_db")
            local_blast    = os.path.join(sample_out_dir, "blast_step5.tsv")
            run_command(
                [tools["makeblastdb"], "-in", local_asm,
                 "-dbtype", "nucl", "-out", local_db],
                slh, description="makeblastdb (local)",
            )
            run_command(
                [tools["blastn"], "-query", global_args.bait,
                 "-db", local_db, "-out", local_blast,
                 "-outfmt", BLAST_FMT_STEP5,
                 "-evalue", "1e-10", "-num_threads", threads_str],
                slh, description="blastn (bait vs local)",
            )

            local_hits = []
            with open(local_blast) as bfh:
                for line in bfh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 13:
                        local_hits.append(rank_local_blast_hit(parts))

            if not local_hits:
                msg = "No BLAST hit in local assembly."
                log.warning(f"  SKIP {sample_id}: {msg}")
                result["status"]  = "SKIP"
                result["message"] = msg
                return result

            candidate_regions = collapse_overlapping_hits(local_hits)
            best_hit = candidate_regions[0]
            second_hit = candidate_regions[1] if len(candidate_regions) > 1 else None
            ambiguity  = summarize_hit_ambiguity(
                best_hit, second_hit, global_args.min_bitscore_margin,
            )

            best_contig = best_hit["sseqid"]
            best_start  = best_hit["sstart"]
            best_end    = best_hit["send"]
            best_pident = best_hit["pident"]
            result["best_identity"]  = best_pident
            result["bitscore_margin"] = ambiguity["bitscore_margin"]

            log.info(
                f"    Best: {best_contig} ({best_start}-{best_end}, "
                f"{best_hit['length']} bp, {best_pident:.1f}%, "
                f"bitscore={best_hit['bitscore']:.1f})"
            )
            if second_hit:
                log.info(
                    f"    2nd : {second_hit['sseqid']} "
                    f"({second_hit['length']} bp, {second_hit['pident']:.1f}%, "
                    f"bitscore={second_hit['bitscore']:.1f}, "
                    f"margin={ambiguity['bitscore_margin']:.1f})"
                )
            log.info(
                f"    Collapsed {len(local_hits)} redundant bait HSP(s) into "
                f"{len(candidate_regions)} distinct candidate region(s)."
            )

            allele_bed = os.path.join(sample_out_dir, "allele_region.bed")
            local_contig_lengths = get_fasta_lengths(local_asm)
            span = describe_full_query_span(
                best_hit, local_contig_lengths[best_contig]
            )
            span["continuation"] = find_span_continuation(
                span, candidate_regions, best_hit
            )
            rs, re_ = span["start"], span["end"]
            result["span_clipped_bp"]       = span["clipped_bp"]
            result["span_clipped_start_bp"] = span["clipped_start_bp"]
            result["span_clipped_end_bp"]   = span["clipped_end_bp"]
            if span["continuation"]:
                result["span_continuation_contig"] = span["continuation"]["sseqid"]
                result["span_continuation_bp"]     = span["continuation"]["recovered_bp"]
                result["span_continuation_overlap_bp"] = (
                    span["continuation"]["overlap_bp"]
                )
            if span["clipped_bp"]:
                log.warning(
                    f"    Projected span clipped by {span['clipped_bp']} bp at a "
                    f"contig boundary: the allele is truncated by that much."
                )
                if span["continuation"]:
                    cont = span["continuation"]
                    log.warning(
                        f"    {cont['recovered_bp']} of the lost bases align to "
                        f"{cont['sseqid']} ({cont['overlap_bp']} bp shared with "
                        f"the reported contig): the locus is split across the "
                        f"local assembly. Not joined -- reported only."
                    )
            if re_ <= rs:
                raise ValueError(
                    f"Invalid projected allele interval {best_contig}:{rs}-{re_}."
                )
            with open(allele_bed, "w") as fh:
                fh.write(f"{best_contig}\t{rs}\t{re_}\n")

            raw_seq = read_fasta_region(local_asm, best_contig, rs, re_)
            log.info(
                f"    Projected full query span: {best_contig}:{rs + 1}-{re_} "
                f"({len(raw_seq)} bp)."
            )
            if best_start > best_end:
                raw_seq = reverse_complement(raw_seq)
                log.info("    Applied reverse-complement (best hit on reverse strand).")

            final_fasta = os.path.join(
                sample_out_dir,
                f"{sample_id}_{global_args.locus}_reconstructed.fasta",
            )
            pubmlst_hdr = (
                f"{sample_id}_{global_args.locus}_reconstructed [sample={sample_id}] "
                f"[length={len(raw_seq)}] [identity={best_pident:.1f}%] "
                f"[strand={'-' if best_start > best_end else '+'}]"
            )
            write_single_fasta(pubmlst_hdr, raw_seq, final_fasta)
            result["allele_file"]   = final_fasta
            result["allele_length"] = len(raw_seq)

            # ── STEP 6: Validate allele ───────────────────────────────────────
            log.info("  STEP 6/7: Validating reconstructed allele...")
            val_out = os.path.join(sample_out_dir, "blast_validation.tsv")
            run_command(
                [tools["blastn"], "-query", final_fasta,
                 "-db", global_args.bait_db,
                 "-out", val_out,
                 "-outfmt",
                 "6 qseqid sseqid length pident qcovs bitscore qlen slen evalue",
                 "-max_target_seqs", "5", "-num_threads", threads_str],
                slh, description="blastn (validation)",
            )

            val_identity = val_qcov = 0.0
            nearest_allele = ""
            validation_hits = []
            with open(val_out) as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 9:
                        validation_hits.append(parts)

            if validation_hits:
                validation_hits.sort(
                    key=lambda parts: (
                        float(parts[5]), float(parts[4]),
                        float(parts[3]), int(parts[2]),
                    ),
                    reverse=True,
                )
                best_validation = validation_hits[0]
                nearest_allele = best_validation[1]
                val_identity = float(best_validation[3])
                val_qcov = float(best_validation[4])
                aligned_length = int(best_validation[2])
                query_length = int(best_validation[6])
                subject_length = int(best_validation[7])
                exact_known = (
                    val_identity == 100.0
                    and val_qcov == 100.0
                    and aligned_length == query_length == subject_length
                )
            else:
                exact_known = False

            result["best_identity"] = val_identity
            result["validation_coverage"] = val_qcov
            result["nearest_allele"] = nearest_allele
            result["exact_known_allele"] = exact_known

            curated_header = (
                f"{sample_id}_{global_args.locus}_reconstructed "
                f"[sample={sample_id}] [locus={global_args.locus}] "
                f"[length={len(raw_seq)}] [nearest={nearest_allele or 'none'}] "
                f"[identity={val_identity:.1f}] [coverage={val_qcov:.1f}] "
                f"[exact_known={'true' if exact_known else 'false'}] "
                f"[strand={'-' if best_start > best_end else '+'}]"
            )
            write_single_fasta(curated_header, raw_seq, final_fasta)

            log.info(
                f"    Validation: nearest={nearest_allele or 'none'}, "
                f"{val_identity:.1f}% id, {val_qcov:.1f}% qcov, "
                f"exact_known={'yes' if exact_known else 'no'}."
            )
            if (val_identity < global_args.min_identity
                    or val_qcov < global_args.min_coverage):
                msg = (
                    "Allele failed validation: "
                    f"identity={val_identity:.1f}% (minimum {global_args.min_identity:.1f}%), "
                    f"coverage={val_qcov:.1f}% (minimum {global_args.min_coverage:.1f}%)."
                )
                log.warning(f"  SKIP {sample_id}: {msg}")
                result["status"]  = "SKIP"
                result["message"] = msg
                return result

            remap = remap_reads_to_final_allele(
                tools, global_args.aligner_name, global_args.threads,
                final_fasta, r1_gz, r2_gz, s_gz,
                sample_out_dir, slh, global_args.use_progress,
                min_base_quality=global_args.min_base_quality,
                min_mapping_quality=global_args.min_mapping_quality,
                mixture_min_fraction=global_args.mixture_min_fraction,
                mixture_min_sites=global_args.mixture_min_sites,
                mixture_min_alt_depth=global_args.mixture_min_alt_depth,
            )
            result["remap_mean_depth"]   = remap["mean_depth"]
            result["remap_breadth_pct"]  = remap["breadth_pct"]
            result["remap_pct_bases_lt5"] = remap["pct_bases_lt5"]
            result["uncertain_base_count"] = remap.get("uncertain_base_count", 0)
            result["uncertain_base_fraction"] = remap.get("uncertain_base_fraction", 0.0)
            result["interior_uncertain_base_count"] = remap.get(
                "interior_uncertain_base_count", 0)
            result["interior_uncertain_base_fraction"] = remap.get(
                "interior_uncertain_base_fraction", 0.0)
            result["internal_zero_depth_positions"] = remap.get(
                "internal_zero_depth_positions", 0)
            result["terminal_margin_bp"] = remap.get("terminal_margin_bp", 0)
            result["mean_base_quality"] = remap["mean_base_quality"]
            result["mean_mapping_quality"] = remap["mean_mapping_quality"]
            result["strand_balance_pct"] = remap["strand_balance_pct"]
            result["candidate_mixed_sites"] = remap["candidate_mixed_sites"]
            result["mixed_site_count"] = remap["mixed_site_count"]
            result["strand_biased_sites"] = remap["strand_biased_sites"]
            result["reference_discordant_sites"] = remap["reference_discordant_sites"]
            result["max_alt_fraction"] = remap["max_alt_fraction"]
            result["median_mixed_fraction"] = remap["median_mixed_fraction"]
            result["mixture_detected"] = remap["mixture_detected"]

            if (remap["breadth_pct"] < global_args.min_remap_breadth
                    or remap["mean_depth"] < global_args.min_remap_depth):
                msg = (
                    f"Final allele failed remap support: "
                    f"breadth={remap['breadth_pct']:.1f}%, "
                    f"mean_depth={remap['mean_depth']:.1f}x"
                )
                log.warning(f"  SKIP {sample_id}: {msg}")
                result["status"]  = "SKIP"
                result["message"] = msg
                return result

            # ── STEP 7: QC scorecard ──────────────────────────────────────────
            log.info("  STEP 7/7: Assessing allele quality...")
            allele_seq = read_single_fasta_sequence(final_fasta)
            qc = assess_allele_quality(
                allele_seq=allele_seq,
                bait_profile=bait_profile,
                validation_identity=val_identity,
                validation_qcov=val_qcov,
                locus=global_args.locus,
                remap_metrics=remap,
                ambiguity_metrics=ambiguity,
                expect_cds=(not global_args.noncoding_locus),
                span_metrics=span,
                exact_known_allele=exact_known,
            )

            # Graph evidence for any interior scaffold placeholder.  The local
            # assembly graph is already on disk from STEP 4, so this costs no
            # extra alignment and turns an opaque held call into a legible one.
            placeholder = assess_placeholder_junctions(
                sequence=allele_seq,
                gfa_path=os.path.join(
                    spades_dir, "assembly_graph_after_simplification.gfa"
                ),
            )
            if placeholder["flags"]:
                qc["flags"] = list(qc["flags"]) + placeholder["flags"]

            result.update({
                "placeholder_runs": placeholder["runs"],
                "placeholder_bp":   placeholder["total_bp"],
                "placeholder_junction_support": placeholder["verdict"],
                "qc_confidence":    qc["confidence"],
                "sequence_confidence": qc["sequence_confidence"],
                "catalogue_status": qc["catalogue_status"],
                "nearest_allele_identity": qc["nearest_allele_identity"],
                "catalogue_review_recommended": qc["catalogue_review_recommended"],
                "qc_length_delta":  qc["length_delta"],
                "qc_length_zscore": qc["length_zscore"],
                "qc_gc_pct":        qc["gc_pct"],
                "qc_gc_deviation":  qc["gc_deviation"],
                "qc_n_count":       qc["n_count"],
                "qc_internal_stops": qc["internal_stops"],
                "qc_coding_frame_used": qc.get("coding_frame_used", ""),
                "qc_internal_stops_placeholder_closed": qc.get(
                    "internal_stops_placeholder_closed", 0
                ),
                "qc_length_mod3":   qc["length_mod3"],
                "qc_length_profile_n": qc.get("length_profile_n", 0),
                "qc_flags":         "; ".join(qc["flags"]) if qc["flags"] else "",
            })

            qc_text = format_qc_report(sample_id, global_args.locus, qc, bait_profile)
            qc_path = os.path.join(
                sample_out_dir,
                f"{sample_id}_{global_args.locus}_qc_report.txt",
            )
            write_qc_report_file(qc_text, qc_path)
            for line in qc_text.splitlines():
                log.info(f"  {line}")
            if not parallel_worker:
                print(qc_text)

        if global_args.cleanup:
            _cleanup_intermediates(sample_out_dir, final_fasta)

        elapsed = time.monotonic() - t_start
        result.update({
            "status":      "SUCCESS",
            "elapsed_sec": round(elapsed, 1),
            "message": (
                f"QC={qc['confidence']} | delta_len={qc['length_delta']:+d} bp | "
                f"nearest={nearest_allele or 'none'} | id={val_identity:.1f}% | "
                f"depth={remap['mean_depth']:.1f}x | "
                f"breadth={remap['breadth_pct']:.1f}% | "
                f"mixture={'yes' if remap['mixture_detected'] else 'no'} | "
                f"{elapsed:.1f}s"
            ),
        })
        log.info(f"  >> SUCCESS: {sample_id} -> {final_fasta} ({elapsed:.1f}s)")
        return result

    except PipelineStepError as e:
        elapsed = time.monotonic() - t_start
        result["elapsed_sec"] = round(elapsed, 1)
        result["message"]     = str(e)
        log.error(f"  >> FAIL {sample_id}: {e}")
        return result
    except Exception as e:
        elapsed = time.monotonic() - t_start
        result["elapsed_sec"] = round(elapsed, 1)
        result["message"]     = f"Unexpected error: {e}"
        log.error(f"  >> CRITICAL ERROR on {sample_id}: {e}")
        return result
    finally:
        if not result.get("elapsed_sec"):
            result["elapsed_sec"] = round(time.monotonic() - t_start, 1)
        result["workflow_status"] = result.get("status", "FAIL")
        result["result_disposition"] = classify_result_disposition(
            result["workflow_status"], result.get("qc_confidence", "")
        )
        if sample_log_handler is not None:
            log.removeHandler(sample_log_handler)
            sample_log_handler.close()


# ==============================================================================
# Intermediate file cleanup
# ==============================================================================


def _cleanup_intermediates(sample_out_dir: str, keep_file: str) -> None:
    """Remove bulky intermediate files while preserving the final allele FASTA,
    per-sample log, QC report, and BED/TSV summaries.

    Args:
        sample_out_dir: Per-sample output directory.
        keep_file:      Path to the final allele FASTA (BWA index files
                        created on this file during remap are also removed).
    """
    log.debug(f"  Cleaning up {sample_out_dir}...")

    for fname in [
        "full_map.bam", "full_map.sorted.bam", "full_map.sorted.bam.bai",
        "target_region.bam", "target_region_rescued.bam", "target_namesorted.bam",
        "target_read_names.txt", "assembly_input.fasta",
        "target_R1.fq.gz", "target_R2.fq.gz", "target_S.fq.gz",
        "allele_remap.bam", "allele_remap.sorted.bam", "allele_remap.sorted.bam.bai",
        "allele_remap.flagstat.txt", "allele_remap_bwa.log",
        "allele_remap_singleton_bwa.log",
    ]:
        p = os.path.join(sample_out_dir, fname)
        if os.path.exists(p):
            os.remove(p)

    for prefix in ("assembly_db", "local_db"):
        for ext in (
            ".nhr", ".nin", ".nsq", ".ndb", ".nog", ".nos", ".not", ".ntf", ".nto",
        ):
            p = os.path.join(sample_out_dir, prefix + ext)
            if os.path.exists(p):
                os.remove(p)

    for ext in (".amb", ".ann", ".bwt", ".pac", ".sa", ".bwt.2bit.64", ".0123"):
        p = os.path.join(sample_out_dir, "assembly_input.fasta" + ext)
        if os.path.exists(p):
            os.remove(p)

    if keep_file:
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa",
                    ".bwt.2bit.64", ".0123", ".fai"):
            p = keep_file + ext
            if os.path.exists(p):
                os.remove(p)
                log.debug(f"    Removed remap index: {os.path.basename(p)}")

    spades_dir = os.path.join(sample_out_dir, "spades_local")
    if os.path.isdir(spades_dir):
        shutil.rmtree(spades_dir)
