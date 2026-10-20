"""
Read mapping helpers: BWA pipeline, samtools utilities, mate rescue,
and remap-based allele validation.
"""

import os
import shutil
import subprocess
import sys
from typing import Optional

from .utils import log, run_command, PipelineStepError
from .io import fastq_has_records, read_single_fasta_sequence
from .support import collect_per_base_support


def run_bwa_to_sorted_bam(
    tools: dict,
    aligner_name: str,
    threads: int,
    assembly: str,
    r1: str,
    r2: Optional[str],
    locus: str,
    bwa_log_path: str,
    sorted_bam_path: str,
    samtools_log_handle,
    use_progress: bool,
) -> None:
    """Execute ``aligner mem | [tqdm] | samtools sort`` in one stream.

    Sorting directly from SAM avoids writing a large unsorted BAM and then
    reading it back, which materially reduces temporary disk I/O for batches.

    Args:
        tools:              Tool path mapping from check_dependencies().
        aligner_name:       'bwa-mem2' or 'bwa'.
        threads:            Thread count.
        assembly:           Reference FASTA path.
        r1:                 Forward reads path.
        r2:                 Reverse reads path.
        locus:              Locus name (used in progress bar label).
        bwa_log_path:       File to capture aligner stderr.
        sorted_bam_path:    Coordinate-sorted output BAM path.
        samtools_log_handle: Open file handle for samtools stderr.
        use_progress:       Whether to pipe through tqdm.
    """
    sort_threads = max(1, threads // 4)
    aligner_threads = max(1, threads - sort_threads)
    bwa_cmd = [
        tools[aligner_name], "mem", "-t", str(aligner_threads), assembly, r1,
    ]
    if r2:
        bwa_cmd.append(r2)
    samtools_cmd = [
        tools["samtools"], "sort", "-@", str(sort_threads),
        "-o", sorted_bam_path, "-",
    ]

    log.info(f"  Mapping reads with {aligner_name}...")
    log.debug(f"  BWA cmd: {' '.join(bwa_cmd)}")

    try:
        with open(bwa_log_path, "w") as bwa_log:
            p_bwa = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=bwa_log)

            if use_progress and "tqdm" in tools:
                tqdm_cmd = [
                    tools["tqdm"], "--bytes", "--desc",
                    f"  Mapping {locus}", "--unit", "B", "--unit_scale",
                ]
                p_tqdm = subprocess.Popen(
                    tqdm_cmd, stdin=p_bwa.stdout,
                    stdout=subprocess.PIPE, stderr=sys.stderr,
                )
                samtools_stdin = p_tqdm.stdout
            else:
                p_tqdm = None
                samtools_stdin = p_bwa.stdout

            p_sam = subprocess.Popen(
                samtools_cmd, stdin=samtools_stdin,
                stderr=samtools_log_handle,
            )
            p_bwa.stdout.close()
            if p_tqdm:
                p_tqdm.stdout.close()

            rc_bwa = p_bwa.wait()
            rc_tqdm = p_tqdm.wait() if p_tqdm else 0
            rc_sam = p_sam.wait()

        if rc_bwa != 0:
            raise PipelineStepError(
                f"{aligner_name} failed (rc={rc_bwa}). See {bwa_log_path}"
            )
        if rc_tqdm != 0:
            raise PipelineStepError(f"tqdm stream failed (rc={rc_tqdm}).")
        if rc_sam != 0:
            raise PipelineStepError(f"samtools sort failed (rc={rc_sam}).")

        log.debug(
            f"  BAM: {sorted_bam_path} ({os.path.getsize(sorted_bam_path):,} bytes)"
        )
        log.info("  Mapping completed.")

    except PipelineStepError:
        raise
    except Exception as e:
        raise PipelineStepError(f"Mapping pipeline error: {e}")


def collect_mapping_stats(
    samtools_path: str,
    bam_path: str,
    stats_path: str,
    sample_log_handle,
    threads: int,
) -> dict:
    """Run samtools flagstat and parse key metrics.

    Args:
        samtools_path:     Path to samtools binary.
        bam_path:          Input BAM.
        stats_path:        File to write raw flagstat output.
        sample_log_handle: Open file handle for stderr.
        threads:           Thread count.

    Returns:
        Dict with keys: total_reads, mapped_reads, mapped_pct, properly_paired.
    """
    log.debug(f"  Running flagstat on {bam_path}...")
    try:
        result = subprocess.run(
            [samtools_path, "flagstat", "-@", str(threads), bam_path],
            capture_output=True, text=True, check=True,
        )
        with open(stats_path, "w") as f:
            f.write(result.stdout)
        import re
        stats: dict = {}
        for line in result.stdout.splitlines():
            if "in total" in line:
                stats["total_reads"] = int(line.split()[0])
            elif "mapped (" in line:
                stats["mapped_reads"] = int(line.split()[0])
                m = re.search(r"\(([\d.]+)%", line)
                stats["mapped_pct"] = float(m.group(1)) if m else 0.0
            elif "properly paired" in line:
                stats["properly_paired"] = int(line.split()[0])
        log.debug(
            f"  Flagstat: {stats.get('mapped_reads','?')}/"
            f"{stats.get('total_reads','?')} mapped "
            f"({stats.get('mapped_pct','?')}%)"
        )
        return stats
    except Exception as e:
        log.warning(f"  Could not collect mapping stats: {e}")
        return {}


def write_read_name_list_from_bam(
    samtools_path: str,
    bam_path: str,
    read_names_path: str,
) -> int:
    """Write unique QNAME entries from a BAM to a text file; return count.

    Streams via subprocess pipe (samtools view | cut | sort -u) to avoid
    loading the full SAM into Python memory for large or paralog-rich loci.

    Args:
        samtools_path:   Path to samtools binary.
        bam_path:        Input BAM.
        read_names_path: Output text file of unique read names.

    Returns:
        Number of unique read names written.
    """
    with open(read_names_path, "w") as fh:
        p_view = subprocess.Popen(
            [samtools_path, "view", bam_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p_cut = subprocess.Popen(
            ["cut", "-f1"], stdin=p_view.stdout, stdout=subprocess.PIPE,
        )
        p_sort = subprocess.Popen(
            ["sort", "-u"], stdin=p_cut.stdout, stdout=fh,
        )
        p_view.stdout.close()
        p_cut.stdout.close()
        rc_sort = p_sort.wait()
        rc_cut = p_cut.wait()
        _, view_stderr = p_view.communicate()

    if p_view.returncode != 0 or rc_cut != 0 or rc_sort != 0:
        detail = view_stderr.decode(errors="replace").strip()
        raise PipelineStepError(
            "Could not extract target read names "
            f"(samtools={p_view.returncode}, cut={rc_cut}, sort={rc_sort}). "
            f"{detail}"
        )

    count = 0
    with open(read_names_path) as fh:
        for _ in fh:
            count += 1
    return count


def select_alignments_by_read_names(
    tools: dict,
    threads: int,
    bam_in: str,
    read_names_path: str,
    bam_out: str,
    log_handle,
) -> None:
    """Extract all alignments whose QNAME appears in read_names_path (mate rescue).

    Args:
        tools:           Tool path mapping.
        threads:         Thread count.
        bam_in:          Sorted, indexed input BAM.
        read_names_path: Text file of read names (one per line).
        bam_out:         Output BAM with rescued mates.
        log_handle:      Open file handle for command output.
    """
    run_command(
        [
            tools["samtools"], "view",
            "-@", str(threads), "-N", read_names_path, "-b", "-o", bam_out, bam_in,
        ],
        log_handle,
        description="samtools view -N (mate rescue)",
    )


def remap_reads_to_final_allele(
    tools: dict,
    aligner_name: str,
    threads: int,
    allele_fasta: str,
    r1: str,
    r2: str,
    singleton: str,
    out_dir: str,
    log_handle,
    use_progress: bool,
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
    mixture_min_fraction: float = 0.05,
    mixture_min_sites: int = 2,
    mixture_min_alt_depth: int = 2,
) -> dict:
    """Remap extracted reads to the reconstructed allele; compute depth/breadth.

    Evidence-based validation layer: the allele is not accepted on BLAST
    similarity alone but must also be supported by remapped reads.

    Args:
        tools:        Tool path mapping.
        aligner_name: 'bwa-mem2' or 'bwa'.
        threads:      Thread count.
        allele_fasta: Final allele FASTA path.
        r1:           Forward reads (extracted from target region).
        r2:           Reverse reads.
        singleton:    Singleton reads.
        out_dir:      Directory for intermediate remap files.
        log_handle:   Open file handle for command output.
        use_progress: Whether to show tqdm progress bar.
        min_base_quality: Minimum base quality included in per-base evidence.
        min_mapping_quality: Minimum mapping quality included in per-base evidence.
        mixture_min_fraction: Minimum non-reference fraction at a mixed site.
        mixture_min_sites: Minimum bidirectionally supported sites for a mixture call.
        mixture_min_alt_depth: Minimum alternative observations at a mixed site.

    Returns:
        Dict with depth/breadth/coverage metrics.
    """
    allele_sorted = os.path.join(out_dir, "allele_remap.sorted.bam")
    bwa_log       = os.path.join(out_dir, "allele_remap_bwa.log")
    support_path  = os.path.join(out_dir, "allele_remap.support.tsv")
    stats_path    = os.path.join(out_dir, "allele_remap.flagstat.txt")

    run_command(
        [tools[aligner_name], "index", allele_fasta],
        log_handle,
        description=f"{aligner_name} index (final allele)",
    )
    mapping_parts = []
    if fastq_has_records(r1) and fastq_has_records(r2):
        paired_bam = os.path.join(out_dir, "allele_remap.paired.bam")
        run_bwa_to_sorted_bam(
            tools, aligner_name, threads, allele_fasta,
            r1, r2, "final-allele-pairs", bwa_log, paired_bam, log_handle,
            use_progress,
        )
        mapping_parts.append(paired_bam)

    if fastq_has_records(singleton):
        singleton_bam = os.path.join(out_dir, "allele_remap.singleton.bam")
        singleton_log = os.path.join(out_dir, "allele_remap_singleton_bwa.log")
        run_bwa_to_sorted_bam(
            tools, aligner_name, threads, allele_fasta,
            singleton, None, "final-allele-singletons", singleton_log,
            singleton_bam, log_handle, use_progress,
        )
        mapping_parts.append(singleton_bam)

    if not mapping_parts:
        raise PipelineStepError("No extracted reads were available for allele remapping.")
    if len(mapping_parts) == 1:
        shutil.move(mapping_parts[0], allele_sorted)
    else:
        run_command(
            [tools["samtools"], "merge", "-@", str(threads), "-f",
             allele_sorted, *mapping_parts],
            log_handle, description="samtools merge (paired + singleton remap)",
        )
        for part in mapping_parts:
            if os.path.exists(part):
                os.remove(part)
    run_command(
        [tools["samtools"], "index", "-@", str(threads), allele_sorted],
        log_handle, description="samtools index (allele remap)",
    )

    allele_len = len(read_single_fasta_sequence(allele_fasta))
    support = collect_per_base_support(
        tools["samtools"], allele_fasta, allele_sorted, support_path,
        allele_length=allele_len,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
        mixture_min_fraction=mixture_min_fraction,
        mixture_min_sites=mixture_min_sites,
        mixture_min_alt_depth=mixture_min_alt_depth,
    )

    flagstat = collect_mapping_stats(
        tools["samtools"], allele_sorted, stats_path, log_handle, threads,
    )

    support.update({
        "mapped_reads":     flagstat.get("mapped_reads", 0),
        "mapped_pct":       flagstat.get("mapped_pct", 0.0),
        "properly_paired":  flagstat.get("properly_paired", 0),
        "singleton_fastq":  singleton,
        "depth_file":       support_path,
        "bam":              allele_sorted,
    })
    return support
