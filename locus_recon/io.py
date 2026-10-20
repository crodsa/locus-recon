"""FASTA, FASTQ, samplesheet, and structured-output helpers."""

import csv
import gzip
import logging
import os
import re
from datetime import datetime
from typing import List, Optional, Tuple

from . import VERSION, PROGRAM

log = logging.getLogger(PROGRAM)

SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


# ==============================================================================
# FASTA utilities
# ==============================================================================


def count_fasta_sequences(fasta_path: str) -> int:
    """Count '>' headers in a FASTA file.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Number of sequences.
    """
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def get_fasta_total_length(fasta_path: str) -> int:
    """Sum all sequence lengths in a FASTA file.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Total number of nucleotide characters.
    """
    total = 0
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def read_fasta_sequences(fasta_path: str) -> List[Tuple[str, str]]:
    """Parse a multi-FASTA into (header, sequence) tuples. Sequences are uppercased.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        List of (header_id, sequence) tuples where header_id is the first word
        after '>'.
    """
    sequences: List[Tuple[str, str]] = []
    header: Optional[str] = None
    seq_parts: List[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    sequences.append((header, "".join(seq_parts).upper()))
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        sequences.append((header, "".join(seq_parts).upper()))
    return sequences


def read_single_fasta_sequence(fasta_path: str) -> str:
    """Read the first sequence from a FASTA file, uppercased.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        First sequence string, or empty string if file is empty.
    """
    seqs = read_fasta_sequences(fasta_path)
    return seqs[0][1] if seqs else ""


def read_fasta_region(
    fasta_path: str,
    sequence_id: str,
    start: int,
    end: int,
) -> str:
    """Read a 0-based, half-open region from a named FASTA sequence.

    The local SPAdes assembly is small, so an in-process lookup avoids a
    separate ``seqtk`` runtime dependency without becoming a memory bottleneck.
    """
    if start < 0 or end < start:
        raise ValueError(f"Invalid FASTA interval: {sequence_id}:{start}-{end}")
    for header, sequence in read_fasta_sequences(fasta_path):
        if header == sequence_id:
            if end > len(sequence):
                raise ValueError(
                    f"FASTA interval {sequence_id}:{start}-{end} exceeds "
                    f"sequence length {len(sequence)}."
                )
            return sequence[start:end]
    raise ValueError(f"Sequence {sequence_id!r} was not found in {fasta_path}.")


def get_fasta_lengths(fasta_path: str) -> dict:
    """Return a mapping of header ID to sequence length for a FASTA file.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Dict of {header_id: sequence_length}.
    """
    return {header: len(seq) for header, seq in read_fasta_sequences(fasta_path)}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (IUPAC aware).

    Args:
        seq: DNA sequence string (any case).

    Returns:
        Reverse-complemented sequence, same case as input.
    """
    table = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn",
    )
    return seq.translate(table)[::-1]


def wrap_fasta_sequence(seq: str, width: int = 80) -> str:
    """Wrap a sequence string to FASTA line width.

    Args:
        seq:   Sequence string.
        width: Line width (default 80).

    Returns:
        Wrapped sequence with newlines.
    """
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def merge_intervals(intervals: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
    """Merge overlapping or adjacent BED-style intervals per contig.

    Args:
        intervals: List of (contig_id, start, end) tuples (0-based, half-open).

    Returns:
        Sorted, merged list of (contig_id, start, end) tuples.
    """
    by_contig: dict = {}
    for contig, start, end in intervals:
        by_contig.setdefault(contig, []).append((start, end))

    merged: List[Tuple[str, int, int]] = []
    for contig, parts in by_contig.items():
        parts.sort()
        cur_start, cur_end = parts[0]
        for start, end in parts[1:]:
            if start <= cur_end + 1:
                cur_end = max(cur_end, end)
            else:
                merged.append((contig, cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((contig, cur_start, cur_end))
    return sorted(merged)


def write_single_fasta(header: str, sequence: str, output_path: str) -> None:
    """Write a single FASTA record to disk.

    Args:
        header:      Full header line (without '>').
        sequence:    Sequence string.
        output_path: Destination file path.
    """
    with open(output_path, "w") as fout:
        fout.write(">" + header + "\n")
        fout.write(wrap_fasta_sequence(sequence) + "\n")


def fastq_has_records(path: str) -> bool:
    """Return whether a plain or gzip-compressed FASTQ contains a record.

    File size is not a valid emptiness check for gzip files because an empty
    gzip stream still has a header and trailer.  Reading one text line is both
    cheap and sufficient for the pipeline's paired/singleton routing.
    """
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    opener = gzip.open if path.lower().endswith(".gz") else open
    try:
        with opener(path, "rt") as handle:
            return bool(handle.readline())
    except (OSError, UnicodeDecodeError):
        return False


# ==============================================================================
# Samplesheet parsing
# ==============================================================================


def parse_samplesheet(path: str) -> List[Tuple[str, str, str, str]]:
    """Read a 4-column TAB-separated samplesheet.

    Expected columns (no header, lines starting with '#' are comments):
        sample_id   assembly_path   r1_path   r2_path

    Args:
        path: Path to the samplesheet TSV.

    Returns:
        List of (sample_id, assembly, r1, r2) tuples.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: On malformed rows, unsafe IDs, or duplicate sample IDs.
    """
    samples: List[Tuple[str, str, str, str]] = []
    seen_ids = set()
    sheet_dir = os.path.dirname(os.path.abspath(path))

    def normalize_input_path(value: str) -> str:
        expanded = os.path.expanduser(value)
        if os.path.isabs(expanded):
            return os.path.normpath(expanded)
        return os.path.normpath(os.path.join(sheet_dir, expanded))

    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for line_no, row in enumerate(reader, 1):
            if not row or row[0].strip().startswith("#"):
                continue
            row = [field.strip() for field in row]
            if not samples and row and row[0].lower() == "sample_id":
                continue
            if len(row) != 4:
                raise ValueError(
                    f"Samplesheet line {line_no}: expected 4 tab-separated "
                    f"columns, found {len(row)}."
                )
            sample_id, assembly, r1, r2 = row
            if not SAFE_SAMPLE_ID.fullmatch(sample_id):
                raise ValueError(
                    f"Samplesheet line {line_no}: unsafe sample ID {sample_id!r}. "
                    "Use letters, numbers, dots, underscores, or hyphens."
                )
            if sample_id in seen_ids:
                raise ValueError(
                    f"Samplesheet line {line_no}: duplicate sample ID {sample_id!r}."
                )
            if not assembly or not r1 or not r2:
                raise ValueError(
                    f"Samplesheet line {line_no}: input paths cannot be empty."
                )
            seen_ids.add(sample_id)
            samples.append((
                sample_id,
                normalize_input_path(assembly),
                normalize_input_path(r1),
                normalize_input_path(r2),
            ))
    return samples


# ==============================================================================
# Structured output
# ==============================================================================


def classify_result_disposition(workflow_status: str, qc_confidence: str) -> str:
    """Map computational completion to a conservative downstream disposition.

    ``workflow_status`` answers whether the reconstruction workflow completed;
    it does not make the resulting sequence safe for allele catalogues.  Only a
    successful HIGH-confidence reconstruction is accepted by default.
    """
    if workflow_status != "SUCCESS":
        return "HOLD"
    return {
        "HIGH": "PASS",
        "MEDIUM": "REVIEW",
        "LOW": "REVIEW",
        "SUSPECT": "HOLD",
    }.get(qc_confidence, "HOLD")


def write_batch_report(
    results: List[dict],
    report_path: str,
    locus: str,
    run_date: Optional[str] = None,
) -> None:
    """Write a comprehensive TSV summary with QC metrics, version, and run date.

    Columns follow a stable order compatible with downstream parsing tools.
    Empty cells are written as '' for SKIP/FAIL rows that have no QC data.

    Args:
        results:     List of per-sample result dicts from process_sample().
        report_path: Output file path.
        locus:       Locus name (e.g. 'aroE').
        run_date:    ISO-format date string; defaults to today.
    """
    run_date = run_date or datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    header = [
        "sample_id", "status", "workflow_status", "result_disposition",
        "locus", "tool_version", "run_date",
        "allele_length", "nearest_allele", "exact_known_allele",
        "best_identity_pct", "validation_coverage_pct",
        "contigs_hit", "target_regions", "extracted_read_names", "spades_mode",
        "mapped_reads",
        "qc_confidence", "sequence_confidence",
        "catalogue_status", "nearest_allele_identity_pct",
        "catalogue_review_recommended",
        "qc_length_delta", "qc_length_zscore",
        "qc_gc_pct", "qc_gc_deviation",
        "qc_n_count", "qc_internal_stops",
        "qc_internal_stops_placeholder_closed", "qc_coding_frame_used",
        "qc_length_mod3", "qc_length_profile_n",
        "placeholder_runs", "placeholder_bp", "placeholder_junction_support",
        "qc_flags",
        "remap_mean_depth", "remap_breadth_pct", "remap_pct_bases_lt5",
        "uncertain_base_count", "uncertain_base_fraction",
        "interior_uncertain_base_count", "interior_uncertain_base_fraction",
        "internal_zero_depth_positions", "terminal_margin_bp",
        "span_clipped_bp", "span_clipped_start_bp", "span_clipped_end_bp",
        "span_continuation_contig", "span_continuation_bp",
        "span_continuation_overlap_bp",
        "mean_base_quality", "mean_mapping_quality", "strand_balance_pct",
        "candidate_mixed_sites", "mixed_site_count", "strand_biased_sites",
        "reference_discordant_sites", "max_alt_fraction", "median_mixed_fraction",
        "mixture_detected",
        "bitscore_margin",
        "elapsed_sec", "allele_file", "message",
    ]

    with open(report_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        for r in results:
            has_qc = bool(r.get("qc_confidence"))
            workflow_status = r.get("workflow_status", r.get("status", "FAIL"))
            disposition = r.get("result_disposition") or classify_result_disposition(
                workflow_status, r.get("qc_confidence", "")
            )
            writer.writerow([
                r["sample_id"],
                workflow_status,  # `status`: same value as `workflow_status`
                workflow_status,
                disposition,
                locus,
                VERSION,
                run_date,
                r.get("allele_length", 0),
                r.get("nearest_allele", ""),
                str(bool(r.get("exact_known_allele", False))).lower() if has_qc else "",
                f"{r['best_identity']:.1f}" if r.get("best_identity") else "",
                f"{r.get('validation_coverage', 0):.1f}" if has_qc else "",
                r.get("contigs_hit", 0),
                r.get("target_regions", 0),
                r.get("extracted_read_names", 0),
                r.get("spades_mode", ""),
                r.get("mapped_reads", 0),
                r.get("qc_confidence", ""),
                r.get("sequence_confidence", r.get("qc_confidence", "")),
                r.get("catalogue_status", "") if has_qc else "",
                f"{r.get('nearest_allele_identity', 0):.1f}" if has_qc else "",
                str(bool(r.get("catalogue_review_recommended", False))).lower()
                if has_qc else "",
                r.get("qc_length_delta", "")    if has_qc else "",
                f"{r.get('qc_length_zscore', 0):.2f}" if has_qc else "",
                f"{r.get('qc_gc_pct', 0):.1f}"        if has_qc else "",
                f"{r.get('qc_gc_deviation', 0):.1f}"  if has_qc else "",
                r.get("qc_n_count", "")         if has_qc else "",
                r.get("qc_internal_stops", "")  if has_qc else "",
                r.get("qc_internal_stops_placeholder_closed", "") if has_qc else "",
                r.get("qc_coding_frame_used", "") if has_qc else "",
                r.get("qc_length_mod3", "")     if has_qc else "",
                r.get("qc_length_profile_n", "") if has_qc else "",
                r.get("placeholder_runs", "")   if has_qc else "",
                r.get("placeholder_bp", "")     if has_qc else "",
                r.get("placeholder_junction_support", "") if has_qc else "",
                r.get("qc_flags", ""),
                f"{r.get('remap_mean_depth', 0):.1f}"     if has_qc else "",
                f"{r.get('remap_breadth_pct', 0):.1f}"    if has_qc else "",
                f"{r.get('remap_pct_bases_lt5', 0):.1f}"  if has_qc else "",
                r.get("uncertain_base_count", "") if has_qc else "",
                f"{r.get('uncertain_base_fraction', 0):.4f}" if has_qc else "",
                r.get("interior_uncertain_base_count", "") if has_qc else "",
                f"{r.get('interior_uncertain_base_fraction', 0):.4f}" if has_qc else "",
                r.get("internal_zero_depth_positions", "") if has_qc else "",
                r.get("terminal_margin_bp", "") if has_qc else "",
                r.get("span_clipped_bp", 0),
                r.get("span_clipped_start_bp", 0),
                r.get("span_clipped_end_bp", 0),
                r.get("span_continuation_contig", ""),
                r.get("span_continuation_bp", 0),
                r.get("span_continuation_overlap_bp", 0),
                f"{r.get('mean_base_quality', 0):.1f}" if has_qc else "",
                f"{r.get('mean_mapping_quality', 0):.1f}" if has_qc else "",
                f"{r.get('strand_balance_pct', 0):.1f}" if has_qc else "",
                r.get("candidate_mixed_sites", "") if has_qc else "",
                r.get("mixed_site_count", "") if has_qc else "",
                r.get("strand_biased_sites", "") if has_qc else "",
                r.get("reference_discordant_sites", "") if has_qc else "",
                f"{r.get('max_alt_fraction', 0):.4f}" if has_qc else "",
                f"{r.get('median_mixed_fraction', 0):.4f}" if has_qc else "",
                str(bool(r.get("mixture_detected", False))).lower() if has_qc else "",
                f"{r.get('bitscore_margin', 0):.1f}"      if has_qc else "",
                r.get("elapsed_sec", 0),
                r.get("allele_file", ""),
                r.get("message", ""),
            ])

    log.info(f"Batch report written to: {report_path}")


def write_allele_catalogs(
    results: List[dict],
    output_dir: str,
    locus: str,
) -> dict:
    """Write accepted, review, hold, and all-candidate FASTA collections.

    The accepted collection holds PASS results only; every technically
    successful reconstruction, whatever its disposition, is in the
    all-candidate collection.

    Args:
        results:    List of per-sample result dicts.
        output_dir: Directory in which to write the catalog.
        locus:      Locus name (used in the output file name).

    Returns a mapping from catalogue role to ``path`` and ``count``.
    """
    paths = {
        "accepted": os.path.join(output_dir, f"{locus}_accepted_alleles.fasta"),
        "review": os.path.join(
            output_dir, f"{locus}_review_required_candidates.fasta"
        ),
        "hold": os.path.join(output_dir, f"{locus}_hold_candidates.fasta"),
        "all_candidates": os.path.join(
            output_dir, f"{locus}_reconstructed_candidates.fasta"
        ),
    }
    contents = {role: [] for role in paths}

    for r in results:
        workflow_status = r.get("workflow_status", r.get("status", "FAIL"))
        if workflow_status != "SUCCESS" or not r.get("allele_file"):
            continue
        allele_file = r["allele_file"]
        if not os.path.exists(allele_file):
            log.warning(f"Candidate file missing for {r['sample_id']}: {allele_file}")
            continue
        with open(allele_file) as handle:
            content = handle.read()
        if content and not content.endswith("\n"):
            content += "\n"

        disposition = r.get("result_disposition") or classify_result_disposition(
            workflow_status, r.get("qc_confidence", "")
        )
        contents["all_candidates"].append(content)
        if disposition == "PASS":
            contents["accepted"].append(content)
        elif disposition == "REVIEW":
            contents["review"].append(content)
        else:
            contents["hold"].append(content)

    output = {}
    for role, path in paths.items():
        with open(path, "w") as handle:
            handle.write("".join(contents[role]))
        output[role] = {"path": path, "count": len(contents[role])}
        log.info(
            "Sequence collection written (%d candidates): %s",
            len(contents[role]), path,
        )
    return output
