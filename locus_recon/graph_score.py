"""Competitive read scoring for enumerated assembly-graph paths.

Enumerating paths through a local assembly graph produces alternative
spellings of an unresolved locus, and their graph depth alone does not say
which spelling the reads support.  This module maps the recruited reads
against every retained path in a single index, so the paths compete for each
read, and scores each path by the sequence that uniquely anchored reads
actually cover.

A rank is not a call.  Paths that share most of their sequence receive
near-identical scores by construction, and a path can only be ranked against
the alternatives that were enumerated.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional

from .mapping import run_bwa_to_sorted_bam
from .utils import log, run_command


def parse_idxstats(text: str) -> Dict[str, dict]:
    """Parse ``samtools idxstats`` output into per-reference counts.

    Args:
        text: Raw idxstats output (reference, length, mapped, unmapped).

    Returns:
        Mapping of reference name to {'length', 'mapped_reads'}; the unmapped
        ``*`` row is dropped.
    """
    stats: Dict[str, dict] = {}
    for line in text.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4 or fields[0] == "*":
            continue
        stats[fields[0]] = {
            "length": int(fields[1]),
            "mapped_reads": int(fields[2]),
        }
    return stats


def aggregate_depth(text: str) -> Dict[str, dict]:
    """Aggregate ``samtools depth -a`` output per reference.

    Args:
        text: Raw depth output (reference, position, depth).

    Returns:
        Mapping of reference name to {'positions', 'covered_positions',
        'total_depth'}.
    """
    aggregate: Dict[str, dict] = {}
    for line in text.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        entry = aggregate.setdefault(
            fields[0],
            {"positions": 0, "covered_positions": 0, "total_depth": 0},
        )
        depth = int(fields[2])
        entry["positions"] += 1
        entry["total_depth"] += depth
        if depth > 0:
            entry["covered_positions"] += 1
    return aggregate


def score_paths(
    idxstats: Dict[str, dict],
    depth_aggregate: Dict[str, dict],
    graph_depths: Optional[Dict[str, float]] = None,
) -> List[dict]:
    """Rank paths by the sequence that uniquely anchored reads cover.

    Reads that map equally well to several paths are given a mapping quality
    of zero by the aligner and are excluded from the depth used here, so the
    score reflects sequence that distinguishes one path from the others rather
    than sequence they share.

    Args:
        idxstats:        Output of parse_idxstats().
        depth_aggregate: Output of aggregate_depth() on a mapping-quality
                         filtered depth run.
        graph_depths:    Optional mean graph depth per candidate, carried
                         through for reporting.

    Returns:
        List of row dicts, best first, each carrying 'rank'.
    """
    rows: List[dict] = []
    for reference, stat in idxstats.items():
        depth = depth_aggregate.get(reference, {})
        positions = depth.get("positions") or stat["length"]
        covered = depth.get("covered_positions", 0)
        total_depth = depth.get("total_depth", 0)
        rows.append({
            "candidate_id": reference,
            "length": stat["length"],
            "mean_graph_depth": (graph_depths or {}).get(reference),
            "mapped_reads": stat["mapped_reads"],
            "unique_mean_depth": total_depth / positions if positions else 0.0,
            "unique_breadth_pct": 100.0 * covered / positions if positions else 0.0,
            "unsupported_bp": max(0, positions - covered),
        })
    rows.sort(key=lambda row: (
        -row["unique_breadth_pct"],
        -row["unique_mean_depth"],
        -row["length"],
        row["candidate_id"],
    ))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
    return rows


def write_scored_summary(rows: List[dict], output_path: str) -> None:
    """Write the ranked path summary as TSV."""
    header = [
        "rank", "candidate_id", "length", "mean_graph_depth", "mapped_reads",
        "unique_mean_depth", "unique_breadth_pct", "unsupported_bp", "nodes",
    ]
    with open(output_path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            graph_depth = row.get("mean_graph_depth")
            handle.write("\t".join([
                str(row["rank"]),
                row["candidate_id"],
                str(row["length"]),
                "" if graph_depth is None else f"{graph_depth:.3f}",
                str(row["mapped_reads"]),
                f"{row['unique_mean_depth']:.3f}",
                f"{row['unique_breadth_pct']:.2f}",
                str(row["unsupported_bp"]),
                row.get("nodes", ""),
            ]) + "\n")


def ranking_is_informative(rows: List[dict]) -> bool:
    """Whether any path carries sequence that uniquely placed reads cover.

    When every path scores zero the enumerated alternatives are
    indistinguishable at this read length: every read fits several of them
    equally well and is discarded by the mapping-quality filter.  The ordering
    is then arbitrary and must not be read as a preference.

    Args:
        rows: Rows from score_paths().

    Returns:
        False when no path has uniquely anchored coverage.
    """
    return any(row["unique_breadth_pct"] > 0 for row in rows)


def competitive_read_scores(
    tools: dict,
    aligner_name: str,
    threads: int,
    paths_fasta: str,
    r1: str,
    r2: Optional[str],
    out_dir: str,
    log_handle,
    use_progress: bool = False,
    min_mapping_quality: int = 20,
    graph_depths: Optional[Dict[str, float]] = None,
) -> List[dict]:
    """Map reads against all paths at once and score each path.

    Args:
        tools:               Tool path mapping from check_dependencies().
        aligner_name:        'bwa-mem2' or 'bwa'.
        threads:             Thread count.
        paths_fasta:         FASTA of the retained graph paths.
        r1:                  Forward reads.
        r2:                  Reverse reads, or None for single-end input.
        out_dir:             Destination for the BAM and text outputs.
        log_handle:          Open handle receiving external tool output.
        use_progress:        Whether tqdm is available for the mapping stream.
        min_mapping_quality: Depth is counted only from alignments at or above
                             this mapping quality, which excludes reads that
                             map equally well to several paths.
        graph_depths:        Optional mean graph depth per candidate.

    Returns:
        Ranked rows from score_paths().
    """
    os.makedirs(out_dir, exist_ok=True)
    sorted_bam = os.path.join(out_dir, "graph_paths.sorted.bam")
    idxstats_path = os.path.join(out_dir, "graph_paths.idxstats.txt")
    depth_path = os.path.join(out_dir, "graph_paths.depth.txt")
    aligner_log = os.path.join(out_dir, "graph_paths_bwa.log")

    run_command(
        [tools[aligner_name], "index", paths_fasta],
        log_handle, description=f"{aligner_name} index (graph paths)",
    )
    run_bwa_to_sorted_bam(
        tools, aligner_name, threads, paths_fasta, r1, r2,
        "graph-paths", aligner_log, sorted_bam, log_handle, use_progress,
    )
    run_command(
        [tools["samtools"], "index", sorted_bam],
        log_handle, description="samtools index (graph paths)",
    )
    with open(idxstats_path, "w") as handle:
        run_command(
            [tools["samtools"], "idxstats", sorted_bam],
            handle, description="samtools idxstats (graph paths)",
        )
    with open(depth_path, "w") as handle:
        run_command(
            [tools["samtools"], "depth", "-a", "-Q", str(min_mapping_quality),
             sorted_bam],
            handle, description="samtools depth (graph paths)",
        )
    rows = score_paths(
        parse_idxstats(open(idxstats_path).read()),
        aggregate_depth(open(depth_path).read()),
        graph_depths=graph_depths,
    )
    log.info(
        f"  Scored {len(rows)} graph path(s) by uniquely anchored coverage "
        f"(MAPQ >= {min_mapping_quality})."
    )
    if rows and not ranking_is_informative(rows):
        log.warning(
            "  No path carries uniquely placed reads: the enumerated paths "
            "cannot be told apart at this read length, and their order is "
            "arbitrary."
        )
    return rows
