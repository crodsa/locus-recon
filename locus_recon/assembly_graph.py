"""Small, dependency-free helpers for resolving paths in SPAdes GFA graphs.

SPAdes can leave a repeated or multi-copy locus as several short FASTA
scaffolds while retaining the complete alternatives in its assembly graph.
This module reads that graph and enumerates terminal, non-cyclic paths without
guessing which path is biologically correct.  Candidate ranking and read-level
validation remain separate steps.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
import re
from typing import Dict, Iterable, List, Tuple

from .io import reverse_complement, wrap_fasta_sequence


OrientedNode = Tuple[str, str]
_CIGAR_OVERLAP = re.compile(r"^(\d+)M$")


@dataclass(frozen=True)
class GraphSegment:
    """A sequence segment and its SPAdes depth annotation."""

    sequence: str
    depth: float = 0.0


@dataclass(frozen=True)
class GraphPath:
    """A fully spelled path through an oriented GFA graph."""

    nodes: Tuple[OrientedNode, ...]
    sequence: str
    mean_depth: float


@dataclass
class AssemblyGraph:
    """Minimal sequence graph representation required for path spelling."""

    segments: Dict[str, GraphSegment]
    adjacency: Dict[OrientedNode, List[Tuple[OrientedNode, int]]]


def _opposite(node: OrientedNode) -> OrientedNode:
    return node[0], "-" if node[1] == "+" else "+"


def _oriented_sequence(segment: GraphSegment, orientation: str) -> str:
    if orientation == "+":
        return segment.sequence
    if orientation == "-":
        return reverse_complement(segment.sequence)
    raise ValueError(f"Invalid GFA orientation: {orientation!r}")


def parse_gfa(path: str) -> AssemblyGraph:
    """Parse ``S`` and ``L`` records from a GFA 1.x file.

    Link overlaps must use the simple ``<number>M`` form emitted by SPAdes.
    Reverse-complement links are added explicitly because GFA links are
    bidirected even though only one orientation is written.
    """

    segments: Dict[str, GraphSegment] = {}
    raw_links: List[Tuple[OrientedNode, OrientedNode, int]] = []

    with open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                if len(fields) < 3:
                    raise ValueError(f"Malformed GFA segment at line {line_no}.")
                segment_id, sequence = fields[1], fields[2].upper()
                if segment_id in segments:
                    raise ValueError(
                        f"Duplicate GFA segment ID {segment_id!r} at line {line_no}."
                    )
                if sequence == "*":
                    raise ValueError(
                        f"GFA segment {segment_id!r} has no embedded sequence."
                    )
                depth = 0.0
                for tag in fields[3:]:
                    if tag.startswith("DP:f:"):
                        depth = float(tag[5:])
                        break
                segments[segment_id] = GraphSegment(sequence, depth)
            elif fields[0] == "L":
                if len(fields) < 6:
                    raise ValueError(f"Malformed GFA link at line {line_no}.")
                match = _CIGAR_OVERLAP.fullmatch(fields[5])
                if not match:
                    raise ValueError(
                        f"Unsupported GFA overlap {fields[5]!r} at line {line_no}; "
                        "expected a simple <number>M overlap."
                    )
                source = (fields[1], fields[2])
                target = (fields[3], fields[4])
                raw_links.append((source, target, int(match.group(1))))

    if not segments:
        raise ValueError(f"GFA contains no sequence segments: {path}")

    adjacency: Dict[OrientedNode, List[Tuple[OrientedNode, int]]] = {
        (segment_id, orientation): []
        for segment_id in segments
        for orientation in ("+", "-")
    }
    for source, target, overlap in raw_links:
        if source[0] not in segments or target[0] not in segments:
            raise ValueError(
                f"GFA link refers to an undefined segment: {source[0]} -> {target[0]}"
            )
        if source[1] not in {"+", "-"} or target[1] not in {"+", "-"}:
            raise ValueError(f"Invalid GFA link orientation: {source} -> {target}")
        adjacency[source].append((target, overlap))
        adjacency[_opposite(target)].append((_opposite(source), overlap))

    for edges in adjacency.values():
        edges.sort(key=lambda item: item[0])
    return AssemblyGraph(segments=segments, adjacency=adjacency)


def spell_path(graph: AssemblyGraph, nodes: Iterable[OrientedNode]) -> GraphPath:
    """Spell a sequence from an oriented path and validate every graph edge."""

    node_tuple = tuple(nodes)
    if not node_tuple:
        raise ValueError("Cannot spell an empty graph path.")

    first = node_tuple[0]
    sequence = _oriented_sequence(graph.segments[first[0]], first[1])
    weighted_depth = graph.segments[first[0]].depth * len(sequence)
    depth_bases = len(sequence)

    for source, target in zip(node_tuple, node_tuple[1:]):
        matching = [
            overlap
            for neighbor, overlap in graph.adjacency.get(source, [])
            if neighbor == target
        ]
        if not matching:
            raise ValueError(f"Nodes are not linked in GFA: {source} -> {target}")
        overlap = matching[0]
        target_sequence = _oriented_sequence(
            graph.segments[target[0]], target[1]
        )
        if overlap > len(sequence) or overlap > len(target_sequence):
            raise ValueError(
                f"Overlap {overlap} exceeds a segment length at {source} -> {target}."
            )
        if overlap and sequence[-overlap:] != target_sequence[:overlap]:
            raise ValueError(
                f"Sequence overlap does not match at {source} -> {target}."
            )
        added = target_sequence[overlap:]
        sequence += added
        weighted_depth += graph.segments[target[0]].depth * len(added)
        depth_bases += len(added)

    return GraphPath(
        nodes=node_tuple,
        sequence=sequence,
        mean_depth=(weighted_depth / depth_bases if depth_bases else 0.0),
    )


def enumerate_terminal_paths(
    graph: AssemblyGraph,
    *,
    max_paths: int = 10000,
    max_nodes: int = 1000,
) -> List[GraphPath]:
    """Enumerate all acyclic source-to-sink paths in a small assembly graph.

    Reverse-complement duplicates are collapsed by sequence.  The function is
    intentionally bounded: complex whole-genome graphs should be reduced to a
    locus component before calling it.
    """

    indegree = {node: 0 for node in graph.adjacency}
    for edges in graph.adjacency.values():
        for target, _overlap in edges:
            indegree[target] = indegree.get(target, 0) + 1
    sources = sorted(
        node
        for node, degree in indegree.items()
        if degree == 0 and graph.adjacency.get(node)
    )
    if not sources:
        raise ValueError("GFA graph has no terminal source; it may be cyclic.")

    node_paths: List[Tuple[OrientedNode, ...]] = []

    def visit(
        node: OrientedNode,
        current: Tuple[OrientedNode, ...],
        visited_segments: frozenset[str],
    ) -> None:
        if len(node_paths) >= max_paths:
            raise ValueError(
                f"GFA path count exceeded the safety limit ({max_paths})."
            )
        if len(current) >= max_nodes:
            raise ValueError(
                f"GFA path length exceeded the safety limit ({max_nodes} nodes)."
            )
        next_path = current + (node,)
        next_visited = visited_segments | {node[0]}
        all_outgoing = graph.adjacency.get(node, [])
        outgoing = [
            target
            for target, _overlap in all_outgoing
            if target[0] not in next_visited
        ]
        if not outgoing:
            # Do not misreport a path that only stopped because it entered a
            # cycle as a biological terminal path.
            if not all_outgoing:
                node_paths.append(next_path)
            return
        for target in outgoing:
            visit(target, next_path, next_visited)

    for source in sources:
        visit(source, (), frozenset())

    if not node_paths:
        raise ValueError("GFA graph has no acyclic source-to-sink path.")

    unique: Dict[str, GraphPath] = {}
    for nodes in node_paths:
        path = spell_path(graph, nodes)
        canonical = min(path.sequence, reverse_complement(path.sequence))
        existing = unique.get(canonical)
        if existing is None or path.nodes < existing.nodes:
            unique[canonical] = path
    return sorted(
        unique.values(),
        key=lambda path: (-len(path.sequence), path.nodes),
    )


def write_paths_fasta(
    paths: Iterable[GraphPath],
    output_path: str,
    *,
    prefix: str = "graph_path",
) -> int:
    """Write graph paths as FASTA records with path and depth metadata."""

    count = 0
    with open(output_path, "w") as handle:
        for count, path in enumerate(paths, 1):
            node_text = ",".join(
                f"{segment_id}{orientation}"
                for segment_id, orientation in path.nodes
            )
            handle.write(
                f">{prefix}_{count} length={len(path.sequence)} "
                f"mean_graph_depth={path.mean_depth:.2f} nodes={node_text}\n"
            )
            handle.write(wrap_fasta_sequence(path.sequence) + "\n")
    return count


# ==============================================================================
# Scaffold-placeholder junctions
# ==============================================================================


def find_interior_n_runs(sequence: str, min_run_bp: int = 1) -> List[Tuple[int, int]]:
    """Locate runs of N that lie strictly inside a sequence.

    A local assembler that scaffolds across a gap it cannot spell writes a run
    of N.  Terminal runs are trimmed elsewhere and carry no junction, so only
    interior runs are returned.

    Args:
        sequence:   Candidate sequence (uppercase).
        min_run_bp: Ignore runs shorter than this.

    Returns:
        List of (start, end) 1-based inclusive coordinates, left to right.
    """
    runs: List[Tuple[int, int]] = []
    for match in re.finditer(r"N+", sequence.upper()):
        start, end = match.start() + 1, match.end()
        if start == 1 or end == len(sequence):
            continue
        if end - start + 1 >= min_run_bp:
            runs.append((start, end))
    return runs


def junction_support(
    path_sequences: Iterable[str],
    left_flank: str,
    right_flank: str,
) -> dict:
    """Test whether two flanks of a placeholder are contiguous in the graph.

    The two anchors are searched in every supplied path spelling and in its
    reverse complement.  A path that carries both anchors with nothing between
    them spells the join that the scaffolder could not spell, which is graph
    evidence for the junction even though no read spans the placeholder.  A
    path that carries both anchors with sequence between them offers a
    different spelling, and the junction is then ambiguous rather than
    supported.

    Args:
        path_sequences: Spelled graph paths (see enumerate_terminal_paths).
        left_flank:     Sequence immediately 5' of the placeholder.
        right_flank:    Sequence immediately 3' of the placeholder.

    Returns:
        Dict with keys: assessed, supported, paths_with_both_flanks, gaps_bp
        (sorted distinct distances between the anchors), min_gap_bp.
    """
    if not left_flank or not right_flank:
        return {
            "assessed": False, "supported": False,
            "paths_with_both_flanks": 0, "gaps_bp": [], "min_gap_bp": None,
        }

    gaps: List[int] = []
    for sequence in path_sequences:
        for spelling in (sequence, reverse_complement(sequence)):
            left_at = spelling.find(left_flank)
            if left_at < 0:
                continue
            right_at = spelling.find(right_flank, left_at + len(left_flank) - 1)
            if right_at < 0:
                continue
            gaps.append(right_at - (left_at + len(left_flank)))
    distinct = sorted(set(gaps))
    return {
        "assessed": True,
        "supported": bool(gaps) and all(gap == 0 for gap in gaps),
        "paths_with_both_flanks": len(gaps),
        "gaps_bp": distinct,
        "min_gap_bp": distinct[0] if distinct else None,
    }


def assess_placeholder_junctions(
    sequence: str,
    gfa_path: str,
    flank_bp: int = 40,
    min_run_bp: int = 10,
    max_paths: int = 10000,
    max_nodes: int = 1000,
) -> dict:
    """Describe every interior placeholder in a candidate against the graph.

    This uses the local assembly graph the pipeline already produced, so it
    adds no alignment work.  The verdict is reported, never scored: the
    confidence tier stays a statement about read support for reported bases,
    and a graph-supported junction is still a junction no read spans.

    Args:
        sequence:   Reconstructed candidate sequence.
        gfa_path:   Local assembly graph (GFA) from the same sample.
        flank_bp:   Anchor length taken either side of each placeholder.
        min_run_bp: Ignore placeholder runs shorter than this.
        max_paths:  Safety limit passed to enumerate_terminal_paths.
        max_nodes:  Safety limit passed to enumerate_terminal_paths.

    Returns:
        Dict with keys: runs, total_bp, verdict
        (GRAPH_SUPPORTED / AMBIGUOUS / NOT_SUPPORTED / NOT_ASSESSED),
        details (per-run dicts) and flags (report strings).
    """
    runs = find_interior_n_runs(sequence, min_run_bp=min_run_bp)
    total_bp = sum(end - start + 1 for start, end in runs)
    summary = {
        "runs": len(runs), "total_bp": total_bp,
        "verdict": "NOT_ASSESSED", "details": [], "flags": [],
    }
    if not runs:
        summary["verdict"] = "NO_PLACEHOLDER"
        return summary
    if not gfa_path or not os.path.isfile(gfa_path):
        summary["flags"].append(
            f"PLACEHOLDER_JUNCTION_NOT_ASSESSED ({len(runs)} interior "
            f"placeholder run(s), {total_bp} bp; no local assembly graph found)"
        )
        return summary

    try:
        graph = parse_gfa(gfa_path)
        paths = [
            path.sequence
            for path in enumerate_terminal_paths(
                graph, max_paths=max_paths, max_nodes=max_nodes
            )
        ]
    except (OSError, ValueError) as exc:
        summary["flags"].append(
            f"PLACEHOLDER_JUNCTION_NOT_ASSESSED ({len(runs)} interior "
            f"placeholder run(s), {total_bp} bp; graph unreadable: {exc})"
        )
        return summary

    verdicts = []
    for start, end in runs:
        left = sequence[max(0, start - 1 - flank_bp) : start - 1]
        right = sequence[end : end + flank_bp]
        support = junction_support(paths, left, right)
        support.update({"start": start, "end": end, "length_bp": end - start + 1})
        if not support["assessed"] or support["paths_with_both_flanks"] == 0:
            support["verdict"] = "NOT_SUPPORTED"
        elif support["supported"]:
            support["verdict"] = "GRAPH_SUPPORTED"
        else:
            support["verdict"] = "AMBIGUOUS"
        verdicts.append(support["verdict"])
        summary["details"].append(support)

    if all(v == "GRAPH_SUPPORTED" for v in verdicts):
        summary["verdict"] = "GRAPH_SUPPORTED"
    elif any(v == "AMBIGUOUS" for v in verdicts):
        summary["verdict"] = "AMBIGUOUS"
    else:
        summary["verdict"] = "NOT_SUPPORTED"

    spans = "; ".join(
        f"{item['start']}-{item['end']} ({item['length_bp']} bp, "
        f"{item['paths_with_both_flanks']} graph path(s) carry both flanks"
        + (
            f", gap {item['min_gap_bp']} bp"
            if item["min_gap_bp"] not in (None, 0) else ""
        )
        + ")"
        for item in summary["details"]
    )
    if summary["verdict"] == "GRAPH_SUPPORTED":
        summary["flags"].append(
            f"PLACEHOLDER_JUNCTION_GRAPH_SUPPORTED ({spans}; the flanks are "
            "contiguous in the local assembly graph, so the placeholder is a "
            "scaffolding gap rather than missing sequence -- no read spans it)"
        )
    elif summary["verdict"] == "AMBIGUOUS":
        summary["flags"].append(
            f"PLACEHOLDER_JUNCTION_AMBIGUOUS ({spans}; the graph offers more "
            "than one spelling across the placeholder)"
        )
    else:
        summary["flags"].append(
            f"PLACEHOLDER_JUNCTION_NOT_SUPPORTED ({spans}; no graph path "
            "carries both flanks)"
        )
    return summary
