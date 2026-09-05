"""Small, dependency-free helpers for resolving paths in SPAdes GFA graphs.

SPAdes can leave a repeated or multi-copy locus as several short FASTA
scaffolds while retaining the complete alternatives in its assembly graph.
This module reads that graph and enumerates terminal, non-cyclic paths without
guessing which path is biologically correct.  Candidate ranking and read-level
validation remain separate steps.
"""

from __future__ import annotations

from dataclasses import dataclass
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
