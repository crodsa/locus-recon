"""Graph/depth consensus for locus copy-number interpretation.

Depth and graph topology answer different questions. The depth estimator
measures continuous pile-up over an assembly representation; a graph exposes
how many distinct flanking contexts the assembler retained. This module keeps
both observations intact and combines them through explicit decision states.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
from pathlib import Path
import subprocess
import tempfile
from typing import Iterable

from .assembly_graph import AssemblyGraph, OrientedNode
from .depth_ratio import DepthRatioResult
from .io import read_fasta_sequences, reverse_complement, wrap_fasta_sequence

__all__ = [
    "GraphTargetHit",
    "GraphTargetDiscoveryResult",
    "GraphContextResult",
    "CopyNumberConsensusResult",
    "infer_graph_contexts",
    "discover_graph_target_hits",
    "reconcile_copy_number",
]


@dataclass(frozen=True)
class GraphTargetHit:
    """One bait-to-GFA-segment alignment retained for context inference."""

    segment_id: str
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    identity: float
    bitscore: float

    @property
    def orientation(self) -> str:
        return "+" if self.subject_start <= self.subject_end else "-"

    @property
    def query_interval(self) -> tuple[int, int]:
        return min(self.query_start, self.query_end), max(self.query_start, self.query_end)


@dataclass
class GraphTargetDiscoveryResult:
    """Bait selection, BLAST parameters, and retained graph target hits."""

    bait_id: str
    bait_length: int
    hits: list[GraphTargetHit]
    min_identity: float
    min_hit_bp: int
    blastn_command: list[str]

    def as_dict(self) -> dict:
        result = asdict(self)
        result["hits"] = [asdict(hit) for hit in self.hits]
        return result


@dataclass
class GraphContextResult:
    """Graph-supported context count and its identifiability state."""

    target_segments: list[str]
    bait_length: int
    bait_covered_bp: int
    bait_coverage_fraction: float
    left_context_count: int
    right_context_count: int
    context_count: int | None
    context_lower_bound: int | None
    status: str
    flags: list[str] = field(default_factory=list)
    left_context_signatures: list[str] = field(default_factory=list)
    right_context_signatures: list[str] = field(default_factory=list)
    min_context_bp: int = 0
    max_context_nodes: int = 0
    max_context_paths: int = 0
    max_context_bp: int = 0

    def as_dict(self) -> dict:
        return asdict(self)

    def as_row(self) -> dict:
        result = self.as_dict()
        for key in (
            "target_segments",
            "flags",
            "left_context_signatures",
            "right_context_signatures",
        ):
            result[key] = ";".join(result[key])
        return result


@dataclass
class CopyNumberConsensusResult:
    """Additive interpretation of unchanged depth and graph evidence."""

    copy_number_call: float | int | None
    copy_number_kind: str
    copy_number_method: str
    copy_number_status: str
    copy_number_lower_bound: float | int | None = None
    copy_number_upper_bound: float | int | None = None
    consensus_flags: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)

    def as_row(self) -> dict:
        result = self.as_dict()
        result["consensus_flags"] = ";".join(self.consensus_flags)
        return result


def _query_union_bp(hits: Iterable[GraphTargetHit], bait_length: int) -> int:
    intervals: list[tuple[int, int]] = []
    for hit in hits:
        start, end = hit.query_interval
        start = max(1, min(start, bait_length))
        end = max(1, min(end, bait_length))
        if start <= end:
            intervals.append((start, end))
    if not intervals:
        return 0
    intervals.sort()
    covered = 0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start + 1
            current_start, current_end = start, end
    return covered + current_end - current_start + 1


def _select_bait(path: str | Path, bait_id: str | None) -> tuple[str, str]:
    records = [
        (header.split()[0], sequence.upper())
        for header, sequence in read_fasta_sequences(str(path))
    ]
    if not records:
        raise ValueError(f"bait FASTA contains no records: {path}")
    if bait_id is None:
        if len(records) != 1:
            available = ", ".join(identifier for identifier, _sequence in records)
            raise ValueError(
                "--bait-id is required when the bait FASTA contains multiple "
                f"records ({available})"
            )
        return records[0]
    matches = [record for record in records if record[0] == bait_id]
    if len(matches) != 1:
        raise ValueError(f"bait identifier was not found exactly once: {bait_id}")
    return matches[0]


def discover_graph_target_hits(
    graph: AssemblyGraph,
    bait_fasta: str | Path,
    *,
    bait_id: str | None = None,
    blastn: str = "blastn",
    min_identity: float = 80.0,
    min_hit_bp: int | None = None,
    min_hit_fraction: float = 0.10,
    temporary_parent: str | Path | None = None,
) -> GraphTargetDiscoveryResult:
    """Align one selected bait to embedded GFA segments with BLASTN."""

    if not 0 < min_identity <= 100:
        raise ValueError("min_identity must be in (0, 100]")
    if not 0 < min_hit_fraction <= 1:
        raise ValueError("min_hit_fraction must be in (0, 1]")
    selected_id, bait_sequence = _select_bait(bait_fasta, bait_id)
    effective_min_hit = (
        min_hit_bp
        if min_hit_bp is not None
        else max(200, int(min_hit_fraction * len(bait_sequence)))
    )
    if effective_min_hit < 1:
        raise ValueError("min_hit_bp must be positive")

    with tempfile.TemporaryDirectory(
        prefix="locus-recon-copy-number-",
        dir=str(temporary_parent) if temporary_parent else None,
    ) as temporary:
        temporary_path = Path(temporary)
        selected_bait = temporary_path / "bait.fasta"
        selected_bait.write_text(
            f">{selected_id}\n{wrap_fasta_sequence(bait_sequence)}\n"
        )
        segments_fasta = temporary_path / "graph_segments.fasta"
        with segments_fasta.open("w") as handle:
            for segment_id in sorted(graph.segments):
                handle.write(f">{segment_id}\n")
                handle.write(wrap_fasta_sequence(graph.segments[segment_id].sequence) + "\n")

        outfmt_fields = "sseqid pident length qstart qend sstart send bitscore"
        command = [
            blastn,
            "-query", str(selected_bait),
            "-subject", str(segments_fasta),
            "-outfmt", f"6 {outfmt_fields}",
            "-perc_identity", str(min_identity),
            "-max_target_seqs", "10000",
        ]
        try:
            completed = subprocess.run(command, capture_output=True, text=True)
        except FileNotFoundError as exc:
            raise RuntimeError(f"BLASTN executable was not found: {blastn}") from exc
        if completed.returncode != 0:
            raise RuntimeError(f"BLASTN graph discovery failed:\n{completed.stderr[:1000]}")

    candidates: list[GraphTargetHit] = []
    for line in completed.stdout.splitlines():
        row_fields = line.rstrip().split("\t")
        if len(row_fields) != 8:
            raise RuntimeError(f"malformed BLASTN graph-discovery row: {line[:200]}")
        segment_id, identity, length, qstart, qend, sstart, send, bitscore = row_fields
        if segment_id not in graph.segments:
            raise RuntimeError(f"BLASTN returned an unknown graph segment: {segment_id}")
        if int(length) < effective_min_hit or float(identity) < min_identity:
            continue
        candidates.append(GraphTargetHit(
            segment_id=segment_id,
            query_start=int(qstart),
            query_end=int(qend),
            subject_start=int(sstart),
            subject_end=int(send),
            identity=float(identity),
            bitscore=float(bitscore),
        ))

    retained: list[GraphTargetHit] = []
    for hit in sorted(
        candidates,
        key=lambda item: (
            -item.bitscore,
            -(item.query_interval[1] - item.query_interval[0] + 1),
            -item.identity,
            item.segment_id,
            item.subject_start,
            item.subject_end,
        ),
    ):
        hit_subject = sorted((hit.subject_start, hit.subject_end))
        overlaps = any(
            old.segment_id == hit.segment_id
            and hit_subject[0] <= max(old.subject_start, old.subject_end)
            and min(old.subject_start, old.subject_end) <= hit_subject[1]
            for old in retained
        )
        if not overlaps:
            retained.append(hit)
    retained.sort(key=lambda item: (item.query_interval, item.segment_id))

    public_command = [
        blastn,
        "-query", "<selected-bait.fasta>",
        "-subject", "<gfa-segments.fasta>",
        "-outfmt", f"6 {outfmt_fields}",
        "-perc_identity", str(min_identity),
        "-max_target_seqs", "10000",
    ]
    return GraphTargetDiscoveryResult(
        bait_id=selected_id,
        bait_length=len(bait_sequence),
        hits=retained,
        min_identity=min_identity,
        min_hit_bp=effective_min_hit,
        blastn_command=public_command,
    )


def _oriented_sequence(graph: AssemblyGraph, node: OrientedNode) -> str:
    sequence = graph.segments[node[0]].sequence
    return sequence if node[1] == "+" else reverse_complement(sequence)


def _canonical_signature(sequence: str) -> str:
    canonical = min(sequence, reverse_complement(sequence))
    return hashlib.sha256(canonical.encode("ascii")).hexdigest()


def _context_anchor(hit: GraphTargetHit, *, query_side: str) -> OrientedNode:
    """Return the oriented graph end pointing away from one bait end."""

    if query_side not in {"left", "right"}:
        raise ValueError(f"invalid query side: {query_side}")
    if query_side == "left":
        orientation = "-" if hit.orientation == "+" else "+"
    else:
        orientation = "+" if hit.orientation == "+" else "-"
    return hit.segment_id, orientation


def _target_flank_context(
    graph: AssemblyGraph,
    hit: GraphTargetHit,
    *,
    query_side: str,
) -> tuple[OrientedNode, str]:
    """Return the outward graph end and sequence between bait and that end."""

    anchor = _context_anchor(hit, query_side=query_side)
    sequence = graph.segments[hit.segment_id].sequence
    subject_low = min(hit.subject_start, hit.subject_end)
    subject_high = max(hit.subject_start, hit.subject_end)
    if subject_low < 1 or subject_high > len(sequence):
        raise ValueError(
            f"target coordinates exceed segment {hit.segment_id}: "
            f"{subject_low}-{subject_high} of {len(sequence)}"
        )
    if anchor[1] == "+":
        flank = sequence[subject_high:]
    else:
        flank = reverse_complement(sequence[:subject_low - 1])
    return anchor, flank


def _collect_contexts(
    graph: AssemblyGraph,
    anchor: OrientedNode,
    *,
    target_segments: frozenset[str],
    min_context_bp: int,
    max_context_nodes: int,
    max_context_paths: int,
    max_context_bp: int,
    target_flank: str = "",
) -> tuple[list[str], set[str]]:
    signatures: set[str] = set()
    flags: set[str] = set()
    stack: list[tuple[OrientedNode, str, frozenset[str], int]] = [
        (anchor, "", frozenset(target_segments), 0)
    ]
    explored_edges = 0

    while stack:
        node, context, visited, nodes_used = stack.pop()
        outgoing = graph.adjacency.get(node, [])
        if not outgoing:
            if len(context) >= min_context_bp:
                signatures.add(_canonical_signature(target_flank + context[:min_context_bp]))
            elif not context and len(target_flank) >= min_context_bp:
                signatures.add(_canonical_signature(target_flank[:min_context_bp]))
            else:
                flags.add("INSUFFICIENT_CONTEXT_SPAN" if context else "NO_CONTEXT_EDGE")
            continue

        for neighbor, overlap in outgoing:
            explored_edges += 1
            if explored_edges > max_context_paths:
                flags.add("MAX_CONTEXT_PATHS_REACHED")
                stack.clear()
                break
            segment_id = neighbor[0]
            if segment_id in visited:
                flags.add("CYCLE_ENCOUNTERED")
                continue
            next_nodes = nodes_used + (0 if segment_id in target_segments else 1)
            if next_nodes > max_context_nodes:
                flags.add("MAX_CONTEXT_NODES_REACHED")
                continue

            added = ""
            if segment_id not in target_segments:
                sequence = _oriented_sequence(graph, neighbor)
                if overlap > len(sequence):
                    flags.add("INVALID_GRAPH_OVERLAP")
                    continue
                added = sequence[overlap:]
            next_context = context + added
            required_context = min(len(next_context), min_context_bp)
            if len(target_flank) + required_context > max_context_bp:
                flags.add("MAX_CONTEXT_BP_REACHED")
                continue
            if len(next_context) >= min_context_bp:
                signatures.add(
                    _canonical_signature(target_flank + next_context[:min_context_bp])
                )
                continue
            stack.append((neighbor, next_context, visited | {segment_id}, next_nodes))

    return sorted(signatures), flags


def _component_max_overlap(
    graph: AssemblyGraph,
    target_segments: Iterable[str],
) -> int:
    """Return the largest overlap in graph components containing target hits."""

    pending = [
        (segment_id, orientation)
        for segment_id in target_segments
        for orientation in ("+", "-")
    ]
    visited: set[OrientedNode] = set()
    maximum = 0
    while pending:
        node = pending.pop()
        if node in visited:
            continue
        visited.add(node)
        for neighbor, overlap in graph.adjacency.get(node, []):
            maximum = max(maximum, overlap)
            if neighbor not in visited:
                pending.append(neighbor)
    return maximum


TRAVERSAL_INCOMPLETENESS_FLAGS = frozenset({
    "MAX_CONTEXT_NODES_REACHED",
    "MAX_CONTEXT_PATHS_REACHED",
    "MAX_CONTEXT_BP_REACHED",
    "INVALID_GRAPH_OVERLAP",
    "CYCLE_ENCOUNTERED",
})


def context_status(
    left_count: int,
    right_count: int,
    flags: Iterable[str],
) -> tuple[str, int | None, int | None]:
    """Classify a pair of side context counts into (status, count, lower bound).

    Every flag in TRAVERSAL_INCOMPLETENESS_FLAGS marks a walk that was pruned
    or aborted before it was exhausted, so the signature sets are subsets of
    the reachable contexts and each side count is a lower bound rather than a
    count. Exactness is claimed only for an exhausted traversal.
    """

    flags = set(flags)
    incomplete = bool(TRAVERSAL_INCOMPLETENESS_FLAGS & flags)

    if left_count and right_count and left_count == right_count and not incomplete:
        status = "MATCHED_CONTEXT_COUNT"
        context_count = left_count
        lower_bound = left_count
    elif incomplete:
        # Matching counts from an incomplete traversal can be an artefact of
        # where the walk stopped rather than a property of the locus: both
        # sides can stop short at the same number. The count is retained as a
        # lower bound -- distal complexity does not erase it -- but exactness
        # is not claimed, and no upper bound is asserted from a subset.
        status = "TRAVERSAL_LIMIT_REACHED"
        context_count = None
        positive = [count for count in (left_count, right_count) if count]
        lower_bound = min(positive) if positive else None
    elif left_count and right_count:
        status = "ASYMMETRIC_CONTEXTS"
        context_count = None
        lower_bound = min(left_count, right_count)
    elif left_count or right_count:
        status = "ONE_SIDED_LOWER_BOUND"
        context_count = None
        lower_bound = max(left_count, right_count)
    else:
        status = "GRAPH_INDETERMINATE"
        context_count = None
        lower_bound = None
    return status, context_count, lower_bound


def infer_graph_contexts(
    graph: AssemblyGraph,
    hits: list[GraphTargetHit],
    *,
    bait_length: int,
    min_query_coverage: float = 0.90,
    min_context_bp: int | None = None,
    max_context_nodes: int = 8,
    max_context_paths: int = 64,
    max_context_bp: int = 5000,
) -> GraphContextResult:
    """Count distinct graph contexts at both query-defined locus ends."""

    if bait_length < 1:
        raise ValueError("bait_length must be positive")
    if not 0 < min_query_coverage <= 1:
        raise ValueError("min_query_coverage must be in (0, 1]")
    for name, value in (
        ("max_context_nodes", max_context_nodes),
        ("max_context_paths", max_context_paths),
        ("max_context_bp", max_context_bp),
    ):
        if value < 1:
            raise ValueError(f"{name} must be positive")
    for hit in hits:
        if hit.segment_id not in graph.segments:
            raise ValueError(f"target segment is absent from graph: {hit.segment_id}")

    target_segments = list(dict.fromkeys(hit.segment_id for hit in hits))
    max_overlap = _component_max_overlap(graph, target_segments)
    effective_min_context = (
        min_context_bp
        if min_context_bp is not None
        else max(100, 2 * max_overlap)
    )
    if effective_min_context < 1:
        raise ValueError("min_context_bp must be positive")
    if effective_min_context > max_context_bp:
        raise ValueError("min_context_bp cannot exceed max_context_bp")

    ordered_hits = sorted(
        hits,
        key=lambda hit: (
            hit.query_interval[0], hit.query_interval[1], hit.segment_id, -hit.bitscore,
        ),
    )
    target_segments = list(dict.fromkeys(hit.segment_id for hit in ordered_hits))
    covered_bp = _query_union_bp(ordered_hits, bait_length)
    coverage = covered_bp / bait_length
    base = {
        "target_segments": target_segments,
        "bait_length": bait_length,
        "bait_covered_bp": covered_bp,
        "bait_coverage_fraction": round(coverage, 6),
        "min_context_bp": effective_min_context,
        "max_context_nodes": max_context_nodes,
        "max_context_paths": max_context_paths,
        "max_context_bp": max_context_bp,
    }
    if not ordered_hits or coverage < min_query_coverage:
        return GraphContextResult(
            **base,
            left_context_count=0,
            right_context_count=0,
            context_count=None,
            context_lower_bound=None,
            status="INSUFFICIENT_TARGET_COVERAGE",
            flags=["BAIT_QUERY_UNION_BELOW_THRESHOLD"],
        )

    left_hit = min(
        ordered_hits,
        key=lambda hit: (hit.query_interval[0], -hit.bitscore, hit.segment_id),
    )
    right_hit = max(
        ordered_hits,
        key=lambda hit: (hit.query_interval[1], hit.bitscore, hit.segment_id),
    )
    target_set = frozenset(target_segments)
    left_anchor, left_flank = _target_flank_context(
        graph, left_hit, query_side="left"
    )
    right_anchor, right_flank = _target_flank_context(
        graph, right_hit, query_side="right"
    )

    def infer_side(anchor: OrientedNode, flank: str) -> tuple[list[str], set[str]]:
        if len(flank) >= max_context_bp:
            return (
                [_canonical_signature(flank[:effective_min_context])],
                {"CONTEXT_RESOLVED_WITHIN_TARGET_SEGMENT"},
            )
        return _collect_contexts(
            graph,
            anchor,
            target_segments=target_set,
            min_context_bp=effective_min_context,
            max_context_nodes=max_context_nodes,
            max_context_paths=max_context_paths,
            max_context_bp=max_context_bp,
            target_flank=flank,
        )

    left_signatures, left_flags = infer_side(
        left_anchor,
        left_flank,
    )
    right_signatures, right_flags = infer_side(
        right_anchor,
        right_flank,
    )
    flags = sorted(left_flags | right_flags)
    left_count, right_count = len(left_signatures), len(right_signatures)
    status, context_count, lower_bound = context_status(left_count, right_count, flags)

    return GraphContextResult(
        **base,
        left_context_count=left_count,
        right_context_count=right_count,
        context_count=context_count,
        context_lower_bound=lower_bound,
        status=status,
        flags=flags,
        left_context_signatures=left_signatures,
        right_context_signatures=right_signatures,
    )


def reconcile_copy_number(
    depth: DepthRatioResult,
    graph: GraphContextResult | None,
) -> CopyNumberConsensusResult:
    """Reconcile graph contexts with the unchanged depth result."""

    flags: list[str] = []
    graph_exact = (
        graph is not None
        and graph.status == "MATCHED_CONTEXT_COUNT"
        and graph.context_count is not None
    )
    graph_lower_bound = (
        graph is not None
        and graph.status in {
            "ONE_SIDED_LOWER_BOUND",
            "ASYMMETRIC_CONTEXTS",
            "TRAVERSAL_LIMIT_REACHED",
        }
        and graph.context_lower_bound is not None
    )
    depth_usable = (
        depth.reliable
        and depth.depth_call in {"MULTICOPY_DEPTH", "SINGLE_COPY_COMPATIBLE"}
    )
    dosage_usable = (
        depth_usable
        and depth.dosage_estimate is not None
        and depth.dosage_status != "NOT_ESTIMATED"
    )
    # Low mapping uniqueness removes depth and never adds it, so a multicopy
    # verdict it accompanies still stands, but its dosage is only a floor.
    depth_floor = (
        not depth.reliable
        and depth.depth_call == "MULTICOPY_DEPTH"
        and depth.dosage_estimate is not None
    )

    if graph_exact:
        assert graph is not None and graph.context_count is not None
        count = graph.context_count
        # A numeric comparison needs an interval that is not biased downward;
        # with low mapping uniqueness the depth interval is only a floor.
        discordant = depth.reliable and not (
            depth.ratio_ci_low <= count <= depth.ratio_ci_high
        )
        if discordant:
            flags.append("DEPTH_GRAPH_NUMERIC_DISCORDANCE")
        # The two lines of evidence agree that the locus is multi-copy but not
        # on the number: naming that a consensus overstates their agreement.
        # The graph count is still returned -- it is the better-resolved of the
        # two -- under a method name that records the disagreement.
        method = (
            "GRAPH_COUNT_OVER_DISCORDANT_DEPTH" if discordant
            else "GRAPH_DEPTH_CONSENSUS"
        )
        if count >= 2 and depth_usable and depth.depth_call == "MULTICOPY_DEPTH":
            return CopyNumberConsensusResult(
                count, "INTEGER_CONTEXT_COUNT", method, "SUPPORTED",
                count, count, flags,
            )
        if count == 1 and depth_usable and depth.depth_call == "SINGLE_COPY_COMPATIBLE":
            return CopyNumberConsensusResult(
                1, "INTEGER_CONTEXT_COUNT", method, "SUPPORTED",
                1, 1, flags,
            )
        if count == 1 and depth.depth_call == "MULTICOPY_DEPTH" and dosage_usable:
            flags.append("GRAPH_SINGLE_CONTEXT_TANDEM_DEPTH")
            return CopyNumberConsensusResult(
                depth.dosage_estimate, "MEAN_DEPTH_DOSAGE",
                "DEPTH_TANDEM_COMPATIBLE", "SUPPORTED", consensus_flags=flags,
            )
        if count >= 2 and depth_usable and depth.depth_call == "SINGLE_COPY_COMPATIBLE":
            flags.append("GRAPH_MULTICOPY_DEPTH_SINGLE_COPY_COMPATIBLE")
            return CopyNumberConsensusResult(
                None, "NOT_ESTIMATED", "NONE", "EVIDENCE_CONFLICT",
                consensus_flags=flags,
            )

    if graph_lower_bound and depth_usable and depth.depth_call == "MULTICOPY_DEPTH":
        assert graph is not None
        return CopyNumberConsensusResult(
            None, "LOWER_BOUND", "GRAPH_LOWER_BOUND", "LOWER_BOUND",
            copy_number_lower_bound=graph.context_lower_bound,
            consensus_flags=flags,
        )

    if depth_floor:
        # Every quantity available is a floor, so the larger one is reported:
        # the depth floor for a single-context tandem array, an integer graph
        # floor when more contexts are resolved than the depth floor implies.
        flags.append("DEPTH_LOWER_BOUND")
        floors = [("DEPTH_ONLY", depth.dosage_estimate)]
        if graph_exact:
            floors.append(("GRAPH_LOWER_BOUND", graph.context_count))
        elif graph_lower_bound:
            floors.append(("GRAPH_LOWER_BOUND", graph.context_lower_bound))
        elif graph is not None:
            flags.append("GRAPH_EVIDENCE_UNRESOLVED")
        method, floor = max(floors, key=lambda item: item[1])
        return CopyNumberConsensusResult(
            None, "LOWER_BOUND", method, "LOWER_BOUND",
            copy_number_lower_bound=floor,
            consensus_flags=flags,
        )

    if dosage_usable:
        if graph is not None:
            flags.append("GRAPH_EVIDENCE_UNRESOLVED")
        return CopyNumberConsensusResult(
            depth.dosage_estimate, "MEAN_DEPTH_DOSAGE", "DEPTH_ONLY", "SUPPORTED",
            consensus_flags=flags,
        )

    if graph is not None and graph.status != "MATCHED_CONTEXT_COUNT":
        flags.append("GRAPH_EVIDENCE_UNRESOLVED")
    if not depth.reliable or depth.depth_call == "DEPTH_INDETERMINATE":
        flags.append("DEPTH_EVIDENCE_UNRELIABLE")
    return CopyNumberConsensusResult(
        None, "NOT_ESTIMATED", "NONE", "INDETERMINATE",
        consensus_flags=flags,
    )
