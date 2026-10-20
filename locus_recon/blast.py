"""BLAST hit parsing, candidate collapsing, and ranking helpers."""

from typing import List, Optional

def parse_discovery_blast_hit(parts: list) -> dict:
    """Parse one target-discovery record in ``BLAST_FMT_STEP1`` order."""
    if len(parts) < 11:
        raise ValueError(f"Expected 11 BLAST discovery fields, found {len(parts)}.")
    return {
        "qseqid": parts[0],
        "sseqid": parts[1],
        "qstart": int(parts[2]),
        "qend": int(parts[3]),
        "sstart": int(parts[4]),
        "send": int(parts[5]),
        "length": int(parts[6]),
        "pident": float(parts[7]),
        "qcov": float(parts[8]),
        "bitscore": float(parts[9]),
        "evalue": float(parts[10]),
    }


def rank_local_blast_hit(parts: list) -> dict:
    """Parse and score one local-assembly BLAST hit (format BLAST_FMT_STEP5).

    Ranking tuple: (bitscore, qcov, pident, length) — bitscore first, so
    paralog or partial-match ambiguity is less likely to win over a true hit.

    Args:
        parts: Tab-split fields from one BLAST output line (13 fields expected).

    Returns:
        Dict with all parsed fields plus a 'rank' tuple for sorting.
    """
    return {
        "qseqid":  parts[0],
        "sseqid":  parts[1],
        "qstart":  int(parts[2]),
        "qend":    int(parts[3]),
        "sstart":  int(parts[4]),
        "send":    int(parts[5]),
        "length":  int(parts[6]),
        "pident":  float(parts[7]),
        "qcov":    float(parts[8]),
        "bitscore": float(parts[9]),
        "evalue":  float(parts[10]),
        "qlen":    int(parts[11]),
        "slen":    int(parts[12]),
        "rank": (float(parts[9]), float(parts[8]), float(parts[7]), int(parts[6])),
    }


def describe_full_query_span(hit: dict, subject_length: int) -> dict:
    """Project the full locus onto a contig and report what the contig cannot hold.

    BLAST reports only the aligned subject interval.  Taking that interval
    verbatim truncates a reconstructed locus whenever one or both query ends
    do not align, so the unaligned query prefix/suffix is projected onto the
    local scaffold.  When the scaffold ends inside the locus that projection
    runs past the contig and must be clipped, and the clipped bases are query
    sequence the reconstruction will silently lack.

    Clipping is 1:1 in projected coordinates, so the number of clipped subject
    bases equals the number of query bases lost, and can never exceed the
    unaligned remainder on that side.

    Args:
        hit:            Parsed BLAST hit (needs qstart/qend/qlen/sstart/send).
        subject_length: Length of the contig the hit lies on.

    Returns:
        Dict with ``start``/``end`` (0-based, half-open subject interval, as
        returned by :func:`estimate_full_query_span`), ``clipped_start_bp`` and
        ``clipped_end_bp`` (subject bases the contig could not supply at each
        of its own ends), ``clipped_bp`` (their sum), and ``unclaimed_query``
        -- the 1-based inclusive query intervals the clipped projection failed
        to cover.
    """
    qstart, qend, qlen = hit["qstart"], hit["qend"], hit["qlen"]
    sstart, send = hit["sstart"], hit["send"]
    reverse = sstart > send

    if not reverse:
        start = sstart - qstart
        end = send + (qlen - qend)
        # Forward strand: the query prefix projects off the contig start and
        # the query suffix off the contig end.
        lost_prefix = max(0, -start)
        lost_suffix = max(0, end - subject_length)
    else:
        start = send - 1 - (qlen - qend)
        end = sstart + (qstart - 1)
        # Reverse strand: orientation flips which query end each contig
        # boundary truncates.
        lost_suffix = max(0, -start)
        lost_prefix = max(0, end - subject_length)

    unclaimed = []
    if lost_prefix:
        unclaimed.append((1, lost_prefix))
    if lost_suffix:
        unclaimed.append((qlen - lost_suffix + 1, qlen))

    return {
        "start": max(0, start),
        "end": min(subject_length, end),
        "clipped_start_bp": max(0, -start),
        "clipped_end_bp": max(0, end - subject_length),
        "clipped_bp": max(0, -start) + max(0, end - subject_length),
        "unclaimed_query": unclaimed,
    }


def estimate_full_query_span(hit: dict, subject_length: int) -> tuple:
    """Estimate full-locus subject bounds from a possibly terminally clipped HSP.

    Thin wrapper over :func:`describe_full_query_span` retained for callers
    that need only the interval.  Use the descriptive form when the fact that
    the projection was clipped matters.

    Returns:
        ``(start, end)`` as a 0-based, half-open subject interval.
    """
    span = describe_full_query_span(hit, subject_length)
    return span["start"], span["end"]


def find_span_continuation(
    span: dict, candidate_regions: List[dict], best_hit: dict
) -> Optional[dict]:
    """Find another candidate contig carrying query sequence the clip discarded.

    A clipped projection means the locus runs off the end of the chosen contig.
    If a different candidate region aligns to the query interval that was lost,
    the locus is split across the local assembly rather than absent from it --
    the discarded sequence is present in the data and was simply not used.
    Those two cases mean different things, so the distinction is reported and
    not acted upon: joining contigs is a scaffolding decision, and this tool
    reports assembly geometry rather than resolving it.

    Args:
        span:              Dict from :func:`describe_full_query_span`.
        candidate_regions: Collapsed candidate hits, best first.
        best_hit:          The candidate the span was projected from.

    Returns:
        ``None`` when no candidate covers the lost interval, else a dict with
        ``sseqid``, ``recovered_bp`` (lost query bases that contig aligns to)
        and ``overlap_bp`` (query overlap shared with the chosen contig, i.e.
        sequence available to anchor a join).
    """
    if not span["unclaimed_query"]:
        return None

    best_id = best_hit["sseqid"]
    best_qstart, best_qend = best_hit["qstart"], best_hit["qend"]
    found = []
    for cand in candidate_regions:
        if cand["sseqid"] == best_id:
            continue
        recovered = 0
        for lo, hi in span["unclaimed_query"]:
            recovered += max(0, min(hi, cand["qend"]) - max(lo, cand["qstart"]) + 1)
        if recovered <= 0:
            continue
        overlap = max(
            0, min(best_qend, cand["qend"]) - max(best_qstart, cand["qstart"]) + 1
        )
        found.append({
            "sseqid": cand["sseqid"],
            "recovered_bp": recovered,
            "overlap_bp": overlap,
        })

    if not found:
        return None
    return max(found, key=lambda d: (d["recovered_bp"], d["overlap_bp"]))


def collapse_overlapping_hits(hits: List[dict], max_gap: int = 20) -> List[dict]:
    """Collapse redundant bait hits that describe the same subject region.

    A multi-allele bait database normally produces many high-scoring hits to
    the *same* reconstructed locus.  Treating those records as independent
    candidates makes the second-best hit test report false paralog warnings.
    This function groups overlapping (or nearly adjacent) subject intervals
    and keeps the best-ranked representative from each distinct region.

    Args:
        hits: Parsed local BLAST hits from :func:`rank_local_blast_hit`.
        max_gap: Maximum gap between intervals that still belong to one region.

    Returns:
        Best representative for each distinct subject region, sorted by rank.
        Representatives also contain ``region_start``, ``region_end``, and
        ``supporting_hits`` fields for provenance.
    """
    if not hits:
        return []

    by_subject: dict = {}
    for hit in hits:
        start = min(hit["sstart"], hit["send"])
        end = max(hit["sstart"], hit["send"])
        by_subject.setdefault(hit["sseqid"], []).append((start, end, hit))

    candidates: List[dict] = []
    for subject, subject_hits in by_subject.items():
        subject_hits.sort(key=lambda item: (item[0], item[1]))
        group: List[dict] = []
        group_start = group_end = 0

        def emit_group() -> None:
            if not group:
                return
            representative = dict(max(group, key=lambda item: item["rank"]))
            representative["region_start"] = group_start
            representative["region_end"] = group_end
            representative["supporting_hits"] = len(group)
            candidates.append(representative)

        for start, end, hit in subject_hits:
            if not group:
                group = [hit]
                group_start, group_end = start, end
            elif start <= group_end + max_gap:
                group.append(hit)
                group_end = max(group_end, end)
            else:
                emit_group()
                group = [hit]
                group_start, group_end = start, end
        emit_group()

    return sorted(candidates, key=lambda item: item["rank"], reverse=True)


def summarize_hit_ambiguity(
    best_hit: dict,
    second_hit: Optional[dict],
    min_bitscore_margin: float = 10.0,
) -> dict:
    """Compare the best and second-best local BLAST hits to flag ambiguity.

    A margin below min_bitscore_margin is flagged as ambiguous, meaning the
    reconstructed allele could plausibly originate from a paralog.

    Args:
        best_hit:            Top-ranked hit dict from rank_local_blast_hit().
        second_hit:          Second-ranked hit, or None if only one hit exists.
        min_bitscore_margin: Minimum bitscore gap to call the best hit unambiguous.

    Returns:
        Dict with keys: second_best_bitscore, second_best_identity,
        bitscore_margin, ambiguous.
    """
    if second_hit is None:
        return {
            "second_best_bitscore":  0.0,
            "second_best_identity":  0.0,
            "bitscore_margin":       999_999.0,
            "ambiguous":             False,
        }
    margin = best_hit["bitscore"] - second_hit["bitscore"]
    return {
        "second_best_bitscore": second_hit["bitscore"],
        "second_best_identity": second_hit["pident"],
        "bitscore_margin":      margin,
        "ambiguous":            margin < min_bitscore_margin,
    }
