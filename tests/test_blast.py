from locus_recon.blast import (
    collapse_overlapping_hits,
    describe_full_query_span,
    estimate_full_query_span,
    find_span_continuation,
    parse_discovery_blast_hit,
    rank_local_blast_hit,
    summarize_hit_ambiguity,
)


def _hit(query, subject, start, end, bitscore, identity=99.0, qcov=100.0):
    return rank_local_blast_hit([
        query,
        subject,
        "1",
        str(abs(end - start) + 1),
        str(start),
        str(end),
        str(abs(end - start) + 1),
        str(identity),
        str(qcov),
        str(bitscore),
        "1e-80",
        str(abs(end - start) + 1),
        "1000",
    ])


def test_discovery_hit_keeps_short_fragment_metadata():
    hit = parse_discovery_blast_hit([
        "allele_1", "contig_7", "1", "95", "401", "495", "95",
        "98.2", "19.0", "160", "1e-30",
    ])
    assert hit["length"] == 95
    assert hit["qcov"] == 19.0
    assert hit["sstart"] == 401


def test_redundant_alleles_do_not_create_false_second_candidate():
    hits = [
        _hit("allele_1", "scaffold_1", 100, 599, 900),
        _hit("allele_2", "scaffold_1", 102, 599, 890),
        _hit("allele_1", "scaffold_9", 50, 549, 700),
    ]
    candidates = collapse_overlapping_hits(hits)
    assert len(candidates) == 2
    assert candidates[0]["sseqid"] == "scaffold_1"
    assert candidates[0]["supporting_hits"] == 2

    ambiguity = summarize_hit_ambiguity(candidates[0], candidates[1], 10.0)
    assert ambiguity["bitscore_margin"] == 200.0
    assert not ambiguity["ambiguous"]


def test_close_distinct_candidates_are_flagged():
    best = _hit("allele_1", "scaffold_1", 1, 500, 900)
    second = _hit("allele_1", "scaffold_2", 1, 500, 895)
    assert summarize_hit_ambiguity(best, second, 10.0)["ambiguous"]


def test_full_query_span_projection_handles_both_strands():
    plus = rank_local_blast_hit([
        "a", "s", "11", "490", "111", "590", "480", "99", "96",
        "800", "1e-50", "500", "1000",
    ])
    reverse = rank_local_blast_hit([
        "a", "s", "11", "490", "590", "111", "480", "99", "96",
        "800", "1e-50", "500", "1000",
    ])
    assert estimate_full_query_span(plus, 1000) == (100, 600)
    assert estimate_full_query_span(reverse, 1000) == (100, 600)


def _span_hit(subject, qstart, qend, qlen, sstart, send, bitscore=1000.0):
    """Build a local hit with explicit query/subject geometry."""
    return rank_local_blast_hit([
        "bait", subject, str(qstart), str(qend), str(sstart), str(send),
        str(abs(send - sstart) + 1), "99.2", "100", str(bitscore), "0.0",
        str(qlen), "100000",
    ])


def test_span_within_contig_reports_no_clip():
    # Locus sits comfortably inside the contig: the projection fits.
    hit = _span_hit("NODE_1", qstart=1, qend=500, qlen=500, sstart=1001, send=1500)
    span = describe_full_query_span(hit, 3000)
    assert span["clipped_bp"] == 0
    assert span["unclaimed_query"] == []
    assert (span["start"], span["end"]) == (1000, 1500)


def test_span_clipped_at_contig_end_reports_lost_query_bases():
    # Hpfe0006 23S, from that run's blast_step5.tsv: the alignment ends on the
    # contig's final base and 317 bp of bait has nowhere to project.
    hit = _span_hit(
        "NODE_1_length_3387_cov_892.092145",
        qstart=1, qend=2570, qlen=2887, sstart=818, send=3387,
    )
    span = describe_full_query_span(hit, 3387)
    assert span["clipped_end_bp"] == 317
    assert span["clipped_start_bp"] == 0
    assert span["clipped_bp"] == 317
    assert span["unclaimed_query"] == [(2571, 2887)]
    # The interval handed downstream is still the clipped one.
    assert (span["start"], span["end"]) == (817, 3387)
    assert estimate_full_query_span(hit, 3387) == (817, 3387)


def test_span_clipped_at_contig_start():
    hit = _span_hit("NODE_2", qstart=101, qend=600, qlen=600, sstart=1, send=500)
    span = describe_full_query_span(hit, 4000)
    assert span["clipped_start_bp"] == 100
    assert span["clipped_end_bp"] == 0
    assert span["unclaimed_query"] == [(1, 100)]
    assert span["start"] == 0


def test_reverse_strand_clip_attributes_the_correct_query_end():
    # Reverse-strand hit running off the contig start: in query orientation the
    # bases lost are the SUFFIX, not the prefix.
    hit = _span_hit("NODE_3", qstart=1, qend=400, qlen=500, sstart=400, send=1)
    span = describe_full_query_span(hit, 4000)
    assert span["clipped_start_bp"] == 100
    assert span["unclaimed_query"] == [(401, 500)]


def test_continuation_found_when_locus_is_split_across_contigs():
    # Hpfe0006 again: the lost 317 bp align to NODE_2, which overlaps the
    # reported contig by 77 bp of bait (q2494-2570).
    best = _span_hit(
        "NODE_1_length_3387_cov_892.092145",
        qstart=1, qend=2570, qlen=2887, sstart=818, send=3387, bitscore=4641.0,
    )
    other = _span_hit(
        "NODE_2_length_784_cov_947.591231",
        qstart=2494, qend=2887, qlen=2887, sstart=1, send=394, bitscore=712.0,
    )
    span = describe_full_query_span(best, 3387)
    cont = find_span_continuation(span, [best, other], best)
    assert cont is not None
    assert cont["sseqid"].startswith("NODE_2")
    assert cont["recovered_bp"] == 317
    assert cont["overlap_bp"] == 77


def test_no_continuation_when_lost_sequence_is_absent_from_the_assembly():
    best = _span_hit("NODE_1", qstart=1, qend=2570, qlen=2887, sstart=818, send=3387)
    # A second candidate that covers only sequence already reported.
    other = _span_hit("NODE_9", qstart=100, qend=900, qlen=2887, sstart=1, send=801)
    span = describe_full_query_span(best, 3387)
    assert find_span_continuation(span, [best, other], best) is None


def test_continuation_is_none_when_nothing_was_clipped():
    hit = _span_hit("NODE_1", qstart=1, qend=500, qlen=500, sstart=1001, send=1500)
    span = describe_full_query_span(hit, 3000)
    assert find_span_continuation(span, [hit], hit) is None
