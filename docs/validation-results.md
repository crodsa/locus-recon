# Validation results of release 1.0.1

This page gives the result of each validation analysis and the folder that holds its inputs, outputs and script.
[validation/README.md](../validation/README.md) explains how to rerun them, and [validation.md](validation.md)
describes the method. No single data set validates every output, so each analysis states what it tests and what
it does not.

## Summary

| Analysis | Folder under `validation/` | What it tests | Result | What it does not show |
|---|---|---|---|---|
| Simulated reconstruction benchmark, 13 samples | `mock_dataset/`, `reference_results/` | Exact recovery across contig breaks; refusal of weak, mixed, paralogous and absent targets | Six of six coded criteria met | Sensitivity or a limit of detection outside this truth model |
| Completeness benchmark, 14 constructed assemblies | `completeness-benchmark/` | Whether the measured loss equals the missing sequence | Equal in all 14 cases, 0 to 822 bp | Discovery on real assemblies; no reads are mapped |
| Closed genomes, five *Helicobacter pylori* isolates | `tier-calibration/` | Sequence and tier against the closed genome, with one bait from strain 26695 | 19 of 25 loci accepted, all 19 identical to the closed genome | Other species; the closed genomes were built in part from short reads of the same isolates |
| Evaluation of the tier, 62 cases | `tier-calibration/` | What an accepted result is worth | 25 accepted, all identical to the truth: 19 of 35 cases with real reads, 3 of 13 with simulated reads, 3 of 14 constructed without reads | The cases are not independent and come from five isolates of one species; the pooled interval that the runner prints is not a property of the tool |
| Checks of single components | `frame-and-graph-evidence/` | Frame search, bait validation, length profile, runs of N, path ranking | Recorded outcome of each check | The checks use constructed inputs |
| Synthetic copy series, 27 data sets | `constructed-copy-series/` | Depth arithmetic and bootstrap for one to three copies | Median absolute error 0.04 copies; all classes correct | Coverage bias of real libraries; reads are error-free and uniform |
| *aphA1* in *Acinetobacter baumannii*, seven runs | `aphA1-copy-number/` | Depth dosage of a real amplification | 10.5 and 77.7 copies against 10 ± 2 and 75 ± 14 by published quantitative PCR | Six runs come from one patient |
| 23S rRNA gene and *gyrB*, five *H. pylori* isolates | `low-copy-23S/` | A small copy number | Depth ratios of 2.55 to 3.22 for two copies; graph count of two in four isolates, lower bound in the fifth | These genomes also guided the design of the graph method |
| *Clostridioides difficile* LIBA-6656, *tcdB* | `liba6656-gfa-rerun/`, `liba6656-depth/`, `depth-copy-number/overlap_null/` | A locus with two divergent copies that collapse in the draft | The workflow held a 1,955-bp fragment; the query-weighted dosage of 2.29 is explained by overlapping alignments | Two full-length sequences were built by hand from the graph. They are hypotheses that await long reads |
| Comparison with ARIBA 2.14.7 and KMA 1.7.1, 39 cases | `comparator-benchmark/` | The sequences that three tools return from the same reads and reference sequences | Exact in 19, 19 and 20 of 25 *H. pylori* loci for Locus-Recon, ARIBA and KMA | Other species, tuned settings, SRST2 or stringMLST |
| External evaluation, 16 species, 128 units of isolate and locus | `external-evaluation/` | Sequence, decision and copy number on public isolates that took no part in development, with a protocol written before the runs | Of 113 units with a truth copy, 98 sequences exact (KMA 99, ARIBA 94); 95 accepted, 89 of them exact; the six accepted and not exact belong to one isolate whose public reads do not match its public genome | Difficult loci; an independent truth; the copy number of amplified genes |
| Targeted tests, 790 simulated samples | `targeted-tests/` | Mixtures with few sites, strand-bias alerts, ends and indels, graph contexts, overlapping spans, bootstrap block size | No inexact sequence accepted; a mixture that differs at one site was accepted in 50 of 50 samples | Platform errors, G+C bias and the repeats of a real genome |
| Bootstrap block size, real data | `block-size-sensitivity/` | Dependence of the interval on the block size | The interval changes up to 5.5-fold with the block size and often excludes the reference value | It is not an interval for the copy number |

## Reconstruction

**Simulated benchmark.** Twelve of 13 samples yielded a candidate and the sample without the target was skipped.
Eleven candidates were identical to the truth. Three were accepted (`HIGH`): the intact draft and the two drafts
broken inside the locus. The depth series was exact at every depth but not accepted: `SUSPECT` at 3× and 5×,
`LOW` at 10× and `MEDIUM` at 20×. Mixtures with 10, 20 and 50% of a second allele had 22, 25 and 25 credible
mixed sites out of 25 differing positions and were held. The 5% mixture had two such sites, the minimum that
defines a mixture. The paralogue control returned a sequence with five wrong bases and was held. This benchmark
runs with two settings that differ from the defaults (`--min-remap-depth 3`, `--min-bitscore-margin 80`).

**Completeness.** The three assemblies without loss were `HIGH`. All nine losses of 10 bp or more raised
`ALLELE_SPAN_CLIPPED_AT_CONTIG_END`. In the two cases with the missing part on a second contig, that contig was
named and the candidate was held. The tiers of this benchmark use fixed read-support values, because no reads are
involved.

**Closed genomes.** Five *gyrB*, five *ureB*, five *glmM* and four 23S rRNA sequences were accepted, and all 19
were identical to the closed genome. The 15 single-copy sequences were 95.7 to 97.5% identical to the bait. No
*cagA* candidate was accepted: the five were 87.9 to 89.5% identical to the bait and matched their genome at every
aligned position, but their length differed from the annotated gene by 6 to 33 bp. The 23S rRNA candidate of
Hpfe0006 was identical to the genome over 2,570 bp and lacked 317 bp; the workflow measured the loss, found the
missing bases on a second scaffold and held the candidate. The 23S rRNA candidate of Hpfe0002 was accepted with
one credible mixed site (position 1,639, alternative fraction 0.47), which may reflect a difference between the
two gene copies or a mixed population.

**What the tiers mean.** Of the 37 cases that were not accepted, 17 were identical to the whole truth sequence,
15 were incomplete but identical to the truth over their own length, two *cagA* candidates overran the gene, two
contained wrong bases and one sample had no candidate. One of the two sequences with a wrong base was `LOW`, that
is, in review (*gyrB* of Hpfe0002 at about 8×, one wrong base in 2,322). A tier below `HIGH` means that evidence
is incomplete, not that the sequence is wrong, and a sequence in review can be wrong.

## External evaluation

The protocol ([`validation/external-evaluation/PROTOCOL.md`](../validation/external-evaluation/PROTOCOL.md)) fixed
the software, the isolates, the loci, the truth and the decision rules for ARIBA and KMA on 18 November 2026,
before any run. Sixteen FDA-ARGOS isolates of 16 species were analysed for eight loci each
(*gyrB*, *rpoB*, *recA*, *dnaK*, *atpD*, *tuf* and the 23S and 16S rRNA genes) with reads subsampled to 100×.
[`validation/external-evaluation/results/README.md`](../validation/external-evaluation/results/README.md) gives
every table, and `DEVIATIONS.md` in that folder every departure from the protocol.

**Units.** The bait rule of the protocol found a bait for 116 of the 128 units. Of these, 81 have one copy in the
complete genome, 6 several identical copies, 26 several different copies and 3 no copy under the rule of the
protocol. The other 12 units were run with a fallback rule and are reported apart.

**One isolate does not match its genome.** The draft of *Neisseria meningitidis* FDAARGOS_209 (run SRR4272082)
aligns to 93.77% of the complete genome GCF_002073575.2 at 97.35% identity. The other 15 drafts align to at least
98.85% of their genome at 99.97% identity or more. This comparison was not in the protocol. The isolate is kept in
every table, and every summary is also given without it; the second set was defined after the results were seen.

**Sequences and decisions (113 units with a truth copy).** Locus-Recon returned a sequence identical to a truth
copy in 98 units, KMA in 99 and ARIBA in 94. Locus-Recon accepted 95 sequences, and 89 were exact. The six that
were accepted and not exact (3 to 86 substitutions) belong to FDAARGOS_209; KMA and ARIBA returned and accepted the
same six sequences. Without that isolate, all 89 accepted sequences were exact; KMA accepted 91 under its primary
rule (one not exact) and ARIBA 83 (one not exact). In the 81 units with one copy the three tools behaved alike.
The decision differed where copies differ (26 units): Locus-Recon accepted 9 sequences, all exact, and held 17, of
which 8 were complete and exact; KMA accepted 15 under its primary rule (2 not exact) and returned an exact
sequence in 19; ARIBA accepted 5 (1 not exact). One exact sequence of a one-copy unit was held with
`STRAND_BIASED_ALT_SUPPORT` as its only flag (1 of 81). Seven of the 95 accepted sequences carried
`MIXED_SITE_BELOW_MIXTURE_CALL`.

**Copy number.** No one-copy unit was called `MULTICOPY_DEPTH` (0 of 81); the depth ratio of these units had a
median of 1.087 (0.787 to 1.589) and rose with the G+C content of the genome. Of the 32 units with two or more
copies, 21 were called `MULTICOPY_DEPTH`. The bootstrap interval of the ratio contained the true copy number in
38 of 113 units. The graph context count equalled the true copy number in 69 units, was a correct lower bound in
32, was below the truth without being marked as a bound in 6, and was absent in 6. The consensus was correct in 83
of 89 units with the status `SUPPORTED`, 6 of 14 with `REVIEW` and 10 of 10 with `LOWER_BOUND`. Five of the six
wrong `SUPPORTED` results are rRNA genes with six or seven copies for which the graph gave two contexts while the
depth ratio was 4.80 to 7.45; the sixth is in FDAARGOS_209. The context count is the number of distinct flanks,
so it is a lower bound when copies share their flanks.

## Targeted tests with simulated reads

[`validation/targeted-tests/README.md`](../validation/targeted-tests/README.md) describes the 790 samples of
section 3 of the protocol and `DEVIATIONS.md` in that folder the departures from it.

**Mixtures with few sites.** A mixture of two alleles that differ at one site was accepted in 50 of 50 samples,
with the major allele as sequence; 41 carried `MIXED_SITE_BELOW_MIXTURE_CALL` and 9 no flag (7 with a minor
fraction of 5%). With two sites, 2, 9 and 10 of 10 samples were held at 5, 10 and 20%; with five sites, 6 and 10 of
10 at 5 and 10%; with ten sites, 9 and 10 of 10. No accepted sequence was inexact.

**Strand-bias alerts.** `STRAND_BIASED_ALT_SUPPORT` held 4 of 160 pure samples whose sequence was exact (3 of 10
at 6,999 bp and 20×, 1 of 10 at 501 bp and 1,000×).

**Ends, indels and pseudogenes.** 133 of 135 sequences were exact; the two that were not (10×) were held. A 12-bp
insertion, a 1-bp deletion and a premature stop codon were recovered exactly and never accepted (0 of 15 each):
the workflow does not accept a correct allele whose length or reading frame departs from the bait.

**Graph contexts and depth.** Identical copies that the assembler collapsed were counted correctly in 40 of 40
samples with two to four copies. Two complete identical copies in the draft gave a depth ratio of 0.000 and no
consensus (20 of 20), a tandem duplicate was not quantified (0 of 10), and copies that differ by 1% gave no graph
count in 69 of 70 samples. No consensus with the status `SUPPORTED` was wrong.

**Overlapping spans.** The dosage separated one copy from two only when the discovery intervals shared 25% of
the query or less. At 65% or more it returned 0.000 as a lower bound for both. Without the graph, no consensus was
`SUPPORTED` when spans overlapped.

**Block size.** In the simulated samples and in the *aphA1* and 23S rRNA data
(`validation/block-size-sensitivity/`), the bootstrap interval changed with the block size and often excluded the
reference value: for *aphA1* it contained it in 3, 2, 2 and 0 of 6 cases with blocks of 100, 250, 500 and
1,000 bp. The interval describes sampling along the locus. It is not an interval for the copy number.

## Comparison with ARIBA and KMA

ARIBA 2.14.7 and KMA 1.7.1 were run with default settings on the reads and reference sequences of three analyses
above. [`validation/comparator-benchmark/README.md`](../validation/comparator-benchmark/README.md) gives the
settings, the scoring and every case.

**_H. pylori_, 25 loci.** The three tools returned the exact sequence of every *gyrB*, *ureB* and *glmM* gene.
For the 23S rRNA gene, KMA was exact in all five isolates, and ARIBA and Locus-Recon in four. ARIBA called T for C
at position 1,639 of Hpfe0002, with its flag for a collapsed repeat and both nucleotides listed (121 reads with
T, 123 with C). Locus-Recon reports this site as mixed and KMA writes it in lower case. ARIBA and KMA returned
the complete gene of Hpfe0006.
Locus-Recon returned 2,570 of the 2,887 bp of Hpfe0006 and held them. No tool returned an exact *cagA* gene.
ARIBA returned fragments of 146 to 1,323 bp in four isolates and nothing in the fifth. KMA returned a consensus
of 3,341 to 3,412 bp at an edit distance of 145 to 239 from the gene, with 13 to 65 bases in lower case. The
Locus-Recon candidates matched the
genome at every aligned base, had ends wrong by 6 to 33 bp, and were not accepted.

**Simulated samples.** All tools were exact for the seven samples with one allele. In the four mixtures KMA
returned both alleles, ARIBA returned one, and Locus-Recon returned the major allele and held it as a mixture.
In the paralogue control ARIBA returned the true allele and the paralogue as two sequences. KMA returned one
sequence with three wrong bases and one extra base, all four in lower case. Locus-Recon returned a sequence with five wrong
bases and held it as ambiguous.

**LIBA-6656.** ARIBA returned 24 fragments of at most 1,421 bp. KMA returned two full-length sequences. The one
on template B7.6 is identical to the chromosome-associated analyst-built candidate. The one on template B7.9 is
the extrachromosomal candidate without four bases (positions 5,616 to 5,619), with the 16 bases around the gap in
lower case; 67 reads span that site and none carries the deletion. The standard Locus-Recon workflow returned one held fragment of 1,955 bp.

**Reading.** When a close reference exists, KMA is the quick route to the sequence, and Locus-Recon does not
return a better one. KMA raises no flag. It writes in lower case a base that is not significantly
overrepresented in the reads, and every KMA sequence that differed from the truth carried that mark at or next
to its errors; the mark does not say what is wrong and cannot show missing bases. ARIBA reported its doubt in
flags and in the read depths of each difference, except for the 5 and 10% mixtures, and in the paralogue control
it alone returned the true allele. Neither tool states a decision. In these 39 cases Locus-Recon accepted no sequence
that differed from the truth; it states why a sequence was not accepted and estimates copy number.

## Copy number

**Synthetic series.** Dosage was 0.94 to 1.10 for one copy, 1.94 to 2.10 for two and 2.92 to 3.14 for three.
None of the nine single-copy data sets was called multicopy, and all 18 others were.

**_aphA1_.** The four isolates with a tobramycin MIC of 0.5 mg/L had dosages of 0.95 to 1.10. MRSN57 had 10.5
copies (10.2 to 10.9) and MRSN58 77.7 (73.0 to 79.0). A derivative selected *in vitro* had 91.4 copies; no
published value exists for comparison.

**23S rRNA gene and _gyrB_.** Depth called every 23S rRNA locus multicopy and every *gyrB* locus single-copy.
The 23S rRNA ratios exceeded two in all five isolates, and the cause was not established. The graph gave two
contexts at both ends in four isolates. In Hpfe0004 the walk met a cycle and the result is a lower bound of two.

**LIBA-6656.** The depth ratio of *tcdB* was 1.21 (1.11 to 1.35), compatible with a single copy, because each
divergent copy lies on its own contigs. The query-weighted dosage was 2.29, but the nine intervals overlap on the
gene and would give 1.84 if every ratio were one. In 87 control windows of the same assembly, dosage rose with
overlap, and the fit predicts 2.21 (1.89 to 2.53) at the overlap of *tcdB*. The dosage is therefore labelled
`REVIEW_OVERLAPPING_SPANS` and is not a count.

## Provenance

Every analysis of this page was run again on 18 and 19 November 2026 with release 1.0.1 and the real external
programs (`validation/regenerate/`), and the records under `validation/` are the output of that run. They carry
the SHA-256 digest of the package sources (`locus_recon/*.py`),
`2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0`, which `validation/common.py` computes.
Against the records of release 1.0.0 (tag `v1.0.0`), 14,383 values were compared: no tier, disposition, sequence,
dosage, ratio or interval changed, and the only new label is `MIXED_SITE_BELOW_MIXTURE_CALL` on the 23S rRNA
sequence of Hpfe0002. ARIBA and KMA were not run again; their outputs of 11 and 13 November 2026 were scored
again, and the rows of Locus-Recon were scored from the regenerated candidates. The older records at the top of
`validation/depth-copy-number/` are those of release 1.0.0. [`validation/INTEGRATION_RUN.md`](../validation/INTEGRATION_RUN.md)
records one run of the complete workflow with the real programs; the unit tests replace those programs with
stubs.
