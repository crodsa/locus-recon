# Methods and validation guide

This document provides a concise methods description and a practical validation design for using Locus-Recon in research. Thresholds must be calibrated for the organism, sequencing protocol, locus family, and allele scheme under study.

## Computational method

For each sample, known alleles for one locus are aligned to the draft assembly with BLASTN. Target discovery retains HSPs meeting minimum identity and aligned-length requirements; it does not require high per-HSP query coverage. HSP coordinates are padded, clipped to contig boundaries, and merged per contig.

Paired reads are aligned to the draft with BWA-MEM2 (or BWA). The aligner stream is coordinate-sorted directly by samtools. Alignments overlapping target intervals seed a read-name set, after which all alignments sharing those names are recovered to rescue mates. Primary paired reads and singletons are converted to FASTQ and locally assembled with SPAdes.

Known alleles are aligned to the local scaffolds. Overlapping hits from multiple bait alleles are collapsed into distinct subject regions before ranking, preventing redundant alleles at one locus from acting as false alternative loci. The best region is ranked by bitscore, query coverage, identity, and aligned length. Its HSP is extended by the unaligned query prefix and suffix, within scaffold bounds, to estimate full locus boundaries. Reverse-strand candidates are reverse-complemented.

The candidate is aligned to a batch-shared BLAST database of known alleles.
Extracted paired reads and singletons are independently remapped to the
candidate. A quality-filtered samtools pileup, with BAQ recalibration and `-A`
so that read pairs not flagged as properly paired against a short candidate
still count, includes zero-depth positions and records base quality, mapping
quality, forward/reverse depth, reference and alternative observations,
alternative fraction, and a two-sided Fisher strand bias test. Minimum identity, query coverage, mean depth, and breadth gates are
applied before the result receives `SUCCESS` status.

A candidate mixed site must exceed the configured alternative fraction and
depth and have alternative observations on both strands. The locus-level
mixture call requires multiple such sites. One-strand alternatives are reported
separately. Quality confidence additionally considers length relative to the
bait distribution, GC deviation, ambiguous bases, reading-frame consistency,
internal stops (the minimum over all six frames), completeness of the projected
locus span, depth patchiness, mean qualities,
strand balance, mixture evidence, reconstructed-reference discordance, and
score separation from the second-best distinct local region. A truncation at
a contig end is scored once: the length checks are applied with the clipped
bases restored. The QC tier and every downgrade reason are emitted in
human-readable and tabular reports.

### Research-mode multi-copy analysis

When the default candidate is truncated or mixed and a duplicated target is
biologically plausible, the retained local SPAdes GFA may contain complete
alternative paths even when `scaffolds.fasta` does not. The
`locus-recon-graph-paths` command exports bounded terminal paths and their
segment/depth metadata without treating them as allele calls.

Candidate copies should be assigned first by independent genomic context
(unique chromosome/plasmid flanks or other copy-specific anchors), then by
complementary graph branches. Reads must be aligned competitively to all
candidates. Reads anchored uniquely to one copy can rescue their mates and
support copy-specific polishing. A final competitive remap should report every
base, mapping/base quality, strand balance, discordance, and residual mixture.
Intervals that contain no copy-specific variant within the library's phasing
distance must be reported as unphaseable unless long-read or other independent
evidence links them.

The repository includes a complete real-data example under
[`validation/liba6656-gfa-rerun/`](../validation/liba6656-gfa-rerun/README.md).
It shows the retained GFA, all exported terminal paths, portable provenance,
candidate-to-path alignments and the distinction between a graph proposal and
a competitively supported, polished candidate.

## Recommended validation panel

Build a locked benchmark before using reconstructed alleles in a study.

The repository includes an executable synthetic panel in
[`validation/`](../validation/README.md). Use it for installation checks and
regression testing, then add real samples representative of the intended study.

| Panel | Construction | Expected outcome |
|---|---|---|
| Intact positive controls | High-quality assemblies with known allele calls and matched reads. | Exact sequence recovery and correct `nearest_allele`. |
| Synthetic fragmentation | Split each known locus at several positions in an otherwise trusted assembly; retain original reads. | Recovery across central and near-terminal breaks. |
| Coverage dilution | Downsample matched reads across a relevant depth series. | Defined sensitivity boundary and appropriate QC downgrades. |
| Divergence series | Use validated alleles spanning the observed within-species diversity. | Calibrated discovery/validation identity thresholds. |
| Paralog controls | Include genomes with a related locus or duplicated target region. | Alternative distinct regions detected or downgraded as ambiguous. |
| Negative controls | Species or samples known to lack the target, plus unrelated reads. | No confident reconstruction. |
| Mixtures/contamination | Mix reads from two known alleles at controlled proportions, while keeping total depth constant. | Explicit mixture detection above a calibrated boundary; no confident chimeric interpretation. |
| Real fragmented isolates | Resolve with an independent method such as long-read/hybrid assembly or targeted sequencing. | Sequence concordance at every reconstructed base. |

## Metrics to report

At minimum, report:

- exact-sequence accuracy, not only locus recovery;
- sensitivity and false-positive rate by fragmentation and depth stratum;
- nearest-allele identity and coverage;
- remap mean depth, breadth, and fraction below 5×;
- mean base and mapping quality and total strand balance;
- mixed-site count, mixture call, alternative-fraction distribution, and
  one-strand alternative count;
- frequency of each QC tier and flag;
- ambiguous-region frequency;
- error modes after manual review;
- Locus-Recon version, full parameters, tool versions, bait checksum, and bait database release.

Do not tune thresholds and estimate final performance on the same isolates. Use a development panel for threshold selection and a separate held-out panel for reporting.

## Manual review checklist

Before accepting a novel or submission-bound sequence:

1. Confirm sample identity and that the reads match the draft assembly.
2. Review contamination and mixed-isolate QC from the upstream workflow.
3. Inspect the local scaffold and both BLAST tables.
4. Inspect `allele_remap.support.tsv`, especially termini, indels,
   low-complexity segments, alternative fractions, and forward/reverse support.
5. Confirm there is no similarly supported distinct region or known paralog.
6. Compare locus length and translation/frame behavior with the scheme definition.
7. Verify sequence orientation and exact scheme boundaries.
8. Confirm with an independent assembly or sequencing method when the result changes a scientific conclusion or represents a new allele.

## Depth-based copy number

Allelic evidence is silent when the copies of a locus are identical, so a
separate module compares filtered depth over the locus with filtered depth over
a single-copy backbone. Its criteria, thresholds and validation record are kept
together in
[`validation/depth-copy-number/PRESPECIFIED_CRITERIA.md`](../validation/depth-copy-number/PRESPECIFIED_CRITERIA.md).

In summary:

- numerator and denominator use identical Q20/MQ20 filters; contigs under 5 kb,
  100 bp at every contig end, and any region passed to `--exclude` are kept out
  of the backbone;
- the confidence interval is a block bootstrap over 500 bp blocks, because
  per-base depth is autocorrelated; the point denominator uses every eligible
  base, while CI resampling uses at most 512 deterministically and evenly
  spaced backbone blocks so chromosome-scale runs remain executable;
- `MULTICOPY_DEPTH` requires the lower confidence bound to exceed 1.5;
- no windowed per-window criterion for local excess is reported: thresholds for
  one would be the P99 of a null distribution over single-copy backbone
  segments, which calibrates the false-positive rate on unamplified sequence
  but carries no measured sensitivity against loci with known partial
  amplifications;
- `ambiguity_index` below 0.70 marks the ratio as a lower bound rather than
  blocking the call, because the MAPQ filter can only remove depth;
- `dosage_estimate` weights each assembly region's ratio by the gene bases it
  resolves, and is left unset only when a multi-region locus lacks gene
  coordinates. Incomplete query coverage qualifies the estimate rather than
  withholding it: `dosage_status` becomes `ESTIMATED_PARTIAL_COVERAGE` and
  `query_coverage_fraction` records the covered proportion. Overlapping spans
  set `REVIEW_OVERLAPPING_SPANS` in precedence over that, because the additive
  value is then not a copy count — see
  `validation/depth-copy-number/overlap_null/` — and `geometry_expected_dosage`
  reports what the same spans return at unit ratios. `copies_estimate` carries
  the same value as `dosage_estimate`.

External validation used seven public Illumina runs of *Acinetobacter baumannii*
strains with independently published *aphA1* copy numbers: four single-copy
strains were called `SINGLE_COPY_COMPATIBLE` and three amplified strains
`MULTICOPY_DEPTH`, with both point estimates that have a published qPCR interval
falling inside it. The module separates one copy from two with margin; it does
not claim to separate two from 2.5. How far the separation extends was measured
afterwards on constructed geometries: an estimate of 2.258 excludes one copy at
p = 0.004 and is compatible with two and with two and a half alike, which is why
no rounding rule is applied to a fractional estimate.

## Graph/depth copy-number consensus

Depth and graph topology are retained as separate estimands. `ratio` and
`dosage_estimate` quantify continuous pile-up over the assembly representation;
`graph_context_count` is an integer count of distinct flanking contexts retained
at both bait-defined locus ends. `locus-recon-copy-number` combines them through
an explicit decision table without modifying the depth algorithm or its
thresholds.

The graph step aligns a selected bait to embedded SPAdes GFA segments, requires
at least 90% union coverage of the bait, and traverses both query-oriented ends
under explicit 5,000 bp, eight-node and 64-path limits. An exact integer is
eligible only when both ends yield the same positive context count. One-sided
or asymmetric results become lower bounds or indeterminate states; cycles and
safety limits remain visible as flags. The context count does not establish
allele-to-copy phasing, replicon closure, physical chromosome location or
expression.

The five-isolate *H. pylori* 23S panel was rerun from the frozen Illumina inputs
and retained SPAdes graphs. All five 23S loci returned two matched contexts and
`copy_number_call=2`; all five single-copy `gyrB` controls returned one. The raw
23S depth ratios remained 2.546–3.221, and the workflow verified byte-equivalent
values for all 23 columns of the depth-only table archived before the graph
analysis. Because these examples informed the graph/depth design, they are a
biological demonstration rather than an independent estimate of exact-copy
sensitivity.

The *aphA1* experiment remains the orthogonal depth-dosage application. A
collapsed tandem array may have one graph context while its read depth supports
many mean copies; this combination returns `DEPTH_TANDEM_COMPATIBLE` and keeps
the continuous dosage. Regression tests freeze the retained seven-run result,
evaluation and summary hashes so the graph logic cannot silently replace the
*aphA1* interpretation with a one-copy call.

## Known limitations

- Short reads may not resolve repeats longer than the library insert size.
- A local assembler can produce a plausible consensus from a mixture or closely
  related paralogs. Mixture evidence is a screen, not haplotype phasing or
  automatic component-allele reconstruction. The graph-path exporter exposes
  alternatives but does not decide which path is biologically correct.
- Fisher strand-bias tests have limited power at low depth and become highly
  sensitive at extreme depth; interpret effect size and read-level evidence as
  well as the p-value.
- Base and mapping qualities are aligner/platform dependent. The bundled
  thresholds were set on a synthetic Illumina-like panel and audited on real
  reads from one *H. pylori* collection sequenced on one platform
  ([`validation/tier-calibration/`](../validation/tier-calibration/README.md)).
- Bait-derived boundary projection assumes the local HSP and target definition are homologous; remap support must cover the projected ends.
- Reading-frame inference from allele fragments can be underdetermined. Start and stop codons are not required because many MLST definitions are internal gene fragments.
- Length and GC distributions describe the supplied bait set and inherit its sampling bias.
- BLAST similarity to a known allele is not proof of orthology.
- A `SUCCESS` status is a minimum computational gate; the QC tier and underlying evidence remain essential.
- The depth ratio measures collapse, not copy number: a locus whose copies the
  assembler separated shows no pile-up. Use `dosage_estimate` when the question
  is total supported dosage, supply the gene coordinates it requires, and
  inspect `dosage_status` plus `dosage_flags` before interpretation.
- Depth-based evidence says nothing about *which* copy carries which allele.
  Assigning a short read to one physical copy of a repeat requires the read to
  overlap a position at which that copy differs from every other copy; whether
  such positions exist is a property of the locus, computable from a closed
  reference before any sample is sequenced.
