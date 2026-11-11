# Validation results of release 1.0.0

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
| Tier calibration, 62 cases | `tier-calibration/` | What an accepted result is worth | 25 accepted, all identical to the truth (nominal exact 95% interval 0.86 to 1.00) | The cases are not independent |
| Checks of single components | `frame-and-graph-evidence/` | Frame search, bait validation, length profile, runs of N, path ranking | Recorded outcome of each check | The checks use constructed inputs |
| Synthetic copy series, 27 data sets | `constructed-copy-series/` | Depth arithmetic and bootstrap for one to three copies | Median absolute error 0.04 copies; all classes correct | Coverage bias of real libraries; reads are error-free and uniform |
| *aphA1* in *Acinetobacter baumannii*, seven runs | `aphA1-copy-number/` | Depth dosage of a real amplification | 10.5 and 77.7 copies against 10 ± 2 and 75 ± 14 by published quantitative PCR | Six runs come from one patient |
| 23S rRNA gene and *gyrB*, five *H. pylori* isolates | `low-copy-23S/` | A small copy number | Depth ratios of 2.55 to 3.22 for two copies; graph count of two in four isolates, lower bound in the fifth | These genomes also guided the design of the graph method |
| *Clostridioides difficile* LIBA-6656, *tcdB* | `liba6656-gfa-rerun/`, `liba6656-depth/`, `depth-copy-number/overlap_null/` | A locus with two divergent copies that collapse in the draft | The workflow held a 1,955-bp fragment; the query-weighted dosage of 2.29 is explained by overlapping alignments | Two full-length sequences were built by hand from the graph. They are hypotheses that await long reads |

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

The outputs were produced during development and their records name the commits at which they were run. The link
to the release is the SHA-256 digest of the package sources (`locus_recon/*.py`),
`3863abc32fb5446a72919c790cb2558217871410f441db1ecf17dde3de55cdcf`, which `validation/common.py` computes.
The report, the scorecard, the checksum list and the verification record of `liba6656-gfa-rerun/` are the
output of the released code; see the [changelog](../CHANGELOG.md).
