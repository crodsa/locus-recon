# Prespecified two-copy 23S rRNA validation

Frozen: 2026-08-27, before running the estimator on this panel.

## Question and evidence boundary

This validation asks whether an Illumina-only analysis can reject a one-copy
model for a locus known from a withheld hybrid-closed chromosome to contain two
full-length copies. Closed chromosomes provide truth only; they are not passed
to Locus-Recon, the draft assembler, or the depth estimator.

Because each hybrid assembly used the same Illumina run together with long
reads, the truth is external to the tested draft/estimator path but is not a
fully independent wet-laboratory copy-number measurement.

The source study is BioProject PRJNA816422. Its publication states that complete
genomes were hybrid assembled with Unicycler v0.4.6 from Illumina and Nanopore
data and that *Helicobacter pylori* has two 23S rRNA copies. Current RefSeq GBFF
structured comments may list only `illumina` as sequencing technology; the
publication is therefore the authority for hybrid-assembly provenance, while
copy count and feature completeness are verified directly from each deposited
closed chromosome.

## Frozen inclusion and selection rules

A sample is eligible only when all conditions hold:

1. one complete circular non-plasmid chromosome is deposited for PRJNA816422
   (additional plasmid replicons are permitted and recorded but excluded from
   chromosome copy-count truth);
2. the study documents hybrid assembly with long and short reads;
3. a paired Illumina WGS run belongs to the identical BioSample;
4. the closed chromosome contains exactly two complete 23S rRNA features;
5. the chromosome contains exactly one complete `gyrB` control CDS;
6. sample identifiers and accession mappings are unambiguous.

Eligible BioSamples are ordered by the numeric strain suffix. The first five
passing isolates are frozen in
`validation/manifests/low_copy_23S_prjna816422.tsv`; no estimator output is used
for selection.

## Tested input and fixed workflow

- Raw paired Illumina reads are downloaded by run accession and checksummed.
- A draft genome is assembled from those reads alone with SPAdes 4.2.0
  `--isolate`; thread and memory values are parameters and are logged.
- A fixed external bait from *H. pylori* strain 26695 is used for 23S rRNA and
  `gyrB`; sample-specific closed truth is never used as a bait.
- Reads are mapped back to the Illumina-only draft with BWA-MEM2 2.2.1 and
  samtools 1.22.1 using recorded parameters.
- Discovery intervals and query spans are obtained without consulting the
  closed-chromosome copy locations.
- Locus depth, dosage geometry, ambiguity, failures, and calibration state are
  retained for every sample and both loci.

## Primary criteria

Thresholds remain the prespecified values; they will not be tuned after
observing this panel.

1. Zero of five single-copy `gyrB` controls may be classified
   `MULTICOPY_DEPTH`.
2. At least four of five two-copy 23S loci (80%) must reject the one-copy model
   as `MULTICOPY_DEPTH`.
3. Every estimate, confidence interval, absolute error, ambiguity metric,
   assembly geometry, dosage status/flags, exclusion, and failure is reported.

Failure of criterion 2 prevents a general claim that two copies are reliably
resolved. It will instead define the observed lower operating boundary.

## Additive graph analysis after depth prespecification

The graph/depth consensus was designed after the first depth-only result had
shown that all five 23S loci rejected one copy but returned continuous ratios of
2.546–3.221. The graph criteria below are therefore a transparent development
acceptance test, not an independent validation added retroactively to the
frozen primary criteria:

1. all five 23S loci must yield matching left/right graph-context counts of two
   and `copy_number_call=2`;
2. all five `gyrB` controls must yield matching counts of one and
   `copy_number_call=1`;
3. all 23 columns in the archived depth-only table must retain identical
   serialized values;
4. graph discovery, traversal limits, flags, consensus states, tool versions,
   parameters and input hashes must be retained in machine-readable outputs.

The result may demonstrate that graph context restores the deposited two-copy
interpretation in these cases. It may not be reported as held-out sensitivity,
universal exact rRNA copy-number recovery, allele-to-copy phasing or physical
replicon closure.

## Secondary constructed series

A deterministic 1/2/3-copy series is evaluated separately across fixed depth
and GC strata with fixed seeds. It is labelled constructed validation and
cannot substitute for the external two-copy panel.

## Source records

- NCBI BioProject PRJNA816422 (accessions and same-BioSample mapping).
- Hu L. et al. *Microbiology Spectrum* 2023;11:e04522-22,
  doi:10.1128/spectrum.04522-22 (hybrid-assembly methods and biological copy
  context).
- The five RefSeq GBFF records in the frozen manifest (direct feature truth).

## Editorial note on the deposited run (added 2026-09-28)

This note follows the run archived in this directory, made with Locus-Recon
1.0.0. The frozen text above is unchanged, and no criterion was relaxed after
the result was seen.

- Primary criteria 1 to 3 are met. All five 23S loci were `MULTICOPY_DEPTH`,
  no `gyrB` control was, and every estimate, interval, error, ambiguity metric,
  geometry, dosage status and flag is reported in `results.tsv`.
- Graph criterion 1 is met in four of five loci. Hpfe0001, Hpfe0002, Hpfe0003
  and Hpfe0006 returned two matched contexts and `copy_number_call=2`. At
  Hpfe0004 both locus ends also yielded two contexts, but the traversal met a
  cycle (`CYCLE_ENCOUNTERED`). Locus-Recon treats a traversal pruned by a cycle
  or a safety limit as incomplete and reports its count as a lower bound, so
  this locus is `TRAVERSAL_LIMIT_REACHED` with a lower bound of two and no exact
  call (`copy_number_kind=LOWER_BOUND`). The bound agrees with the two deposited
  copies but is not the exact call the criterion requires.
- Graph criterion 2 is met: all five `gyrB` controls returned one context and
  `copy_number_call=1`.
- Graph criterion 3 names 23 columns. Two of them, the dosage-profile columns
  `profile_call` and `profile_calibrated`, are not reported by the released
  estimator; the other 21 retain identical serialized values.
- Graph criterion 4 is met by `results.tsv`, `closed_truth_audit.json` and
  `workflow_provenance.json`.

`evaluation.json` combines primary criteria 1 and 2 with graph criteria 1 to 3
in `all_primary_criteria_passed`, which is therefore `false`, and the workflow
exits with status 1 when it reproduces this result.
