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
