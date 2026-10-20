# Prespecified constructed 1/2/3-copy series

Frozen: 2026-08-27, before generating or analysing the series.

> **Editorial note.** The windowed profile named in the design and in
> criterion 4 is not part of Locus-Recon 1.0.0 (see the editorial note in
> [`../depth-copy-number/PRESPECIFIED_CRITERIA.md`](../depth-copy-number/PRESPECIFIED_CRITERIA.md)),
> so the committed table carries no profile state. The record below is kept as
> written; criteria 1 to 3 are evaluated exactly as prespecified.

This secondary validation isolates depth-estimator behaviour under exact,
known geometry. It does not represent independent biological validation.

## Design

- truth copies: 1, 2, and 3 identical tandem copies;
- nominal haploid depth: 20x, 50x, and 100x;
- locus GC: 0.30, 0.50, and 0.70;
- one deterministic replicate per 27-cell factorial combination;
- 100,000 bp single-copy backbone, one 4,000 bp representative in the mapping
  reference, 150 bp paired reads, and 350 bp fragments;
- fixed base seed 20260817, error-free reads, BWA-MEM2/BWA mapping and samtools
  sorting;
- Q20/MQ20 filters, 500 bp bootstrap blocks, 1,000 resamples, seed
  20260817, one exact 4 kb/20-window calibrated profile;
- no threshold fitting or post-result sample exclusion.

The simulated source carries the stated tandem copy count; reads from it are
mapped to the one-copy representative, modelling repeat collapse. Uniform
fragment sampling deliberately omits library-specific GC bias and sequencing
error. The GC strata therefore test sequence-composition/mapping behaviour,
not a calibrated wet-lab GC-bias model.

## Criteria

1. No one-copy case is classified `MULTICOPY_DEPTH`.
2. Every two- and three-copy case is classified `MULTICOPY_DEPTH`.
3. Median absolute dosage error across all 27 cases is at most 0.25 copies.
4. Every case retains its point estimate, confidence interval, ambiguity,
   profile state, tool versions, seed, parameters, and input checksums.

Any failure remains in the result table and prevents a passing constructed
series; thresholds are not changed.
