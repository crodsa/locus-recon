# Changelog

## 1.0.0 — 2026-09-05

Initial public release. Locus-Recon reconstructs a target locus from short
reads when the assembly has broken it across contigs or collapsed
near-identical copies into one consensus, reports how much depth that locus
carries relative to a single-copy backbone, and states how far each result is
supported by the reads.

### Locus reconstruction

- Target discovery from a curated bait set using partial HSPs with length and
  identity criteria, so a locus split across contigs still seeds a
  reconstruction. Redundant bait hits at one region are collapsed before
  paralogue ambiguity is assessed, and unaligned bait termini are projected
  onto the local scaffolds rather than truncating a candidate to a single BLAST
  HSP.
- Local read recruitment, mate rescue and per-region SPAdes assembly, with an
  automatic single-cell-mode retry for low or non-uniform coverage. Samples
  that cannot be recovered are reported as `SKIP` data conditions rather than
  software failures.
- Independent remap validation of every candidate: quality-filtered per-base
  evidence carrying base quality, mapping quality, forward and reverse depth,
  strand balance, allele fractions and Fisher strand-bias tests, including
  rescued singleton reads, written to a per-candidate
  `allele_remap.support.tsv` for inspection.
- Explicit mixed-allele detection requiring alternative support on both strands
  at several sites, with configurable quality, fraction, site-count and
  alternative-depth thresholds. Mixture, strand-bias and
  reconstructed-reference-discordance flags feed the QC scorecard and the batch
  report.
- Nearest-known-allele and exact-match fields, enforced final query coverage as
  well as identity, and coding-length checks against the bait mode rather than
  an assumed multiple of three.

### Sequence confidence and catalogue relationship

- `sequence_confidence` (aliased `qc_confidence`) scores read support only:
  per-base certainty, depth and breadth of competitive remapping, mixture
  evidence, length plausibility and span completeness. The per-base certainty
  verdict carries explicit reason codes and is evaluated over the allele
  interior, outside a terminal margin measured per sample from the alignments.
- The relationship to the bait catalogue is reported on its own axis by
  `classify_catalogue_status()`: `EXACT_MATCH`, `NONEXACT_MATCH` (at or above
  97% identity to the nearest allele), `DIVERGENT_FROM_REFERENCE` (85% to 97%)
  or `HIGHLY_DIVERGENT_FROM_REFERENCE` (below 85%), with the identity to the
  nearest allele and an advisory review recommendation. Catalogue distance does
  not constrain the confidence tier: the catalogue records what has been
  sequenced before, not what is real, so a well-supported reconstruction of an
  under-represented lineage is `PASS`. Three descriptive flags name the
  distance — `CATALOGUE_DIVERGENCE_BELOW_HIGH`, `CATALOGUE_DIVERGENT` and
  `CATALOGUE_HIGHLY_DIVERGENT`.
- `result_disposition` follows the sequence tier; "accepted" is reserved for
  `PASS`.

### Completeness of the reported span

- `describe_full_query_span()` returns the projected locus span together with
  the number of bases each contig boundary could not supply and the query
  intervals the clip failed to cover, so a reconstruction that was silently
  truncated at a contig end cannot be reported as `HIGH` on a flawless per-base
  profile.
- `find_span_continuation()` searches the remaining candidate regions for a
  different contig aligning to the lost interval and reports how many lost
  bases it carries and how much query sequence it shares with the reported
  contig. A clip with no continuation means the sequence is absent from the
  local assembly; a clip with one means the locus is split across contigs and
  the sequence is present but unused.
- `ALLELE_SPAN_CLIPPED_AT_CONTIG_END` reports per-end and total truncation and
  blocks `HIGH`; `LOCUS_SPLIT_ACROSS_CONTIGS` is emitted in addition when a
  continuation exists, names that contig and the shared sequence, and is treated
  as severe, because sequence demonstrably present in the data was dropped.
  Only clips at or above `SPAN_CLIP_THRESHOLDS["min_flag_bp"]` (10 bp) are
  treated as evidence of truncation, since terminal alignments routinely stop a
  base or two short for reasons that carry no completeness information.
- Six report columns: `span_clipped_bp`, `span_clipped_start_bp`,
  `span_clipped_end_bp`, `span_continuation_contig`, `span_continuation_bp`,
  `span_continuation_overlap_bp`.

Locus-Recon does not join the split contigs. Scaffolding across an assembly gap
is a decision about which of several possible layouts is correct; the tool
reports the geometry and leaves the join to the analyst.

### Locus depth and copy number

- `locus-recon-depth-ratio` computes a normalised locus-to-backbone depth ratio
  under the same Q20/MQ20 filters as the per-base evidence table, with a block
  bootstrap confidence interval (500 bp blocks, 1000 resamples, at most 512
  evenly spaced backbone blocks so chromosome-scale runs stay executable) and an
  `ambiguity_index` that marks a ratio as a lower bound when the assembler kept
  copies apart and mapping quality collapsed. The locus GC fraction is reported
  as a description of the locus; no GC-matched denominator is applied.
- `copies_estimate`, a length-weighted copy number combining what the assembler
  resolved with what it collapsed,
  `copies = SUM_r (aligned_gene_bp_r * ratio_r) / gene_length`, reachable with
  `--gene-span` and `--gene-length`. It reduces exactly to the ratio for a
  single-fragment locus and is left at `None`, not guessed, when a locus spans
  several fragments without gene coordinates.
- `MULTICOPY_DEPTH` requires the lower bootstrap bound to clear 1.5, a value
  asserted rather than fitted: over 300 single-copy 7,104 bp backbone segments
  from two real alignments the rule rejected single copy 0 times, with a largest
  observed lower bound of 1.365. Runs record which threshold was used and
  whether it was asserted or calibrated, and a calibration can only make
  rejection harder.
- `locus-recon-copy-number` preserves every depth field and optionally combines
  it with bounded, bait-anchored SPAdes graph-context inference, distinguishing
  integer context counts, continuous mean-depth dosage, lower bounds, evidence
  conflicts and indeterminate states, with deterministic graph/depth
  reconciliation and provenance-rich JSON and TSV output.

### Structural evidence

- The bounded `locus-recon-graph-paths` research command exports alternative
  paths retained in a local SPAdes GFA graph when FASTA scaffolds collapse or
  truncate a duplicated locus, with reverse-complement link handling, path
  spelling under overlap validation, cycle and path-count bounds, length
  filters, graph-depth metadata and a machine-readable path summary.
- A documented conservative multi-copy workflow built on compartment and flank
  anchors, complementary graph branches, competitive read mapping, per-base
  polishing and explicit reporting of unphaseable intervals.

### Reading frame

- Internal stop counts are minimised over all six reading frames, not the
  forward three, and the frame that achieved the minimum is reported in
  `qc_coding_frame_used` (`+1`..`+3`, `-1`..`-3`). A reconstruction inherits the
  orientation of the bait that recruited it, so a bait catalogue supplied
  antisense to the coding strand would otherwise produce stop codons the
  candidate does not have. The start and stop codon checks read the candidate
  in the same frame.
- A length that is not a multiple of three in a candidate whose projected span
  is clipped by a contig boundary is reported as the descriptive
  `FRAME_LENGTH_SHIFT_TRUNCATED`. The truncation itself is scored by
  `ALLELE_SPAN_CLIPPED_AT_CONTIG_END` and measured by `span_clipped_bp`, so the
  same fact does not constrain the tier twice.
- Stop codons produced by the frameshift of a scaffold placeholder are
  attributed to it: `qc_internal_stops_placeholder_closed` gives the count with
  interior placeholders excised, and `PLACEHOLDER_FRAMESHIFT_EXPLAINS_STOPS`
  reports both numbers when it falls. `INTERNAL_STOPS` is still raised, because
  the delivered sequence does contain them.

### Bait database validation

- Record identifiers longer than the 50-character `makeblastdb -parse_seqids`
  limit are rejected by name, before any external program runs, instead of
  surfacing as a raw `makeblastdb` failure.
- Alignment gap symbols are named as such, with the advice to remove the gaps,
  rather than only as invalid DNA symbols.
- Records sharing no 25-mer with the full-length records of the set are
  reported as probable paralogues or unrelated fragments. This is a warning,
  not an error, and the comparison is against the records at least half the
  length of the longest one so that two fragments cannot vouch for each other.
- Length flags carry the number of bait alleles the median and IQR were
  computed from (`qc_length_profile_n`), and a profile below five alleles adds
  the descriptive `CATALOGUE_LENGTH_PROFILE_UNDERPOWERED`.

### Scaffold placeholders and graph evidence

- A run of N inside a candidate comes from the local assembler scaffolding
  across a gap it could not spell. The flanks of each interior run are searched
  in the local assembly graph the same run produced, so the check costs no
  extra alignment: `placeholder_runs`, `placeholder_bp` and
  `placeholder_junction_support` (`GRAPH_SUPPORTED`, `AMBIGUOUS`,
  `NOT_SUPPORTED`, `NOT_ASSESSED`, `NO_PLACEHOLDER`) with a descriptive
  `PLACEHOLDER_JUNCTION_*` flag. The verdict is reported and never scored: the
  tier remains a statement about read support for reported bases, and a
  graph-supported junction is still a junction no read spans.

### Graph path scoring

- `locus-recon-graph-paths --reads-r1` (optionally `--reads-r2`) indexes all
  retained paths together so the reads compete for placement, and scores each
  path by the sequence uniquely anchored reads cover: `mapped_reads`,
  `unique_mean_depth`, `unique_breadth_pct`, `unsupported_bp` and `rank`, with
  `--threads`, `--min-mapping-quality` and `--score-dir`. Without reads the
  command behaves as before; if no aligner is available the unranked summary is
  kept and the reason reported.
- When no path carries uniquely placed reads, the command states that the
  enumerated paths cannot be told apart at that read length and that their
  order is arbitrary, rather than presenting a tie as a preference.

### Reproducibility and validation

- A deterministic 13-sample validation dataset covering exact long-read truth,
  central and near-terminal fragmentation, 3x to 20x depth, 5% to 50% allele
  mixtures, a closely related paralogue and a negative control, with an
  end-to-end benchmark runner, a machine-readable evaluator, regression tests
  and compact reference results.
- A tier calibration over every case in the repository with known truth, giving
  the acceptance rate per prespecified case class and the precision of the top
  tier: all 11 single-copy intact cases accepted and exact, no false accepts in
  37 cases, and 21 of 26 withheld cases exact over the span they reported. The
  closed-genome stage runs end to end on public reads with a fixed external
  bait.
- A deposit reusing the curated `tcdB` alleles and local assembly graph already
  in the repository to validate orientation-independent frame metrics, the bait
  database errors, length-profile reporting, the three placeholder junction
  verdicts, and competitive path scoring both against constructed truth and on
  the published read pair, where the enumerated paths cannot be separated.
- A 14-case deterministic completeness benchmark in which the withheld locus
  bases are known by construction, covering intact, flush, one-sided and
  two-sided truncation, reverse orientation and three split-contig layouts.
- Prespecified criteria, null distributions, estimator predictions, post hoc
  resolution distributions at one, two and 2.5 copies, and full
  external-validation output for the depth module under
  `validation/depth-copy-number/`, each with the script that regenerates it.
- Per-run reproducibility manifests with checksums, parameters, inputs and tool
  versions; deterministic report ordering for parallel runs; dry-run
  validation; bounded numeric command-line arguments; memory-per-sample
  control; safe sample and locus names; duplicate detection; optional
  samplesheet headers; samplesheet-relative paths.
- Top-level Conda dependencies pinned to exact versions in `environment.yml`:
  Python 3.11.13, pip 25.2, tqdm 4.67.1, BLAST+ 2.16.0, BWA-MEM2 2.2.1,
  samtools 1.22.1 and SPAdes 4.2.0. The `test` extra is self-contained and
  selects Biopython 1.85 on Python 3.9 or 1.87 on Python 3.10 and newer,
  matching the declared Python 3.9, 3.11 and 3.12 CI matrix.
