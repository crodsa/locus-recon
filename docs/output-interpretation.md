# Output and interpretation reference

This page documents the files, report fields, confidence system, and review
logic implemented by Locus-Recon 1.0. Start with the
[README result guide](../README.md#read-the-results) if you want a shorter
introduction.

## The four decisions are separate

Locus-Recon reports four related but independent decisions.

1. **Reconstruction status** asks whether the sample completed the minimum
   workflow gates. It is `SUCCESS`, `SKIP`, or `FAIL`.
2. **Sequence confidence** asks whether the reported bases are supported by the
   read evidence, whether the sequence is internally coherent, and whether the
   locus span is complete. It is `HIGH`, `MEDIUM`, `LOW`, or `SUSPECT`, is
   assigned only after `SUCCESS`, and is reported by `sequence_confidence`
   (with `qc_confidence` retained as a compatibility alias).
3. **Result disposition** maps the preceding evidence to an action: `PASS`,
   `REVIEW`, or `HOLD`.
4. **Catalogue status** asks how the reconstructed sequence relates to the
   supplied reference catalogue: `EXACT_MATCH`, `NONEXACT_MATCH`,
   `DIVERGENT_FROM_REFERENCE`, or `HIGHLY_DIVERGENT_FROM_REFERENCE`. It is
   reported by `catalogue_status`, `nearest_allele_identity_pct`,
   `catalogue_review_recommended`, `exact_known_allele`, and `nearest_allele`.

These decisions must not be collapsed into one allele call. A sample can be
`SUCCESS` but `SUSPECT`, or `SUCCESS` and `HIGH` while still not matching a
known allele exactly.

**Sequence confidence and catalogue status are independent by design.** Once a
candidate has cleared the discovery requirement of sufficient similarity to the
bait to be the intended locus, divergence from the nearest curated allele no
longer influences sequence confidence. A reconstruction can be an accurate,
fully supported sequence of a genuinely novel or divergent allele; in that case
it is reported as `HIGH` sequence confidence with a catalogue status indicating
divergence, and the catalogue relationship is what a curator reviews. Sequence
confidence is reduced only by read-evidence, sequence-integrity or
span-completeness observations: unresolved bases, thin or uneven depth,
low base or mapping quality, strand imbalance, mixture, reference discordance,
ambiguous best hits, interior zero-depth positions, coding-frame defects,
excessive ambiguous bases, GC anomalies, length deviation, and truncation of the
projected locus span.

## Output layout

A run using `--locus target_locus` creates a structure like this:

```text
results/target_locus/
├── run_manifest.json
├── locus_recon_batch_target_locus.log
├── locus_recon_report_target_locus.tsv
├── target_locus_accepted_alleles.fasta
├── target_locus_review_required_candidates.fasta
├── target_locus_hold_candidates.fasta
├── target_locus_reconstructed_candidates.fasta
├── target_locus_reconstructed_alleles.fasta  # legacy PASS-only alias
├── bait_database/
│   ├── target_locus.*
│   └── makeblastdb.log
└── isolate_001/
    ├── isolate_001_target_locus_reconstructed.fasta
    ├── isolate_001_target_locus_qc_report.txt
    ├── isolate_001_target_locus_recon.log
    ├── blast_step1.tsv
    ├── target_regions.bed
    ├── mapping_flagstat.txt
    ├── blast_step5.tsv
    ├── allele_region.bed
    ├── blast_validation.tsv
    ├── allele_remap.support.tsv
    └── ... optional retained intermediate files
```

The exact intermediate set depends on where a sample stops and whether
`--cleanup` is enabled. A candidate FASTA can exist for a later `SKIP` because
the sequence is written before the final similarity and remap gates. Only the
batch report determines workflow status. Candidate collections are partitioned
by `result_disposition`; file existence is never a verdict.

## Standalone graph/depth copy-number output

`locus-recon-copy-number` writes a provenance-rich JSON record and optionally a
one-row TSV. It does not add columns to the reconstruction batch report. The
record retains every `DepthRatioResult` field and adds three evidence layers:

| Layer | Principal fields | Interpretation |
|---|---|---|
| Depth | `ratio`, `ratio_ci_low`, `ratio_ci_high`, `depth_call`, `dosage_estimate`, `dosage_status`, `dosage_flags` | Continuous pile-up over the assembly representation. These fields are unchanged from `locus-recon-depth-ratio`. |
| Graph | `left_context_count`, `right_context_count`, `graph_context_count`, `graph_context_lower_bound`, `graph_context_status`, `graph_context_flags` | Bait-anchored flanking contexts retained by the GFA. An exact count requires matching positive counts at both ends. |
| Consensus | `copy_number_call`, `copy_number_kind`, `copy_number_method`, `copy_number_status`, bounds, `consensus_flags` | Deterministic reconciliation that preserves the estimand and uncertainty state. |

Never interpret `copy_number_call` without `copy_number_kind`. An
`INTEGER_CONTEXT_COUNT` describes graph contexts; a `MEAN_DEPTH_DOSAGE`
describes a continuous depth-supported mean. `DEPTH_TANDEM_COMPATIBLE` is the
expected combination when a tandem array has one flanking context but elevated
depth. `DEPTH_GRAPH_NUMERIC_DISCORDANCE` retains disagreement between an exact
graph count and the depth interval without overwriting either value.

The JSON `provenance` object records input paths and SHA-256 hashes, effective
depth and graph parameters, external-tool versions, and the retained BLASTN
target hits. TSV list fields are semicolon-delimited; empty cells and JSON nulls
mean not identified or not applicable, not zero.

## Batch-level files

| File | Purpose | Retain? |
|---|---|---|
| `locus_recon_report_<locus>.tsv` | Stable, tab-separated sample summary and principal interpretation table. | Yes; this is the first file to analyze. |
| `<locus>_accepted_alleles.fasta` | `PASS` candidates only (`SUCCESS` + `HIGH`). | Yes; retain the TSV and provenance too. |
| `<locus>_review_required_candidates.fasta` | `REVIEW` candidates (`MEDIUM` or `LOW`). | Yes while adjudication is ongoing. |
| `<locus>_hold_candidates.fasta` | Reconstructed `HOLD` candidates (`SUSPECT`). | Yes for audit; do not treat as accepted alleles. |
| `<locus>_reconstructed_candidates.fasta` | Every successful reconstruction, irrespective of disposition. | Yes for a complete audit trail. |
| `<locus>_reconstructed_alleles.fasta` | Compatibility alias containing the same PASS-only records as `accepted_alleles`. | Only for legacy consumers. |
| `run_manifest.json` | Program/Python versions, complete command, parsed parameters, platform, bait and samplesheet SHA-256 values, bait profile, resolved sample paths and file metadata, executable paths, and detected tool versions. | Yes; primary provenance record. |
| `locus_recon_batch_<locus>.log` | Main-process setup, progress, dependency, and batch-summary messages. | Yes. |
| `bait_database/` | BLAST nucleotide database built once for the whole batch from the supplied bait FASTA. | Retain when preserving a complete run; it can be regenerated from the checksummed bait. |
| `bait_database/makeblastdb.log` | Output from construction of the shared bait database. | Retain for troubleshooting. |

## Per-sample final and evidence files

| File | What it contains | Interpretation |
|---|---|---|
| `<sample>_<locus>_reconstructed.fasta` | One reconstructed candidate, oriented to the best bait hit. The header records sample, locus, length, nearest bait record, identity, coverage, exact-known status, and strand. | Use only with the batch row. Its existence alone does not mean `SUCCESS`. |
| `<sample>_<locus>_qc_report.txt` | Human-readable scorecard for length, similarity, remap support, per-base evidence, composition, coding-frame behavior, flags, and confidence. | Written for successful samples after QC classification. |
| `<sample>_<locus>_recon.log` | Commands, progress, warnings, and errors for that sample. | First diagnostic file for `SKIP`, `FAIL`, or an unexpected metric. |
| `blast_step1.tsv` | All bait-to-draft BLAST HSPs before the implemented discovery identity/length filtering. | Shows whether fragments were present but below discovery thresholds. |
| `target_regions.bed` | Padded, contig-clipped, merged target intervals in 0-based, half-open BED coordinates. | Defines where read-name recruitment was seeded. |
| `mapping_flagstat.txt` | `samtools flagstat` summary for mapping all input reads to the draft assembly. | Helps detect an assembly/read mismatch or globally poor mapping. |
| `blast_step5.tsv` | Bait-to-local-assembly BLAST hits used for candidate-region ranking and ambiguity assessment. | Inspect for alternative local regions or truncated hits. |
| `allele_region.bed` | Projected candidate interval on the selected local SPAdes scaffold. | Connects the final FASTA to its local-assembly coordinates. |
| `blast_validation.tsv` | Candidate-to-known-allele BLAST results, up to five targets, used to choose the nearest bait record and enforce final similarity/coverage gates. | Distinguishes exact known matches from closest-known similarity. |
| `allele_remap.support.tsv` | One quality-filtered evidence row per candidate position, including zero-depth positions. | Principal file for manual base-level review. Preserved by `--cleanup`. |

## Optional intermediates retained without `--cleanup`

The following files support alignment-level or assembler-level investigation:

| File or directory | Role |
|---|---|
| `assembly_input.fasta` and BWA/BLAST indexes | Per-sample copy of the supplied assembly and its temporary indexes. |
| `full_map.sorted.bam` and `.bai` | All paired reads mapped to the draft assembly. |
| `bwa_pipe.log` | Aligner and sorting stream messages for full-read mapping. |
| `target_region.bam` | Alignments overlapping the selected BED regions. |
| `target_read_names.txt` | Unique read names used to seed mate rescue. |
| `target_region_rescued.bam` | All alignments in the full BAM sharing a selected read name. |
| `target_namesorted.bam` | Name-sorted rescued alignments used to regenerate FASTQ. |
| `target_R1.fq.gz`, `target_R2.fq.gz`, `target_S.fq.gz` | Recovered paired reads and singletons supplied to local assembly and final remapping. |
| `spades_local/` | Complete local SPAdes working directory; `scaffolds.fasta` is the candidate-search assembly. |
| `local_db.*` | Temporary BLAST database for the local SPAdes scaffolds. |
| `allele_remap.sorted.bam` and `.bai` | Paired and singleton selected reads remapped to the reconstructed candidate. |
| `allele_remap.flagstat.txt` | `samtools flagstat` output for final-candidate remapping. |
| `allele_remap_bwa.log`, `allele_remap_singleton_bwa.log` | Aligner messages for paired and singleton remapping. |

`allele_remap.paired.bam` and `allele_remap.singleton.bam` are transient remap
parts. In a normal completed remap they are moved or merged into
`allele_remap.sorted.bam` and removed; one can remain after an interrupted step.

See [cleanup and reruns](cli-reference.md#cleanup-and-reruns) for the exact
files removed after a successful sample.

## Batch report columns

Columns are emitted in the stable order below. Cells without applicable QC data
are empty for `SKIP` and `FAIL` rows.

### Identity and run context

| Column | Meaning |
|---|---|
| `sample_id` | Unique samplesheet identifier and per-sample directory name. |
| `status` | Compatibility alias identical to `workflow_status`. |
| `workflow_status` | `SUCCESS`, `SKIP`, or `FAIL`. See [status semantics](#status-semantics). |
| `result_disposition` | `PASS`, `REVIEW`, or `HOLD`; controls the candidate collection into which a successful reconstruction is written. |
| `locus` | Value passed to `--locus`. |
| `tool_version` | Locus-Recon package version. |
| `run_date` | Batch start date/time recorded by the CLI. |

### Reconstructed sequence and nearest known allele

| Column | Meaning |
|---|---|
| `allele_length` | Candidate length in base pairs; may be present for a later `SKIP`. |
| `nearest_allele` | Subject identifier of the best validation hit in the supplied bait database. It is a nearest reference, not automatically an assigned allele. |
| `exact_known_allele` | `true` only when identity and query coverage are both 100% and aligned, query, and subject lengths are all equal. |
| `best_identity_pct` | Percent identity to the selected nearest bait allele during final validation. |
| `nearest_allele_identity_pct` | The same percent identity, reported in the catalogue-status block. It does not constrain sequence confidence. |
| `catalogue_status` | `EXACT_MATCH`, `NONEXACT_MATCH` (≥97% but not exact), `DIVERGENT_FROM_REFERENCE` (85–97%), or `HIGHLY_DIVERGENT_FROM_REFERENCE` (<85%). |
| `catalogue_review_recommended` | `true` whenever the sequence is not an exact catalogue match, so a curator decides whether it is a new allele. |
| `sequence_confidence` | Sequence confidence tier; identical to `qc_confidence`, named explicitly to separate it from catalogue status. |
| `validation_coverage_pct` | Percentage of the reconstructed query covered by the selected validation hit. |

### Discovery, recruitment, and assembly

| Column | Meaning |
|---|---|
| `contigs_hit` | Number of draft contigs with discovery HSPs that passed identity and aligned-length filters. |
| `target_regions` | Number of padded target intervals after overlapping/adjacent intervals were merged per contig. |
| `extracted_read_names` | Number of unique read names recruited from target-overlapping alignments before all alignments with those names were rescued. |
| `spades_mode` | `standard` or `single-cell-retry`; empty if local assembly was not reached or completed. |
| `mapped_reads` | Number of reads reported as mapped to the full draft assembly. |

### Overall QC and sequence profile

| Column | Meaning |
|---|---|
| `qc_confidence` | `HIGH`, `MEDIUM`, `LOW`, or `SUSPECT`; populated only for `SUCCESS`. |
| `qc_length_delta` | Candidate length minus the bait-set median length, in bp. |
| `qc_length_zscore` | Candidate length relative to bait mean and standard deviation; when bait length variance is effectively zero, a difference greater than one bp is represented by `999.0`. |
| `qc_gc_pct` | Candidate GC percentage. |
| `qc_gc_deviation` | Absolute percentage-point difference between candidate GC and bait-set mean GC. |
| `qc_n_count` | Number of `N` bases in the candidate. |
| `qc_internal_stops` | Minimum internal-stop count across all six reading frames; zero when coding checks are disabled. Evaluating both strands keeps an antisense bait catalogue from producing stop codons the candidate does not have. |
| `qc_coding_frame_used` | The frame that achieved that minimum (`+1`..`+3`, `-1`..`-3`), and the frame in which the start and stop codon checks were read. |
| `qc_length_mod3` | Candidate length modulo three for coding-locus checks. |
| `qc_length_profile_n` | Number of bait alleles the length median and IQR were computed from. Below five, length flags are accompanied by `CATALOGUE_LENGTH_PROFILE_UNDERPOWERED` and should be read as weakly determined. |
| `placeholder_runs`, `placeholder_bp` | Interior runs of N in the candidate, and their total length. These come from the local assembler scaffolding across a gap it could not spell. |
| `placeholder_junction_support` | Whether the flanks of those runs are contiguous on a path through the local assembly graph: `GRAPH_SUPPORTED`, `AMBIGUOUS` (the graph spells extra sequence across the gap), `NOT_SUPPORTED` (no path carries both flanks), `NOT_ASSESSED` (no readable graph), `NO_PLACEHOLDER`. Reported, never scored: a graph-supported junction is still a junction no read spans. |
| `qc_flags` | Semicolon-separated reasons for confidence downgrading. |

### Final read support and possible mixture

| Column | Meaning |
|---|---|
| `remap_mean_depth` | Mean quality-filtered depth across the complete candidate length, including uncovered positions. |
| `remap_breadth_pct` | Percentage of candidate positions with at least one quality-filtered observation. |
| `remap_pct_bases_lt5` | Percentage of positions below 5× quality-filtered depth. Descriptive only; does not affect the confidence tier. |
| `uncertain_base_count`, `uncertain_base_fraction` | Positions across the full candidate length at which the reported base is not adequately supported. |
| `interior_uncertain_base_count`, `interior_uncertain_base_fraction` | The same count restricted to positions more than `terminal_margin_bp` from either end. This fraction is the tier-affecting statistic. |
| `internal_zero_depth_positions` | Interior positions with no quality-filtered observations. |
| `terminal_margin_bp` | Mean aligned read length measured from the remap BAM, used as the interior/edge boundary. Measured per sample, not configured. |
| `span_clipped_bp` | Bases the projected locus span lost to a contig boundary; the amount by which the reported allele is truncated. Split per end as `span_clipped_start_bp` and `span_clipped_end_bp`. |
| `span_continuation_contig`, `span_continuation_bp`, `span_continuation_overlap_bp` | If another candidate contig carries the lost sequence: its name, how many of the lost bases it carries, and the query sequence it shares with the reported contig. Empty when no continuation was found. |
| `mean_base_quality` | Mean Phred base quality across included pileup observations. |
| `mean_mapping_quality` | Mean Phred mapping quality across included pileup observations. |
| `strand_balance_pct` | Overall forward/reverse observation balance; 100 is equal and 0 is entirely one-stranded. |
| `candidate_mixed_sites` | Positions meeting configured alternative-depth and alternative-fraction requirements before strand filtering. |
| `mixed_site_count` | Candidate mixed positions with alternative observations on both strands and without significant implemented strand bias. |
| `strand_biased_sites` | Candidate mixed positions with one-strand alternative support or Fisher strand-bias `p < 0.01`. |
| `reference_discordant_sites` | Positions at which the dominant quality-filtered base nearly fixes an alternative to the reconstructed reference under the configured mixture-fraction boundary. |
| `max_alt_fraction` | Largest alternative fraction observed at any candidate position. |
| `median_mixed_fraction` | Median alternative fraction across credible bidirectional mixed sites; descriptive, not a calibrated organism-abundance estimate. |
| `mixture_detected` | `true` when `mixed_site_count` reaches `--mixture-min-sites`. |

### Ambiguity and bookkeeping

| Column | Meaning |
|---|---|
| `bitscore_margin` | Best local-region BLAST bitscore minus the second-best distinct-region bitscore. A value below `--min-bitscore-margin` is ambiguous. `999999.0` means no second distinct region was available for comparison. |
| `elapsed_sec` | Per-sample wall-clock time in seconds. |
| `allele_file` | Path to the per-sample candidate FASTA. Its presence is not equivalent to `SUCCESS`. |
| `message` | Concise success summary, skip reason, or failure message. |

## Status semantics

### `SUCCESS`

The candidate passed all four configurable minimum gates:

- final identity ≥ `--min-identity` (default 85%);
- final query coverage ≥ `--min-coverage` (default 80%);
- remap mean depth ≥ `--min-remap-depth` (default 5×); and
- remap breadth ≥ `--min-remap-breadth` (default 90%).

`SUCCESS` does not require `HIGH` QC. Mixture, ambiguity, structural, quality,
and per-base certainty flags are evaluated after these gates and can make a
successful candidate `SUSPECT`. Coverage-patchiness flags are descriptive and
cannot by themselves change the tier.

### `SKIP`

`SKIP` is an expected no-result condition for that sample rather than an
external-command or unhandled software failure. Implemented examples include:

- no draft contig passed bait-fragment discovery;
- no reads overlapped selected target regions;
- no reads remained after mate rescue;
- standard and single-cell SPAdes could not model the specific low/non-uniform
  coverage condition handled by the pipeline;
- SPAdes produced no scaffolds;
- no bait hit was found in the local assembly;
- the candidate failed final identity or query-coverage gates; or
- the candidate failed final remap mean-depth or breadth gates.

Because a `SKIP` can reflect data limitations, divergence, or real absence, it
must be interpreted with its `message`, log, and upstream assembly/read QC.

### `FAIL`

`FAIL` represents a missing/empty/unreadable input recorded at batch validation,
a failed external command within a sample, a worker exception, or another
unexpected error. Batch-wide dependency and shared-database failures terminate
the command rather than generating normal per-sample `FAIL` rows.

## QC confidence logic

The bait set defines expected length and GC behavior. Length tolerance is the
bait interquartile range (IQR) when nonzero; otherwise it is the greater of the
bait length standard deviation and 3 bp. Coding checks infer plausible forward
frames by finding the frame or frames with the fewest internal stops across the
bait records. Start and stop codons are not required because many typing loci
are internal gene fragments.

Fixed flag boundaries are summarized below. These confidence thresholds
are distinct from the configurable minimum `SUCCESS` gates.

| Evidence | No high-tier flag | Medium boundary | Low boundary; worse is severe |
|---|---:|---:|---:|
| Identity to nearest bait | ≥97% | ≥93% | ≥85% |
| Validation query coverage | ≥95% | ≥90% | ≥80% |
| Fraction of `N` bases | ≤0.1% | ≤0.5% | ≤1.0% |
| GC deviation from bait mean | ≤3 percentage points | ≤5 pp | ≤8 pp |
| Remap mean depth | ≥20× | ≥10× | ≥5× |
| Remap breadth | ≥99% | ≥97% | ≥90% |
| Positions below 5× (descriptive; see below) | ≤1% | ≤5% | ≤15% |
| Interior uncertain bases | ≤1% | ≤5% | ≤15% |
| Mean base quality | ≥Q30 | ≥Q25 | ≥Q20 |
| Mean mapping quality | ≥Q40 | ≥Q30 | ≥Q20 |
| Overall strand balance | ≥40% | ≥25% | ≥10% |

Additional severe evidence includes an internal stop, length-modulo frame
shift, ambiguous distinct region, detected allele mixture, strand-biased
alternative support, near-fixed reference discordance, or an interior
position with no read support at all (`INTERNAL_ZERO_DEPTH_GAP`), which can
indicate a chimeric join.

### Completeness is checked separately from per-base support

Every statistic in the table above describes the bases that were reported. None
of them can describe bases that were never reported, so a truncated allele can
carry a flawless per-base profile. Completeness is therefore asserted at the
point where the locus is projected onto its contig: the full bait span is
projected from the alignment, and if that projection runs past either end of
the contig, the overhang is sequence the reconstruction cannot contain.

The overhang is reported as `span_clipped_bp` and, above 10 bp, raises
`ALLELE_SPAN_CLIPPED_AT_CONTIG_END`, which blocks `HIGH`. The 10 bp floor
exists because terminal alignments routinely stop a base or two short of the
query end for reasons that carry no completeness information.

Locus-Recon then asks whether the lost sequence exists elsewhere in the local
assembly. If another candidate contig aligns to the lost query interval, the
locus is split across contigs rather than absent from the data, and
`LOCUS_SPLIT_ACROSS_CONTIGS` names that contig (`span_continuation_contig`),
how many of the lost bases it carries (`span_continuation_bp`) and how much
query sequence it shares with the reported contig
(`span_continuation_overlap_bp`) — the overlap available to anchor a join.
That case is treated as severe, because sequence demonstrably present in the
data was dropped from the result.

**Locus-Recon does not join the contigs.** Choosing a layout across an
assembly gap is a scaffolding decision, and the tool reports assembly geometry
rather than resolving it. The split is surfaced with the evidence needed to
make the join, and the join is left to the operator.

Positions below 5× are reported but **descriptive**: they do not enter any
tier decision. On a short extracted-and-remapped locus this statistic is
largely fixed by read length, mean depth and locus length rather than by
whether the reconstructed bases are correct — in the mock benchmark every
sub-5× position in every successful pure control lies at a contig terminus.
The tier-affecting evidence statistic is instead the **interior uncertain-base
fraction**: the fraction of positions more than one mean read length from
either end at which the reported base is not adequately supported (see
*Per-base certainty* below). Positions within one read length of an end are
structurally depth-limited, because no read can begin outside the
reconstructed allele.

Classification is rule-based and explainable:

- no tier-affecting flags gives `HIGH`; a `HIGH` result may still carry
  flags marked `descriptive`;
- one modest flag with acceptable length behavior generally gives `MEDIUM`;
- up to two modest flags with broader acceptable length behavior can give
  `LOW`;
- any severe flag, or a larger combination of deviations, gives `SUSPECT`;
- an isolated length deviation with no severe flag is assigned `MEDIUM`,
  independently of catalogue identity, and is additionally marked
  `POSSIBLE_NOVEL_ALLELE` when identity is at least 97%; and
- catalogue-relationship flags (`CATALOGUE_*`, `POSSIBLE_NOVEL_ALLELE`) are
  descriptive and never enter the tier decision.

These labels prioritize review. They are not calibrated probabilities.

## QC flag reference

Flag text contains the observed value in parentheses. The table groups exact
implemented flag names by the evidence they represent.

| Evidence family | Implemented flags | Meaning |
|---|---|---|
| Length | `LENGTH_MARGINAL`, `LENGTH_DEVIANT`, `LENGTH_ANOMALOUS` | Candidate length is outside 1×, 2×, or 3× the effective bait tolerance around the median. |
| Possible length variant | `POSSIBLE_NOVEL_ALLELE` | Length differs, identity is at least 97%, and no severe flag was present. Descriptive: a review hypothesis, not a confirmed new allele, and not a tier determinant. |
| Catalogue relationship | `CATALOGUE_DIVERGENCE_BELOW_HIGH`, `CATALOGUE_DIVERGENT`, `CATALOGUE_HIGHLY_DIVERGENT` | Identity to the nearest catalogue allele is below 97%, 93%, or 85%, respectively. Descriptive: these annotate novelty and never constrain sequence confidence. |
| Validation coverage | `COVERAGE_BELOW_HIGH`, `REDUCED_COVERAGE`, `LOW_COVERAGE` | Query coverage is below 95%, 90%, or 80%, respectively. A `LOW_COVERAGE` successful row normally cannot occur with default minimum gates but can occur after threshold changes. |
| Ambiguous bases | `SOME_N_CONTENT`, `ELEVATED_N_CONTENT`, `HIGH_N_CONTENT` | `N` fraction exceeds 0.1%, 0.5%, or 1.0%. |
| GC composition | `GC_BELOW_HIGH`, `GC_SHIFT`, `GC_ANOMALY` | Absolute GC deviation exceeds 3, 5, or 8 percentage points. |
| Coding structure | `INTERNAL_STOPS`, `FRAME_LENGTH_SHIFT` | The best of all six reading frames contains internal stops, or candidate length modulo three differs from the bait mode. Disabled by `--noncoding-locus`. |
| Coding structure, descriptive | `FRAME_LENGTH_SHIFT_TRUNCATED` | The length is not a multiple of three in a candidate already known to be clipped by a contig boundary, where the arithmetic follows from the truncation. Does not constrain the tier; the truncation itself is scored by `ALLELE_SPAN_CLIPPED_AT_CONTIG_END`. |
| Scaffold placeholder, descriptive | `PLACEHOLDER_JUNCTION_GRAPH_SUPPORTED`, `PLACEHOLDER_JUNCTION_AMBIGUOUS`, `PLACEHOLDER_JUNCTION_NOT_SUPPORTED`, `PLACEHOLDER_JUNCTION_NOT_ASSESSED` | Graph evidence about the join across an interior run of N, with the coordinates of the run and how many graph paths carry both flanks. Does not constrain the tier. |
| Mean remap depth | `REMAP_DEPTH_BELOW_HIGH`, `REDUCED_REMAP_DEPTH`, `LOW_REMAP_DEPTH` | Mean depth is below 20×, 10×, or 5×. |
| Remap breadth | `REMAP_BREADTH_BELOW_HIGH`, `REDUCED_REMAP_BREADTH`, `LOW_REMAP_BREADTH` | Breadth is below 99%, 97%, or 90%. |
| Coverage patchiness (**descriptive**) | `REMAP_PATCHINESS_BELOW_HIGH`, `MODEST_REMAP_PATCHINESS`, `PATCHY_REMAP_SUPPORT` | More than 1%, 5%, or 15% of candidate positions are below 5×. Reported with a `descriptive` annotation and excluded from every tier decision. |
| Per-base certainty | `UNCERTAIN_BASES_BELOW_HIGH`, `ELEVATED_UNCERTAIN_BASES`, `WIDESPREAD_UNCERTAIN_BASES` | More than 1%, 5%, or 15% of interior positions are uncertain. Tier-affecting. |
| Interior support gap | `INTERNAL_ZERO_DEPTH_GAP` | One or more interior positions have no quality-filtered observations. Severe. |
| Allele truncated by a contig boundary | `ALLELE_SPAN_CLIPPED_AT_CONTIG_END` | The projected locus span ran off the end of the contig, so the reported allele is short by the stated number of bases. Blocks `HIGH`. |
| Locus split across contigs | `LOCUS_SPLIT_ACROSS_CONTIGS` | The clipped bases align to another candidate contig, which is named together with the query sequence it shares with the reported contig. Severe. |
| Distinct-region ambiguity | `AMBIGUOUS_BEST_HIT` | The best and second-best distinct local regions are separated by less than `--min-bitscore-margin`. Redundant bait hits to one region are collapsed first. |
| Base quality | `BASE_QUALITY_BELOW_HIGH`, `REDUCED_BASE_QUALITY`, `LOW_BASE_QUALITY` | Mean base quality is below Q30, Q25, or Q20. |
| Mapping quality | `MAPPING_QUALITY_BELOW_HIGH`, `REDUCED_MAPPING_QUALITY`, `LOW_MAPPING_QUALITY` | Mean mapping quality is below Q40, Q30, or Q20. |
| Overall strand balance | `STRAND_BALANCE_BELOW_HIGH`, `STRAND_IMBALANCE`, `SEVERE_STRAND_IMBALANCE` | Strand balance is below 40%, 25%, or 10%. |
| Mixed reads | `ALLELE_MIXTURE` | The configured number of credible bidirectional mixed positions was reached. |
| Strand artifact risk | `STRAND_BIASED_ALT_SUPPORT` | One or more candidate alternative sites were one-stranded or had Fisher strand-bias `p < 0.01`. |
| Candidate/reference conflict | `REFERENCE_DISCORDANCE` | One or more positions were nearly fixed for an alternative base relative to the reconstructed candidate. |

## Per-base support columns

`allele_remap.support.tsv` is generated by `samtools mpileup -aa` after applying
`--min-base-quality` and `--min-mapping-quality`. It therefore includes
zero-depth positions and has one row per reconstructed coordinate.

| Column | Meaning |
|---|---|
| `contig` | Reconstructed FASTA sequence identifier. |
| `position` | 1-based position reported by mpileup. |
| `reference_base` | Base in the reconstructed candidate. |
| `depth` | Number of parsed quality-filtered nucleotide observations. |
| `forward_depth`, `reverse_depth` | Observations on each alignment strand. |
| `strand_balance_pct` | Per-position forward/reverse balance, from 0 to 100. |
| `ref_forward`, `ref_reverse` | Reconstructed-reference observations by strand. |
| `alt_base` | Most abundant non-reference base or indel representation. |
| `alt_forward`, `alt_reverse` | Alternative observations by strand. |
| `alt_fraction` | `alt_depth / depth`. |
| `mean_base_quality` | Mean base quality at the position. |
| `mean_mapping_quality` | Mean read-mapping quality at the position. |
| `strand_bias_pvalue` | Two-sided Fisher exact p-value comparing reference and alternative strand counts. |
| `candidate_mixed_site` | Meets configured depth and alternative-fraction bounds before strand filtering. |
| `credible_mixed_site` | Candidate site with alternative support on both strands and no implemented significant strand bias. |
| `strand_biased_site` | Candidate site with one-strand alternative support or Fisher `p < 0.01`. |
| `reference_discordant_site` | Alternative fraction is greater than `1 - mixture_min_fraction`. |
| `is_low_coverage` | Quality-filtered depth is below 5×. A coverage statement only. |
| `is_uncertain` | The reported base is not adequately supported at this position. |
| `uncertainty_reason` | Semicolon-separated reason codes when `is_uncertain` is true; empty otherwise. |

### Per-base certainty

`is_low_coverage` and `is_uncertain` are distinct and are never conflated. Low
filtered depth is compatible with a certain base: at 4× with concordant,
adequately scored observations seen from both strands the base is resolved.
Conversely, a position at 20× split evenly between two bases is uncertain
despite ample coverage. The certainty model is therefore not a depth cutoff.

| Reason code | Condition |
|---|---|
| `ZERO_DEPTH` | No quality-filtered observations. |
| `INSUFFICIENT_FILTERED_DEPTH` | Fewer than two quality-filtered observations, so no observation can be corroborated. |
| `NO_BIDIRECTIONAL_SUPPORT` | Below the bidirectional-support depth threshold, the dominant base was seen on one strand only. |
| `EXCESS_ALTERNATIVE_SUPPORT` | Alternative fraction is at or above the maximum compatible with a single resolved base. |
| `CREDIBLE_MIXED_SITE` | The position is a credible bidirectional mixed site, so a single base cannot be asserted. |
| `LOW_BASE_QUALITY` | Mean base quality at the position is below `--min-base-quality`. |
| `LOW_MAPPING_QUALITY` | Mean mapping quality at the position is below `--min-mapping-quality`. |

The default candidate mixed-site requirements are alternative fraction from
0.05 through 0.95 and at least two alternative observations. The default
locus-level call requires at least two credible mixed sites. These values must
be calibrated for the library protocol, depth, organism, and purpose.

## Practical review paths

### Exact known result

1. Confirm `status=SUCCESS` and inspect `qc_confidence` and `qc_flags`.
2. Confirm `exact_known_allele=true` and the expected `nearest_allele`.
3. Check remap breadth, low-depth termini, mixture fields, and upstream sample
   identity/QC.
4. Apply the typing scheme's own acceptance rules.

An exact bait match can still be `LOW` or `SUSPECT` if read support or other
evidence is weak.

### Novel-looking result

`exact_known_allele=false` means only that the sequence did not exactly match a
bait record. It does not prove novelty, and on its own it does not reduce
sequence confidence: a divergent sequence with complete, unambiguous read
support is reported as `HIGH` with a divergent `catalogue_status`, which routes
the decision to a curator rather than treating novelty as a defect.

1. Inspect `catalogue_status`, `nearest_allele_identity_pct`, query coverage,
   length, and `POSSIBLE_NOVEL_ALLELE` or other flags.
2. Verify the supplied bait database is current and has the correct boundaries.
3. Review `blast_validation.tsv`, `allele_region.bed`, the local scaffold, and
   every discordant/low-depth position in `allele_remap.support.tsv`.
4. Exclude paralogy, mixture, contamination, and sample/assembly mismatch.
5. Confirm with an independent assembly or sequencing method when the result is
   submission-bound or changes a scientific conclusion.

### Ambiguous distinct regions

`AMBIGUOUS_BEST_HIT` compares different local-assembly regions after multiple
bait alleles hitting the same region have been collapsed. Inspect the distinct
scaffold regions in `blast_step5.tsv`; consider paralogy, repeated sequence,
contamination, or competing assemblies. Do not select a region solely because
it appears first in the report.

### Weak or patchy remap support

Filter `allele_remap.support.tsv` on `is_uncertain=true` first, and read
`uncertainty_reason` — this separates positions that are merely shallow from
positions where the base is genuinely unresolved. Compare the result against
`terminal_margin_bp`: uncertainty confined to within one read length of either
end is the expected terminal coverage ramp and does not indicate a
reconstruction problem. Interior uncertainty does. Interior troughs, and in
particular any `internal_zero_depth_positions`, can indicate a repeat, indel,
chimeric join, misassembly, or mixture; inspect the BAM and the SPAdes
graph/files when intermediates were retained.

### Mixture or strand-biased alternatives

For `ALLELE_MIXTURE`, inspect all `credible_mixed_site=true` rows and their
fractions. For `STRAND_BIASED_ALT_SUPPORT`, inspect read ends, alignment
positions, mapping quality, and strand counts. Run upstream whole-genome
contamination or mixed-isolate QC. Locus-Recon does not phase the component
alleles or estimate their organism-level abundance.

## Interpretation boundaries

- BLAST similarity to a known allele is not proof of orthology.
- A local assembler can make a plausible consensus from a mixture or paralog.
- Fisher strand tests have limited power at low depth and high sensitivity at
  extreme depth; inspect counts and effect sizes as well as p-values.
- Bait-derived boundary projection assumes homologous target definitions;
  verify read support across projected ends.
- Bait length and GC distributions inherit the content and sampling bias of the
  supplied allele set.
- No-detection remains uncertain until divergence, read support, assembly
  quality, and biological absence have been distinguished.

For study-level validation and manual-review requirements, continue to the
[methods and validation guide](validation.md).
