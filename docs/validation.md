# Methods and Validation Guide

This document describes the computational methods implemented in `Locus-Recon` 1.0.0 and provides a practical protocol for validating and calibrating the software on new organisms, sequencing protocols, or allele schemes.

## 1. Computational method

### Step 1: Fragment-tolerant target discovery
For each sample, curated reference alleles (`--bait`) for the target locus are aligned against the draft assembly using `blastn`. Discovery retains High-Scoring Segment Pairs (HSPs) meeting `--min-discovery-identity` (default `80.0%`) and `--min-hsp-length` (default `50` bp). Because discovery gates on aligned fragment length rather than whole-query coverage, a locus fractured across multiple contig ends still seeds read recruitment. Passing HSP coordinates are padded by `--region-padding` (default `300` bp), clipped to contig boundaries, and merged per contig (`target_regions.bed`).

### Step 2: Paired-end read mapping and mate rescue
Paired-end reads are mapped to the draft assembly with `bwa-mem2` (or `bwa`) and coordinate-sorted by `samtools` (`full_map.sorted.bam`). Alignments overlapping `target_regions.bed` define a target read-name set (`target_read_names.txt`), after which `samtools view -N` extracts every alignment in the full BAM sharing those read names (`target_region_rescued.bam`). This rescues unmapped or distally mapped mates whose pairs span assembly gaps or contig ends. Rescued alignments are name-sorted and converted into paired (`target_R1.fq.gz`, `target_R2.fq.gz`) and singleton (`target_S.fq.gz`) FASTQ files.

### Step 3: Local *de novo* assembly
Recruited paired and singleton reads are assembled locally with `spades.py` in default mode (`--phred-offset 33`, using `--threads` and `--memory-per-sample`). If standard SPAdes fails with `Invalid kmer coverage histogram` (common in very small or non-uniform read subsets), `Locus-Recon` automatically cleans the working directory and retries with `--sc` (`single-cell-retry`).

### Step 4: Distinct-region ranking and bait-projected boundary extraction
Reference bait alleles are aligned to the local SPAdes scaffolds (`blast_step5.tsv`). Overlapping hits from multiple bait alleles to the same scaffold interval are collapsed into distinct subject regions before ranking, preventing redundant alleles in a large typing catalogue from masquerading as paralogous loci. Regions are ranked by bitscore, query coverage, percent identity, and aligned length, and the bitscore margin between the best and second-best distinct region is compared against `--min-bitscore-margin` (default `10.0`).

To extract the complete locus rather than a single partial HSP, `describe_full_query_span()` projects the unaligned 5' prefix and 3' suffix of the best-matching bait allele onto the selected local scaffold (`allele_region.bed`) and reverse-complements reverse-strand hits. If the projected locus span runs past either end of the local contig by at least 10 bp, the missing overhang is recorded in `span_clipped_bp` (`span_clipped_start_bp`, `span_clipped_end_bp`) and flagged (`ALLELE_SPAN_CLIPPED_AT_CONTIG_END`). `find_span_continuation()` then checks whether any other candidate contig in the local assembly aligns to the missing query interval; if so, the locus is flagged as split across local contigs (`LOCUS_SPLIT_ACROSS_CONTIGS`, reporting `span_continuation_contig`, `span_continuation_bp`, and `span_continuation_overlap_bp`). `Locus-Recon` never joins separate contigs across an assembly gap.

If the extracted sequence contains interior runs of `N` bases introduced by scaffold gap-filling, `assess_placeholder_junctions()` inspects `spades_local/assembly_graph_after_simplification.gfa` to determine whether the flanking sequences are connected on a continuous graph path (`GRAPH_SUPPORTED`, `AMBIGUOUS`, `NOT_SUPPORTED`, `NOT_ASSESSED`, or `NO_PLACEHOLDER`).

### Step 5: Catalogue comparison and independent read remapping
The candidate sequence is aligned against the shared BLAST database of the reference bait (`blast_validation.tsv`) to identify the `nearest_allele`, `best_identity_pct`, `validation_coverage_pct`, and `exact_known_allele`. Independently, the recruited paired and singleton reads are remapped directly to the reconstructed candidate sequence (`allele_remap.sorted.bam`). A quality-filtered pileup (`samtools mpileup -A -aa`, applying BAQ recalibration, `--min-base-quality`, and `--min-mapping-quality`) evaluates every position of the candidate (`allele_remap.support.tsv`), including zero-depth positions and read pairs not marked as properly paired against a short reference (`-A`).

### Step 6: Orthogonal quality classification
A candidate reaches `workflow_status = SUCCESS` only if it clears all four configurable minimum gates:
- `--min-identity` (default `85.0%`)
- `--min-coverage` (default `80.0%`)
- `--min-remap-depth` (default `5.0×`)
- `--min-remap-breadth` (default `90.0%`)

For `SUCCESS` samples, `Locus-Recon` reports two orthogonal axes:
1. **`sequence_confidence` (`HIGH`, `MEDIUM`, `LOW`, `SUSPECT`) and `result_disposition` (`PASS`, `REVIEW`, `HOLD`):** Determined strictly by intrinsic read support, structural coherence, and span completeness (remap depth and breadth, interior uncertain-base fraction outside `terminal_margin_bp`, interior zero-depth gaps, bidirectional mixed sites, strand-biased alternative sites, reference discordance, distinct-region bitscore margin, N content, GC deviation, reading-frame integrity across all six frames, length relative to the bait IQR/stdev tolerance, and contig-end truncation).
2. **`catalogue_status` (`EXACT_MATCH`, `NONEXACT_MATCH`, `DIVERGENT_FROM_REFERENCE`, `HIGHLY_DIVERGENT_FROM_REFERENCE`):** Describes sequence identity to the nearest allele in the supplied bait database without penalizing `sequence_confidence`. Consequently, a genuinely novel or divergent allele with complete, unambiguous read support achieves `HIGH` confidence (`PASS`) while setting `catalogue_review_recommended = true`.

---

## 2. Depth-based copy number and dosage (`locus-recon-depth-ratio`)

Allelic mixture evidence is silent when duplicated copies of a locus are identical in sequence (for example, a recent tandem amplification of an antimicrobial resistance gene or conserved rRNA operons). When identical copies collapse into a single representative in the draft assembly, reads from all copies pile onto that interval.

`locus-recon-depth-ratio` quantifies this pile-up directly from `<sample>/full_map.sorted.bam`:

- **Matched quality filtering:** Both numerator (target locus intervals) and denominator (single-copy chromosomal backbone) are extracted with `samtools depth -a` using identical base-quality (`--min-bq 20`), mapping-quality (`--min-mq 20`), and SAM flag (`UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY`) filters.
- **Backbone definition:** Contigs shorter than 5,000 bp, 100 bp at each contig terminus, the target `--locus` intervals, and any user-specified `--exclude` intervals (such as plasmids, rRNA operons, or mobile elements) are excluded from the denominator.
- **Block-bootstrap confidence interval:** Because depth along a chromosome is spatially autocorrelated, uncertainty is estimated by resampling contiguous 500 bp blocks (`--n-boot 1000`, `--seed 20260817`), capping backbone resampling at `--max-backbone-bootstrap-blocks 512` evenly spaced blocks while keeping every eligible backbone base in the point median.
- **Multi-copy decision boundary:** `depth_call = MULTICOPY_DEPTH` is assigned when `ratio_ci_low > 1.5` (`--call-threshold 1.5`). Because within-locus block bootstrapping measures local precision rather than genome-wide coverage wave variance, a threshold of `1.5` prevents false rejections of single-copy loci caused by routine chromosomal coverage fluctuations.
- **Mapping uniqueness (`ambiguity_index`):** Computed as filtered median depth divided by unfiltered (`Q0/MQ0`) median depth across the locus. When the draft assembly retains multiple nearly identical copies, short reads map with `MAPQ = 0` and filtered depth drops. If `ambiguity_index < 0.70`, `reliable` is set to `false` and `dosage_status` becomes `LOWER_BOUND` (`LOW_MAPPING_UNIQUENESS`). Because MAPQ filtering can only remove depth, a locus whose filtered lower bound still exceeds `1.5` is reported as `MULTICOPY_DEPTH` with its ratio interpreted as a lower bound.
- **Per-query-base dosage (`dosage_estimate`):** When a locus is fragmented across multiple draft intervals and `--gene-span` plus `--gene-length` are supplied, `query_copy_profile()` computes mean copy multiplicity across covered query bases:
  - Complete, non-overlapping spans give `dosage_status = ESTIMATED`.
  - Incomplete query coverage gives `ESTIMATED_PARTIAL_COVERAGE` along with deterministic `[dosage_lower, dosage_upper]` bounds.
  - Overlapping query spans set `dosage_status = REVIEW_OVERLAPPING_SPANS` (`OVERLAPPING_QUERY_SPANS`) and report `geometry_expected_dosage` (the dosage the same span geometry would yield if every fragment had a depth ratio of `1.0`), preventing overlapping blast hits on a fragmented single-copy gene from being misread as a multi-copy count.

---

## 3. Graph/depth copy-number reconciliation (`locus-recon-copy-number`)

While depth measures read pile-up onto collapsed sequence, the draft assembly graph (`assembly_graph_after_simplification.gfa`) records how many distinct flanking contexts branch from the ends of a collapsed repeat. `locus-recon-copy-number` combines `DepthRatioResult` with bait-anchored GFA traversal (`infer_graph_contexts()` and `reconcile_copy_number()`):

- **Bait-anchored GFA discovery:** Aligns the target bait against GFA segment sequences (`--graph-min-identity 80.0`, `--graph-min-hit-bp max(200, 0.10 * bait_length)`), requiring at least 90% union bait coverage (`--graph-min-query-coverage 0.90`) for an exact integer graph count.
- **Bounded flanking traversal:** Traverses outward from both the 5' (`left`) and 3' (`right`) bait ends across non-target graph links up to `--graph-max-context-bp` (5,000 bp), `--graph-max-context-nodes` (8 nodes), and `--graph-max-context-paths` (64 paths), requiring at least `max(100, 2 * max_link_overlap)` bp of post-branch sequence to register a distinct context signature.
- **Consensus reconciliation:** Preserves the distinction between an integer count of flanking contexts (`INTEGER_CONTEXT_COUNT`) and continuous read-depth dosage (`MEAN_DEPTH_DOSAGE`):
  - When left and right context counts match (`N = L = R >= 1`) without hitting cycles or traversal bounds, and agrees with the depth class/interval, `copy_number_method = GRAPH_DEPTH_CONSENSUS`.
  - When an exact multi-copy graph count (`N >= 2`) falls outside the bootstrap depth interval while depth also calls `MULTICOPY_DEPTH`, the integer context count is reported under `GRAPH_COUNT_OVER_DISCORDANT_DEPTH` with flag `DEPTH_GRAPH_NUMERIC_DISCORDANCE`.
  - When the graph shows a single flanking context (`N = 1`) but depth calls `MULTICOPY_DEPTH` (characteristic of a tandem array collapsed into a single graph loop), `Locus-Recon` preserves the continuous depth dosage under `DEPTH_TANDEM_COMPATIBLE`.
  - When traversal encounters a cycle (`CYCLE_ENCOUNTERED`), exceeds safety bounds (`PATH_LIMIT_REACHED`), or has asymmetric left/right counts (`ASYMMETRIC_ORIENTATION_COUNTS`), or when `ambiguity_index < 0.70`, the result is reported as a `LOWER_BOUND` (`copy_number_lower_bound`) rather than an unsupported point estimate.

---

## 4. Local graph path enumeration and competitive scoring (`locus-recon-graph-paths`)

When the standard `locus-recon` workflow flags a locus as `SUSPECT` due to `ALLELE_MIXTURE`, `AMBIGUOUS_BEST_HIT`, or contig-end truncation, the local SPAdes graph (`<sample>/spades_local/assembly_graph_after_simplification.gfa`) may contain alternative source-to-sink paths that were collapsed or broken in `scaffolds.fasta`.

`locus-recon-graph-paths` enumerates terminal paths through the local GFA (`--min-length`, `--max-length`, `--max-paths`, `--max-nodes`) and, when given `--reads-r1` (and `--reads-r2`), indexes all retained paths together as a single multi-FASTA reference for competitive read alignment:
- Reads that map equally well to shared sequence across multiple paths receive `MAPQ = 0` and are excluded by `--min-mapping-quality` (default `20`).
- Only reads uniquely anchored to path-specific sequence contribute to `unique_mean_depth`, `unique_breadth_pct`, and `unsupported_bp`.
- If no path carries uniquely anchored reads (`unique_breadth_pct == 0` across all paths), the command explicitly reports that the alternative paths cannot be distinguished at the library's read length/insert size and require long-read or hybrid sequencing to phase.

---

## 5. Designing a study-specific validation panel

Before deploying `Locus-Recon` on a new species, typing scheme, or sequencing workflow, run the bundled synthetic benchmark ([`validation/README.md`](../validation/README.md)) to confirm installation integrity, and then evaluate a study-specific control panel covering:

| Control type | Construction | Target behavior |
|---|---|---|
| **Intact positive controls** | Closed or high-quality assemblies with known allele calls and matched paired-end reads | `SUCCESS`, `HIGH` (`PASS`), 100% concordance with reference truth |
| **Synthetic fragmentation** | Split known target loci at central and near-terminal positions in trusted assemblies; supply original paired reads | Complete sequence recovery across contig breaks |
| **Coverage dilution series** | Subsample matched reads across relevant depths (e.g., 40×, 20×, 10×, 5×) | Characterize the sensitivity floor and verify progressive QC tier downgrades (`HIGH` → `MEDIUM` → `LOW` → `SUSPECT`) |
| **Divergent / novel alleles** | Hold out known divergent alleles from the `--bait` file or test lineages at the edge of within-species diversity | Confirm `HIGH` sequence confidence (`PASS`) paired with `NONEXACT_MATCH` or `DIVERGENT_FROM_REFERENCE` catalogue status |
| **Paralog & multi-copy controls** | Genomes carrying a known paralog, 2-copy repeat (e.g., rRNA), or tandem amplification | Verify `AMBIGUOUS_BEST_HIT` / `ALLELE_MIXTURE` gating and evaluate copy number with `locus-recon-copy-number` |
| **Two-strain mixtures** | Blend reads from two isolates with distinct known alleles at 5%, 10%, 20%, and 50% minor fractions | Verify `ALLELE_MIXTURE` detection (`SUSPECT` / `HOLD`) without silent chimerism |
| **Target-negative controls** | Isolates or related species confirmed to lack the target locus | Verify clean `SKIP` at discovery without false-positive reconstruction |

---

## 6. Pre-submission manual review checklist

Before submitting a novel allele to a public nomenclature repository (e.g., PubMLST or BIGSdb) or reporting a clinical/epidemiological conclusion from a `REVIEW` or `NONEXACT_MATCH` call:

1. **Verify sample and read provenance:** Confirm in `run_manifest.json` and `mapping_flagstat.txt` that R1/R2 reads belong to the draft assembly and that global mapping rate is normal.
2. **Check whole-genome contamination:** Confirm upstream species identification and contamination metrics for the isolate.
3. **Inspect `allele_remap.support.tsv`:** Verify that `interior_uncertain_base_count == 0`, `internal_zero_depth_positions == 0`, and every novel polymorphic site has bidirectional read support without strand bias (`strand_biased_site == false`).
4. **Verify locus completeness:** Confirm `span_clipped_bp == 0` and `placeholder_runs == 0`.
5. **Exclude paralog interference:** Inspect `blast_step5.tsv` and `bitscore_margin` to verify that no competing local scaffold region exists.
6. **Confirm reading-frame integrity:** For coding loci, verify `qc_internal_stops == 0` and check that length differences preserve the reading frame (`qc_length_mod3`).
7. **Seek orthogonal confirmation when required:** Use long-read sequencing, hybrid assembly, or targeted PCR/Sanger sequencing when a locus remains `SUSPECT` (`HOLD`) or when short reads cannot phase multi-copy repeats.
