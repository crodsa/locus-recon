# Changelog

## 1.0.0 — 2026-09-28

Initial public release of `Locus-Recon`, a reference-guided, assembly-aware toolkit for targeted bacterial locus reconstruction, orthogonal per-base quality adjudication, and read-depth/graph copy-number estimation from paired-end short reads.

### Targeted locus reconstruction (`locus-recon`)

- **Fragment-tolerant target discovery:** Identifies target regions in fragmented draft assemblies from partial BLASTN HSPs (`--min-discovery-identity 80.0`, `--min-hsp-length 50`) rather than requiring full-length contig matches, padding and merging candidate intervals (`--region-padding 300`) before read recruitment.
- **Read recruitment and mate-pair rescue:** Maps paired reads to the draft assembly with `bwa-mem2` (or `bwa` fallback), extracts reads overlapping target intervals, and rescues all alignments sharing those read names via `samtools view -N` (`samtools >= 1.12`) to recover unmapped or distally mapped mates across contig breaks.
- **Local *de novo* SPAdes assembly:** Assembles recruited paired and singleton reads locally in default SPAdes mode with automatic `--sc` (`single-cell-retry`) fallback when low or uneven coverage triggers a k-mer histogram failure.
- **Distinct-region ranking and bait-projected boundaries:** Collapses redundant hits from multi-allele bait catalogues into distinct scaffold intervals before ranking, flags competing regions via `--min-bitscore-margin`, and projects unaligned 5' and 3' bait termini onto the selected scaffold (`describe_full_query_span()`) to recover full-length loci.
- **Independent competitive remapping and BAQ pileup:** Remaps recruited paired and singleton reads to the reconstructed candidate and evaluates every base (`samtools mpileup -A -aa`) for filtered depth, strand balance, base and mapping quality, Fisher's exact strand bias, and bidirectional mixed sites.

### Four-axis orthogonal quality control and disposition routing

- **Separation of workflow completion, sequence support, and catalogue novelty:**
  1. `workflow_status` (`SUCCESS`, `SKIP`, `FAIL`) gates on minimum identity, query coverage, remap depth, and remap breadth.
  2. `sequence_confidence` (`HIGH`, `MEDIUM`, `LOW`, `SUSPECT`) and `result_disposition` (`PASS`, `REVIEW`, `HOLD`) evaluate intrinsic read support, structural coherence, and span completeness without penalizing sequence divergence from the reference catalogue.
  3. `catalogue_status` (`EXACT_MATCH`, `NONEXACT_MATCH`, `DIVERGENT_FROM_REFERENCE`, `HIGHLY_DIVERGENT_FROM_REFERENCE`) and `catalogue_review_recommended` classify novelty relative to the supplied bait alleles.
  4. **Bait-projected span completeness** (`span_clipped_bp`, `span_clipped_start_bp`, `span_clipped_end_bp`, `ALLELE_SPAN_CLIPPED_AT_CONTIG_END`) and **local continuation detection** (`find_span_continuation()`, `LOCUS_SPLIT_ACROSS_CONTIGS`, `span_continuation_contig`, `span_continuation_bp`, `span_continuation_overlap_bp`) prevent truncated or split reconstructions from receiving `HIGH` confidence (`PASS`) despite high depth over the retained bases.
- **Disposition-partitioned FASTA outputs:** Automatically separates reconstructed sequences into `<locus>_accepted_alleles.fasta` (`PASS`), `<locus>_review_required_candidates.fasta` (`REVIEW`), `<locus>_hold_candidates.fasta` (`HOLD`), and `<locus>_reconstructed_candidates.fasta` (all `SUCCESS` candidates).
- **Orientation-invariant reading-frame and placeholder diagnostics:** Evaluates internal stop codons across all six reading frames (`qc_coding_frame_used`), attributes stops caused by interior `N`-run frameshifts (`qc_internal_stops_placeholder_closed`, `PLACEHOLDER_FRAMESHIFT_EXPLAINS_STOPS`), and inspects the local SPAdes GFA (`placeholder_junction_support`) to check whether interior `N` runs are bridged by a graph path.
- **Pre-flight bait database validation:** Rejects FASTA headers exceeding the 50-character `makeblastdb -parse_seqids` limit, flags alignment gap characters (`-`, `.`), warns on probable off-target/paralog records lacking shared 25-mers with full-length alleles, and flags underpowered length profiles (`<5` alleles, `CATALOGUE_LENGTH_PROFILE_UNDERPOWERED`).

### Copy-number, dosage, and assembly-graph commands

- **`locus-recon-depth-ratio`:** Computes quality-matched (`Q20/MQ20`) locus-to-backbone depth ratios from `<sample>/full_map.sorted.bam`, with 500 bp block-bootstrap confidence intervals, mapping uniqueness diagnostics (`ambiguity_index`), multi-copy calling (`ratio_ci_low > 1.5`), and per-query-base dosage estimation (`dosage_estimate`, `dosage_status`, `geometry_expected_dosage`) across multi-region loci.
- **`locus-recon-copy-number`:** Reconciles continuous depth dosage with bounded, bait-anchored flanking-context counts in draft SPAdes GFA 1.x graphs (`INTEGER_CONTEXT_COUNT`, `MEAN_DEPTH_DOSAGE`, `LOWER_BOUND`, `INDETERMINATE`), preserving continuous dosage for tandem arrays (`DEPTH_TANDEM_COMPATIBLE`) and explicitly flagging numeric discordance (`GRAPH_COUNT_OVER_DISCORDANT_DEPTH`) or cycle/limit-truncated graph traversals.
- **`locus-recon-graph-paths`:** Enumerates bounded source-to-sink paths from a local SPAdes GFA graph and optionally ranks candidate paths by uniquely anchored read breadth (`--reads-r1`, `--reads-r2`, `--min-mapping-quality 20`), explicitly reporting when short reads cannot phase alternative paths.

### Reproducibility and installation verification

- **Run provenance:** Writes `run_manifest.json` recording package and Python versions, CLI arguments, SHA-256 checksums of the bait and samplesheet, resolved sample file metadata, and external binary versions.
- **Deterministic 13-sample synthetic benchmark (`validation/`):** Bundled dataset generator, end-to-end runner (`validation/run_mock_validation.py`), evaluator, and reference results covering intact, centrally and terminally fragmented, depth-diluted (3×–20×), mixed-allele (5%–50%), paralogous, and negative controls.
