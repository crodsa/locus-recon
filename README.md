# Locus-Recon

**Reference-guided reconstruction, orthogonal quality adjudication, and copy-number estimation for fragmented bacterial loci from short reads**

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](CHANGELOG.md)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue.svg)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

`Locus-Recon` is a command-line toolkit for bacterial genomics and molecular typing (MLST, cgMLST, wgMLST, virulence factors, and antimicrobial resistance determinants). When a target locus is broken across contig boundaries, truncated at a scaffold end, or collapsed across repeated copies in a short-read draft assembly, standard assembly-only BLAST searches return missing calls (`LNF`), partial fragments, or silent chimeras.

Given a draft assembly, the paired-end reads used to build it, and a multi-FASTA catalogue of known alleles for the target locus, `Locus-Recon` recruits locus-associated read pairs (including unmapped and distally mapped mates), performs targeted local *de novo* assembly, projects full locus boundaries, and validates the reconstructed candidate by independent competitive read remapping.

---

## Table of Contents

1. [Overview and Command Suite](#overview-and-command-suite)
2. [How Locus-Recon Works](#how-locus-recon-works)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Preparing Inputs](#preparing-inputs)
6. [Running `locus-recon`](#running-locus-recon)
7. [Read the Results](#read-the-results)
8. [Copy-Number, Dosage, and Assembly-Graph Commands](#copy-number-dosage-and-assembly-graph-commands)
9. [Documentation Index](#documentation-index)
10. [Citation and License](#citation-and-license)

---

## Overview and Command Suite

`Locus-Recon` provides four command-line entry points that operate on a shared set of alignments and assembly artifacts:

| Command | Purpose | Primary Inputs | Primary Outputs |
|---|---|---|---|
| **`locus-recon`** | Batch targeted locus reconstruction, per-base read-support auditing, and disposition routing (`PASS` / `REVIEW` / `HOLD`). | `--samplesheet`, `--bait`, `--locus`, `--main-output-dir` | `locus_recon_report_<locus>.tsv`, disposition-partitioned FASTAs, per-sample `allele_remap.support.tsv`, `run_manifest.json` |
| **`locus-recon-depth-ratio`** | Quality-matched (`Q20/MQ20`) locus-to-backbone read-depth ratio and per-query-base dosage estimation with 500 bp block-bootstrap confidence intervals. | `<sample>/full_map.sorted.bam`, `--locus` intervals (optional `--gene-span`, `--gene-length`) | Console summary and optional one-row TSV (`--tsv`) |
| **`locus-recon-copy-number`** | Reconciles continuous read-depth dosage with bait-anchored flanking context counts in draft SPAdes assembly graphs (`GFA 1.x`). | `<sample>/full_map.sorted.bam`, `--locus`, optional draft `--gfa` + `--bait` | Provenance-rich JSON (`--output-json`) and optional TSV (`--output-tsv`) |
| **`locus-recon-graph-paths`** | Enumerates bounded source-to-sink paths from a local SPAdes GFA graph and optionally scores them by competitive read mapping (`MAPQ >= 20`). | `<sample>/spades_local/assembly_graph_after_simplification.gfa` (optional `--reads-r1`, `--reads-r2`) | Candidate paths FASTA (`--output`) and ranked/unranked summary TSV (`--summary`) |

---

## How Locus-Recon Works

### Seven-step targeted reconstruction workflow (`locus-recon`)

```text
Draft assembly (FASTA) + Curated bait alleles (FASTA) + Paired-end reads (FASTQ)
  │
  ├─► 1. Fragment-tolerant discovery (BLASTN: identity >= 80%, HSP length >= 50 bp; padded +/- 300 bp)
  ├─► 2. Full-assembly read mapping (BWA-MEM2 / BWA) & coordinate sorting (full_map.sorted.bam)
  ├─► 3. Target read recruitment & mate-pair rescue across all contigs (samtools view -N)
  ├─► 4. Targeted local de novo assembly (SPAdes default mode, automatic --sc fallback)
  ├─► 5. Distinct-region ranking, bait-projected boundary extraction & span-clip / continuation checks
  ├─► 6. Catalogue comparison (BLASTN) & independent competitive read remapping (samtools mpileup -A -aa)
  └─► 7. Four-axis orthogonal QC adjudication & disposition-partitioned FASTA routing
```

### Four-axis orthogonal adjudication

`Locus-Recon` separates four questions that single-score typing pipelines conflate:

1. **Workflow status (`workflow_status`):** Did the sample clear the minimum computational gates (`SUCCESS`, `SKIP`, or `FAIL`)?
2. **Sequence confidence (`sequence_confidence` / `qc_confidence`) and Disposition (`result_disposition`):** Are the reported bases supported by the reads, internally coherent, and complete across the projected locus span?
   - `HIGH` → **`PASS`** (written to `<locus>_accepted_alleles.fasta`)
   - `MEDIUM` or `LOW` → **`REVIEW`** (written to `<locus>_review_required_candidates.fasta`)
   - `SUSPECT` → **`HOLD`** (written to `<locus>_hold_candidates.fasta`)
3. **Catalogue status (`catalogue_status`):** How closely does the reconstructed sequence match the supplied reference catalogue (`EXACT_MATCH`, `NONEXACT_MATCH` at $\ge 97\%$ identity, `DIVERGENT_FROM_REFERENCE` at $85\text{--}97\%$, or `HIGHLY_DIVERGENT_FROM_REFERENCE` at $<85\%$)?
   - **Sequence confidence is independent of catalogue distance:** A genuinely novel allele with complete, unambiguous read support is classified as `HIGH` (`PASS`) while `catalogue_review_recommended=true` alerts the curator to its novelty.
4. **Projected span completeness (`span_clipped_bp`):** Because per-base pileup metrics evaluate only the bases present in a contig, a truncated allele can exhibit flawless read depth over its retained fragment. `Locus-Recon` projects the full bait span onto the local contig:
   - If a contig boundary clips $\ge 10\text{ bp}$ of the projected locus, `ALLELE_SPAN_CLIPPED_AT_CONTIG_END` blocks `HIGH` confidence (`PASS`).
   - If another local contig carries the missing sequence, `LOCUS_SPLIT_ACROSS_CONTIGS` downgrades the call to `SUSPECT` (`HOLD`) and reports the continuation contig, recovered bases (`span_continuation_bp`), and shared overlap (`span_continuation_overlap_bp`). `Locus-Recon` never silently concatenates separate contigs across an assembly gap.

---

## Installation

### System requirements

- **Operating System:** Linux (x86_64) or macOS; on Windows, use WSL2.
- **Conda / Mamba:** [Miniforge](https://github.com/conda-forge/miniforge) is recommended to install the pinned external bioinformatics binaries alongside the Python package.

### Install from source with Conda/Mamba

```bash
git clone https://github.com/crodsa/locus-recon.git
cd locus-recon

mamba env create -f environment.yml
mamba activate locus-recon
python -m pip install .

locus-recon --version
```

*(If `mamba` is not installed, substitute `conda env create -f environment.yml` and `conda activate locus-recon`.)*

### Pinned environment components (`environment.yml`)

| Component | Pinned Version | Role in `Locus-Recon` |
|---|---:|---|
| Python | `3.11.13` | Package runtime (compatible with Python $\ge 3.9$). |
| pip | `25.2` | Installs the local `locus-recon` Python package. |
| BLAST+ | `2.16.0` | `blastn` and `makeblastdb` for fragment discovery, scaffold ranking, and validation. |
| BWA-MEM2 | `2.2.1` | Primary short-read aligner (`bwa` is automatically used as a fallback if `bwa-mem2` is absent). |
| samtools | `1.22.1` | Coordinate/name sorting, mate rescue (`samtools view -N`, requires $\ge 1.12$), BAQ pileup, and depth extraction. |
| SPAdes | `4.2.0` | Local *de novo* assembly of recruited reads. |
| numpy | `1.26.4` | Optional acceleration for 500 bp block-bootstrap confidence intervals in `depth_ratio`. |
| biopython | `1.87` | Sequence parsing utilities for auxiliary validation workflows. |
| tqdm | `4.67.1` | Optional progress display during alignment streaming. |

---

## Quick Start

### 1. Pre-flight validation (`--dry-run`)

Verify your installation, external executables, and the bundled 13-sample synthetic dataset without running read alignment or assembly:

```bash
locus-recon \
  --samplesheet validation/mock_dataset/samples.tsv \
  --main-output-dir validation/dry_run_check \
  --bait validation/mock_dataset/mock_locus_alleles.fasta \
  --locus mockLocus \
  --dry-run
```

### 2. Run the bundled 13-sample end-to-end benchmark

Run the complete synthetic benchmark (intact, centrally and terminally fragmented, 3×–20× depth series, 5%–50% two-allele mixtures, close paralog, and negative control) and verify all six core evaluation criteria:

```bash
python validation/run_mock_validation.py \
  --force \
  --threads 2 \
  --memory-gb 4
```

This writes the full reconstruction output to `validation/mock_results/` and the machine-readable evaluation to `validation/mock_results/validation/validation_summary.json` (see [`validation/README.md`](validation/README.md)).

---

## Preparing Inputs

### 1. Samplesheet (`--samplesheet`)

Provide a 4-column tab-separated file (`sample_id`, `assembly_path`, `r1_path`, `r2_path`). A header row starting with `sample_id` is optional, lines starting with `#` are ignored, and relative file paths are resolved relative to the directory containing the samplesheet:

```tsv
sample_id	assembly_path	r1_path	r2_path
isolate_001	assemblies/isolate_001.fasta	reads/isolate_001_R1.fastq.gz	reads/isolate_001_R2.fastq.gz
isolate_002	assemblies/isolate_002.fasta	reads/isolate_002_R1.fastq.gz	reads/isolate_002_R2.fastq.gz
```

- **Sample IDs** must be unique, start with a letter or digit, and contain only `A-Z`, `a-z`, `0-9`, `.`, `_`, or `-`.
- **Reads (`r1_path`, `r2_path`)** must be the paired-end short reads corresponding to each isolate's draft assembly (plain `.fastq`/`.fq` or `.fastq.gz`/`.fq.gz`).

### 2. Curated Bait FASTA (`--bait`)

Provide a multi-FASTA file containing known reference alleles for **one target locus** (e.g., downloaded from PubMLST, BIGSdb, Pasteur, or a curated virulence/AMR database).

At startup, `Locus-Recon` profiles the bait file and enforces pre-flight quality checks:
- **Unique, BLAST-safe identifiers:** The first whitespace-delimited token after `>` must be unique and $\le 50$ characters (enforced before `makeblastdb -parse_seqids` runs).
- **Unaligned DNA sequences:** Sequences must contain only valid IUPAC nucleotide symbols (`ACGTRYSWKMBDHVN`). Alignment gap characters (`-` or `.`) are rejected with an explicit diagnostic message.
- **Off-target / fragment screening:** Any record sharing no 25-mer with the full-length alleles in the file triggers a pre-flight warning for potential paralogs or unrelated sequences.
- **Empirical length and GC baseline:** `Locus-Recon` computes the allele count (`qc_length_profile_n`), length range, median, mode, interquartile range (`IQR`), standard deviation, and mean GC content. Effective length tolerance for `HIGH` confidence is `IQR` when `IQR > 0`, or `max(stdev, 3 bp)` otherwise (if fewer than 5 alleles are supplied, length flags are accompanied by `CATALOGUE_LENGTH_PROFILE_UNDERPOWERED`).

---

## Running `locus-recon`

### Standard batch run

```bash
locus-recon \
  --samplesheet study/samplesheet.tsv \
  --main-output-dir results/aroE \
  --bait refs/aroE_alleles.fasta \
  --locus aroE \
  --threads 8 \
  --parallel-samples 4 \
  --memory-per-sample 16
```

### Resource planning and intermediate cleanup

- **CPU allocation:** Maximum logical CPUs $\approx$ `--threads` $\times$ `--parallel-samples` (in the example above, $8 \times 4 = 32$ threads).
- **Memory allocation:** Maximum SPAdes RAM $\approx$ `--memory-per-sample` $\times$ `--parallel-samples` (in the example above, $16 \times 4 = 64\text{ GB}$).
- **Non-coding loci (`--noncoding-locus`):** Pass `--noncoding-locus` when targeting rRNA genes (e.g., 23S rRNA), intergenic regions, or promoters so six-frame stop-codon and length-modulo-3 checks are disabled.
- **Intermediate cleanup (`--cleanup`):** Removes bulky per-sample BAMs, temporary FASTQs, and `spades_local/` directories after a sample reaches `SUCCESS`, while preserving all FASTAs, BLAST tables, BED files, logs, QC reports, and `allele_remap.support.tsv`.
  > **Important:** Do **not** pass `--cleanup` if you plan to run `locus-recon-depth-ratio`, `locus-recon-copy-number`, or `locus-recon-graph-paths` on those samples, as those companion commands read `<sample>/full_map.sorted.bam`, `<sample>/spades_local/assembly_graph_after_simplification.gfa`, and `<sample>/target_R1.fq.gz` / `target_R2.fq.gz`.

### Key `locus-recon` options

| Option | Default | Description |
|---|---:|---|
| `--samplesheet TSV` | *required* | 4-column TSV (`sample_id`, `assembly_path`, `r1_path`, `r2_path`). |
| `--main-output-dir DIR` | *required* | Output directory for batch reports, catalogues, and sample folders. |
| `-b`, `--bait FASTA` | *required* | Multi-FASTA of known reference alleles for the target locus. |
| `-l`, `--locus NAME` | *required* | Target locus identifier (`[A-Za-z0-9][A-Za-z0-9._-]*`). |
| `-t`, `--threads N` | `8` | Threads per sample. |
| `-j`, `--parallel-samples N` | `1` | Number of samples processed concurrently. |
| `--memory-per-sample GB` | `16` | Maximum memory (GB) passed to each SPAdes process. |
| `--min-discovery-identity PCT` | `80.0` | Minimum BLASTN percent identity for Step 1 draft discovery HSPs. |
| `--min-hsp-length BP` | `50` | Minimum aligned length (bp) for Step 1 draft discovery HSPs. |
| `--region-padding BP` | `300` | Flanking bases added around discovery HSPs before read recruitment. |
| `--min-identity PCT` | `85.0` | Minimum BLASTN identity to nearest bait allele for `SUCCESS`. |
| `--min-coverage PCT` | `80.0` | Minimum query coverage of the reconstructed candidate for `SUCCESS`. |
| `--min-bitscore-margin SCORE` | `10.0` | Minimum bitscore separation between best and second-best distinct local regions. |
| `--min-remap-depth X` | `5.0` | Minimum mean filtered remap depth across the candidate for `SUCCESS`. |
| `--min-remap-breadth PCT` | `90.0` | Minimum percentage of candidate bases covered at $\ge 1\times$ for `SUCCESS`. |
| `--min-base-quality PHRED` | `20` | Minimum base quality (`-q`) in `samtools mpileup`. |
| `--min-mapping-quality PHRED` | `20` | Minimum mapping quality (`-Q`) in `samtools mpileup`. |
| `--mixture-min-fraction FRAC` | `0.05` | Minor-allele fraction threshold $(0, 0.5)$ for candidate mixed sites. |
| `--mixture-min-sites N` | `2` | Minimum credible bidirectional mixed sites required to call `ALLELE_MIXTURE`. |
| `--mixture-min-alt-depth N` | `2` | Minimum quality-filtered alternative reads required at a candidate mixed site. |
| `--noncoding-locus` | `off` | Disable coding-frame stop-codon and modulo-3 length checks. |
| `--cleanup` | `off` | Delete bulky BAM/FASTQ/SPAdes intermediates after sample `SUCCESS`. |
| `--dry-run` | `off` | Validate inputs, bait profile, dependencies, and manifest without processing reads. |

See the [Command-Line Reference](docs/cli-reference.md) for full details and exit codes.

---

## Read the Results

### Output directory structure

```text
results/aroE/
├── run_manifest.json                             # Full provenance, parameters, tool versions & SHA-256 hashes
├── locus_recon_batch_aroE.log                    # Batch execution log
├── locus_recon_report_aroE.tsv                   # Primary 55-column batch summary table
├── aroE_accepted_alleles.fasta                   # PASS candidates only (SUCCESS + HIGH confidence)
├── aroE_review_required_candidates.fasta         # REVIEW candidates (SUCCESS + MEDIUM or LOW confidence)
├── aroE_hold_candidates.fasta                    # HOLD candidates (SUCCESS + SUSPECT confidence)
├── aroE_reconstructed_candidates.fasta           # All SUCCESS reconstructions (audit collection)
├── bait_database/                                # Shared BLAST database built from --bait
└── isolate_001/
    ├── isolate_001_aroE_reconstructed.fasta      # Candidate sequence oriented to best bait match
    ├── isolate_001_aroE_qc_report.txt            # Human-readable QC scorecard and flag breakdown
    ├── isolate_001_aroE_recon.log                # Sample-level execution log
    ├── allele_remap.support.tsv                  # Per-base Q20/MQ20 pileup & mixture evidence table
    ├── blast_step1.tsv                           # Draft assembly discovery HSPs
    ├── target_regions.bed                        # Padded draft intervals used for read recruitment
    ├── blast_step5.tsv                           # Local SPAdes scaffold hits
    ├── allele_region.bed                         # Projected candidate coordinates on local scaffold
    └── blast_validation.tsv                      # Candidate vs. reference bait BLASTN results
```

### How to interpret the batch report (`locus_recon_report_<locus>.tsv`)

Always start with these columns in `locus_recon_report_<locus>.tsv`:

| Column | Values | Action / Meaning |
|---|---|---|
| `workflow_status` (`status`) | `SUCCESS`, `SKIP`, `FAIL` | Whether a candidate passed minimum identity, coverage, depth, and breadth gates. `SKIP` indicates an expected data condition (e.g., target absent or unbridgeable gap); `FAIL` indicates a file or execution error. |
| `sequence_confidence` (`qc_confidence`) | `HIGH`, `MEDIUM`, `LOW`, `SUSPECT` | Confidence in the reconstructed nucleotide sequence based on read support, frame integrity, and span completeness. |
| `result_disposition` | `PASS`, `REVIEW`, `HOLD` | Operational routing: accept `PASS` (`HIGH`), inspect `REVIEW` (`MEDIUM`/`LOW`), and withhold `HOLD` (`SUSPECT`). |
| `catalogue_status` | `EXACT_MATCH`, `NONEXACT_MATCH`, `DIVERGENT_FROM_REFERENCE`, `HIGHLY_DIVERGENT_FROM_REFERENCE` | Relationship to the reference bait catalogue (`nearest_allele`, `nearest_allele_identity_pct`, `exact_known_allele`). |
| `catalogue_review_recommended` | `true`, `false` | `true` whenever `exact_known_allele=false` so novel allele candidates are routed for nomenclature curation. |
| `span_clipped_bp` | integer $\ge 0$ | Bases of the projected bait span lost to a local contig end ($\ge 10\text{ bp}$ raises `ALLELE_SPAN_CLIPPED_AT_CONTIG_END` and blocks `HIGH`). |
| `span_continuation_contig` | contig ID or empty | Populated when another local contig carries the clipped sequence (`LOCUS_SPLIT_ACROSS_CONTIGS`, downgraded to `SUSPECT`). |
| `mixture_detected` / `mixed_site_count` | `true`/`false`, integer | Whether $\ge$ `--mixture-min-sites` positions carry bidirectional minor-allele support without Fisher strand bias ($p \ge 0.01$). |
| `qc_flags` | semicolon-delimited list | Exact reasons for any confidence downgrade (see [`docs/output-interpretation.md`](docs/output-interpretation.md#qc-flag-reference)). |

---

## Copy-Number, Dosage, and Assembly-Graph Commands

When a target locus is duplicated or amplified in tandem, identical copies produce no mixed polymorphic sites. `Locus-Recon` provides three specialized commands to quantify copy number and explore alternative graph topologies.

### 1. Normalised depth ratio and dosage (`locus-recon-depth-ratio`)

Computes the quality-filtered (`Q20/MQ20`) median read depth across one or more target discovery intervals divided by the median depth of the single-copy chromosomal backbone ($\ge 5\text{ kb}$ contigs, $100\text{ bp}$ end-trimmed, excluding target and `--exclude` intervals) using `<sample>/full_map.sorted.bam`:

```bash
# Single-interval locus
locus-recon-depth-ratio \
  --bam results/aphA1/isolate_001/full_map.sorted.bam \
  --locus NODE_7:250422-251237 \
  --exclude NODE_12 \
  --tsv results/aphA1/isolate_001_depth.tsv

# Multi-fragment locus with BLAST query coordinates (--gene-span and --gene-length)
locus-recon-depth-ratio \
  --bam results/target/isolate_002/full_map.sorted.bam \
  --locus NODE_3:100-950 \
  --locus NODE_8:1-620 \
  --gene-span 1-851 \
  --gene-span 852-1471 \
  --gene-length 1471 \
  --tsv results/target/isolate_002_depth.tsv
```

- **`depth_call`:** `MULTICOPY_DEPTH` when the 95% block-bootstrap lower bound exceeds `1.5` (`ratio_ci_low > 1.5`); otherwise `SINGLE_COPY_COMPATIBLE`.
- **`ambiguity_index`:** Ratio of filtered (`Q20/MQ20`) to unfiltered (`Q0/MQ0`) locus depth. Values $< 0.70$ indicate multi-mapping across uncollapsed assembly copies (`LOW_MAPPING_UNIQUENESS`) and mark the dosage as a `LOWER_BOUND`.
- **`dosage_estimate` and `dosage_status`:** Combines assembly fragment coverage with per-region depth ratios (`ESTIMATED`, `ESTIMATED_PARTIAL_COVERAGE`, `REVIEW_OVERLAPPING_SPANS` with `geometry_expected_dosage`, `LOWER_BOUND`, or `NOT_ESTIMATED`).

### 2. Joint graph/depth copy-number consensus (`locus-recon-copy-number`)

Combines `locus-recon-depth-ratio` with bounded, bait-anchored traversal of the draft SPAdes assembly graph (`assembly_graph_after_simplification.gfa`) to reconcile discrete flanking-context counts with continuous read-depth dosage:

```bash
locus-recon-copy-number \
  --bam results/23S/isolate_001/full_map.sorted.bam \
  --locus NODE_19:901-3789 \
  --gfa draft_assemblies/isolate_001/assembly_graph_after_simplification.gfa \
  --bait refs/23S_reference.fasta \
  --bait-id 23S_rRNA \
  --output-json results/23S/isolate_001.copy-number.json \
  --output-tsv results/23S/isolate_001.copy-number.tsv
```

Always interpret `copy_number_call` together with `copy_number_kind`, `copy_number_method`, and `copy_number_status`:

| `copy_number_kind` | `copy_number_method` | Meaning |
|---|---|---|
| `INTEGER_CONTEXT_COUNT` | `GRAPH_DEPTH_CONSENSUS` | Both 5' and 3' bait ends have matching flanking context counts ($L = R \ge 1$) with no cycle/limit truncation, concordant with read depth. |
| `INTEGER_CONTEXT_COUNT` | `GRAPH_COUNT_OVER_DISCORDANT_DEPTH` | Exact multi-copy graph context count ($L = R \ge 2$) when depth also calls `MULTICOPY_DEPTH` but the bootstrap depth interval does not contain the integer count (`DEPTH_GRAPH_NUMERIC_DISCORDANCE`). |
| `MEAN_DEPTH_DOSAGE` | `DEPTH_TANDEM_COMPATIBLE` | Graph shows 1 flanking context ($L = R = 1$) while depth calls `MULTICOPY_DEPTH` (characteristic of a collapsed tandem repeat array); preserves continuous `dosage_estimate`. |
| `MEAN_DEPTH_DOSAGE` | `DEPTH_ONLY` | No `--gfa` was supplied, or graph topology was unresolved while depth yielded a reliable estimate. |
| `LOWER_BOUND` | `GRAPH_DEPTH_LOWER_BOUND` / `DEPTH_LOWER_BOUND` | Graph traversal hit a cycle (`CYCLE_ENCOUNTERED`), path/node cap, or asymmetric ends, or `ambiguity_index < 0.70`; reports the conservative floor in `copy_number_lower_bound`. |

### 3. Local graph path enumeration and competitive read scoring (`locus-recon-graph-paths`)

Exports bounded terminal (source-to-sink) paths from a local SPAdes GFA graph (`<sample>/spades_local/assembly_graph_after_simplification.gfa`) and optionally ranks them by competitive read placement (`MAPQ >= 20`):

```bash
locus-recon-graph-paths \
  --gfa results/tcdB/isolate_001/spades_local/assembly_graph_after_simplification.gfa \
  --output results/tcdB/isolate_001.graph_paths.fasta \
  --summary results/tcdB/isolate_001.graph_paths.tsv \
  --min-length 6000 \
  --max-length 8000 \
  --reads-r1 results/tcdB/isolate_001/target_R1.fq.gz \
  --reads-r2 results/tcdB/isolate_001/target_R2.fq.gz \
  --threads 4
```

When `--reads-r1` is provided, all enumerated paths are indexed together so reads compete for alignment. The summary TSV reports `rank`, `candidate_id`, `length`, `mean_graph_depth`, `mapped_reads`, `unique_mean_depth`, `unique_breadth_pct`, `unsupported_bp`, and `nodes`. If every read maps equally well to multiple paths (`unique_breadth_pct == 0`), the command states explicitly that the candidate paths cannot be phased at the library's read length.

---

## Documentation Index

- **[Command-Line Reference](docs/cli-reference.md):** Complete option tables, argument validation rules, resource planning, `--cleanup` behavior, `run_manifest.json` schema, and exit codes for all four CLI commands.
- **[Output and Interpretation Reference](docs/output-interpretation.md):** Detailed guide to all output files, the 55 columns of `locus_recon_report_<locus>.tsv`, the 20 columns of `allele_remap.support.tsv`, QC confidence thresholds, and every QC flag.
- **[Troubleshooting Guide](docs/troubleshooting.md):** Step-by-step diagnostics for installation issues, `SKIP` / `FAIL` states, and `SUSPECT` / `REVIEW` quality flags.
- **[Methods and Validation Guide](docs/validation.md):** Algorithmic details, copy-number and graph reconciliation rules, study-level validation panel design, and pre-submission curation checklists.
- **[Synthetic Benchmark Guide](validation/README.md):** Instructions and reference results for the bundled 13-sample end-to-end verification suite.

---

## Citation and License

If you use `Locus-Recon` in published work, please cite the software release ([`CITATION.cff`](CITATION.cff)) and the version/release date of the reference allele catalogue used as `--bait`.

- **Software license:** MIT License ([`LICENSE`](LICENSE)).
- **Bundled synthetic benchmark and documentation data:** Creative Commons Attribution 4.0 International (`CC BY 4.0`, [`LICENSE-DATA.md`](LICENSE-DATA.md)).
