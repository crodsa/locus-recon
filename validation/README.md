# Synthetic Installation and Regression Benchmark

This directory contains a deterministic, self-contained 13-sample benchmark for verifying a `Locus-Recon` installation end to end. It tests the complete reconstruction, read-recruitment, local assembly, remapping, mixture-detection, paralog-ambiguity, and quality-control pipeline without requiring external network access.

## Benchmark design

The benchmark dataset (`mock_dataset/`) is generated deterministically from a fixed pseudo-random seed (`SEED = 731921`) and models a 507 bp coding target locus (`mockLocus`) across 13 controlled scenarios:

- **Ground truth:** Complete synthetic chromosome FASTAs (`complete_genome_A.fasta`, `complete_genome_B.fasta`, `complete_genome_paralog.fasta`, `complete_genome_negative.fasta`) and locus-spanning mock long reads (`mock_long_reads_A.fastq.gz`).
- **Intact baseline (`baseline_intact_A`):** Intact draft chromosome with 40× paired-end short reads (150 bp reads, 350 bp insert size, 0.001 substitution error rate).
- **Fragmented drafts (`fragmented_mid_A`, `fragmented_terminal_A`):** Target locus split across two contigs at 50% (central break) or 12% (near-terminal break) of the locus length, paired with 40× short reads.
- **Coverage dilution series (`depth_03x_A`, `depth_05x_A`, `depth_10x_A`, `depth_20x_A`):** Central-split draft assembly paired with 3×, 5×, 10×, and 20× short reads.
- **Two-allele mixtures (`mixture_05pct_B`, `mixture_10pct_B`, `mixture_20pct_B`, `mixture_50pct_B`):** Central-split draft assembly paired with 60× total reads containing 5%, 10%, 20%, and 50% minor-allele (`mockLocus_B`) reads.
- **Close paralog (`paralog_control`):** Draft assembly and 60× reads carrying both `mockLocus_A` and a paralogous copy differing at 14 positions in a separate genomic region.
- **Target-negative control (`negative_control`):** Draft assembly and 40× reads lacking the target locus.

```text
complete genome + spanning long-read truth
                 |
                 +-- intact draft + 40x short reads
                 +-- central split + 40x short reads
                 +-- near-terminal split + 40x short reads
                 +-- fragmented draft + 3x / 5x / 10x / 20x reads
                 +-- allele A:B mixtures at 95:5 / 90:10 / 80:20 / 50:50
                 +-- target plus close paralog
                 `-- target-negative genome
```

Every setting, mutated codon position, break coordinate, and SHA-256 file checksum is recorded in `mock_dataset/dataset_manifest.json`.

## Running the benchmark

Activate the `locus-recon` Conda environment and install the package:

```bash
mamba env create -f environment.yml
mamba activate locus-recon
python -m pip install .
```

### 1. Fast pre-flight check (`--dry-run`)

Validate the benchmark inputs, bait database, external tool dependencies, and run manifest without running read alignment or local assembly:

```bash
locus-recon \
  --samplesheet validation/mock_dataset/samples.tsv \
  --main-output-dir validation/dry_run_check \
  --bait validation/mock_dataset/mock_locus_alleles.fasta \
  --locus mockLocus \
  --dry-run
```

### 2. Full end-to-end benchmark

Run all 13 samples through `locus-recon` and evaluate the outputs against the ground-truth manifest:

```bash
python validation/run_mock_validation.py \
  --force \
  --threads 2 \
  --memory-gb 4
```

Options for `run_mock_validation.py`:
- `--dataset-dir PATH`: Benchmark dataset directory (default: `validation/mock_dataset`).
- `--results-dir PATH`: Output directory for the run (default: `validation/mock_results`).
- `--threads INT`: Threads per sample (default: `2`).
- `--memory-gb INT`: SPAdes memory cap in GB (default: `4`).
- `--regenerate`: Rebuild `mock_dataset/` from scratch before running.
- `--force`: Overwrite an existing `--results-dir`.

### 3. Regenerating the synthetic dataset only

To verify deterministic dataset generation without running the pipeline:

```bash
python validation/generate_mock_dataset.py --force
```

## Evaluation criteria and reference results

`run_mock_validation.py` writes its evaluation to:

```text
validation/mock_results/validation/
├── validation_results.tsv
├── validation_summary.json
└── workflow_provenance.json
```

Six core criteria must pass (`all_core_criteria_passed: true`):

1. `baseline_exact`: Exact sequence reconstruction for the intact 40× control.
2. `fragmented_exact`: Exact sequence reconstruction for both central and near-terminal fragmented controls.
3. `no_false_mixture_in_successful_pure_controls`: Zero false-positive mixture calls across all pure baseline, fragmented, and depth-series controls.
4. `mixtures_10pct_and_above_detected`: Explicit `ALLELE_MIXTURE` detection in the 10%, 20%, and 50% mixed-allele controls.
5. `paralog_ambiguity_flagged`: `AMBIGUOUS_BEST_HIT` flagged in the paralog control.
6. `negative_not_reconstructed`: No false reconstruction (`SKIP`) in the target-negative control.

Compact reference outputs from Linux x86_64 (`Locus-Recon 1.0.0`, `BWA-MEM2 2.2.1`, `samtools 1.22.1`, `SPAdes 4.2.0`, `BLAST+ 2.16.0`) are provided in [`reference_results/`](reference_results/BENCHMARK_REPRODUCTION.md):

| Scenario | Observed workflow status | Sequence accuracy | QC confidence (`result_disposition`) | Key evidence |
|---|---|---|---|---|
| Intact 40× control | `SUCCESS` | Exact (507 bp, 100%) | `HIGH` (`PASS`) | Clean per-base support |
| Central & terminal splits (40×) | `SUCCESS` | Exact (507 bp, 100%) | `HIGH` (`PASS`) | Assembled across break |
| Depth series (20× / 10× / 5× / 3×) | `SUCCESS` | Exact (507 bp, 100%) | `MEDIUM` / `LOW` / `SUSPECT` / `SUSPECT` | Confidence scales with depth |
| Mixtures (5% / 10% / 20% / 50%) | `SUCCESS` | 100% to major allele | `SUSPECT` (`HOLD`) | `ALLELE_MIXTURE` (2, 22, 25, 25 sites) |
| Close paralog control | `SUCCESS` | 99.01% (5 mismatches) | `SUSPECT` (`HOLD`) | `AMBIGUOUS_BEST_HIT` + `ALLELE_MIXTURE` (14 sites) |
| Target-negative control | `SKIP` | Not reconstructed | — | Stopped at Step 1 (`No contigs matched bait`) |

See [`reference_results/validation_summary.json`](reference_results/validation_summary.json), [`reference_results/validation_results.tsv`](reference_results/validation_results.tsv), and [`reference_results/locus_recon_report_mockLocus.tsv`](reference_results/locus_recon_report_mockLocus.tsv) for the complete reference tables.
