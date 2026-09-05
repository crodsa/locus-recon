# Mock truth-set validation

This directory contains a deterministic, distributable benchmark for the complete Locus-Recon workflow. It is designed to test scientific behavior rather than only Python functions.

## Benchmark components

The generated dataset contains:

- a complete synthetic chromosome and perfect locus truth;
- locus-spanning mock long reads representing independent long-read/hybrid evidence;
- an intact short-read baseline;
- central and near-terminal locus fragmentation;
- 3×, 5×, 10×, and 20× short-read depth samples;
- controlled 5%, 10%, 20%, and 50% two-allele mixtures;
- a genome containing a closely related paralog in a distinct genomic region;
- a negative genome with no target locus.

```text
complete genome + spanning long-read truth
                 |
                 +-- intact draft + 40x short reads
                 +-- central split + 40x short reads
                 +-- near-terminal split + 40x short reads
                 +-- intact draft + 3x / 5x / 10x / 20x reads
                 +-- allele A:B mixtures at 95:5 / 90:10 / 80:20 / 50:50
                 +-- target plus close paralog
                 `-- target-negative genome
```

The generator uses a fixed seed, records every setting and changed allele position in `dataset_manifest.json`, and writes gzip files with deterministic timestamps.

## Generate only

```bash
python validation/generate_mock_dataset.py
```

The committed `mock_dataset/` can be recreated exactly:

```bash
python validation/generate_mock_dataset.py --force
```

## Run the complete benchmark

Activate the environment and install the package first:

```bash
mamba env create -f environment.yml
mamba activate locus-recon
python -m pip install -e '.[test]'

python validation/run_mock_validation.py --force
```

Useful runner options are `--threads`, `--memory-gb`, `--dataset-dir`, and
`--results-dir`. Existing generated results are refused unless `--force` is
given, which helps prevent accidental mixing of benchmark runs.

The runner performs reconstruction and then creates:

```text
validation/mock_results/validation/
├── validation_results.tsv
└── validation_summary.json
```

Core pass criteria are:

1. exact baseline reconstruction;
2. exact reconstruction for both fragmented controls;
3. no mixture false positives among successful pure samples;
4. detection of all mixtures at 10% or above;
5. ambiguity flagging in the paralog control;
6. no successful reconstruction in the negative control.

The 5% mixture and low-depth series are reported as empirical detection/sensitivity boundaries but are not required core passes. This prevents a small synthetic benchmark from overstating universal limits.

## Bundled reference result

The compact files in `reference_results/` were regenerated with the final
Locus-Recon 1.0 release, BWA-MEM2 2.2.1, samtools 1.22.1, SPAdes 4.2.0, and
BLAST+ 2.16.0+. The execution record, including the environment and the outcome
of each criterion, is in
[`BENCHMARK_REPRODUCTION.md`](reference_results/BENCHMARK_REPRODUCTION.md).

| Test | Observed result |
|---|---|
| Core criteria | 6/6 passed |
| Intact 40x control | Exact 507 bp allele |
| Central and terminal fragmentation | Exact 507 bp allele in both |
| Depth series | 3x and 5x skipped; exact allele at 10x and 20x |
| Controlled mixtures | Detected at 5%, 10%, 20%, and 50% |
| Pure successful controls | No mixture false positives |
| Close paralog | Ambiguity flagged |
| Negative control | Not reconstructed |

View [the JSON summary](reference_results/validation_summary.json),
[per-sample evaluation](reference_results/validation_results.tsv), or
[full batch report](reference_results/locus_recon_report_mockLocus.tsv). The
reference output is intentionally small; bulky BAM and SPAdes intermediates are
generated locally and excluded from the distribution.

## Real-data GFA rerun

The recovered `ERR467623` paired reads made it possible to regenerate a
contemporary LIBA-6656 local assembly graph with the final Locus-Recon 1.0
toolchain. The complete compact deposit is under
[`liba6656-gfa-rerun/`](liba6656-gfa-rerun/README.md). It contains the GFA,
the exact ungapped bait FASTA, all 512 exported terminal paths, portable
parameters and checksums, and the alignment of the two deposited reconstructed
candidates against those paths.

The rerun restores an executable graph/path record, but it is not described as
the exact historical GFA. The 512 paths are combinatorial proposals through
nine variable regions; they are not 512 accepted alleles and still require the
competitive-read, per-base, polishing and context evidence deposited with the
publication package.

## Characteristic failure modes

Two error modes are visible in the bundled record, and both are worth reading
before a status field is trusted on its own.

The paralog control returns a sequence that is not exactly right. At 83.4x mean
remap depth and 99.8% breadth the run is reported as `SUCCESS`, yet the returned
allele is 99.01% identical to the designated truth, which over 507 bp is five
mismatches. The target and its paralog differ at fourteen positions and at no
other position, so those five errors sit on discriminating sites: the consensus
took the paralog base. The evidence carries the warning that the status does
not, since the run reports fourteen candidate sites, eight of them mixed at a
median minor fraction of 0.458, six with strand-biased support, the ambiguity
flag raised and a QC verdict of `SUSPECT`. A locus with a near-identical
paralog is read from the evidence table, never from `status` alone.

The low-evidence modes fail in the safer direction. At nominal 3x and 5x no
allele is returned, and the negative control yields no candidate contig, so
scarcity of evidence produces silence rather than a fabricated allele.

## Truth model

The complete chromosome FASTAs are the exact sequence truth. Noisy mock long
reads independently span the target and represent the evidence from which a
complete long-read or hybrid assembly would be obtained. Locus-Recon does not
consume those reads; its short-read reconstruction is compared directly with
the complete-genome locus they support.

Each sample's expected category, depth, mixture fraction, truth allele, and
read-pair count is stored in `truth/expected_results.tsv`. The generator seed,
sequence changes, break coordinates, insert model, error rate, and SHA-256 hash
of every dataset file are recorded in `dataset_manifest.json`.

This benchmark is a regression and method-development control. It does not replace validation with real isolates, real platform-specific error profiles, contamination controls, and independently generated hybrid assemblies.

## Reproducible real-data workflows

Publication-facing analyses are now executable from accession-level inputs
instead of being represented only by deposited result files. Shared helpers in
`validation/common.py` fail on missing inputs or non-zero commands and record
SHA-256 checksums, the exact Git commit, parameters, and tool versions.

The seven-run *aphA1* analysis is retained in full:

```bash
python validation/workflows/run_aphA1_validation.py \
  --reference work/aphA1/MRSN56.fasta \
  --reads-dir work/aphA1/reads \
  --output-dir work/aphA1/results
```

Its frozen manifest is `manifests/aphA1_copy_number.tsv`. The 816 bp gene
interval supplies the aggregate depth and dosage estimate; the historical,
separate 4 kb interval supplies the calibration-eligible 20-window profile.
The workflow requires all seven runs and evaluates the two published qPCR
intervals without using them to tune a threshold.

The external low-copy validation is prespecified in
[`low-copy-23S/PRESPECIFICATION.md`](low-copy-23S/PRESPECIFICATION.md) and can be
run with:

```bash
python -m pip install -e '.[validation]'
python validation/workflows/run_low_copy_23S.py \
  --work-dir work/low-copy-23S --threads 4 --memory-gb 8
```

The runner downloads and verifies the five frozen paired Illumina datasets and
closed truth records, generates SPAdes drafts from Illumina only, maps the reads
back to those drafts, discovers 23S rRNA and `gyrB` from a fixed external
`NC_000915.1` bait, evaluates the released depth rule unchanged, and counts bait-
anchored contexts in the corresponding SPAdes GFA. The archived rerun reports
two contexts for every 23S locus and one for every `gyrB` control while the
workflow verifies that all 23 historical depth fields remain unchanged. Closed
sample-specific chromosomes are read only by the truth-audit stage. Circular
plasmid replicons are recorded separately and do not count as chromosomes or
contribute to chromosome feature truth. Compact completed outputs are archived
under `validation/low-copy-23S/`. These five cases informed the graph/depth
design and are therefore a development/biological demonstration, not an
independent exact-copy performance estimate.
