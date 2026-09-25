# Troubleshooting guide

This guide connects Locus-Recon 1.0.0 messages and unexpected results to the
evidence files that can explain them. It does not treat every negative result as
a software error or every reconstruction as biologically correct.

For output definitions, see the
[output and interpretation reference](output-interpretation.md). For option
defaults and ranges, see the [command-line reference](cli-reference.md).

## Diagnose in this order

1. Read the terminal's final batch summary.
2. Open `locus_recon_report_<locus>.tsv` and find the sample's `status` and
   `message`.
3. Open `<sample>/<sample>_<locus>_recon.log` if the per-sample directory exists.
4. Inspect the evidence file named in the relevant section below.
5. Verify sample identity, assembly/read pairing, upstream contamination QC,
   bait scheme/release, and the exact command in `run_manifest.json` before
   changing thresholds.

A missing file can be expected when the sample stopped before that workflow
stage. `--cleanup` removes bulky files only after `SUCCESS`; principal tables,
BED intervals, logs, and per-base support remain, but the inputs of the
companion commands do not.

## Installation and startup

### `locus-recon: command not found`

**Likely meaning**

- the `locus-recon` Conda environment is not active;
- `python -m pip install .` was not run inside the source checkout; or
- the terminal has not loaded the Conda shell initialization.

**Check**

```bash
conda env list
conda activate locus-recon
python -m pip show locus-recon
which locus-recon
locus-recon --version
```

On Windows under WSL, run these commands in the Linux shell, not PowerShell.
If the package is absent, return to the repository directory and run
`python -m pip install .` while the environment is active.

### `Missing required tools` or `No supported aligner found`

**Likely meaning**

The active environment does not expose one or more of `blastn`, `makeblastdb`,
`samtools`, `spades.py`/`spades`, and either `bwa-mem2` or `bwa` on `PATH`.

**Check**

```bash
which blastn makeblastdb samtools bwa-mem2 spades.py
samtools --version
spades.py --version
```

**Corrective action**

Create or reactivate the supplied environment:

```bash
mamba env create -f environment.yml
mamba activate locus-recon
python -m pip install .
```

If the environment already exists but has drifted, compare its package list
with `environment.yml` or recreate it under a new name. Do not satisfy the check
with arbitrary executables of unknown provenance.

### `samtools ... >= 1.12 is required`

Mate rescue uses `samtools view -N`, so older samtools versions are rejected.
Activate or recreate the supplied environment. Record the resolved executable
path and version from the startup log or `run_manifest.json`; another samtools
earlier on `PATH` may be shadowing the intended environment.

### Bait FASTA validation error

Implemented bait errors include an empty file, duplicate sequence ID, empty
record, or a non-IUPAC DNA symbol.

**Inspect**

- the bait file named in the command;
- the record identifier immediately following each `>`; and
- whitespace, amino-acid symbols, punctuation, and gap characters in sequences.

Locus-Recon uses only the first whitespace-delimited word of each FASTA header
as its identifier. Two descriptive headers can therefore still collide. Fix the
source FASTA rather than removing valid biological diversity.

### Samplesheet parsing error

Common messages report the line number, the wrong number of columns, an unsafe
or duplicate sample ID, or an empty path.

**Inspect**

- open the file in a text editor that can display tabs;
- confirm exactly four tab-separated fields per data row;
- confirm a header, if present, begins with `sample_id`;
- confirm all IDs are unique and contain only permitted characters; and
- remember that relative paths are resolved from the samplesheet directory.

Saving a spreadsheet as CSV creates comma-separated data and is not valid.
Export as tab-separated text/TSV.

### `Bait file not found`, `Samplesheet not found`, or `Missing input files`

Check `pwd`, the command stored in `run_manifest.json` when available, and the
resolved paths in the error. Assembly/read paths in the samplesheet are relative
to the samplesheet, while the `--bait`, `--samplesheet`, and
`--main-output-dir` arguments are resolved from the shell invocation.

Do not replace a missing read file with reads from another isolate merely to
complete a row.

### `One or both FASTQ files are empty or unreadable`

Locus-Recon checks that both R1 and R2 contain at least one readable text line;
`.gz` paths are opened as gzip streams.

**Check**

```bash
gzip -t reads/isolate_001_R1.fastq.gz
gzip -t reads/isolate_001_R2.fastq.gz
```

Also confirm file permissions, transfer completeness, and that the files are
actually FASTQ data. This startup check is not a full FASTQ-format or read-pair
consistency validator; run the study's normal read QC as well.

### Dry run creates files

This is expected. `--dry-run` creates the main output directory, initializes
the batch log, builds the shared bait BLAST database, and writes
`run_manifest.json`. It does not map or assemble sample reads. Use a separate
dry-run output directory if those setup artifacts should not share a real-run
directory.

## Fragment discovery

### `No contigs matched bait`

**Meaning**

BLAST produced no draft HSP that met both `--min-discovery-identity` and
`--min-hsp-length`. The message reports raw and passing hit counts.

**Inspect**

- `blast_step1.tsv` for raw HSP identity, aligned length, query position, draft
  contig, and score;
- the bait scheme, locus, sequence orientation, and database release;
- assembly FASTA validity and sample identity; and
- upstream assembly completeness, contamination, and coverage.

**Possible actions**

- correct an incorrect bait or assembly;
- add validated divergent alleles to the curated bait set;
- cautiously lower discovery identity or HSP length only after benchmarking
  off-target and paralog behavior; or
- use another evidence source if the assembly is too incomplete.

**Do not conclude**

This message alone does not distinguish true biological absence from sequence
divergence, an assembly gap that removed nearly the entire locus, or unsuitable
bait data.

### Many contigs or target regions pass discovery

This can reflect a fragmented locus, homologous/paralogous regions, repetitive
sequence, contamination, or permissive thresholds.

Inspect `blast_step1.tsv` and `target_regions.bed`. Compare contig identity,
aligned length, coverage, annotation, and upstream contamination evidence.
Reducing `--region-padding` changes read recruitment but does not remove the
underlying multiple-hit evidence. Raising discovery thresholds can hide true
divergent fragments, so validate any change with positive and negative controls.

## Read mapping and recruitment

### Very low mapping in `mapping_flagstat.txt`

The R1/R2 files may not match the supplied assembly, the assembly may represent
a different processing state, or read/assembly quality may be poor.

Compare sample identifiers and checksums from the upstream workflow. Inspect
`bwa_pipe.log` and the complete flagstat. Do not continue biological
interpretation until an assembly/read swap has been excluded.

### `No reads mapped to target regions`

**Meaning**

The full assembly mapping completed, but no alignment overlapped the padded BED
intervals in `target_regions.bed`.

**Inspect**

- `target_regions.bed` and `blast_step1.tsv`;
- `mapping_flagstat.txt` and `full_map.sorted.bam` when retained;
- whether R1/R2 are the reads used for the assembly; and
- contig identifiers and interval bounds.

If overall mapping is healthy but target intervals lack reads, the bait HSP can
be an unsupported assembly artifact or the locus can be at extremely low depth.
Increasing `--region-padding` only helps when informative reads lie just outside
the current intervals; it cannot repair a sample mismatch.

### `No reads were extracted after mate rescue`

Target-overlapping read names were found, but no usable paired or singleton
FASTQ records remained after all matching alignments were selected, name-sorted,
and converted.

Inspect `target_region.bam`, `target_read_names.txt`,
`target_region_rescued.bam`, `target_namesorted.bam`, and the three extracted
FASTQs. Check read-name conventions, secondary/supplementary alignments,
truncation, and BAM-to-FASTQ messages in the per-sample log.

Do not substitute independently processed reads without confirming they still
match the assembly and sample.

## Local assembly

### SPAdes reports an invalid k-mer coverage histogram

Locus-Recon handles one specific SPAdes failure automatically. When the SPAdes
log contains `Invalid kmer coverage histogram`, the workflow removes the failed
working directory and retries with `--sc`. The batch report records
`single-cell-retry` if that retry succeeds.

If both standard and `--sc` modes encounter the same condition, the sample is
`SKIP` with a low/non-uniform-coverage message. Inspect recruited read counts,
FASTQs, read QC, and locus depth. This expected skip should not be relabelled as
a software failure.

### `SPAdes produced no scaffolds`

**Inspect**

- `<sample>_<locus>_recon.log`;
- `spades_local/spades.log` and other SPAdes diagnostics;
- `target_R1.fq.gz`, `target_R2.fq.gz`, and `target_S.fq.gz`; and
- the number of recruited read names in the batch report.

Very low depth, short or poor-quality reads, a non-uniform selected subset,
contamination, or insufficient bridging pairs can prevent assembly. Increasing
padding may recruit more flanking pairs, but it also increases off-target data.
Benchmark the change rather than assuming more reads are always better.

### SPAdes command failure other than the handled histogram condition

The sample returns `FAIL`. Inspect the exact command and SPAdes log for memory,
storage, permission, version, or malformed-read errors. Confirm available disk
space and that `--memory-per-sample × --parallel-samples` fits the node. Do not
classify an execution failure as locus absence.

### `No BLAST hit in local assembly`

SPAdes produced scaffolds, but none aligned to any bait allele under the BLAST
search itself.

Inspect `spades_local/scaffolds.fasta`, `blast_step5.tsv`, recruited reads, and
discovery evidence. The selected reads may not bridge the target, the local
assembly may be off-target, or the locus may be too divergent/repetitive for the
current design. Lowering final gates cannot help when no local BLAST hit exists.

## Candidate validation

### `Allele failed validation`

The message reports final identity and query coverage and their configured
minimums. A per-sample candidate FASTA and `blast_validation.tsv` can already
exist even though the status is `SKIP`.

**Inspect**

- all validation hits in `blast_validation.tsv`;
- candidate length and projected bounds in `allele_region.bed`;
- `blast_step5.tsv` and the selected local scaffold; and
- bait correctness and completeness.

Low query coverage can indicate a truncated scaffold, incorrect boundary
projection, a structural difference, or a non-target homolog. Low identity can
indicate divergence, a paralog, contamination, or an incorrect local assembly.
Do not promote the existing candidate FASTA as a successful allele merely
because it was written.

### `Final allele failed remap support`

The candidate passed similarity/coverage gates but failed the configured mean
depth or breadth requirement. The message reports both values.

Inspect `allele_remap.support.tsv`, remap BAM and flagstat when retained,
extracted FASTQs, and the candidate ends. Zero-depth termini can reflect
projected boundaries unsupported by reads; internal gaps can indicate a repeat,
indel, misassembly, or mixture.

Lowering the remap gate changes which candidates receive `SUCCESS`; it does not
create missing evidence. Calibrate a change with independently known truth and
report it explicitly.

## Successful results with QC flags

### `SUCCESS` but `MEDIUM`, `LOW`, or `SUSPECT`

This is expected behavior, not a contradiction. `SUCCESS` applies only the
minimum identity, query-coverage, mean-depth, and breadth gates. QC confidence
then considers stricter similarity/support thresholds plus length, GC,
ambiguous bases, frame consistency, patchiness, qualities, strands, mixture,
reference discordance, and distinct-region ambiguity.

Read `qc_flags`, the human-readable QC report, and the exact flag definitions in
the [output reference](output-interpretation.md#qc-flag-reference). Only
`PASS` records enter `accepted_alleles.fasta`; inspect the explicit review and
hold collections rather than treating every reconstructed candidate as an
accepted allele.

### Low remap breadth or patchy support

Open `allele_remap.support.tsv` and inspect or plot:

- `position` versus `depth`;
- `alt_fraction`;
- `forward_depth` and `reverse_depth`;
- base and mapping qualities; and
- credible, strand-biased, and reference-discordant markers.

```text
Healthy support:  25 28 31 30 29 27 26 24 25 28  ... across the locus
Patchy support:   22 21  0  0  1 18 20  0  0 19  ... investigate
```

Retain intermediates when investigating a new study so BAM alignments and the
SPAdes local assembly are available.

### `AMBIGUOUS_BEST_HIT`

The best and second-best **distinct local regions** had a bitscore separation
below `--min-bitscore-margin`. Multiple known bait alleles hitting the same
region are collapsed first and do not by themselves create this flag.

Inspect separate regions in `blast_step5.tsv` and `spades_local/scaffolds.fasta`.
Consider a paralog, duplicated/repetitive sequence, contamination, or competing
local assemblies. The program still reports the top-ranked candidate, but the
flag is severe and the confidence is `SUSPECT`.

### `ALLELE_MIXTURE`

At least `--mixture-min-sites` positions met the configured alternative depth
and fraction rules, had alternative support on both strands, and were not
significantly strand-biased by the implemented Fisher test.

Inspect rows with `credible_mixed_site=true` in
`allele_remap.support.tsv`. Compare their fractions, qualities, read positions,
and consistency across the locus. Run whole-genome contamination/mixed-isolate
QC and check sample provenance. Locus-Recon does not determine the component
haplotypes or prove that the signal reflects two viable isolates.

### `STRAND_BIASED_ALT_SUPPORT`

One or more candidate alternative positions were supported on only one strand
or had Fisher strand-bias `p < 0.01`. Inspect rows with
`strand_biased_site=true`, especially read ends, low-complexity sequence,
indels, mapping quality, and strand counts. This can reflect an artifact, but it
is treated as severe until resolved.

### `REFERENCE_DISCORDANCE`

At one or more positions, the quality-filtered observations were nearly fixed
for a base different from the reconstructed candidate. Inspect
`reference_discordant_site=true` rows, the candidate FASTA, remap BAM, and local
scaffold. Possible explanations include a consensus/local-assembly mismatch,
alignment issue, indel representation, or mixed sequence. Do not edit the FASTA
manually without reconstructable evidence and provenance.

### `INTERNAL_STOPS` or `FRAME_LENGTH_SHIFT`

Confirm that the target is a coding locus and that bait records use consistent
biological boundaries and orientation. Remember that typing definitions can be
internal gene fragments, so the workflow does not require start or terminal
stop codons. Inspect the sequence in the bait-supported frame and consider
indels, ambiguous bases, a paralog, or assembly error.

Use `--noncoding-locus` only when coding-frame checks are biologically
inappropriate, not merely to suppress an inconvenient flag.

### `POSSIBLE_NOVEL_ALLELE`

This flag means the candidate has an isolated length deviation, identity of at
least 97%, and no severe flag. It is a review hypothesis. Confirm the bait
database is current, verify exact scheme boundaries, inspect every changed or
inserted/deleted position, exclude paralogy/mixture, and seek independent
confirmation before claiming or submitting novelty.

The flag is descriptive. An isolated length deviation with no severe flag is
assigned `MEDIUM` whether or not identity reaches 97%, so the label records a
hypothesis about the catalogue rather than setting the tier.

### `CATALOGUE_DIVERGENCE_BELOW_HIGH`, `CATALOGUE_DIVERGENT`, `CATALOGUE_HIGHLY_DIVERGENT`

The reconstruction is further from the nearest catalogue allele than 97%, 93% or
85% identity respectively. These flags are descriptive: they do not reduce
sequence confidence, because distance from a curated catalogue describes the
catalogue, not the quality of the evidence supporting the reported bases. A
sequence carrying one of them can still be `HIGH`/`PASS`. Treat it as a
candidate new or divergent allele: confirm the bait database is the intended
scheme and boundaries, check `nearest_allele`, and apply the typing scheme's own
rules for accepting a new allele. If the sequence is instead expected to match a
known allele, the divergence usually points to the wrong bait database, the
wrong locus boundaries, a paralogue, or a mixed sample — check
`AMBIGUOUS_BEST_HIT`, `ALLELE_MIXTURE` and `REFERENCE_DISCORDANCE` first.

## Batch and rerun behavior

### A mixed batch exits `0` although some samples failed

Exit `0` means that at least one real sample reached `SUCCESS` or a dry run
passed. It does not mean every sample succeeded. Use the TSV `status` column and
count `SUCCESS`, `SKIP`, and `FAIL` rows explicitly.

### The command exits `1` although useful diagnostic files exist

Exit `1` is expected when no real sample reaches `SUCCESS`, a dry run has bad
sample inputs, or a batch-wide input/dependency/shared-database step fails.
Retained logs and intermediate evidence can still explain the failure.

### A rerun already has an output directory

Locus-Recon reuses the directory and overwrites or rebuilds its named log,
manifest, report, catalog, databases, per-sample log, and intermediates. The
SPAdes working directory is explicitly removed before local assembly.

Use a new main output directory for every bait release or parameterization that
must remain traceable. Compare manifests rather than assuming an old directory
still represents one coherent run.

### A companion command cannot find its input

`locus-recon-depth-ratio` and `locus-recon-copy-number` read
`<sample>/full_map.sorted.bam`; `locus-recon-graph-paths` reads the local graph
in `<sample>/spades_local/` and the recruited reads `target_R1.fq.gz` and
`target_R2.fq.gz`. `--cleanup` deletes all of them after a successful sample.
Rerun that sample without `--cleanup` into a new output directory.
`locus-recon-copy-number` also needs the assembly graph of the
draft itself, which comes from the SPAdes run that built the draft, not from
Locus-Recon.

### Disk usage remains high after `--cleanup`

Cleanup applies only to samples that reached `SUCCESS`; skipped and failed
samples retain available intermediates for diagnosis. Batch-level bait
databases and retained evidence tables/logs also remain. Review exact contents
against the [cleanup list](cli-reference.md#cleanup-and-reruns) before removing
anything manually, and preserve the run's provenance and evidence requirements.

## When to stop tuning and seek independent evidence

Do not keep relaxing thresholds solely until a candidate becomes `SUCCESS`.
Use an independent long-read/hybrid assembly, targeted sequencing, or another
validated method when:

- alternative distinct regions remain unresolved;
- the candidate changes a scientific or epidemiological conclusion;
- a novel allele is proposed for database submission;
- projected ends lack read support;
- mixed or reference-discordant positions remain unexplained; or
- repeated parameter changes produce incompatible candidates.

The [methods and validation guide](validation.md) provides a full validation
panel and manual-review checklist.
