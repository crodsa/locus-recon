# Command-line reference

This page documents installation, all public command-line options, resource
planning, cleanup, provenance, and exit behavior for Locus-Recon 1.0.0. For a
first run, begin with the [README quick start](../README.md#quick-start).

## Installation and dependencies

### Recommended source installation

Locus-Recon is installed from its GitHub source checkout. The provided Conda
environment supplies the external bioinformatics programs and pins the runtime
Python version.

```bash
git clone https://github.com/crodsa/locus-recon.git
cd locus-recon

mamba env create -f environment.yml
mamba activate locus-recon
python -m pip install .

locus-recon --version
```

Use `conda env create -f environment.yml` if Conda is available but Mamba is
not. [Miniforge](https://github.com/conda-forge/miniforge) provides both.

The environment currently pins:

| Component | Environment constraint | Runtime behavior |
|---|---:|---|
| Python | `3.11.13` | Package metadata supports Python ≥3.9. CI tests Python 3.9, 3.11, 3.12, and 3.13 on Linux. |
| pip | `25.2` | Installs the local Python package. |
| tqdm | `4.67.1` | Provides optional byte-stream progress display. |
| BLAST+ | `2.16.0` | `blastn` and `makeblastdb` are required. |
| BWA-MEM2 | `2.2.1` | Preferred short-read aligner. The code accepts `bwa` as a fallback if `bwa-mem2` is absent. |
| samtools | `1.22.1` | The code enforces samtools ≥1.12 because mate rescue uses `samtools view -N`. |
| SPAdes | `4.2.0` | Local assembly runs in SPAdes' default mode; `--careful` is not used. |

At startup, Locus-Recon resolves executable paths, checks the samtools version,
and records detected paths and versions in `run_manifest.json`. A missing
required dependency terminates the batch before sample processing.

The project verifies Python/package behavior on Linux in CI. The bundled
end-to-end reference run was generated on Linux x86_64. Native Windows
execution is not validated by the project; use a Linux environment such as WSL
when working from Windows.

## Command synopsis

```text
locus-recon \
  --samplesheet TSV \
  --main-output-dir DIR \
  --bait FASTA \
  --locus NAME \
  [runtime options] \
  [discovery and validation thresholds] \
  [per-base evidence options] \
  [--noncoding-locus]
```

Run `locus-recon --help` to print the installed version's own synopsis.

## Required arguments

| Option | Value | Meaning and validation |
|---|---|---|
| `--samplesheet` | TSV path | Four tab-separated columns: `sample_id`, `assembly_path`, `r1_path`, `r2_path`. A header is optional. Relative input paths are resolved from the samplesheet directory. |
| `--main-output-dir`, `--main_output_dir` | directory | Root directory for batch and per-sample outputs. It is created if absent. Existing named outputs can be overwritten on rerun. |
| `-b`, `--bait` | FASTA path | Multi-FASTA of known alleles for the one requested locus. Sequence IDs must be unique; records must be non-empty and contain IUPAC DNA symbols. |
| `-l`, `--locus` | name | Name used in output paths and reports. It must start with a letter or number and contain only letters, numbers, dots, underscores, or hyphens. |

### Samplesheet rules

```tsv
sample_id	assembly_path	r1_path	r2_path
isolate_001	data/isolate_001.fasta	reads/isolate_001_R1.fastq.gz	reads/isolate_001_R2.fastq.gz
```

- Exactly four columns are required.
- Blank lines and lines whose first field begins with `#` are ignored.
- The first data-like line is treated as a header when its first field is
  `sample_id`, case-insensitively.
- Sample IDs follow the same safe-character rule as the locus name and must be
  unique.
- Input paths cannot be empty. `~` is expanded; relative paths are resolved
  from the samplesheet directory.
- Assemblies must be non-empty. R1 and R2 may be plain or gzip-compressed, but
  both must contain readable records.

## Runtime options

| Option | Default | Accepted value | Behavior |
|---|---:|---|---|
| `-t`, `--threads` | `8` | integer ≥1 | Threads available to each sample. Mapping internally allocates approximately three quarters to BWA and one quarter to samtools sorting; the value is also passed to SPAdes and other samtools operations. |
| `-j`, `--parallel-samples`, `--parallel_samples` | `1` | integer ≥1 | Maximum samples processed concurrently. Worker count is capped by the number of valid samples. |
| `--memory-per-sample` | `16` | integer ≥1 GB | Maximum memory passed to each SPAdes process with `-m`. It is not a batch-wide memory cap. |
| `-v`, `--verbose` | off | flag | Prints debug messages to standard error in addition to writing debug logs. |
| `--no-progress`, `--no_progress` | off | flag | Disables the byte-stream progress display during mapping. Useful for non-interactive logs. |
| `--cleanup` | off | flag | Removes bulky intermediates only after a sample reaches `SUCCESS`. Final sequence and principal summaries remain. |
| `--dry-run` | off | flag | Validates bait, samplesheet, file availability, executable dependencies, shared bait database, and manifest without mapping or assembling sample reads. |

## Discovery and final validation options

| Option | Default | Accepted value | What it controls |
|---|---:|---|---|
| `--min-discovery-identity` | `80.0` | 0–100% | Minimum identity for a bait-to-draft HSP to seed a target region. |
| `--min-hsp-length` | `50` | integer ≥1 bp | Minimum aligned length for a discovery HSP. Discovery intentionally does not gate on HSP query coverage, so partial fragments can seed reconstruction. |
| `--region-padding`, `--region_padding` | `300` | integer ≥0 bp | Sequence added on both sides of each passing draft HSP before read selection, clipped at contig boundaries. Overlapping or adjacent intervals are merged. |
| `--min-identity`, `--min_identity` | `85.0` | 0–100% | Minimum candidate identity to the nearest bait allele during final validation. |
| `--min-coverage`, `--min_coverage` | `80.0` | 0–100% | Minimum final BLAST query coverage of the reconstructed candidate. |
| `--min-bitscore-margin`, `--min_bitscore_margin` | `10.0` | number ≥0 | Minimum score separation between the best and second-best distinct local regions. A smaller margin creates an ambiguity flag but does not itself prevent `SUCCESS`. |
| `--min-remap-depth`, `--min_remap_depth` | `5.0` | number ≥0 × | Minimum mean quality-filtered depth after selected reads are remapped to the candidate. |
| `--min-remap-breadth`, `--min_remap_breadth` | `90.0` | 0–100% | Minimum percentage of candidate positions with at least one quality-filtered remap observation. |

The final identity, coverage, mean-depth, and breadth values are minimum
`SUCCESS` gates. The fixed QC classifier uses stricter boundaries for `HIGH`,
`MEDIUM`, `LOW`, and `SUSPECT`; see the
[confidence reference](output-interpretation.md#qc-confidence-logic).

### Changing discovery thresholds

Lowering `--min-discovery-identity` or `--min-hsp-length` can recover more
divergent or shorter fragments, but it can also recruit unrelated homologs or
paralogs. Increasing `--region-padding` can retrieve more bridging mates, but it
also expands the recruited neighborhood. Calibrate these choices with positive,
fragmented, paralog, mixture, and negative controls rather than changing them
only to force a desired sample to reconstruct.

## Per-base evidence and mixture options

| Option | Default | Accepted value | What it controls |
|---|---:|---|---|
| `--min-base-quality` | `20` | integer ≥0 Phred | Minimum base quality supplied to `samtools mpileup`; lower-quality bases do not contribute to depth or alternative evidence. Base qualities are BAQ-recalibrated first. |
| `--min-mapping-quality` | `20` | integer ≥0 Phred | Minimum alignment mapping quality supplied to mpileup. |

The pileup is run with `-A`, so a read whose pair the aligner could not mark as
properly paired, which is common against a candidate of a few hundred bases,
still counts on its own strand. Without it, those reads would be dropped
selectively and could manufacture strand imbalance at alternative sites.
| `--mixture-min-fraction` | `0.05` | number >0 and <0.5 | Lower alternative-fraction boundary for a candidate mixed site. The upper boundary is `1 - value`; values above that are reference-discordant. |
| `--mixture-min-sites` | `2` | integer ≥1 | Number of credible bidirectional mixed positions required for `mixture_detected=true`. |
| `--mixture-min-alt-depth` | `2` | integer ≥1 | Minimum quality-filtered alternative observations required at a candidate mixed site. |

A candidate mixed position must meet fraction and alternative-depth rules.
Alternative observations must be present on both strands and must not have a
two-sided Fisher strand-bias `p < 0.01` to become a credible mixed position.
That p-value boundary is fixed and not configurable. The locus-level mixture call is not a
whole-genome contamination screen and does not phase component alleles.

## Locus options

| Option | Default | Behavior |
|---|---:|---|
| `--noncoding-locus`, `--noncoding_locus` | off | Disables coding-frame integrity checks. Use only for a genuinely non-coding target or a definition for which coding-frame checks are inappropriate. |

When coding checks are active, the bait set is evaluated in all six reading
frames, and the frame or tied frames with the fewest internal stops are
reported as plausible. For a candidate, `internal_stops` is the minimum over
all six frames and `qc_coding_frame_used` names the frame that achieved it.
Candidate length modulo three is compared with the bait length mode. Start and
terminal stop codons are not required.

## Utility options

| Option | Behavior |
|---|---|
| `-h`, `--help` | Print help and exit. |
| `--version` | Print the installed Locus-Recon version and exit. |

Where a table lists an underscore spelling beside the hyphenated one, both are
accepted; the hyphenated form is used throughout this documentation.

## Resource planning

### CPU

Approximate maximum logical CPUs requested by sample processing are:

```text
threads per sample × concurrent samples
```

For example:

```bash
locus-recon \
  --samplesheet study/samples.tsv \
  --main-output-dir results/target_locus \
  --bait study/target_locus_alleles.fasta \
  --locus target_locus \
  --threads 8 \
  --parallel-samples 4 \
  --memory-per-sample 12
```

This permits approximately 32 logical CPUs. Sample-level parallelism often
helps more than assigning very high thread counts to one short local assembly.

### Memory

Approximate maximum SPAdes allowance is:

```text
memory per sample × concurrent samples
```

The example permits up to 48 GB across four concurrent SPAdes jobs. Reserve
additional memory for the operating system, process overhead, mapping, and BAM
sorting. `--memory-per-sample` is forwarded to each SPAdes invocation; it does
not reserve or enforce the total node memory.

### Storage and I/O

All input reads are mapped to their draft assembly, so read access and sorted
BAM writing commonly dominate I/O. Use a suitable high-throughput filesystem
for a large batch. `--cleanup` reduces retained storage only after successful
samples; failed and skipped samples retain available intermediates for diagnosis.

Independent loci can run as independent jobs when total CPU, memory, and I/O
allow. Use a separate main output directory per locus and parameterization.

## Dry-run behavior

`--dry-run` performs the following work:

1. creates the main output directory and initializes the batch log;
2. profiles and validates the bait FASTA;
3. parses the samplesheet and validates assembly/R1/R2 file presence and basic
   readability;
4. resolves required executables and checks the samtools version;
5. builds the shared bait BLAST database;
6. writes `run_manifest.json`; and
7. exits without mapping reads or creating normal per-sample reconstruction
   directories.

The dry run exits nonzero if any sample input is missing/empty/unreadable or a
batch-wide validation/dependency/database step fails.

## Local assembly behavior

SPAdes receives the rescued paired reads, recovered singletons, the per-sample
thread count, a Phred offset of 33, and the per-sample memory allowance, and
runs in its default mode; `--careful` is not used.

If standard SPAdes fails and its log contains `Invalid kmer coverage histogram`,
Locus-Recon removes the failed SPAdes directory and retries with `--sc`. A
successful retry is recorded as `single-cell-retry`. If the retry
encounters the same condition, the sample is `SKIP`. Other SPAdes command
failures remain execution failures.

## Cleanup and reruns

### What `--cleanup` removes

Cleanup runs only after a sample reaches `SUCCESS`. It removes:

- the copied input assembly and its BWA indexes;
- per-sample draft/local BLAST database index files;
- full-map, target-region, rescued, name-sorted, and final-remap BAMs and indexes;
- the target read-name list;
- extracted R1, R2, and singleton FASTQs;
- final-allele BWA/FASTA index sidecars;
- final-remap flagstat and aligner logs; and
- the complete `spades_local/` directory.

It preserves, among other summaries:

- the final reconstructed FASTA;
- per-sample log and QC report;
- `blast_step1.tsv`, `blast_step5.tsv`, and `blast_validation.tsv`;
- `target_regions.bed` and `allele_region.bed`;
- `allele_remap.support.tsv`;
- `mapping_flagstat.txt` and the full-map aligner log; and
- all batch-level report, catalog, log, database, and manifest files.

Omit `--cleanup` when validating a new study, investigating a problematic
sample, or retaining BAM/SPAdes evidence for audit. It also removes the inputs
of the companion commands: `full_map.sorted.bam`, which `locus-recon-depth-ratio`
and `locus-recon-copy-number` read, and the local graph in `spades_local/` and
the recruited reads `target_R1.fq.gz` and `target_R2.fq.gz`, which
`locus-recon-graph-paths` reads. Leave it off for any sample you may analyse
with those commands.

### Rerunning into an existing directory

The batch log, report, catalog, manifest, bait database, per-sample log, and
named intermediates are rebuilt or overwritten. A pre-existing per-sample
SPAdes directory is deleted before local assembly. Old unrelated files are not
systematically versioned or removed.

Use a new main output directory whenever parameterizations or bait releases
must remain independently auditable. Do not rely on an existing directory to
preserve the previous run.

## Run manifest

`run_manifest.json` uses schema version `1.0` and records:

- program and package version;
- start date/time and complete `sys.argv` command list;
- Python version and platform string;
- all parsed arguments, including resolved output/bait/samplesheet paths and
  derived runtime fields;
- bait path, size, modification time, SHA-256 checksum, and profile;
- samplesheet path, size, modification time, and SHA-256 checksum;
- each sample's resolved assembly/R1/R2 path, size, and modification time; and
- resolved executable paths and their detected version strings.

The manifest checksums the bait and samplesheet, not every assembly and FASTQ.
If full input-content integrity is required, create and retain a separate
checksum manifest for all study inputs.

## Exit codes

| Exit code | Meaning |
|---:|---|
| `0` | Help/version completed, a dry run passed, or a real batch produced at least one `SUCCESS` sample. A mixed batch can still contain `SKIP` or `FAIL` rows. |
| `1` | A required file/dependency/shared step failed, a dry run found invalid sample files, or a real batch produced no `SUCCESS` sample. |
| `2` | Argument parsing or structured-input validation failed, including malformed samplesheets, invalid numeric ranges, unsafe locus names, or invalid bait content. |

Always inspect the batch TSV rather than using the process exit code to infer
that every sample succeeded.

## Bundled end-to-end benchmark

After installing the complete Conda environment, run the deterministic
13-sample synthetic workflow with:

```bash
python validation/run_mock_validation.py \
  --force \
  --threads 2 \
  --memory-gb 4
```

The runner regenerates the mock dataset, performs reconstruction, and writes a
machine-readable evaluation under `validation/mock_results/validation/`.
Existing results are refused unless `--force` is supplied. Additional runner
options include `--dataset-dir` and `--results-dir`.

See the [benchmark manual](../validation/README.md), bundled
[reference summary](../validation/reference_results/validation_summary.json),
and [methods guide](validation.md). The synthetic benchmark is an installation
and regression control, not a substitute for validation on real isolates and
platform-specific libraries.

## Development installation

This section is for contributors, not routine users.

```bash
mamba env create -f environment.yml
mamba activate locus-recon
python -m pip install -e '.[test]'

python -m pytest -q
python -m build
```

The build command requires the `build` package if it is not already installed.
CI runs the Python tests on Python 3.9, 3.11, 3.12, and 3.13, runs `ruff check`,
and separately verifies that the source distribution and wheel can be built.

## Auxiliary commands

Three commands complement the main `locus-recon` batch workflow. They act on
files the workflow or the draft assembly already produced; none realigns reads
or re-runs discovery. A sample run with `--cleanup` no longer has their inputs.

### `locus-recon-graph-paths`

Exports bounded source-to-sink paths retained in a local SPAdes graph, for loci
where the FASTA scaffold collapsed or truncated a repeat. Research mode: it
assigns no truth label, and every exported path still requires competitive read
validation.

```bash
locus-recon-graph-paths --gfa SAMPLE/spades_local/assembly_graph_after_simplification.gfa --output SAMPLE.graph_paths.fasta --summary SAMPLE.graph_paths.tsv --min-length 6000 --max-length 8000
```

| option | meaning |
|---|---|
| `--gfa` | `assembly_graph_after_simplification.gfa` from the local assembly, in `<sample>/spades_local/` |
| `--output` | FASTA destination for the retained paths |
| `--summary` | path metadata TSV (default `OUTPUT.paths.tsv`) |
| `--min-length`, `--max-length` | discard paths outside this length range (`--max-length 0` disables) |
| `--max-paths`, `--max-nodes` | safety bounds on enumeration |
| `--prefix` | identifier prefix for exported candidates |
| `--reads-r1`, `--reads-r2` | recruited reads for competitive scoring (`<sample>/target_R1.fq.gz`, `<sample>/target_R2.fq.gz`); all retained paths are indexed together so the reads compete for placement |
| `--threads` | threads for competitive scoring (default 4) |
| `--min-mapping-quality` | depth is counted only from alignments at or above this quality, excluding reads that fit several paths equally well (default 20) |
| `--score-dir` | destination for the scoring BAM and text outputs (default: `graph_path_scores/` beside `--output`) |

Without `--reads-r1` the summary columns are `candidate_id`, `length`,
`mean_graph_depth`, `nodes`. With reads they become `rank`, `candidate_id`,
`length`, `mean_graph_depth`, `mapped_reads`, `unique_mean_depth`,
`unique_breadth_pct`, `unsupported_bp`, `nodes`, ordered by uniquely anchored
breadth. A rank orders the enumerated alternatives; it is not an allele call,
and paths sharing most of their sequence score alike by construction. If no
aligner is available the unranked summary is kept and the reason is printed.

### `locus-recon-depth-ratio`

Normalised locus depth ratio, for loci whose copies are identical and therefore
invisible to allelic evidence. Takes `<sample>/full_map.sorted.bam`, the reads
mapped to the draft assembly during read recruitment.

```bash
locus-recon-depth-ratio --bam SAMPLE/full_map.sorted.bam --locus contig_7:250422-251237 --exclude plasmid_1 --exclude plasmid_2 --tsv SAMPLE.locus_depth.tsv
```

| option | meaning |
|---|---|
| `--bam` | reads mapped to the draft assembly |
| `--locus` | `contig:start-end`; repeat once per discovery interval |
| `--gene-span` | query coordinates of the bait covered by the matching `--locus` interval; repeat in the same order |
| `--gene-length` | length of one copy of the gene, in bp; required with `--gene-span` |
| `--exclude` | region kept out of the single-copy backbone; repeatable |
| `--reference` | draft assembly FASTA; reports the locus GC content as a descriptive field |
| `--min-bq`, `--min-mq` | quality filters, applied identically to locus and backbone (default 20/20) |
| `--call-threshold` | lower confidence bound above which single copy is rejected (default 1.5) |
| `--seed`, `--n-boot` | block-bootstrap seed and resamples (default 20260817, 1000) |
| `--max-backbone-bootstrap-blocks` | deterministic maximum number of evenly spaced 500 bp backbone blocks used for CI resampling (default 512); the point denominator still uses all eligible bases |
| `--tsv` | write the full result as a one-row TSV |

`--gene-span` and `--gene-length` are what turn the ratio into
`dosage_estimate` for a locus split across several assembly regions;
`copies_estimate` carries the same value under a shorter name. Incomplete query
coverage qualifies the dosage rather than withholding it
(`ESTIMATED_PARTIAL_COVERAGE`, with `query_coverage_fraction`); overlapping
query spans set `REVIEW_OVERLAPPING_SPANS` and report overlap metrics.
Thresholds and their null distributions are documented in
[`validation/depth-copy-number/PRESPECIFIED_CRITERIA.md`](../validation/depth-copy-number/PRESPECIFIED_CRITERIA.md).

### `locus-recon-copy-number`

Add optional SPAdes graph-context evidence to the unchanged depth result. JSON
is the required machine-readable output; TSV is optional. The graph is the one
the SPAdes run that produced the draft assembly wrote, whose contigs the BAM and
`--locus` coordinates refer to, not the local graph in `<sample>/spades_local/`.

```bash
locus-recon-copy-number \
  --bam SAMPLE/full_map.sorted.bam \
  --locus NODE_19:901-3789 \
  --gfa draft_assembly/SAMPLE/assembly_graph_after_simplification.gfa \
  --bait loci.fasta --bait-id 23S \
  --output-json SAMPLE.23S.copy-number.json \
  --output-tsv SAMPLE.23S.copy-number.tsv
```

All depth options (`--locus`, `--gene-span`, `--gene-length`, filters,
bootstrap, exclusions and the call threshold) retain the meanings documented
above. Additional options are:

| option | meaning |
|---|---|
| `--gfa` | SPAdes GFA 1.x of the draft assembly, with embedded segment sequences. Must be supplied together with `--bait`. |
| `--bait` | FASTA containing the target sequence used to locate the locus in GFA segments. |
| `--bait-id` | Exact first-token FASTA identifier; required when `--bait` has more than one record. |
| `--blastn` | BLASTN executable (default `blastn`). |
| `--graph-min-identity` | Minimum retained bait-to-segment HSP identity (default 80%). |
| `--graph-min-hit-bp` | Minimum HSP length; default is the greater of 200 bp or 10% of bait length. |
| `--graph-min-query-coverage` | Minimum union coverage of the bait before an integer graph count is eligible (default 0.90). |
| `--graph-min-context-bp` | Minimum post-branch context sequence; default is `max(100, 2 × largest GFA overlap)`. |
| `--graph-max-context-nodes` | Maximum traversed non-target nodes per path (default 8). |
| `--graph-max-context-paths` | Maximum explored edges per bait side (default 64). |
| `--graph-max-context-bp` | Maximum distance from a bait end attributed to locus context (default 5,000 bp). |
| `--output-json` | Required output containing arrays, nulls, depth, graph, consensus and provenance fields. |
| `--output-tsv` | Optional one-row TSV with semicolon-delimited list fields. |

`copy_number_call` is not self-interpreting. Read it with
`copy_number_kind`, `copy_number_method`, `copy_number_status`, bounds and flags.
An `INTEGER_CONTEXT_COUNT` is distinct from a `MEAN_DEPTH_DOSAGE`. A matched
single graph context plus `MULTICOPY_DEPTH` returns the continuous dosage as
`DEPTH_TANDEM_COMPATIBLE`; this preserves tandem arrays such as *aphA1* rather
than incorrectly replacing their dosage with one. An integer count outside the
depth interval is returned under `GRAPH_COUNT_OVER_DISCORDANT_DEPTH` with
`DEPTH_GRAPH_NUMERIC_DISCORDANCE`. When `ambiguity_index` is below 0.70 and
depth still rejects one copy, the dosage is only a floor: the result is a
`LOWER_BOUND` flagged `DEPTH_LOWER_BOUND`, carrying the larger of the depth and
graph floors. Supplying invalid graph inputs is an error; valid but unresolved
topology is retained as an explicit uncertainty state.
