# Comparison with ARIBA and KMA

This folder holds a comparison of Locus-Recon 1.0.0 with two established tools that also return a locus
sequence from short reads: ARIBA 2.14.7 and KMA 1.7.1. The three tools received the same reads and the same
reference sequences, and the sequences they returned were scored against the same truth.

The script is [`../run_comparator_benchmark.py`](../run_comparator_benchmark.py). It runs ARIBA and KMA. The rows
of Locus-Recon are read from the tables that are already in this repository; Locus-Recon is not run again.

## What was compared

| Data set | Cases | Reads | Reference sequences given to each tool | Truth |
|---|---|---|---|---|
| Simulated benchmark | 13 samples, one locus | `../mock_dataset/reads/` | Alleles A, B and C of `mock_locus_alleles.fasta` | The simulated allele; both alleles in a mixture |
| *Helicobacter pylori* | 5 isolates × 5 loci = 25 | Public Illumina runs of `../tier-calibration/closed_genome_manifest.tsv` | One bait per locus from strain 26695 | The annotated gene in the closed genome of the same isolate |
| *Clostridioides difficile* LIBA-6656 | *tcdB* | ERR467623 | The 111 *tcdB* alleles of `../liba6656-gfa-rerun/` | None. The two analyst-built candidates serve as the comparison |

The two candidates of LIBA-6656 (`liba6656_analyst_built_candidates.fasta`) are hypotheses, not a truth. Agreement
with them shows that two methods reach the same sequence from the same reads. It does not confirm the sequence.

## Settings

- ARIBA: `ariba prepareref` and `ariba run` with default settings. A second configuration adds
  `--nucmer_min_id 80`, because the default of 90% is close to the identity of *cagA* to its bait.
- KMA: `kma index` and `kma -ipe` with default settings.
- Coding loci were declared as genes to ARIBA and the 23S rRNA gene as non-coding. `ariba prepareref` rejected
  the three simulated alleles as genes, so they were declared non-coding
  (`catalogue_sequences_declared_noncoding_for_ariba` in `PROVENANCE.json`).
- 16 threads on one Linux server. The Conda environment is in `conda_environment.yml`. ARIBA 2.14.7 needs
  `setuptools` below 81.
- The simulated benchmark of Locus-Recon ran with two settings that differ from its defaults
  (`--min-remap-depth 3`, `--min-bitscore-margin 80`); all other Locus-Recon results use the defaults.

## How a returned sequence is scored

Each returned sequence is compared with every truth sequence of its locus, on both strands.

| Relation | Meaning |
|---|---|
| identical | The same sequence |
| uncalled bases, otherwise identical | Same length; differs only at N |
| incomplete, identical over its own length | An exact part of the truth |
| truth contained, extra bases | The whole truth plus bases at an end |
| identical over the overlap, different boundaries | The aligned part is identical; at least one sequence continues beyond it |
| differs | At least one aligned base or gap column differs |

Only `identical` counts as exact. `identity_pct`, `aligned_bp`, `mismatches` and `gap_columns` describe the longest
local alignment (BLAST+ 2.17.0). `edit_distance` is the global edit distance between the whole returned sequence
and the whole truth sequence, so it also counts missing and extra bases. `truth_sequences_identical` gives how
many distinct truth sequences of the case were returned exactly, for example `2 of 2` when both alleles of a
mixture are returned. When a tool returns several sequences, the row describes the best one and
`all_returned_sequences` lists them all.

## Results

Exact sequences among the 25 *H. pylori* cases: Locus-Recon 19, ARIBA 19 (both
configurations) and KMA 20.

| Cases | Locus-Recon | ARIBA | KMA |
|---|---|---|---|
| *gyrB*, *ureB*, *glmM* (15) | 15 exact, all accepted | 15 exact | 15 exact |
| 23S rRNA gene (5) | 4 exact and accepted; Hpfe0002 with one mixed site reported at position 1,639. Hpfe0006: 2,570 of 2,887 bp, held, with the missing 317 bp located on a second scaffold | 4 exact. Hpfe0002: T for C at position 1,639, with the flag `variants_suggest_collapsed_repeat` | 5 exact. Hpfe0002: position 1,639 in lower case |
| *cagA* (5) | No exact sequence and none accepted: candidates identical to the genome over the overlap, with ends wrong by 6 to 33 bp | No exact sequence: fragments of 372 to 1,323 bp in four isolates, none in the fifth (624 to 2,457 bp with `--nucmer_min_id 80`) | No exact sequence: consensus of 3,341 to 3,412 bp at an edit distance of 145 to 239 from the gene, without a warning |
| Simulated, one allele (7) | 7 exact; 3 accepted, 4 not accepted at 3× to 20× | 7 exact | 7 exact; 273 of 507 bases in lower case at 3× |
| Simulated mixtures (4) | Major allele exact; a mixture is reported and the sample is held in all 4 | One allele, exact: the major one at 5, 10 and 20%, the second one at 50%. A flag for mixed bases at 20 and 50% only | Both alleles exact in all 4, with depths that follow the mixing fraction |
| Paralogue control (1) | A sequence with wrong bases, held as ambiguous | The true allele, exact, and the paralogue as a second sequence (`region_assembled_twice`) | A sequence with three wrong bases and one gap; no warning |
| Sample without the locus (1) | Nothing returned | Nothing returned | Nothing returned |
| LIBA-6656 *tcdB* | One fragment of 1,955 bp, held: 23 mixed sites, both ends missing, locus split | 24 fragments of at most 1,421 bp (25 with `--nucmer_min_id 80`) | Two full-length sequences and a 584-bp fragment |

The two full-length sequences that KMA returns for LIBA-6656 are close to the two analyst-built candidates. The
consensus on template B7.6 (7,104 bp) is identical to the chromosome-associated candidate. The consensus on
template B7.9 (7,100 bp) is the extrachromosomal candidate without four bases, AGGA at positions 5,616 to 5,619,
and has no mismatch. The reads do not support that deletion: in the competitive alignment of Supplementary Data 1,
67 reads span positions 5,606 to 5,629 and none carries a deletion (`kma_B7_9_gap_read_check.tsv`).

KMA took 5.5 to 9.1 s per *H. pylori* isolate for the five loci and ARIBA 17.4 to 26.8 s. Locus-Recon was not
timed here; it needs a draft assembly and runs a local assembly and two read mappings for each locus.

## What the comparison shows

- For a single-copy locus with a close reference, the three tools return the same sequence. Locus-Recon does not
  return better sequences than KMA or ARIBA in these data.
- KMA is the quickest route to a sequence when a close template exists. It recovered the 23S rRNA gene that
  Locus-Recon left incomplete, and it separated the two *tcdB* copies of LIBA-6656 because the catalogue holds a
  close allele of each.
- A consensus does not say when it is wrong. KMA returned a sequence that differs from the truth for every
  *cagA* gene and for the paralogue control, and one that lacks four bases for one *tcdB* copy, each time without
  a warning.
- Locus-Recon accepted no wrong sequence. It also reports why a sequence was not accepted and estimates how
  many copies the locus has. The comparators report depth but no copy number.

The comparison has limits. It covers one species with real truth, one bait strain and default settings. SRST2 and
stringMLST were not run. The Locus-Recon
candidates of the *H. pylori* and simulated cases are not stored in this repository, so their relation to the
truth is derived from the deposited comparison fields.

## Files

| File | Content |
|---|---|
| `comparator_benchmark_results.tsv` | One row per case and tool configuration: 39 cases × 4 = 156 rows |
| `comparator_benchmark_summary.tsv` | Counts of each relation per group of cases and tool configuration |
| `comparator_returned_sequences.fasta` | Every sequence returned by ARIBA and KMA |
| `comparator_truth_sequences.fasta` | The truth sequences, and the two LIBA-6656 candidates |
| `*_catalogue.fasta`, `*_catalogue_ariba_metadata.tsv` | The reference sequences as given to the tools |
| `liba6656_analyst_built_candidates.fasta` | Input: the two analyst-built *tcdB* candidates |
| `kma_B7_9_gap_read_check.tsv` | Read evidence at the four bases that the KMA consensus lacks |
| `raw_tool_outputs.tar.gz` | Reports, sequences, logs and timing records of every ARIBA and KMA run |
| `PROVENANCE.json`, `TOOL_VERSIONS.txt`, `conda_environment.yml`, `conda_list.txt` | Software, inputs and parameters of the run |

## How to rerun

```bash
conda create --name lr-comparators --channel conda-forge --channel bioconda \
    python=3.12 ariba=2.14.7 kma=1.7.1 blast=2.17.0 "setuptools<81"
conda activate lr-comparators
python -m pip install .

# full run: downloads about 3 GB of public reads into the work directory
python validation/run_comparator_benchmark.py --workdir /path/to/work --threads 16

# scoring only, from the tool outputs kept here (needs BLAST+ and network access to NCBI for the annotations)
mkdir work && tar -xzf validation/comparator-benchmark/raw_tool_outputs.tar.gz -C work
python validation/run_comparator_benchmark.py --workdir work --outdir work/out --score-only
```

Without `--outdir`, the full run writes into this folder.
