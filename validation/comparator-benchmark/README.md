# Comparison with ARIBA and KMA

This folder holds a comparison of Locus-Recon with two established tools that also return a locus
sequence from short reads: ARIBA 2.14.7 and KMA 1.7.1. The three tools received the same reads and the same
reference sequences, and the sequences they returned were scored against the same truth. The sequences of
ARIBA and KMA are those of their runs of 11 and 13 November 2026. The sequences of Locus-Recon are those of
the regeneration of the validation analyses with release 1.0.1, of 18 November 2026.

The script is [`../run_comparator_benchmark.py`](../run_comparator_benchmark.py). It runs ARIBA and KMA. With the
option `--locus-recon-evidence` it scores the candidate sequences that the regeneration kept
(`evidence/` of the regeneration package, `validation_evidence/release_1.0.1/regenerated-1.0.1/` of Supplementary Data 1 of the article) with the function that scores the sequences of the other two tools.
Locus-Recon itself is not run by this script. Without that option the rows of Locus-Recon are read from
tables of this repository, as in release 1.0.0.

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
- The tools were run on 11 November 2026 (*H. pylori*, LIBA-6656, and KMA on the simulated samples) and on
  13 November 2026 (ARIBA on the simulated samples), from source distributions of release 1.0.0 with the same
  package code. `PROVENANCE.json` and the first scoring date from 13 November 2026.
- On 18 November 2026 the outputs of ARIBA and KMA kept in `raw_tool_outputs.tar.gz` were scored again, without
  running the two tools, together with the Locus-Recon candidates of release 1.0.1
  (`LOCUS_RECON_ROWS_PROVENANCE.json`). Every field of the ARIBA and KMA rows was reproduced.
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
| 23S rRNA gene (5) | 4 exact and accepted; Hpfe0002 with one mixed site reported at position 1,639. Hpfe0006: 2,570 of 2,887 bp, held, with the missing 317 bp located on a second scaffold | 4 exact. Hpfe0002: T for C at position 1,639, with the flag `variants_suggest_collapsed_repeat` and both nucleotides listed (121 reads with T, 123 with C) | 5 exact. Hpfe0002: position 1,639 in lower case |
| *cagA* (5) | No exact sequence and none accepted: candidates identical to the genome over the overlap, with ends wrong by 6 to 33 bp | No exact sequence: fragments of 146 to 1,323 bp in four isolates, none in the fifth (at most 2,457 bp with `--nucmer_min_id 80`) | No exact sequence: consensus of 3,341 to 3,412 bp at an edit distance of 145 to 239 from the gene; 13 to 65 bases in lower case |
| Simulated, one allele (7) | 7 exact; 3 accepted, 4 not accepted at 3× to 20× | 7 exact | 7 exact; 273 of 507 bases in lower case at 3× |
| Simulated mixtures (4) | Allele A exact (the major allele at 5, 10 and 20%); a mixture is reported and the sample is held in all 4 | One allele, exact: the major one at 5, 10 and 20%, the second one at 50%. A flag for mixed bases at 20 and 50% only | Both alleles exact in all 4, with depths that follow the mixing fraction |
| Paralogue control (1) | A sequence with wrong bases, held as ambiguous | The true allele, exact, and the paralogue as a second sequence (`region_assembled_twice`) | A sequence with three wrong bases and one extra base, all four in lower case |
| Sample without the locus (1) | Nothing returned | Nothing returned | Nothing returned |
| LIBA-6656 *tcdB* | One fragment of 1,955 bp, held: 23 mixed sites, both ends missing, locus split | 24 fragments of at most 1,421 bp (25 with `--nucmer_min_id 80`), with the flags `region_assembled_twice` and `variants_suggest_collapsed_repeat` | Two full-length sequences and a 584-bp fragment |

The two full-length sequences that KMA returns for LIBA-6656 are close to the two analyst-built candidates. The
consensus on template B7.6 (7,104 bp) is identical to the chromosome-associated candidate. The consensus on
template B7.9 (7,100 bp) is the extrachromosomal candidate without four bases, AGGA at positions 5,616 to 5,619,
and has no mismatch; KMA wrote the 16 bases around the gap in lower case. The reads do not support that deletion: in the competitive alignment of Supplementary Data 1,
67 reads span positions 5,606 to 5,629 and none carries a deletion (`kma_B7_9_gap_read_check.tsv`).

KMA took 5.5 to 9.1 s per *H. pylori* isolate for the five loci and ARIBA 17.4 to 26.8 s. Locus-Recon was not
timed on these reads; it needs a draft assembly and runs a local assembly and two read mappings for each locus.

## What the comparison shows

- For a single-copy locus with a close reference, the three tools return the same sequence. Locus-Recon does not
  return better sequences than KMA or ARIBA in these data.
- KMA is the quickest route to a sequence when a close template exists. KMA and ARIBA recovered the 23S rRNA
  gene that Locus-Recon left incomplete. KMA separated the two *tcdB* copies of LIBA-6656, except for four bases
  of one, because the catalogue holds a close allele of each.
- KMA returned a sequence that differs from the truth for every *cagA* gene and for the paralogue control, and
  one that lacks four bases for one *tcdB* copy. It raises no flag. It writes in lower case a base that is not
  significantly overrepresented in the reads, and each of these sequences had lower-case bases at or next to
  its errors. That mark does not say what is wrong and cannot show bases that are missing. The `.res` report
  gave a template coverage of 93.8 to 95.8% for the *cagA* sequences, but the true genes are also shorter than
  the bait. Correct sequences can carry lower-case bases too (273 of 507 bases at 3×, the mixed site of the 23S
  rRNA gene of Hpfe0002). The column `lowercase_bases` of the result table gives the count for the sequence
  that represents each case; the other sequences are in `comparator_returned_sequences.fasta`.
- ARIBA reported its doubt in flags and in the nucleotides and read depths that it lists for each difference
  from the reference. It did so for the wrong base of the 23S rRNA gene, for the fragmented *cagA* and *tcdB*
  assemblies, for the paralogue and for two of the four mixtures; it did not flag the 5 and 10% mixtures. In
  the paralogue control it was the only tool that returned the true allele exactly.
- Neither comparator states a decision about the sequence as a whole.
- Locus-Recon accepted no wrong sequence. It turns the evidence into a decision, states why a sequence was not
  accepted, measures how much of the locus is missing and estimates how many copies the locus has. The
  comparators report depth but no copy number.

The comparison has limits. It covers one species with real truth, one bait strain and default settings. SRST2 and
stringMLST were not run. In release 1.0.0 the Locus-Recon
candidates of the *H. pylori* and simulated cases had not been kept, and their relation to the truth was
derived from deposited comparison fields. They are now scored from the sequences, which are in
`comparator_returned_sequences.fasta`. `old_vs_new.tsv` of the regeneration package lists every field of the earlier
and of the new rows.

## Files

| File | Content |
|---|---|
| `comparator_benchmark_results.tsv` | One row per case and tool configuration: 39 cases × 4 = 156 rows |
| `comparator_benchmark_summary.tsv` | Counts of each relation per group of cases and tool configuration |
| `comparator_returned_sequences.fasta` | Every sequence returned by ARIBA and KMA, and every candidate of Locus-Recon |
| `LOCUS_RECON_ROWS_PROVENANCE.json` | The Locus-Recon reports and candidates that were scored, with their checksums |
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

# the same, with the Locus-Recon rows scored from the candidates of the regeneration
python validation/run_comparator_benchmark.py --workdir work --outdir work/out --score-only \
    --locus-recon-evidence /path/to/regenerated-1.0.1/evidence
```

Without `--outdir`, the full run writes into this folder.
