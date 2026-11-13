# Locus-Recon

**Recover the sequence of one bacterial locus, and evidence on its copy number, from the short reads you already have.**

Draft assemblies often fail at the loci that matter most. A typing gene is split between contigs, two related
copies collapse into one consensus, or an amplified resistance gene keeps its allele while its dosage changes.
Locus-Recon goes back to the paired-end reads. It recruits the reads of one locus, assembles them locally, and
reports the candidate sequence together with the evidence for it.

The report keeps four questions apart:

| Question | Field | Values |
|---|---|---|
| Did the workflow finish? | `workflow_status` | `SUCCESS`, `SKIP` (no gradable candidate), `FAIL` (input or tool error) |
| Do the reads support each reported base? | `sequence_confidence` | `HIGH`, `MEDIUM`, `LOW`, `SUSPECT` |
| Is the reported span the whole locus? | `span_clipped_bp`, `span_continuation_contig` | missing bases, and the contig that carries them |
| How does the sequence relate to your allele catalogue? | `catalogue_status` | exact, non-exact, divergent, highly divergent |

Identity to the catalogue is reported but does not lower the confidence tier, so a well-supported new allele is
accepted and marked for curation. Three companion commands estimate copy number from read depth, count flanking
contexts in the assembly graph, and list the alternative paths through a local graph.

Release 1.0.0. Python 3.9 or later; Linux or macOS; MIT licence.

## Contents

- [Is it suitable for my data?](#is-it-suitable-for-my-data)
- [Install](#install)
- [Check the installation](#check-the-installation)
- [Run it on your data](#run-it-on-your-data)
- [Read the results](#read-the-results)
- [Copy number](#copy-number)
- [What has been tested](#what-has-been-tested)
- [Limits](#limits)
- [Documentation](#documentation)
- [Citation and licence](#citation-and-licence)

## Is it suitable for my data?

You need, for each isolate:

- a draft assembly (FASTA);
- the paired-end Illumina reads of the same isolate (R1 and R2), ideally those used for the assembly;
- a FASTA file of known alleles of **one** locus, the *baits*.

One run reconstructs one locus in one or many isolates. Run the command again for another locus.

If you only need the sequence of a single-copy gene and your catalogue holds a close allele, a read-based
consensus tool such as [KMA](https://bitbucket.org/genomicepidemiology/kma) gives it in seconds, without an
assembly. In our [comparison](validation/comparator-benchmark/README.md), KMA, ARIBA and Locus-Recon returned the
same sequence for such loci. Use Locus-Recon when you need to know whether that sequence can be trusted (weak or
mixed read support, a gene cut at a contig end, a paralogue, a locus whose length varies) or how many copies
the locus has.

Locus-Recon does not use long reads, unpaired reads or metagenomes, and it does not test the whole genome for
contamination. It ships no allele catalogue: you supply the baits, for example from [PubMLST](https://pubmlst.org/).
A `SKIP` result does not prove that the locus is absent.

## Install

```bash
git clone https://github.com/crodsa/locus-recon.git
cd locus-recon
conda env create --file environment.yml     # or: mamba env create --file environment.yml
conda activate locus-recon
python -m pip install .
```

The environment pins Python 3.11.13, BLAST+ 2.16.0, BWA-MEM2 2.2.1, samtools 1.22.1 and SPAdes 4.2.0.
If you do not have Conda, install [Miniforge](https://github.com/conda-forge/miniforge) first.
Activate the environment in every new terminal with `conda activate locus-recon`.

To install exactly the released version, use the tag: `git clone --branch v1.0.0 https://github.com/crodsa/locus-recon.git`.

## Check the installation

```bash
locus-recon --version          # locus-recon 1.0.0

# checks inputs, tools and the bait database; maps and assembles nothing
locus-recon \
  --samplesheet validation/mock_dataset/samples.tsv \
  --main-output-dir validation/mock_dry_run \
  --bait validation/mock_dataset/mock_locus_alleles.fasta \
  --locus mockLocus --threads 2 --dry-run

# optional: the full 13-sample simulated benchmark (a few seconds per sample)
python validation/run_mock_validation.py --force --threads 2 --memory-gb 4
```

The dry run ends with `Dry run passed`. The benchmark writes its evaluation under `validation/mock_results/validation/`.

## Run it on your data

**1. Baits.** One FASTA file with the known alleles of one locus, all with the same start, end and orientation.
Record identifiers must be unique and at most 50 characters long. Remove alignment gaps. Do not mix loci or
paralogues. The file is checked before anything runs, and each problem is reported with the record that causes it.
For an rRNA gene or another non-coding target add `--noncoding-locus`.

**2. Sample sheet.** A tab-separated file with four columns. Relative paths are read from the folder of the sheet.
A template is in [`examples/samplesheet.tsv`](examples/samplesheet.tsv).

```tsv
sample_id	assembly_path	r1_path	r2_path
isolate_001	data/isolate_001.fasta	reads/isolate_001_R1.fastq.gz	reads/isolate_001_R2.fastq.gz
isolate_002	data/isolate_002.fasta	reads/isolate_002_R1.fastq.gz	reads/isolate_002_R2.fastq.gz
```

**3. Run.** Check the set-up with `--dry-run`, then run the batch.

```bash
locus-recon \
  --samplesheet study/samples.tsv \
  --main-output-dir results/target_locus \
  --bait study/target_locus_alleles.fasta \
  --locus target_locus \
  --threads 8 --parallel-samples 1 --memory-per-sample 16
```

CPU use is about `--threads × --parallel-samples`, and SPAdes memory about `--memory-per-sample × --parallel-samples`.
`--cleanup` deletes the per-sample alignments and local assemblies after a successful sample. Leave it off if you
intend to use the copy-number commands, which read those files.

## Read the results

Open the files in this order:

1. `locus_recon_report_<locus>.tsv`: one row per sample;
2. `<locus>_accepted_alleles.fasta`: the accepted sequences. Sequences that need review, held sequences and all
   candidates are written to three further FASTA files;
3. `<sample>/<sample>_<locus>_qc_report.txt`: the scorecard of one sample;
4. `<sample>/allele_remap.support.tsv`: one row of read evidence per base;
5. `run_manifest.json`: command, parameters, tool versions and checksums.

| Result | `result_disposition` | What it means | What to do |
|---|---|---|---|
| `SUCCESS` + `HIGH` | `PASS` | No flag was raised | Use it; check provenance and scheme rules before submitting a new allele |
| `SUCCESS` + `MEDIUM` or `LOW` | `REVIEW` | One or two moderate deviations | Read `qc_flags` and the per-base table before using the sequence. A sequence in review can be wrong |
| `SUCCESS` + `SUSPECT` | `HOLD` | A severe flag, or three or more deviations | Do not use it as an allele call without independent confirmation |
| `SKIP` | `HOLD` | No gradable candidate | Read `message`; this is not proof of absence |
| `FAIL` | `HOLD` | An input or an external tool failed | Fix the problem and run the sample again |

Three situations deserve attention:

- **Mixed sites.** `mixture_detected` and `ALLELE_MIXTURE` mean that the reads support two sequences at two or more
  positions. The cause may be two divergent copies, a paralogue or a mixed culture.
- **Locus cut by a contig end.** `ALLELE_SPAN_CLIPPED_AT_CONTIG_END` gives the number of missing bases.
  `LOCUS_SPLIT_ACROSS_CONTIGS` names the contig that carries them. High depth on the reported part says nothing
  about the missing part.
- **Divergent from the catalogue.** `catalogue_status` recommends curation. It does not lower the tier.

Every column and flag is defined in [docs/output-interpretation.md](docs/output-interpretation.md).

## Copy number

Identical copies collapse into one consensus and leave no mixed site. Depth and the assembly graph are the only
evidence that short reads keep.

```bash
# depth of the locus relative to the single-copy backbone, with a bootstrap interval
locus-recon-depth-ratio \
  --bam SAMPLE/full_map.sorted.bam \
  --locus contig_7:250422-251237 \
  --exclude plasmid_1 --tsv locus_depth.tsv

# the same estimate, plus the number of distinct flanking contexts in the graph of the draft assembly
locus-recon-copy-number \
  --bam SAMPLE/full_map.sorted.bam --locus NODE_19:901-3789 \
  --gfa draft_assembly/SAMPLE/assembly_graph_after_simplification.gfa \
  --bait loci.fasta --bait-id 23S \
  --output-tsv SAMPLE.23S.copy-number.tsv

# alternative paths through a local graph, when the local assembly collapsed divergent copies
locus-recon-graph-paths \
  --gfa SAMPLE/spades_local/assembly_graph_after_simplification.gfa \
  --output SAMPLE.graph_paths.fasta --summary SAMPLE.graph_paths.tsv \
  --min-length 6000 --max-length 8000
```

`SAMPLE/full_map.sorted.bam` is written by a reconstruction run without `--cleanup`. How to read the output:

- `MULTICOPY_DEPTH` means that the lower limit of the interval is above 1.5. `SINGLE_COPY_COMPATIBLE` is
  compatibility, not proof of one copy.
- Depth measures large amplifications well. For a small count, such as a two-copy rRNA gene, depth overestimated
  the number in our tests and the graph count was the better guide.
- A locus broken across contigs needs `--gene-span` and `--gene-length`. If the spans overlap on the gene, the
  dosage is labelled `REVIEW_OVERLAPPING_SPANS`. Compare it with `geometry_expected_dosage`: overlap alone can
  produce a value near two.
- When `ambiguity_index` is below 0.70 the estimate is a lower bound.
- A path listed by `locus-recon-graph-paths` is a proposal, not an allele call.

Details and worked examples are in the [user guide](docs/user-guide.md#when-every-copy-is-identical-the-depth-ratio).

## What has been tested

All inputs, outputs and scripts are under [`validation/`](validation/README.md). The numbers below are those of
release 1.0.0; [docs/validation-results.md](docs/validation-results.md) gives the detail.

| Test | Data | Result |
|---|---|---|
| Accepted sequences | 62 cases with known truth: 13 simulated, 14 constructed, 35 from real reads of five closed *Helicobacter pylori* genomes | 25 accepted, all 25 identical to the truth. The cases are not independent |
| New alleles | Real reads, one bait from another strain | Accepted sequences were 2.5 to 4.3% divergent from the bait. *cagA*, about 11% divergent and variable in length, was not accepted in any of five isolates |
| Completeness | 14 constructed assemblies | Missing bases measured exactly, from 0 to 822 bp |
| Mixtures and paralogue | Simulated reads | Mixtures with 10 to 50% of a second allele and a paralogue control were held |
| Depth, synthetic | 27 datasets with one to three copies | Median absolute error 0.04 copies |
| Depth, real amplification | *aphA1* in seven *Acinetobacter baumannii* runs | 10.5 and 77.7 copies, against 10 ± 2 and 75 ± 14 by published quantitative PCR |
| Low copy number | 23S rRNA gene (two copies) in five *H. pylori* isolates | Depth ratios of 2.55 to 3.22; the graph gave two copies in four isolates and a lower bound of two in the fifth |
| Collapsed divergent copies | *tcdB* of *Clostridioides difficile* LIBA-6656 | The workflow held a well-covered 1,955-bp fragment with 23 mixed sites. A dosage of 2.29 was shown to arise from overlapping alignments |
| Comparison with ARIBA 2.14.7 and KMA 1.7.1 | The same reads and reference sequences: 25 *H. pylori* loci, 13 simulated samples, LIBA-6656 | Exact sequences in 19 (Locus-Recon), 19 (ARIBA) and 20 (KMA) of the 25 *H. pylori* loci. KMA returned sequences that differ from the truth for *cagA* and for the paralogue control, marked only by lower-case bases; Locus-Recon accepted none of these and stated why |

## Limits

- Accuracy on real reads was tested in one species, with baits from one strain.
- The thresholds of the confidence tier are fixed constants. They were chosen during development and were not
  estimated from data. Validate them on a panel of your own organism and library type.
- A locus whose length varies among strains needs baits that represent that variation, or its candidates will not
  be accepted.
- The 1.5 threshold of the depth call was set by the author and was not fitted.
- Locus-Recon did not return better sequences than ARIBA or KMA when a close reference existed, and it is slower
  and needs a draft assembly. KMA is the quick route to the sequence. Locus-Recon adds the evidence that says
  when the sequence can be trusted and how many copies there are. SRST2 and stringMLST were not compared.
- When the copies of a locus are divergent and the catalogue holds a close allele of each, KMA can separate them
  and the standard Locus-Recon workflow cannot: for *tcdB* of LIBA-6656 it returned one held fragment.
- Short reads cannot close a replicon, assign distant variants to a copy or replace long reads when a result
  changes a conclusion.

## Documentation

| Document | Content |
|---|---|
| [User guide](docs/user-guide.md) | The method in plain language, inputs, a full run, results, copy number, incomplete and collapsed loci |
| [Command-line reference](docs/cli-reference.md) | Every option and default, resources, clean-up, exit codes |
| [Outputs and interpretation](docs/output-interpretation.md) | Every file, column, tier rule and flag |
| [Troubleshooting](docs/troubleshooting.md) | Errors and unexpected results, step by step |
| [Methods](docs/validation.md) | The computational method and how to design your own validation |
| [Validation results](docs/validation-results.md) | Results of each validation analysis and where its files are |
| [Validation manual](validation/README.md) | How to rerun the bundled analyses |
| [Changelog](CHANGELOG.md) | What each version contains |

Questions and bug reports: <https://github.com/crodsa/locus-recon/issues>.

## Citation and licence

A manuscript that describes Locus-Recon 1.0.0 is in preparation. Until it is published, please cite the software
name and version, this repository, and the name and release of the allele database that you used as baits.
[`CITATION.cff`](CITATION.cff) holds the citation metadata, and the DOI of the archived release will be added there.

Software: [MIT licence](LICENSE). Synthetic data, validation records and figures: [CC BY 4.0](LICENSE-DATA.md).
Public sequencing runs keep the terms of their sources.

Release 1.0.0 is the tag `v1.0.0`. Its package code is identified by the SHA-256 digest
`3863abc32fb5446a72919c790cb2558217871410f441db1ecf17dde3de55cdcf` of `locus_recon/*.py`, and the
[changelog](CHANGELOG.md) lists what the release contains.
