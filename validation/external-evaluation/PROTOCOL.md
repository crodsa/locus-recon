# Protocol of the external evaluation and of the targeted tests

Written on 18 November 2026, before any of the analyses below was run. It answers a critical review of the
manuscript, which found that the evaluation of release 1.0.0 was small, partly used during development and
concentrated in *Helicobacter pylori*. Any departure from this protocol will be listed in `DEVIATIONS.md` in
this folder, with its reason. Results will be reported whether or not they favour Locus-Recon.

## 1 What is frozen

- Locus-Recon 1.0.1, commit `00ef038` of the branch `release-1.0.1`. Package code digest (SHA-256 of
  `locus_recon/*.py`, computed by `validation/common.py`): `2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0`.
- No threshold, default or line of package code will change after this commit. If a defect is found, it will be
  reported, and results obtained before and after any correction will be kept apart.
- All commands use default options, except `--noncoding-locus` for rRNA genes and the options that set threads
  and memory. External tools are those of `environment.yml`. ARIBA 2.14.7 and KMA 1.7.1 are run with default
  options, as in `validation/comparator-benchmark/`.

## 2 External evaluation

### 2.1 Isolates

Source: FDA-ARGOS (BioProject PRJNA231221), which holds complete genomes with Illumina reads of the same
isolate. None of these isolates, species or libraries was used to develop or test Locus-Recon.
*H. pylori*, *Acinetobacter baumannii* and *Clostridioides difficile* are excluded because they were.

Sixteen species were listed from the species names of the project, before any run or sequence was examined.
They span 30 to 68% G+C and five phyla or classes. For
each species the isolates with a current complete RefSeq assembly were ordered by FDA-ARGOS strain number, and
the first that had an Illumina paired-end whole-genome run of at least 40× in the European Nucleotide Archive
was taken. When an isolate has several such runs, the one with most bases is used. `isolate_selection_log.tsv`
lists every record examined. In all sixteen species the first isolate met the rule.

| Species | Isolate | Assembly (truth) | Illumina run | Instrument | Depth | G+C (%) | Bait reference |
|---|---|---|---|---|---|---|---|
| *Staphylococcus aureus* | FDAARGOS_1 | GCF_001019255.2 | SRR1955492 | HiSeq 2500 | 534× | 33.0 | NCTC 8325 (GCF_000013425.1) |
| *Campylobacter jejuni* | FDAARGOS_262 | GCF_002209005.1 | SRR4125014 | HiSeq 2500 | 640× | 30.0 | NCTC 11168 (GCF_000009085.1) |
| *Listeria monocytogenes* | FDAARGOS_57 | GCF_001188045.3 | SRR1950322 | HiSeq 2000 | 438× | 38.0 | EGD-e (GCF_000196035.1) |
| *Enterococcus faecalis* | FDAARGOS_324 | GCF_003030425.1 | SRR5419826 | HiSeq 2500 | 397× | 37.5 | T5 (GCF_000393015.1) |
| *Streptococcus pyogenes* | FDAARGOS_149 | GCF_001559175.2 | SRR2820674 | HiSeq 2500 | 680× | 38.5 | NCTC12064 (GCF_900475035.1) |
| *Haemophilus influenzae* | FDAARGOS_199 | GCF_002073475.2 | SRR4123928 | HiSeq 2500 | 675× | 38.0 | FDAARGOS_1560 (GCF_020736045.1) |
| *Bacillus cereus* | FDAARGOS_780 | GCF_013267475.1 | SRR9641652 | HiSeq 4000 | 561× | 35.5 | ATCC 14579 (GCF_006094295.1) |
| *Vibrio parahaemolyticus* | FDAARGOS_51 | GCF_001188185.2 | SRR1784833 | HiSeq 2000 | 221× | 45.5 | RIMD 2210633 substr. RIMD 2210633 (GCF_000196095.1) |
| *Escherichia coli* | FDAARGOS_144 | GCF_002944935.1 | SRR2820658 | HiSeq 2500 | 219× | 50.5 | K-12 substr. MG1655 (GCF_000005845.2) |
| *Salmonella enterica* | FDAARGOS_54 | GCF_001188125.2 | SRR1950246 | HiSeq 2000 | 207× | 52.0 | LT2 (GCF_000006945.2) |
| *Neisseria meningitidis* | FDAARGOS_209 | GCF_002073575.2 | SRR4272082 | HiSeq 2000 | 640× | 52.0 | PartJ-Nmeningitidis-RM8376 (GCF_022869645.1) |
| *Corynebacterium diphtheriae* | FDAARGOS_197 | GCF_002073375.2 | SRR4271515 | HiSeq 2000 | 452× | 53.5 | NCTC11397 (GCF_001457455.1) |
| *Klebsiella pneumoniae* | FDAARGOS_156 | GCF_002951635.1 | SRR2821776 | HiSeq 2500 | 185× | 57.0 | HS11286 (GCF_000240185.1) |
| *Mycobacterium tuberculosis* | FDAARGOS_750 | GCF_009730235.1 | SRR9176616 | HiSeq 4000 | 415× | 65.5 | H37Rv (GCF_000195955.2) |
| *Pseudomonas aeruginosa* | FDAARGOS_501 | GCF_003812905.1 | SRR8178086 | HiSeq 4000 | 339× | 66.0 | PAO1 (GCF_000006765.1) |
| *Bordetella pertussis* | FDAARGOS_1365 | GCF_016904355.1 | SRR13674456 | NovaSeq 6000 | 610× | 67.5 | H640 (GCF_004008975.1) |

Reads are subsampled to 100× of the genome length with a fixed seed (20261118) when the run holds more. The
draft assembly is made with SPAdes in isolate mode, as for the full-depth *H. pylori* libraries.

### 2.2 Loci and baits

Eight loci per isolate, fixed here: five protein-coding genes expected in one copy (*gyrB*, *rpoB*, *recA*,
*dnaK*, *atpD*) and three whose copy number varies among species (23S rRNA gene, 16S rRNA gene, *tuf*).
This gives 128 units of isolate and locus.

The bait of each locus is one sequence, taken from the RefSeq reference genome of the species (last column of
the table), which is never the evaluated isolate. A protein-coding locus is found in the reference by its gene
symbol in the RefSeq annotation, and an rRNA gene by its product name; when several copies exist, the first by
coordinate is the bait. If the reference annotation lacks the gene, the unit is recorded as not available and
is not replaced.

### 2.3 Truth

The bait is aligned to the complete genome of the isolate with BLAST+. Every alignment of at least 80%
identity that covers at least 90% of the bait defines one copy, and the copy is extended to the full span of
the bait. The truth of a unit is the set of these copies, and its copy number is their count. A locus expected
in one copy that has more is analysed in the stratum of its observed copy number. A unit without any copy is a
negative control.

The truth is not fully independent of the reads: FDA-ARGOS assemblies combine PacBio and Illumina data of the
same isolate. This limit will be stated with the results.

### 2.4 Outcomes of reconstruction

For each unit and each tool: whether a sequence was returned; whether it is identical to a copy of the truth
over its whole length; if not, whether it is incomplete but identical over its length, or carries wrong, extra
or missing bases; and, when the truth copies differ from each other, which copy it matches. The scoring code
is that of `validation/run_comparator_benchmark.py`, extended to read this manifest.

For Locus-Recon the disposition (`PASS`, `REVIEW`, `HOLD`), the tier and every flag are recorded. Reported:

1. the accepted sequences that are not identical to the truth, each listed by isolate and locus;
2. the proportion of units accepted, by stratum (one copy; several identical copies; several different copies);
3. among the units not accepted, how many were complete and exact, incomplete, or wrong;
4. accepted sequences that carry `MIXED_SITE_BELOW_MIXTURE_CALL`, and whether the truth copies differ at that site.

The unit of generalisation is the isolate. Counts are given per isolate and per locus. Any interval is obtained
by resampling isolates, and no interval will treat the 128 units as independent.

### 2.5 Decision rules for KMA and ARIBA, fixed here

Neither tool states a decision, so a rule is needed to compare decisions and not only sequences.

- **KMA, primary rule:** accept the consensus when it has no lower-case base and no `n`, and the template
  coverage in the `.res` file is 100%. **Secondary rule:** the same without the condition on coverage.
- **ARIBA rule:** accept when one sequence is returned for the cluster and the flag has `assembled`,
  `assembled_into_one_contig` and `unique_contig` set, has `region_assembled_twice`, `scaffold_graph_bad`,
  `assembly_fail` and `variants_suggest_collapsed_repeat` unset, and, for a protein-coding locus, has
  `complete_gene` set.

For each tool the units are cross-tabulated by decision (accepted or not) and result (exact or not). The value
that the decision of Locus-Recon adds is read from these tables, and so is any value that it does not add.

### 2.6 Outcomes of copy number

For every unit: depth ratio and its interval, depth call, dosage and its status, graph context count and
status, and the consensus with its status. Reported against the copy number of the truth:

1. for one-copy units, the proportion called `MULTICOPY_DEPTH` (false alerts of the 1.5 threshold);
2. for units with two or more copies, the proportion called `MULTICOPY_DEPTH`;
3. the error of the ratio and of the dosage, by copy number, G+C content of the locus and of the genome,
   locus length and number of discovery intervals;
4. how often the bootstrap interval contains the truth;
5. how often the graph count equals the truth, is a correct lower bound, or is wrong;
6. the consensus value and status against the truth, with every case of `SUPPORTED` that is wrong listed.

## 3 Targeted tests

Simulated reads carry substitution errors, as in `validation/generate_mock_dataset.py`, with fixed seeds. Each
scenario states its truth, and every replicate is reported.

| Test | Design | Reported |
|---|---|---|
| Mixtures with few sites | Two alleles that differ at 1, 2, 5 or 10 sites; minor fraction 5, 10, 20, 35 or 50%; 60×; 10 replicates each (200 runs), and 20 pure controls | Per cell: held as mixture, accepted, accepted with the mixed-site flag. Detection limit by number of sites |
| Strand-bias alerts | Pure samples: loci of 500, 1,500, 3,000 and 7,000 bp at 20, 60, 200 and 1,000×, 10 replicates; and every one-copy unit of the external evaluation | Rate of `STRAND_BIASED_ALT_SUPPORT` by length and depth |
| Ends, indels and pseudogenes | Sample alleles that differ from the bait by substitutions within 5, 20 or 50 bp of an end; a 3-bp deletion; a 12-bp insertion; a 1-bp deletion; two compensating indels 30 bp apart; a premature stop codon. Depths 10, 20 and 60×, 5 replicates | Exact or not; tier; uncertain bases in the interior and in the whole sequence; frame flags against the true frame |
| Graph contexts | Loci in 2, 3 and 4 copies with different flanks; 2 copies that share one flank for 500 or 2,000 bp; a tandem duplicate; 1 copy with a 1-kb repeat in a flank; 1 copy in a mixture of two strains that differ at one site of the flank; 2 copies in such a mixture; copies identical or 1% divergent. 60×, 5 replicates | Context counts, status and consensus against the true copy number |
| Block size of the bootstrap | Blocks of 100, 250, 500 (default) and 1,000 bp for the *aphA1* runs, the 23S rRNA gene of *H. pylori* and the external units, through the library function | Interval width and whether it contains the reference value |
| Overlapping spans | A locus in 1 or 2 copies, in drafts cut so that the discovery intervals overlap on 0, 10, 25, 45, 65, 83 or 90% of the query. 60×, 5 replicates | Dosage, `geometry_expected_dosage` and status for both copy numbers at each overlap |

## 4 Existing analyses

The analyses of release 1.0.0 (simulated benchmark, *H. pylori*, *aphA1*, LIBA-6656) will be run again with
1.0.1 where their inputs are available. Sequences, per-base tables, logs, parameters and checksums will be
kept, so that each statement of exactness can be checked against a sequence. Records made by an earlier version
will be kept apart and labelled.

## 5 What this protocol does not provide

- A truth that is independent of short reads.
- Libraries from more than one consortium: all sixteen come from FDA-ARGOS.
- Loci chosen for difficulty. The eight loci were chosen because every species has them; long, variable or
  repeated virulence genes such as *cagA* are not represented.
- An amplified resistance gene with a reference value. The copy numbers evaluated are those of rRNA genes
  and *tuf*.
