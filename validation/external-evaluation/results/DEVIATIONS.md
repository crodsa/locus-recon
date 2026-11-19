# Departures from the protocol of the external evaluation

Protocol: `protocol_PROTOCOL.md` (section 2, and the two targeted tests of section 3 that use the external
units). Each entry states what the protocol says, what was done, why, and where the effect can be seen. No
threshold, default or line of `locus_recon/` was changed. Entries 1 to 3 change which units enter which
summary; the others do not change a result.

## 1 Twelve units have no bait under the rule of protocol 2.2

**Protocol:** a protein-coding bait is found by its gene symbol in the RefSeq annotation of the reference; if
the annotation lacks the gene, the unit is not available and is not replaced.

**Found:** 12 of the 128 units lack the exact symbol. The annotation of *S. aureus* NCTC 8325
(GCF_000013425.1) gives no symbol to *gyrB*, *atpD* and *tuf*; that of *K. pneumoniae* HS11286
(GCF_000240185.1) gives no gene symbol at all; and *tuf* is annotated `tufA` and `tufB` in *E. coli* K-12,
*S. enterica* LT2 and *P. aeruginosa* PAO1 (`logs/bait_symbol_probe_2026-11-18.txt`).

**Done:** the protocol set is the 116 units that the rule finds, and every headline number is given for it.
The 12 others were run in the same way with a bait found by a fallback rule, as a supplementary set
that is labelled in every table (`unit_set`) and summarised apart. The fallback was written after the
reference annotations had been read and before any tool was run on any reads: first the symbol followed by a
paralogue letter (`tufA`, `tufB`; first by coordinate), then the product name of the CDS (`DNA gyrase subunit
B`; `DNA-directed RNA polymerase subunit beta`; `recombinase RecA` or `recombinase A`; `molecular chaperone
DnaK`; `F0F1 ATP synthase subunit beta`; `elongation factor Tu`).

- Symbol with a paralogue letter: FDAARGOS_144 tuf, FDAARGOS_54 tuf, FDAARGOS_501 tuf.
- Product name: FDAARGOS_1 gyrB, FDAARGOS_1 atpD, FDAARGOS_1 tuf, FDAARGOS_156 gyrB, FDAARGOS_156 rpoB, FDAARGOS_156 recA, FDAARGOS_156 dnaK, FDAARGOS_156 atpD, FDAARGOS_156 tuf.

**Why:** the genes are present in these references and only their symbols are missing, so leaving the units
out loses information; pooling them would break the rule that was registered.

## 2 Three units have no truth copy although the locus is present

**Protocol:** a unit without any copy (no alignment of at least 80% identity over at least 90% of the bait) is
a negative control.

**Found:** in 3 units of the protocol set no alignment covers 90% of the bait, yet alignments of at
least 80% identity cover 92.5% to 100.0% of it in pieces
(`baits.tsv`, `bait_fraction_aligned`): FDAARGOS_199 23S, FDAARGOS_780 recA, FDAARGOS_54 23S. The bait
of the 23S rRNA gene of *H. influenzae* and of *S. enterica* aligns to each copy of the isolate in two or more
pieces, as intervening sequences that differ between reference and isolate would produce, and the two halves
of the *recA* bait align 328 bp apart in the genome of *B. cereus* FDAARGOS_780
(`per_unit/<isolate>/<locus>/bait_vs_complete_genome.blastn.tsv`).

**Done:** these units are not counted as negative controls. They form a stratum of their own ("no copy under
the protocol rule, bait aligned in parts"), in which exactness and copy number are not defined; what each tool
returned and decided is reported. The limit that separates this stratum from a negative control (alignments
covering at least 50% of the bait) was set when the bait alignments were read, before any tool result existed.
No unit fell below it, so the evaluation has no negative control.

## 3 One isolate has reads that do not match its complete genome

**Protocol:** the truth of an isolate is its complete genome. The protocol foresaw that the truth is not
independent of the reads; it did not foresee the opposite.

**Found:** after scoring, the six sequences that Locus-Recon accepted and that differ from the truth all
belonged to *N. meningitidis* FDAARGOS_209, and KMA and ARIBA returned the same sequences. Each draft was then
compared with its complete genome with MUMmer dnadiff, which the protocol does not include. The draft of
FDAARGOS_209 aligns to 93.77% of the genome at 97.35% identity with 45,838 substitutions;
the other 15 align to at least 98.85% at 99.97% or more (`isolates_and_drafts.tsv`,
`concordance_reports/`). In this isolate 11 of the 15 truth copies are identical to the bait, which comes from
GCF_022869645.1 (PartJ-Nmeningitidis-RM8376); the sequences that Locus-Recon returned for the eight loci of this isolate differ from the truth at 0.2 to 4.0% of their bases (`external_units.tsv`).

**Done:** nothing was removed. FDAARGOS_209 is in every table and in every summary of the protocol set. An
additional set, "protocol, isolates whose draft agrees with the complete genome", repeats every summary
without it. Its limits (99.9% identity of the one-to-one alignments and 98.0% of the genome aligned) were chosen
after the results were seen; the 15 other isolates are far from them. Both sets of numbers must be read
together, and the second is not a registered analysis.

**Not examined:** why run SRR4272082 differs from assembly GCF_002073575.2 (another strain, a mixture, or a
mislabelled record). The draft is 2,940,763 bp in 3,193 contigs for a genome of 2,181,232 bp.

## 4 Reading of the product name of rRNA genes

The annotation of *S. aureus* NCTC 8325 writes `23S Ribosomal RNA` and `16S Ribosomal RNA`. Product names of
rRNA genes are compared without regard to case, and these two units are in the protocol set.

## 5 Additions that are not in the protocol

- Column `returned_occurs_in_complete_genome` of `external_units.tsv`: whether the returned sequence occurs,
  base for base, in the complete genome. It is descriptive and enters no decision table.
- The dnadiff comparison of entry 3.
- `rerun_check_FDAARGOS_262.tsv`: a second run of one isolate with the final runner.
- `counts_per_isolate.tsv`, `counts_per_locus.tsv`, `copy_number_error_by_factor.tsv` and
  `external_headline.tsv` are further arrangements of the registered outcomes.

## 6 Choices that the protocol leaves open

They are listed in `README.md` ("Choices that the protocol leaves open") and were fixed before any tool was
run on the reads: the subsampling command, `contigs.fasta` as the draft, the BLAST task of the truth and the
merging of overlapping alignments, the selection of intervals from `blast_step1.tsv`, and the scoring of a
continuous dosage. The comparator benchmark runs ARIBA a second time with `--nucmer_min_id 80`; that
configuration is not in the protocol and was not run.

## 7 How the run was carried out

- **Downloads.** One connection to the ENA file server gave about 0.3 to 1 MB per second and often closed
  early without an error. A first trial with one isolate stopped at the MD5 check of a truncated file, and no
  result came from it. The runner then received a download that resumes and tests the size that ENA reports.
  Every file used has the MD5 that ENA publishes (`input_reads.tsv`).
- **Versions of the runner.** The tool stages ran with the runner whose SHA-256 begins `cfb1367846c5`
  (assemble, reconstruct, copy-number, block-size), the comparators and the second subsampling with `6c319149b31d`,
  and the dnadiff and scoring stages with the final file, `ef61c35b0b82` (`PROVENANCE.json`, `stages`). The changes
  between them are the fifth stratum of entry 2, the dnadiff stage and the additional set of entry 3, and
  report tables. The functions that build the commands of SPAdes, Locus-Recon, ARIBA and KMA did not change
  after the first tool run, and the second run of FDAARGOS_262 with the final file reproduced all 584 values.
- **Order and threads.** Isolates were processed as their reads arrived, one at a time, with 8 threads (SPAdes
  with 8, two Locus-Recon runs of 4 threads at a time); the comparators ran three isolates at a time with 4
  threads each. The stage records in `PROVENANCE.json` therefore name the last isolate processed.
- **Bootstrap blocks.** `estimate_locus_copy_number()` keeps its default cap of 512 backbone blocks at every
  block size, so the length of backbone that is resampled grows with the block.
