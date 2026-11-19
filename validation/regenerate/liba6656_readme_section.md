## Verification with the release code

On 18 November 2026 the two commands above were rerun with Locus-Recon 1.0.1
(source digest `2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0`)
from inputs whose checksums equal those in `INPUT_CHECKSUMS.tsv`, once with four
threads and once with six, with a 12 GB SPAdes memory cap, under Python
3.11.13 and the tool versions listed above. `verify_release_rerun.py` compares
the reruns with this deposit and writes `RELEASE_VERIFICATION.json`. The same
verification was made with release 1.0.0 on 28 September 2026, with four and
eight threads; its record is in the history of the repository (tag `v1.0.0`),
and `../regenerated-1.0.1/old_vs_new.tsv` lists both records value by value.

- Both reruns recruited the same 6,195 read names and wrote byte-identical
  graphs and reported sequences.
- Their graph is the deposited graph up to numbering: the same 26 segment
  sequences and the same 32 links once segments are matched by sequence.
  Segment identifiers differ, and two segments differ in k-mer coverage, by 2
  of 27,945 and by 8 of 69,108 k-mers.
- The 512 terminal paths exported from it are identical to the deposited paths
  in sequence and numbering, so paths 210 and 315 and the candidate alignments
  above are unchanged. In the path summary only the node identifiers differ,
  and the mean graph depth by at most 0.001.
- The standard workflow reported one 1,955 bp `tcdB` sequence at `SUSPECT`,
  disposition `HOLD`. It carried `ALLELE_MIXTURE` (23 bidirectional sites,
  median alternative fraction 0.404) and `ALLELE_SPAN_CLIPPED_AT_CONTIG_END`
  (1,345 bp at the contig start and 3,804 bp at its end, 5,149 bp in all);
  `LOCUS_SPLIT_ACROSS_CONTIGS` located 1,345 of the missing bases on a second
  contig, and `FRAME_LENGTH_SHIFT_TRUNCATED` attributes its length, which is
  not a multiple of three, to that truncation. Mixed sites near an allele
  fraction of 0.4 are consistent with two divergent copies collapsed into one
  contig, the situation the two deposited candidates resolve; the workflow
  withholds the consensus rather than reporting it as an allele.
- The reported sequence and the QC report have the same bytes as those that
  release 1.0.0 wrote. The report row differs from that of 1.0.0 in its
  version, run date and run time only.

Three kinds of file are in this directory. The graph, the terminal paths, the
two selected paths, the candidate alignments and `GFA_RERUN_PROVENANCE.json`
are those of 24 August 2026, made with the tools of release 1.0.0; they are the
reference of the verification and were not replaced. `RELEASE_VERIFICATION.json`
and the three `standard_workflow_*` files are the output of release 1.0.1. The
baits, the script and this text are inputs or documentation.
