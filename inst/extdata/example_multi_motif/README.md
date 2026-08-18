# Example 1: multi-motif locus

This directory contains only small, locus-restricted inputs needed for the
package demonstration. The source archive remains read-only.

## Coordinates

- Genome: mouse mm10.
- Window: `chr1:176724198-176724298`, one-based and closed.
- `motif_hits.bed` is BED6 and therefore zero-based, half-open on disk.
- For example, BED `chr1 176724243 176724254` imports as the eleven-base
  `GRanges` interval `chr1:176724244-176724254`.
- The BED `score` field is the relative FIMO score multiplied by 1000 and
  rounded to an integer, as required by BED6. Divide it by 1000 for plotting.
- The five real strict calls have displayed relative scores from 0.944 to
  0.971. Four are minus-strand and the ERR1 call is plus-strand.

## Files

- `atac_seq_signal.bw`: pooled fragment-CPM ATAC-seq signal.
- `chrombpnet_count_contribution.bw`: signed, per-base ChromBPNet counts-head
  contribution scores.
- `motif_hits.bed`: five strict FIMO-derived HOCOMOCO-v14 motif matches in
  BED6 format; `name` is the MEME motif ID.
- `motifs.meme`: the five corresponding HOCOMOCO-v14 probability matrices.
- `locus_sequence.fa`: the 101-bp mm10 reference sequence; its header records
  the genomic interval because this is a locus fragment rather than a whole
  chromosome FASTA.
- `provenance.tsv`: archived source paths, sizes, checksums, extraction tools,
  genome build, and interval.

BigWig files retain one-based `GRanges` alignment when imported through
`rtracklayer`; BigWig's internal coordinate representation is handled by the
library. The FASTA sequence begins at genomic position 176724198, as stated in
its header and enforced by the demo. No preparation scripts are retained.
