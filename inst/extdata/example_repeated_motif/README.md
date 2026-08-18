# Example 2: repeated-motif locus

This compact real-data fixture tests the package tracks at a locus with two
separate matches of the same CTCF PWM. All archived sources remained read-only.

## Locus and selection

- Genome/window: mouse mm10 `chr11:103102780-103102890`, one-based and closed
  (111 bp).
- The broader source peak is `chr11:103102664-103103302`; the archived v0.0.35
  analysis classifies it as a promoter-like enhancer cCRE (`pELS`).
- Oct4_500 ChromBPNet contribution reaches 0.049902. Its two principal positive
  features align with the two displayed CTCF FIMO sites.
- The CTCF-like TF-MoDISco `pattern_0.174` seqlet has the maximum exported
  seqlet score (1000). Its corresponding `pattern_0` probability matrix is
  aligned to the displayed reference and shown with the standard nucleotide
  palette.
- Two spatially distinct strict calls of the exact same
  `CTCF.H14CORE.0.P.B` PWM have relative FIMO scores 0.940 and 0.884. This makes
  the local score-to-fill mapping visible without conflating related PWMs.

## Coordinates and formats

- `atac_seq_signal.bw`: Oct4_500 pooled fragment-CPM ATAC-seq signal.
- `chrombpnet_count_contribution.bw`: signed per-base Oct4_500 ChromBPNet
  counts-head contribution scores.
- `locus_sequence.fa`: the 111-bp mm10 reference sequence.
- `tfmodisco_seqlet.bed`: BED6 placement of TF-MoDISco `pattern_0`.
- `tfmodisco_pattern.meme`: the corresponding 50-position TF-MoDISco PPM.
- `fimo_hits.bed`: BED6 strict CTCF calls. BED coordinates are zero-based,
  half-open; `score` stores relative FIMO score multiplied by 1000 and rounded.
- `fimo_motifs.meme`: the matching HOCOMOCO-v14 CTCF PPM.
- `provenance.tsv`: archived paths, byte sizes, SHA-256 checksums, extraction
  coordinates, and tool versions.

For example, BED `chr11 103102811 103102831` imports as the one-based closed
`chr11:103102812-103102831` CTCF hit. The archived TF-MoDISco seqlet BED reports
`103102789-103102838`, a 49-bp half-open interval for a 50-position matrix. A
bidirectional PWM scan of the 111-bp fixture sequence places the aggregate
50-position `pattern_0` matrix at BED `chr11 103102792 103102842` on the minus
strand, which imports as `chr11:103102793-103102842`. This sequence-validated
placement is explicitly tested and recorded in `provenance.tsv`.
