# GvizRegulatoryTracks 0.1.0

- `ScoreSequenceTrack()` and its shorter alias `DynSeqTrack()` draw signed,
  score-scaled nucleotide sequences with standard Gviz quantitative axes.
- Score tracks accept numeric vectors, scored `GRanges`, BigWig, and bedGraph;
  reference bases can come from FASTA, 2bit, BSgenome, DNAStringSet, or a Gviz
  `SequenceTrack`.
- Wide score views switch automatically to resolution-aware Gviz signal
  rendering and support shared symmetric scales across tracks.
- `MotifLogoTrack()` positions strand-aware information-content logos from
  BED or `GRanges` hits and MEME, matrix, or `universalmotif` motifs.
- Overlapping motif matches use deterministic lanes and fixed physical logo
  heights; wide views switch automatically to genomic ranges.
- Motifs can be colored by nucleotide, motif identity, quantitative score, or
  combined identity and score encodings using fill, brightness, opacity, and
  borders.
- Quantitative motif encodings include configurable continuous legends and
  visible-range, global, or fixed score limits.
- Minus-strand motifs are reverse-complemented into reference-genome
  orientation, including reversed Gviz axes.
- Track boundaries, nucleotide palettes, score normalization, aggregation,
  axes, labels, motif dimensions, and legends are configurable.
- Two reproducible real-data examples use compact BigWig, BED, MEME, and FASTA
  fixtures with explicit coordinate and provenance documentation.
