# GvizRegulatoryTracks 0.1.0

- Hardened score parsing, coordinate and display-parameter validation, and
  range-native overlap detection ahead of the initial release.
- Numeric factor and character metadata columns are accepted when every
  non-missing value parses as a finite number.
- `ScoreSequenceTrack()` and its shorter alias `DynSeqTrack()` draw signed,
  score-scaled nucleotide sequences with standard Gviz quantitative axes.
- Score tracks accept numeric vectors, scored `GRanges`, BigWig, and bedGraph;
  reference bases can come from FASTA, 2bit, BSgenome, DNAStringSet, or a Gviz
  `SequenceTrack`.
- Wide score views switch automatically to resolution-aware Gviz signal
  rendering and support shared symmetric scales across tracks.
- Internal signal tracks disable redundant Gviz range collapsing, avoiding a
  matrix-slot failure when only a subset of adjacent score ranges collapses.
- Motif overview ranges now resolve their score-controlled fill, opacity, and
  border through one tested path shared with the per-hit score styles.
- `MotifLogoTrack()` positions strand-aware information-content logos from
  BED or `GRanges` hits and MEME, matrix, or `universalmotif` motifs.
- Overlapping motif matches use deterministic lanes and fixed physical logo
  heights; wide views switch automatically to genomic ranges.
- Motifs can be colored by nucleotide, motif identity, quantitative score, or
  combined identity and score encodings using fill, brightness, opacity, and
  borders.
- Quantitative motif encodings include configurable continuous legends and
  visible-range, global, or fixed score limits.
- Continuous legends now reserve an internal title band, preventing multiline
  titles from crossing the legend border in compact tracks.
- Added `DirectionalRangeTrack()` for deterministic stranded range lanes with
  shallow, fixed-physical-width tips that never exceed the rectangle height.
- Directional ranges can alternatively use rectangular bodies with subtle
  internal strand glyphs and above-range labels.
- Directional range bodies now use the full closed `GRanges` interval, so
  width-one ranges remain visible, and strand glyphs follow reversed axes.
- `MotifLogoTrack()` supports an independent display-label column, and both
  motif and directional labels are kept inside visible viewport edges on
  forward and reversed axes.
- Visible-range motif score scaling uses parsed numeric values for factor and
  character metadata rather than factor level codes.
- Minus-strand motifs are reverse-complemented into reference-genome
  orientation, including reversed Gviz axes.
- Track boundaries, nucleotide palettes, score normalization, aggregation,
  axes, labels, motif dimensions, and legends are configurable.
- Two reproducible real-data examples use compact BigWig, BED, MEME, and FASTA
  fixtures with explicit coordinate and provenance documentation.
