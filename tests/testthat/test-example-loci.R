test_that("repeated-motif fixture preserves sequence, score, and motif coordinates", {
    skip_if_not_installed("Biostrings")
    skip_if_not_installed("rtracklayer")
    fixture_dir <- system.file(
        "extdata", "example_repeated_motif",
        package = "GvizRegulatoryTracks"
    )
    expect_true(nzchar(fixture_dir))

    region <- GenomicRanges::GRanges(
        "chr11", IRanges::IRanges(103102780L, 103102890L)
    )
    sequence <- as.character(Biostrings::readDNAStringSet(
        file.path(fixture_dir, "locus_sequence.fa")
    )[[1]])
    contribution <- rtracklayer::import(
        file.path(fixture_dir, "chrombpnet_count_contribution.bw"),
        which = region
    )
    atac <- rtracklayer::import(
        file.path(fixture_dir, "atac_seq_signal.bw"), which = region
    )
    hits <- rtracklayer::import(file.path(fixture_dir, "fimo_hits.bed"))
    tfmodisco <- rtracklayer::import(
        file.path(fixture_dir, "tfmodisco_seqlet.bed")
    )

    expect_equal(nchar(sequence), 111L)
    expect_identical(IRanges::start(contribution), 103102780L:103102890L)
    expect_true(any(contribution$score < 0) && any(contribution$score > 0))
    expect_equal(max(contribution$score), 0.049902)
    expect_equal(sum(IRanges::width(atac)), 111L)
    expect_identical(as.character(hits$name), rep("CTCF.H14CORE.0.P.B", 2L))
    expect_identical(IRanges::start(hits), c(103102812L, 103102851L))
    expect_identical(IRanges::end(hits), c(103102831L, 103102870L))
    expect_identical(as.character(GenomicRanges::strand(hits)), c("+", "+"))
    expect_equal(hits$score / 1000, c(0.940, 0.884))
    expect_identical(as.character(tfmodisco$name), "pattern_0")
    expect_identical(IRanges::start(tfmodisco), 103102793L)
    expect_identical(IRanges::end(tfmodisco), 103102842L)
    expect_identical(IRanges::width(tfmodisco), 50L)
    expect_identical(as.character(GenomicRanges::strand(tfmodisco)), "-")
})

test_that("TF-MoDISco pattern placement matches the reference on the minus strand", {
    skip_if_not_installed("Biostrings")
    skip_if_not_installed("universalmotif")
    fixture_dir <- system.file(
        "extdata", "example_repeated_motif",
        package = "GvizRegulatoryTracks"
    )
    sequence <- as.character(Biostrings::readDNAStringSet(
        file.path(fixture_dir, "locus_sequence.fa")
    )[[1]])
    motif_object <- universalmotif::read_meme(
        file.path(fixture_dir, "tfmodisco_pattern.meme")
    )
    if (is.list(motif_object)) motif_object <- motif_object[[1]]
    motif <- methods::slot(motif_object, "motif")
    rownames(motif) <- c("A", "C", "G", "T")
    pwm <- log2((motif + 1e-4) / 0.25)
    score_sequence <- function(x, matrix) {
        bases <- strsplit(x, "", fixed = TRUE)[[1]]
        sum(matrix[cbind(match(bases, rownames(matrix)), seq_along(bases))])
    }
    offset <- 103102793L - 103102780L + 1L
    matched <- substr(sequence, offset, offset + ncol(pwm) - 1L)
    reverse_pwm <- GvizRegulatoryTracks:::.reverse_complement_ppm(motif)
    reverse_pwm <- log2((reverse_pwm + 1e-4) / 0.25)
    starts <- seq_len(nchar(sequence) - ncol(pwm) + 1L)
    plus_scores <- vapply(starts, function(i) {
        score_sequence(substr(sequence, i, i + ncol(pwm) - 1L), pwm)
    }, numeric(1L))
    minus_scores <- vapply(starts, function(i) {
        score_sequence(substr(sequence, i, i + ncol(pwm) - 1L), reverse_pwm)
    }, numeric(1L))

    expect_gt(score_sequence(matched, reverse_pwm), score_sequence(matched, pwm))
    expect_identical(which.max(minus_scores), offset)
    expect_gt(minus_scores[offset], max(plus_scores))
})

test_that("all track encodings render at the repeated-motif locus", {
    skip_if_not_installed("Biostrings")
    skip_if_not_installed("rtracklayer")
    fixture_dir <- system.file(
        "extdata", "example_repeated_motif",
        package = "GvizRegulatoryTracks"
    )
    region <- GenomicRanges::GRanges(
        "chr11", IRanges::IRanges(103102780L, 103102890L)
    )
    sequence <- as.character(Biostrings::readDNAStringSet(
        file.path(fixture_dir, "locus_sequence.fa")
    )[[1]])
    hits <- rtracklayer::import(file.path(fixture_dir, "fimo_hits.bed"))
    hits$relative_fimo_score <- hits$score / 1000
    meme <- file.path(fixture_dir, "fimo_motifs.meme")
    tfmodisco_hits <- rtracklayer::import(
        file.path(fixture_dir, "tfmodisco_seqlet.bed")
    )
    tfmodisco_meme <- file.path(fixture_dir, "tfmodisco_pattern.meme")
    colors <- c(
        A = "#2E8B57", C = "#2878B5", G = "#E69F00", T = "#D64B40"
    )
    motif_palette <- c(`CTCF.H14CORE.0.P.B` = "#7E57C2")

    dynseq <- DynSeqTrack(
        sequence,
        file.path(fixture_dir, "chrombpnet_count_contribution.bw"),
        region = region, genome = "mm10", colors = colors,
        yTicksAt = c(-0.04, 0, 0.04), col.axis = "#333333"
    )
    nucleotide <- MotifLogoTrack(
        tfmodisco_hits, tfmodisco_meme,
        motifColorBy = "nucleotide", colors = colors,
        scoreAesthetic = "none", showAxis = FALSE
    )
    identity <- MotifLogoTrack(
        hits, meme, motifColorBy = "motif_id", motifPalette = motif_palette,
        scoreAesthetic = "none", showAxis = FALSE
    )
    score_fill <- MotifLogoTrack(
        hits, meme, scoreColumn = "relative_fimo_score",
        scoreAesthetic = "fill", scoreLimits = "view", showAxis = FALSE
    )
    combined <- MotifLogoTrack(
        hits, meme, motifColorBy = "motif_id", motifPalette = motif_palette,
        scoreColumn = "relative_fimo_score",
        scoreAesthetic = c("brightness", "opacity", "border"),
        scoreLimits = "view", showAxis = FALSE
    )

    expect_s4_class(dynseq, "CustomTrack")
    expect_s4_class(dynseq, "NumericTrack")
    expect_equal(Gviz::displayPars(dynseq)$yTicksAt, c(-0.04, 0, 0.04))
    expect_identical(
        S4Vectors::mcols(methods::slot(nucleotide, "variables")$hits)$lane,
        1L
    )
    expect_false(methods::slot(nucleotide, "variables")$showScoreLegend)
    expect_false(methods::slot(identity, "variables")$showScoreLegend)
    expect_true(methods::slot(score_fill, "variables")$showScoreLegend)
    expect_true(methods::slot(combined, "variables")$showScoreLegend)

    atac <- Gviz::DataTrack(
        rtracklayer::import(
            file.path(fixture_dir, "atac_seq_signal.bw"), which = region
        ),
        genome = "mm10", chromosome = "chr11", type = "histogram"
    )
    tracks <- list(atac, dynseq, nucleotide, identity, score_fill, combined)
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file, width = 10, height = 10, useDingbats = FALSE)
    expect_no_error(Gviz::plotTracks(
        tracks, chromosome = "chr11", from = 103102780L, to = 103102890L
    ))
    expect_no_error(Gviz::plotTracks(
        tracks, chromosome = "chr11", from = 103102780L, to = 103102890L,
        reverseStrand = TRUE
    ))
    expect_no_error(Gviz::plotTracks(
        list(dynseq, nucleotide), chromosome = "chr11",
        from = 103102200L, to = 103103300L
    ))
    grDevices::dev.off()
    expect_gt(file.info(file)$size, 1000)
})
