test_that("contribution limits are shared and symmetric", {
    expect_equal(contributionYlim(c(-2, 1), c(3, 0), padding = 0), c(-3, 3))
    expect_equal(contributionYlim(c(0, 0), padding = 0), c(-1, 1))
})

test_that("contribution limits parse metadata and validate padding", {
    factor_scores <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(1:3, width = 1),
        score = factor(c("0.2", "0.8", NA))
    )
    character_scores <- factor_scores
    S4Vectors::mcols(character_scores)$score <- c("-0.5", "0.25", NA)

    expect_equal(contributionYlim(factor_scores, padding = 0), c(-0.8, 0.8))
    expect_equal(contributionYlim(character_scores, padding = 0), c(-0.5, 0.5))
    expect_equal(contributionYlim(c(1, NA, Inf), padding = 0), c(-1, 1))

    S4Vectors::mcols(factor_scores)$score <- factor(c("0.2", "bad", NA))
    expect_error(contributionYlim(factor_scores), "non-numeric")
    S4Vectors::mcols(factor_scores)$score <- c("0.2", "Inf", NA)
    expect_error(contributionYlim(factor_scores), "non-finite text")
    expect_error(contributionYlim(1, padding = -0.1), "padding")
    expect_error(contributionYlim(1, padding = Inf), "padding")
})

test_that("numeric contribution data align to one-based coordinates", {
    track <- ScoreSequenceTrack("ACGT", c(1, -2, 0.5, 0),
                                chromosome = "chr1", start = 101)
    expect_s4_class(track, "CustomTrack")
    expect_s4_class(track, "NumericTrack")
    dat <- methods::slot(track, "variables")$data
    expect_equal(dat$position, 101:104)
    expect_equal(dat$base, c("A", "C", "G", "T"))
    expect_equal(dat$score, c(1, -2, 0.5, 0))
    expect_error(
        ScoreSequenceTrack("A", 1, chromosome = "chr1", start = 0),
        "one-based"
    )
    expect_error(
        ScoreSequenceTrack("A", 1, chromosome = "chr1", start = 1.5),
        "one-based"
    )
    expect_error(
        ScoreSequenceTrack(scores = c(1, 2), chromosome = "chr1", start = 10,
                           end = 12),
        "start \\+ length"
    )
    exact_region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10, 11))
    expect_no_error(ScoreSequenceTrack(scores = c(1, 2), region = exact_region))
    expect_error(ScoreSequenceTrack(scores = 1, region = exact_region), "region.*width")
    expect_error(
        ScoreSequenceTrack(scores = c(1, 2), region = exact_region,
                           chromosome = "chr2"),
        "chromosome.*region"
    )
})

test_that("ScoreSequenceTrack and DynSeqTrack are primary compatible constructors", {
    score <- ScoreSequenceTrack(
        "ACGT", c(1, 2, 3, 4), chromosome = "chr1", start = 10
    )
    dynseq <- DynSeqTrack(
        "ACGT", c(1, 2, 3, 4), chromosome = "chr1", start = 10
    )
    expect_s4_class(score, "CustomTrack")
    expect_s4_class(score, "NumericTrack")
    expect_s4_class(dynseq, "CustomTrack")
    expect_s4_class(dynseq, "NumericTrack")
    expect_s4_class(methods::slot(score, "variables")$signalTrack, "DataTrack")
    expect_true(methods::slot(score, "variables")$showBoundary)
    expect_identical(methods::slot(score, "variables")$boundaryColor, "#222222")
    expect_equal(scoreYlim(c(-2, 1), padding = 0), c(-2, 2))
})

test_that("score-track boundaries are configurable and validated", {
    track <- ScoreSequenceTrack(
        "AC", c(-1, 1), chromosome = "chr1", start = 10,
        showBoundary = FALSE, boundaryColor = "grey40", boundaryWidth = 1.2,
        boundaryInset = 0.1
    )
    vars <- methods::slot(track, "variables")
    expect_false(vars$showBoundary)
    expect_identical(vars$boundaryColor, "grey40")
    expect_equal(vars$boundaryWidth, 1.2)
    expect_equal(vars$boundaryInset, 0.1)
    expect_error(
        ScoreSequenceTrack("A", 1, chromosome = "chr1", start = 1,
                           showBoundary = NA),
        "showBoundary"
    )
    expect_error(
        ScoreSequenceTrack("A", 1, chromosome = "chr1", start = 1,
                           boundaryColor = "not-a-color"),
        "boundaryColor"
    )
    expect_error(
        ScoreSequenceTrack("A", 1, chromosome = "chr1", start = 1,
                           boundaryInset = 0.5),
        "boundaryInset"
    )
})

test_that("DynSeq axes use Gviz NumericTrack display parameters", {
    track <- DynSeqTrack(
        "ACGT", c(-1, -0.25, 0.5, 1), chromosome = "chr1", start = 10,
        ylim = c(-1, 1), yTicksAt = c(-1, 0, 1),
        col.axis = "#123456", cex.axis = 0.7
    )
    pars <- Gviz::displayPars(track)

    expect_s4_class(track, "NumericTrack")
    expect_equal(pars$ylim, c(-1, 1))
    expect_equal(pars$yTicksAt, c(-1, 0, 1))
    expect_identical(pars$col.axis, "#123456")
    expect_equal(pars$cex.axis, 0.7)
    expect_true(pars$showAxis)

    Gviz::displayPars(track) <- list(col.axis = "#654321", cex.axis = 0.9)
    expect_identical(Gviz::displayPars(track)$col.axis, "#654321")
    expect_equal(Gviz::displayPars(track)$cex.axis, 0.9)

    hidden <- DynSeqTrack(
        "AC", c(-1, 1), chromosome = "chr1", start = 10,
        showAxis = FALSE
    )
    expect_false(Gviz::displayPars(hidden)$showAxis)
})

test_that("reference sequence is retrieved through Gviz SequenceTrack inputs", {
    reference <- Biostrings::DNAStringSet(c(chr1 = paste(rep("ACGT", 50), collapse = "")))
    track <- ScoreSequenceTrack(
        scores = c(1, 2, 3, 4), reference = reference,
        genome = "test", chromosome = "chr1", start = 101
    )
    dat <- methods::slot(track, "variables")$data
    expect_equal(dat$position, 101:104)
    expect_equal(dat$base, c("A", "C", "G", "T"))

    sequence_track <- Gviz::SequenceTrack(reference, chromosome = "chr1", genome = "test")
    track2 <- ScoreSequenceTrack(
        scores = c(1, 2, 3, 4), reference = sequence_track,
        chromosome = "chr1", start = 105
    )
    expect_equal(methods::slot(track2, "variables")$data$base,
                 c("A", "C", "G", "T"))
})

test_that("GRanges scores expand runs and reject gaps or overlaps", {
    gr <- GenomicRanges::GRanges("chr2", IRanges::IRanges(c(10, 12), c(11, 13)), score = c(1, -1))
    track <- ScoreSequenceTrack("ACGT", gr)
    expect_equal(methods::slot(track, "variables")$data$score, c(1, 1, -1, -1))

    gap <- GenomicRanges::GRanges("chr2", IRanges::IRanges(c(10, 13), c(11, 13)), score = c(1, 2))
    expect_error(ScoreSequenceTrack("ACGT", gap), "does not cover")
    expect_no_error(ScoreSequenceTrack("ACGT", gap, missing = "zero"))

    overlap <- GenomicRanges::GRanges("chr2", IRanges::IRanges(c(10, 11), c(12, 13)), score = c(1, 2))
    expect_error(ScoreSequenceTrack("ACGT", overlap), "overlapping")

    zero_start <- GenomicRanges::GRanges("chr2", IRanges::IRanges(0, 0), score = 1)
    expect_error(ScoreSequenceTrack("A", zero_start), "one-based")
})

test_that("GRanges score metadata are parsed strictly", {
    factor_scores <- GenomicRanges::GRanges(
        "chr2", IRanges::IRanges(10:12, width = 1),
        score = factor(c("1.5", "-0.5", "0"))
    )
    track <- ScoreSequenceTrack("ACG", factor_scores)
    expect_equal(methods::slot(track, "variables")$data$score, c(1.5, -0.5, 0))

    character_scores <- factor_scores
    S4Vectors::mcols(character_scores)$score <- c("2", "0.25", "-1")
    track <- ScoreSequenceTrack("ACG", character_scores)
    expect_equal(methods::slot(track, "variables")$data$score, c(2, 0.25, -1))

    S4Vectors::mcols(character_scores)$score <- c("2", "bad", "-1")
    expect_error(ScoreSequenceTrack("ACG", character_scores), "non-numeric")
    S4Vectors::mcols(character_scores)$score <- c(2, NA, -1)
    expect_error(ScoreSequenceTrack("ACG", character_scores), "finite, non-missing")
    S4Vectors::mcols(character_scores)$score <- c(2, Inf, -1)
    expect_error(ScoreSequenceTrack("ACG", character_scores), "finite, non-missing")
})

test_that("score ranges use range-native overlap and width validation", {
    huge_overlap <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 50000000L), c(100000000L, 150000000L)),
        score = c(1, 2)
    )
    expect_error(
        GvizRegulatoryTracks:::.expand_granges_scores(
            huge_overlap, "score", "error"
        ),
        "overlapping"
    )

    target_scores <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 200000000L), c(100000000L, 300000000L)),
        score = c(1, 2)
    )
    target <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10, 20))
    expanded <- GvizRegulatoryTracks:::.expand_granges_scores(
        target_scores, "score", "error", targetRange = target
    )
    expect_equal(expanded$positions, 10:20)
    expect_equal(expanded$scores, rep(1, 11))

    zero_width <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(start = 10, width = 0), score = 1
    )
    expect_error(ScoreSequenceTrack(scores = zero_width), "positive-width")
})

test_that("explicit GRanges coordinates are validated without truncation", {
    scores <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(10:11, width = 1), score = c(1, 2)
    )
    expect_no_error(ScoreSequenceTrack(scores = scores, start = 10))
    expect_error(ScoreSequenceTrack(scores = scores, start = 10.9), "one-based")
})

test_that("score-track display parameters fail early when invalid", {
    make_track <- function(...) {
        ScoreSequenceTrack("A", 1, chromosome = "chr1", start = 1, ...)
    }
    expect_error(make_track(maxBases = 0), "maxBases")
    expect_error(make_track(maxBases = 1.5), "maxBases")
    expect_error(make_track(maxBases = NA_real_), "maxBases")
    expect_error(make_track(maxBases = .Machine$integer.max + 1), "maxBases")
    expect_error(make_track(showAxis = NA), "showAxis")
    expect_error(make_track(showAxis = 1), "showAxis")
    expect_error(make_track(negativeAlpha = -0.1), "negativeAlpha")
    expect_error(make_track(negativeAlpha = 1.1), "negativeAlpha")
    expect_error(make_track(signalWindowSize = 0), "signalWindowSize")
    expect_error(make_track(baselineColor = "not-a-color"), "baselineColor")
    expect_error(
        make_track(colors = c(A = "bad", C = "blue", G = "gold", T = "red")),
        "colors"
    )
})

test_that("normalization and complements are correct", {
    track <- ScoreSequenceTrack("ACGTN", c(-2, -1, 0, 1, 2),
                                chromosome = "chr1", start = 1,
                                normalization = "maxabs")
    expect_equal(methods::slot(track, "variables")$data$score,
                 c(-1, -0.5, 0, 0.5, 1))
    expect_equal(GvizRegulatoryTracks:::.complement_bases(c("A", "C", "G", "T", "N")),
                 c("T", "G", "C", "A", "N"))
})

test_that("automatic contribution display switches from letters to signal", {
    expect_equal(
        GvizRegulatoryTracks:::.contribution_display_mode("auto", 100, 250, TRUE),
        "letters"
    )
    expect_equal(
        GvizRegulatoryTracks:::.contribution_display_mode("auto", 500, 250, TRUE),
        "signal"
    )
    expect_equal(
        GvizRegulatoryTracks:::.contribution_display_mode("auto", 10, 250, FALSE),
        "signal"
    )
    signal_only <- ScoreSequenceTrack(
        scores = c(1, -1), chromosome = "chr1", start = 10
    )
    expect_false(methods::slot(signal_only, "variables")$hasSequence)
    pars <- Gviz::displayPars(methods::slot(signal_only, "variables")$signalTrack)
    expect_equal(pars$aggregation, "mean")
    expect_equal(pars$window, "auto")
    expect_false(pars$collapse)
})

test_that("score files are imported only for the selected region", {
    skip_if_not_installed("rtracklayer")
    bedgraph <- tempfile(fileext = ".bedGraph")
    writeLines(c("chr1\t100\t102\t1.5", "chr1\t102\t104\t-0.5"), bedgraph)
    region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(101, 104))
    track <- ScoreSequenceTrack(scores = bedgraph, region = region)
    dat <- methods::slot(track, "variables")$data
    expect_equal(dat$position, 101:104)
    expect_equal(dat$score, c(1.5, 1.5, -0.5, -0.5))

    custom <- ScoreSequenceTrack(
        scores = bedgraph, region = region,
        importFunction = function(file, selection) {
            GenomicRanges::GRanges(
                "chr1", IRanges::IRanges(101, 104), score = 2
            )
        }
    )
    expect_equal(methods::slot(custom, "variables")$data$score, rep(2, 4))
})

test_that("letter outlines come from a standard ggseqlogo font", {
    glyphs <- GvizRegulatoryTracks:::.load_glyph_font("roboto_medium")
    expect_true(all(c("A", "C", "G", "T", "N") %in% names(glyphs)))
    expect_true(all(vapply(glyphs[c("A", "C", "G", "T")], nrow, integer(1)) > 5L))
})

test_that("real demo FASTA and BigWigs align at all 101 genomic bases", {
    skip_if_not_installed("Biostrings")
    skip_if_not_installed("rtracklayer")
    fixture_dir <- system.file(
        "extdata", "example_multi_motif", package = "GvizRegulatoryTracks"
    )
    expect_true(nzchar(fixture_dir))
    sequence <- as.character(Biostrings::readDNAStringSet(
        file.path(fixture_dir, "locus_sequence.fa")
    )[[1]])
    region <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(176724198L, 176724298L)
    )
    track <- DynSeqTrack(
        sequence = sequence,
        scores = file.path(fixture_dir, "chrombpnet_count_contribution.bw"),
        region = region
    )
    dat <- methods::slot(track, "variables")$data
    expect_equal(nchar(sequence), 101L)
    expect_equal(dat$position, 176724198L:176724298L)
    expect_equal(paste0(dat$base, collapse = ""), sequence)
    expect_true(any(dat$score < 0) && any(dat$score > 0))

    atac <- rtracklayer::import(
        file.path(fixture_dir, "atac_seq_signal.bw"), which = region
    )
    expect_equal(sum(IRanges::width(atac)), 101L)
    expect_false(any(GenomicRanges::countOverlaps(atac, atac,
                                                  ignore.strand = TRUE) > 1L))
})

test_that("contribution tracks draw in Gviz in both orientations", {
    track <- ScoreSequenceTrack("ACGTACGT", seq(-1, 1, length.out = 8),
                                chromosome = "chr1", start = 100)
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file, width = 7, height = 2)
    expect_no_error(Gviz::plotTracks(track, chromosome = "chr1", from = 100, to = 107))
    expect_no_error(Gviz::plotTracks(track, chromosome = "chr1", from = 100, to = 107,
                                     reverseStrand = TRUE))
    expect_no_error(Gviz::plotTracks(track, chromosome = "chr1", from = 1, to = 1000))
    grDevices::dev.off()
    expect_gt(file.info(file)$size, 0)
})
