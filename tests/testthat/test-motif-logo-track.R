test_that("PPMs are normalized and reverse-complemented", {
    ppm <- rbind(A = c(1, 0), C = c(0, 1), G = c(0, 0), T = c(0, 0))
    rc <- GvizRegulatoryTracks:::.reverse_complement_ppm(ppm)
    expect_equal(rownames(rc), c("A", "C", "G", "T"))
    expect_equal(max.col(t(rc)), c(3, 4))

    counts <- ppm * 10
    expect_equal(GvizRegulatoryTracks:::.as_ppm_matrix(counts, "PCM"), ppm)

    oriented <- GvizRegulatoryTracks:::.oriented_motif_ppm
    expect_equal(oriented(ppm, "-", FALSE, TRUE), rc)
    expect_equal(oriented(ppm, "-", FALSE, FALSE), ppm)
    expect_equal(oriented(ppm, "-", TRUE, TRUE), ppm)
    expect_equal(oriented(ppm, "+", TRUE, TRUE), rc)
})

test_that("motif widths and IDs are validated", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hit <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10, 12), strand = "+",
                                  motif_id = "M1")
    expect_s4_class(MotifLogoTrack(hit, list(M1 = motif)), "CustomTrack")
    expect_error(MotifLogoTrack(GenomicRanges::resize(hit, 2), list(M1 = motif)), "widths")
    expect_error(MotifLogoTrack(hit, list(M2 = motif)), "Missing motif")

    zero_width <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(start = 10, width = 0), motif_id = "M1"
    )
    expect_error(MotifLogoTrack(zero_width, list(M1 = motif)), "positive-width")
})

test_that("overlapping motifs receive deterministic lanes", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hits <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 12, 20), width = 3),
                                   motif_id = c("M1", "M1", "M1"))
    track <- MotifLogoTrack(hits, list(M1 = motif))
    lanes <- S4Vectors::mcols(methods::slot(track, "variables")$hits)$lane
    expect_equal(lanes, c(1L, 2L, 1L))
})

test_that("motif lane geometry uses fixed physical heights without overlap", {
    layout_one <- GvizRegulatoryTracks:::.motif_lane_layout(
        1L, trackHeight = 20, motifHeight = 10, laneGap = 0.30,
        showLabels = TRUE
    )
    layout_tall <- GvizRegulatoryTracks:::.motif_lane_layout(
        1L, trackHeight = 60, motifHeight = 10, laneGap = 0.30,
        showLabels = TRUE
    )
    expect_equal(layout_one$motifHeight, 10)
    expect_equal(layout_tall$motifHeight, 10)
    expect_gt(layout_tall$baselines, layout_one$baselines)

    layout_lanes <- GvizRegulatoryTracks:::.motif_lane_layout(
        3L, trackHeight = 60, motifHeight = 10, laneGap = 0.30,
        showLabels = TRUE
    )
    lane_tops <- layout_lanes$baselines + layout_lanes$motifHeight +
        layout_lanes$labelHeight
    expect_true(all(lane_tops[-length(lane_tops)] < layout_lanes$baselines[-1L]))

    layout_short <- GvizRegulatoryTracks:::.motif_lane_layout(
        3L, trackHeight = 18, motifHeight = 10, laneGap = 0.30,
        showLabels = TRUE
    )
    expect_lt(layout_short$scale, 1)
    expect_lt(layout_short$motifHeight, 10)
    short_tops <- layout_short$baselines + layout_short$motifHeight +
        layout_short$labelHeight
    expect_true(all(short_tops[-length(short_tops)] < layout_short$baselines[-1L]))
})

test_that("motif height is stored and validated", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hit <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(10, 12), motif_id = "M1"
    )
    track <- MotifLogoTrack(hit, list(M1 = motif), motifHeight = 8)
    expect_equal(methods::slot(track, "variables")$motifHeight, 8)
    expect_error(MotifLogoTrack(hit, list(M1 = motif), motifHeight = 0),
                 "positive")
})

test_that("motif-track boundaries are configurable and validated", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hit <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(10, 12), motif_id = "M1"
    )
    track <- MotifLogoTrack(
        hit, list(M1 = motif), showBoundary = FALSE,
        boundaryColor = "grey40", boundaryWidth = 1.2, boundaryInset = 0.1
    )
    vars <- methods::slot(track, "variables")
    expect_false(vars$showBoundary)
    expect_identical(vars$boundaryColor, "grey40")
    expect_equal(vars$boundaryWidth, 1.2)
    expect_equal(vars$boundaryInset, 0.1)
    expect_error(
        MotifLogoTrack(hit, list(M1 = motif), boundaryWidth = 0),
        "boundaryWidth"
    )
})

test_that("score legend width expands to contain its title", {
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file, width = 7, height = 2)
    on.exit(grDevices::dev.off(), add = TRUE)
    short <- GvizRegulatoryTracks:::.score_legend_width("Score", 16)
    long <- GvizRegulatoryTracks:::.score_legend_width(
        "Relative FIMO score", 16
    )
    expect_gte(short, 16)
    expect_gt(long, 16)
})

test_that("automatic motif display switches from logos to ranges", {
    expect_equal(GvizRegulatoryTracks:::.motif_display_mode("auto", 100, 1000), "logo")
    expect_equal(GvizRegulatoryTracks:::.motif_display_mode("auto", 2000, 1000), "ranges")
})

test_that("motif scores map independently to fill, border, and opacity", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 20), width = 3),
        motif_id = c("M1", "M1"), relative_score = c(0.2, 0.8)
    )
    fill <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "fill", scoreLimits = c(0, 1)
    )
    fill_style <- methods::slot(fill, "variables")$hitStyles
    fill_vars <- methods::slot(fill, "variables")
    expect_false(anyNA(fill_style$fill))
    expect_false(identical(fill_style$fill[1], fill_style$fill[2]))
    expect_true(all(is.na(fill_style$border)))
    expect_true(fill_vars$showScoreLegend)
    expect_equal(fill_vars$scoreLimits, c(0, 1))

    combined <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = c("border", "opacity"), scoreLimits = c(0, 1),
        scoreOpacity = c(0.3, 1)
    )
    combined_style <- methods::slot(combined, "variables")$hitStyles
    expect_false(anyNA(combined_style$border))
    expect_lt(combined_style$alpha[1], combined_style$alpha[2])
    expect_error(
        MotifLogoTrack(hits, list(M1 = motif), scoreAesthetic = "fill"),
        "scoreColumn"
    )
})

test_that("motif score metadata are parsed strictly", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 20), width = 3),
        motif_id = "M1", relative_score = factor(c("0.2", "0.8"))
    )
    factor_track <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "fill", scoreLimits = c(0, 1)
    )
    expect_equal(
        methods::slot(factor_track, "variables")$hitStyles$score,
        c(0.2, 0.8)
    )

    S4Vectors::mcols(hits)$relative_score <- c("0.1", "0.9")
    character_track <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "fill", scoreLimits = c(0, 1)
    )
    expect_equal(
        methods::slot(character_track, "variables")$hitStyles$score,
        c(0.1, 0.9)
    )

    S4Vectors::mcols(hits)$relative_score <- c("0.1", "bad")
    expect_error(
        MotifLogoTrack(hits, list(M1 = motif), scoreColumn = "relative_score",
                       scoreAesthetic = "fill"),
        "non-numeric"
    )
    S4Vectors::mcols(hits)$relative_score <- c(0.1, NA)
    expect_error(
        MotifLogoTrack(hits, list(M1 = motif), scoreColumn = "relative_score",
                       scoreAesthetic = "fill"),
        "finite, non-missing"
    )
    S4Vectors::mcols(hits)$relative_score <- c(0.1, Inf)
    expect_error(
        MotifLogoTrack(hits, list(M1 = motif), scoreColumn = "relative_score",
                       scoreAesthetic = "fill"),
        "finite, non-missing"
    )
})

test_that("motif-track display parameters fail early when invalid", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hit <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(10, 12), motif_id = "M1", score = 0.5
    )
    make_track <- function(...) MotifLogoTrack(hit, list(M1 = motif), ...)

    expect_error(make_track(maxBases = 0), "maxBases")
    expect_error(make_track(maxBases = 1.5), "maxBases")
    expect_error(make_track(maxRangeLabels = -1), "maxRangeLabels")
    expect_error(make_track(maxRangeLabels = 1.5), "maxRangeLabels")
    expect_error(make_track(rangeAlpha = -0.1), "rangeAlpha")
    expect_error(make_track(rangeAlpha = 1.1), "rangeAlpha")
    expect_error(make_track(scoreBorderWidth = -0.1), "scoreBorderWidth")
    expect_error(make_track(showLabels = NA), "showLabels")
    expect_error(make_track(showStrand = 1), "showStrand")
    expect_error(make_track(showScoreLegend = NA), "showScoreLegend")
    expect_error(make_track(scoreLegendDigits = 1.5), "scoreLegendDigits")
    expect_error(make_track(rangeFill = "not-a-color"), "rangeFill")
    expect_error(make_track(scorePalette = c("white", "not-a-color")),
                 "scorePalette")
    expect_error(make_track(scoreLimits = c(1, 0)), "scoreLimits")
    expect_error(
        make_track(colors = c(A = "bad", C = "blue", G = "gold", T = "red")),
        "colors"
    )
})

test_that("nucleotide logo colors can be shared without score encoding", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hit <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(10, 12), motif_id = "M1"
    )
    palette <- c(
        A = "#2E8B57", C = "#2878B5", G = "#E69F00", T = "#D64B40"
    )
    track <- MotifLogoTrack(
        hit, list(M1 = motif), motifColorBy = "nucleotide",
        colors = palette, scoreAesthetic = "none"
    )
    vars <- methods::slot(track, "variables")

    expect_identical(unname(vars$colors[names(palette)]), unname(palette))
    expect_true(all(is.na(vars$hitStyles$fill)))
    expect_false(vars$showScoreLegend)
})

test_that("score legends represent mapped limits and can be disabled", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 20), width = 3), strand = c("+", "-"),
        motif_id = "M1", relative_score = c(0.2, 0.8)
    )
    track <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "opacity", scoreOpacity = c(0.2, 1),
        showLabels = TRUE, showStrand = TRUE
    )
    vars <- methods::slot(track, "variables")
    legend_styles <- GvizRegulatoryTracks:::.score_legend_styles(vars, n = 5L)
    expect_true(vars$showScoreLegend)
    expect_equal(vars$scoreLimitMode, "view")
    expect_null(vars$scoreLimits)
    expect_equal(vars$globalScoreLimits, c(0.2, 0.8))
    expect_equal(vars$scoreLegendHeight, 23)
    expect_equal(range(legend_styles$alpha), c(0.2, 1))
    expect_length(unique(legend_styles$fill), 1L)
    expect_equal(
        GvizRegulatoryTracks:::.motif_hit_labels(vars$hits, TRUE),
        c("M1 (+)", "M1 (-)")
    )

    hidden <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "fill", showScoreLegend = FALSE
    )
    expect_false(methods::slot(hidden, "variables")$showScoreLegend)

    unscored <- MotifLogoTrack(hits, list(M1 = motif), scoreAesthetic = "none")
    expect_false(methods::slot(unscored, "variables")$showScoreLegend)
})

test_that("view, global, and fixed score limits are resolved predictably", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 20, 30), width = 3),
        motif_id = "M1", relative_score = c(0.2, 0.5, 0.8)
    )
    view_track <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "fill"
    )
    view_vars <- methods::slot(view_track, "variables")
    visible <- GvizRegulatoryTracks:::.visible_motif_styles(
        view_vars, c(FALSE, TRUE, TRUE)
    )
    expect_equal(visible$limits, c(0.5, 0.8))

    one <- GvizRegulatoryTracks:::.visible_motif_styles(
        view_vars, c(TRUE, FALSE, FALSE)
    )
    expect_equal(one$limits, c(0.15, 0.25))

    global_track <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "fill", scoreLimits = "global"
    )
    global_vars <- methods::slot(global_track, "variables")
    global <- GvizRegulatoryTracks:::.visible_motif_styles(
        global_vars, c(FALSE, TRUE, TRUE)
    )
    expect_equal(global$limits, c(0.2, 0.8))

    fixed_track <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "fill", scoreLimits = c(0, 1)
    )
    fixed_vars <- methods::slot(fixed_track, "variables")
    fixed <- GvizRegulatoryTracks:::.visible_motif_styles(
        fixed_vars, c(FALSE, TRUE, TRUE)
    )
    expect_equal(fixed$limits, c(0, 1))

    factor_hits <- hits
    S4Vectors::mcols(factor_hits)$relative_score <- factor(
        c("0.2", "0.5", "0.8")
    )
    factor_track <- MotifLogoTrack(
        factor_hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "fill", scoreLimits = "view"
    )
    factor_vars <- methods::slot(factor_track, "variables")
    factor_view <- GvizRegulatoryTracks:::.visible_motif_styles(
        factor_vars, c(TRUE, FALSE, TRUE)
    )
    expect_equal(factor_view$limits, c(0.2, 0.8))
    expect_equal(factor_view$styles$scaledScore, c(0, 1))
})

test_that("multiple motif IDs receive stable categorical colors", {
    m1 <- matrix(0.25, nrow = 4, ncol = 3,
                 dimnames = list(c("A", "C", "G", "T"), NULL))
    m2 <- matrix(0.25, nrow = 4, ncol = 4,
                 dimnames = list(c("A", "C", "G", "T"), NULL))
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 20, 30), width = c(3, 4, 3)),
        motif_id = c("M1", "M2", "M1"), relative_score = c(0.2, 0.5, 0.8)
    )
    palette <- c(M1 = "#0072B2", M2 = "#D55E00")
    track <- MotifLogoTrack(
        hits, list(M1 = m1, M2 = m2), motifColorBy = "motif_id",
        motifPalette = palette, showLabels = TRUE
    )
    vars <- methods::slot(track, "variables")
    expect_equal(vars$motifColors, palette)
    expect_equal(vars$hitStyles$fill, unname(palette[c("M1", "M2", "M1")]))
    expect_identical(vars$hitStyles$fill[1], vars$hitStyles$fill[3])

    automatic <- MotifLogoTrack(
        hits[c(2, 1, 3)], list(M1 = m1, M2 = m2), motifColorBy = "motif_id"
    )
    expect_equal(names(methods::slot(automatic, "variables")$motifColors), c("M1", "M2"))
})

test_that("motif identity fill combines with score brightness opacity and border", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 20), width = 3),
        motif_id = c("M1", "M1"), relative_score = c(0.2, 0.8)
    )
    track <- MotifLogoTrack(
        hits, list(M1 = motif), motifColorBy = "motif_id",
        motifPalette = c(M1 = "#0072B2"), scoreColumn = "relative_score",
        scoreAesthetic = c("brightness", "opacity", "border"),
        scoreLimits = c(0, 1), scoreColor = "black"
    )
    style <- methods::slot(track, "variables")$hitStyles
    expect_identical(style$fill[1], style$fill[2])
    expect_gt(style$brightness[1], style$brightness[2])
    expect_lt(style$alpha[1], style$alpha[2])
    expect_false(identical(style$border[1], style$border[2]))
    expect_error(
        MotifLogoTrack(
            hits, list(M1 = motif), motifColorBy = "motif_id",
            scoreColumn = "relative_score", scoreAesthetic = "fill"
        ),
        "mutually|cannot be combined"
    )
})

test_that("overview ranges inherit motif score aesthetics", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 20), width = 3),
        motif_id = c("M1", "M1"), relative_score = c(0, 1)
    )
    track <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = c("brightness", "opacity", "border"),
        scoreLimits = c(0, 1), display = "ranges",
        scoreBrightness = c(0.97, 0), scoreOpacity = c(0.18, 1)
    )
    vars <- methods::slot(track, "variables")
    range_styles <- GvizRegulatoryTracks:::.resolve_motif_range_styles(
        vars$hitStyles, vars
    )
    expect_false(identical(range_styles$fill[1], range_styles$fill[2]))
    expect_lt(range_styles$alpha[1], range_styles$alpha[2])
    expect_false(identical(range_styles$border[1], range_styles$border[2]))
})

test_that("score legend reserves space for multiline titles", {
    one_line <- GvizRegulatoryTracks:::.score_legend_layout("PWM score")
    two_line <- GvizRegulatoryTracks:::.score_legend_layout("Relative\nPWM score")
    expect_equal(two_line$title_just, c("center", "top"))
    expect_lt(two_line$gradient_max, one_line$gradient_max)
    expect_gt(two_line$gradient_max, two_line$gradient_min)
})

test_that("a single score color generates a light-to-base fill gradient", {
    motif <- matrix(0.25, nrow = 4, ncol = 3,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 20), width = 3),
        motif_id = "M1", relative_score = c(0, 1)
    )
    track <- MotifLogoTrack(
        hits, list(M1 = motif), scoreColumn = "relative_score",
        scoreAesthetic = "fill", scoreColor = "#111111",
        scoreBrightness = c(0.85, 0)
    )
    fill <- methods::slot(track, "variables")$hitStyles$fill
    expect_gt(mean(grDevices::col2rgb(fill[1])), mean(grDevices::col2rgb(fill[2])))
    expect_equal(toupper(fill[2]), "#111111")
})

test_that("BED hits and MEME motif files are accepted directly", {
    skip_if_not_installed("rtracklayer")
    skip_if_not_installed("universalmotif")
    motif <- universalmotif::create_motif("ACG", name = "M1")
    meme <- tempfile(fileext = ".meme")
    universalmotif::write_meme(list(motif), meme)
    bed <- tempfile(fileext = ".bed")
    writeLines("chr1\t9\t12\tM1\t10\t-", bed)
    track <- MotifLogoTrack(bed, meme)
    vars <- methods::slot(track, "variables")
    expect_equal(as.character(S4Vectors::mcols(vars$hits)$motif_id), "M1")
    expect_equal(IRanges::start(vars$hits), 10L)
    expect_equal(IRanges::end(vars$hits), 12L)

    one_base_motif <- universalmotif::create_motif("A", name = "M0")
    one_base_meme <- tempfile(fileext = ".meme")
    universalmotif::write_meme(list(one_base_motif), one_base_meme)
    one_base_bed <- tempfile(fileext = ".bed")
    writeLines("chr1\t0\t1\tM0\t10\t+", one_base_bed)
    boundary <- MotifLogoTrack(one_base_bed, one_base_meme)
    boundary_hits <- methods::slot(boundary, "variables")$hits
    expect_equal(IRanges::start(boundary_hits), 1L)
    expect_equal(IRanges::end(boundary_hits), 1L)
    expect_equal(IRanges::width(boundary_hits), 1L)
})

test_that("in-memory motif hits require positive one-based coordinates", {
    motif <- matrix(0.25, nrow = 4, ncol = 1,
                    dimnames = list(c("A", "C", "G", "T"), NULL))
    zero_start <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(0, 0), motif_id = "M1"
    )
    expect_error(MotifLogoTrack(zero_start, list(M1 = motif)), "one-based")
})

test_that("real multi-motif fixtures preserve FIMO coordinates, strands, and scores", {
    skip_if_not_installed("rtracklayer")
    skip_if_not_installed("universalmotif")
    fixture_dir <- system.file(
        "extdata", "example_multi_motif", package = "GvizRegulatoryTracks"
    )
    expect_true(nzchar(fixture_dir))
    bed <- file.path(fixture_dir, "motif_hits.bed")
    meme <- file.path(fixture_dir, "motifs.meme")
    hits <- rtracklayer::import(bed)
    hits$relative_fimo_score <- hits$score / 1000
    ids <- as.character(hits$name)
    expect_equal(IRanges::start(hits),
                 c(176724244L, 176724245L, 176724245L, 176724250L, 176724221L))
    expect_equal(IRanges::end(hits),
                 c(176724254L, 176724261L, 176724261L, 176724261L, 176724231L))
    expect_equal(IRanges::start(hits),
                 c(176724243L, 176724244L, 176724244L, 176724249L,
                   176724220L) + 1L)
    prepared <- MotifLogoTrack(hits, meme)
    motifs <- methods::slot(prepared, "variables")$ppms
    expect_equal(
        unname(IRanges::width(hits)),
        unname(vapply(ids, function(id) ncol(motifs[[id]]), integer(1)))
    )
    expect_setequal(as.character(GenomicRanges::strand(hits)), c("+", "-"))
    expect_equal(range(hits$relative_fimo_score), c(0.944, 0.971))
    expect_true(all(vapply(motifs, function(x) {
        isTRUE(all.equal(unname(colSums(x)), rep(1, ncol(x)), tolerance = 1e-8))
    }, logical(1))))

    palette <- c(
        `PO5F1.H14CORE.1.P.B` = "#0072B2",
        `PO5F1.H14CORE.0.P.B` = "#CC79A7",
        `SOX2.H14CORE.0.P.B` = "#D55E00",
        `NANOG.H14CORE.0.P.B` = "#009E73",
        `ERR1.H14CORE.0.PSM.A` = "#E69F00"
    )
    track <- MotifLogoTrack(
        hits, meme, motifColorBy = "motif_id", motifPalette = palette,
        scoreColumn = "relative_fimo_score",
        scoreAesthetic = c("brightness", "opacity", "border"),
        scoreLimits = "view", showLabels = TRUE, showStrand = TRUE
    )
    vars <- methods::slot(track, "variables")
    expect_equal(length(unique(vars$hitStyles$fill)), 5L)
    expect_true(all(diff(vars$hitStyles$scaledScore[
        order(hits$relative_fimo_score)
    ]) > 0))
    expect_s4_class(track, "CustomTrack")
})

test_that("motif tracks draw in Gviz in both orientations", {
    motif <- rbind(A = c(0.9, 0.1, 0.1, 0.1),
                   C = c(0.05, 0.8, 0.1, 0.1),
                   G = c(0.03, 0.05, 0.7, 0.1),
                   T = c(0.02, 0.05, 0.1, 0.7))
    hit <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 103), strand = "-",
                                  motif_id = "M1")
    track <- MotifLogoTrack(hit, list(M1 = motif), showLabels = TRUE)
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file, width = 7, height = 2)
    expect_no_error(Gviz::plotTracks(track, chromosome = "chr1", from = 95, to = 108))
    expect_no_error(Gviz::plotTracks(track, chromosome = "chr1", from = 95, to = 108,
                                     reverseStrand = TRUE))
    expect_no_error(Gviz::plotTracks(track, chromosome = "chr1", from = 1, to = 2500))
    grDevices::dev.off()
    expect_gt(file.info(file)$size, 0)
})

test_that("motif tracks support an independent display-label column", {
    motif <- rbind(
        A = c(0.8, 0.1), C = c(0.1, 0.1),
        G = c(0.05, 0.7), T = c(0.05, 0.1)
    )
    hit <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(10, 11), strand = "+",
        motif_id = "pattern_8", tomtom_label = "SOX2 / p8"
    )
    track <- MotifLogoTrack(
        hit, list(pattern_8 = motif), labelColumn = "tomtom_label",
        showLabels = TRUE
    )
    prepared_hits <- methods::slot(track, "variables")$hits
    expect_equal(S4Vectors::mcols(prepared_hits)$label, "SOX2 / p8")
    expect_equal(
        GvizRegulatoryTracks:::.motif_hit_labels(prepared_hits, TRUE),
        "SOX2 / p8 (+)"
    )
    expect_error(
        MotifLogoTrack(hit, list(pattern_8 = motif), labelColumn = "missing"),
        "lacks label column"
    )
    layout <- GvizRegulatoryTracks:::.motif_label_layout(
        c(hit, GenomicRanges::shift(hit, 40), GenomicRanges::shift(hit, 88)),
        c(1, 100)
    )
    expect_equal(layout$just, c("left", "center", "right"))
    expect_true(all(layout$x > 1 & layout$x < 100))
    reverse_layout <- GvizRegulatoryTracks:::.motif_label_layout(
        c(hit, GenomicRanges::shift(hit, 40), GenomicRanges::shift(hit, 88)),
        c(100, 1)
    )
    expect_equal(reverse_layout$just, c("right", "center", "left"))
    expect_true(all(reverse_layout$x > 1 & reverse_layout$x < 100))
})
