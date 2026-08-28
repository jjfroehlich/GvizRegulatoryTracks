test_that("directional geometry uses a fixed, full-height tip", {
    long_widths <- GvizRegulatoryTracks:::.directional_tip_width(
        c(100, 500), 12
    )
    expect_equal(long_widths, c(12, 12))
    expect_equal(
        GvizRegulatoryTracks:::.directional_tip_width(10, 12), 4.5
    )

    plus <- GvizRegulatoryTracks:::.directional_polygon(
        10, 50, "+", 1, 0.6, 8
    )
    minus <- GvizRegulatoryTracks:::.directional_polygon(
        10, 50, "-", 1, 0.6, 8
    )
    unstranded <- GvizRegulatoryTracks:::.directional_polygon(
        10, 50, "*", 1, 0.6, 8
    )
    expect_equal(range(plus$y), c(0.7, 1.3))
    expect_equal(range(minus$y), c(0.7, 1.3))
    expect_equal(plus$x[3], 50.5)
    expect_equal(minus$x[3], 9.5)
    expect_length(unstranded$x, 4)

    one_base <- GvizRegulatoryTracks:::.directional_polygon(
        10, 10, "+", 1, 0.6, 0.45
    )
    one_base_rectangle <- GvizRegulatoryTracks:::.directional_rectangle(
        10, 10, 1, 0.6
    )
    expect_equal(diff(range(one_base$x)), 1)
    expect_equal(diff(range(one_base_rectangle$x)), 1)
})

test_that("directional tracks assign stable lanes and feature colors", {
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 15, 40), width = 10),
        strand = c("+", "-", "*"),
        feature = c("A", "B", "A"),
        label = c("one", "two", "three")
    )
    track <- DirectionalRangeTrack(
        hits, featureColumn = "feature", labelColumn = "label",
        fill = c(A = "#112233", B = "#445566"), showLabels = FALSE
    )
    vars <- methods::slot(track, "variables")
    expect_equal(S4Vectors::mcols(vars$hits)$lane, c(1L, 2L, 1L))
    expect_equal(vars$fills, c("#112233", "#445566", "#112233"))
    expect_false(vars$showLabels)

    rectangular <- DirectionalRangeTrack(
        hits, geometry = "rectangle", directionIndicator = "arrow",
        labelPosition = "above"
    )
    rectangular_vars <- methods::slot(rectangular, "variables")
    expect_equal(rectangular_vars$geometry, "rectangle")
    expect_equal(rectangular_vars$directionIndicator, "arrow")
    expect_equal(rectangular_vars$labelPosition, "above")
    rectangle <- GvizRegulatoryTracks:::.directional_rectangle(
        10, 50, 1, 0.6
    )
    expect_equal(rectangle$x, c(9.5, 50.5, 50.5, 9.5))
    expect_equal(range(rectangle$y), c(0.7, 1.3))

    dense <- DirectionalRangeTrack(hits, stacking = "dense")
    expect_equal(
        S4Vectors::mcols(methods::slot(dense, "variables")$hits)$lane,
        rep(1L, 3)
    )
})

test_that("directional track validates inputs and mappings", {
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(10, width = 10), feature = "A"
    )
    expect_error(DirectionalRangeTrack(GenomicRanges::GRanges()), "non-empty")
    expect_error(
        DirectionalRangeTrack(hits, featureColumn = "missing"),
        "lacks feature"
    )
    expect_error(
        DirectionalRangeTrack(hits, featureColumn = c("feature", "other")),
        "one non-empty column name"
    )
    expect_error(
        DirectionalRangeTrack(hits, labelColumn = ""),
        "one non-empty column name"
    )
    bad_labels <- hits
    S4Vectors::mcols(bad_labels)$label <- NA_character_
    expect_error(
        DirectionalRangeTrack(bad_labels, labelColumn = "label"),
        "labels must be non-missing"
    )
    expect_error(
        DirectionalRangeTrack(
            GenomicRanges::GRanges(
                c("chr1", "chr2"),
                IRanges::IRanges(c(10, 20), width = 10), feature = "A"
            )
        ),
        "one chromosome"
    )
    expect_error(
        DirectionalRangeTrack(
            c(hits, GenomicRanges::shift(hits, 20)),
            fill = c(B = "#112233", C = "#445566")
        ),
        "lacks colors"
    )
    expect_error(
        DirectionalRangeTrack(
            c(hits, GenomicRanges::shift(hits, 20)),
            fill = c(A = "#112233", A = "#445566")
        ),
        "unique, non-empty names"
    )
})

test_that("directional glyphs and labels follow the displayed axis", {
    expect_equal(
        GvizRegulatoryTracks:::.directional_glyphs(c("+", "-", "*")),
        c(">", "<", "")
    )
    expect_equal(
        GvizRegulatoryTracks:::.directional_glyphs(
            c("+", "-", "*"), reverseAxis = TRUE
        ),
        c("<", ">", "")
    )
    forward <- GvizRegulatoryTracks:::.viewport_label_layout(
        c(1, 50, 100), c(1, 100)
    )
    reverse <- GvizRegulatoryTracks:::.viewport_label_layout(
        c(1, 50, 100), c(100, 1)
    )
    expect_equal(forward$just, c("left", "center", "right"))
    expect_equal(reverse$just, c("right", "center", "left"))
    expect_true(all(reverse$x > 1 & reverse$x < 100))
})

test_that("directional tracks draw pointed and rectangular ranges both ways", {
    hits <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(10, 25, 40), width = c(1, 10, 15)),
        strand = c("+", "-", "*"), feature = c("A", "B", "A"),
        label = c("one", "two", "three")
    )
    pointed <- DirectionalRangeTrack(
        hits, labelColumn = "label", showLabels = TRUE,
        directionIndicator = "arrow"
    )
    rectangular <- DirectionalRangeTrack(
        hits, geometry = "rectangle", directionIndicator = "arrow"
    )
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file, width = 7, height = 3)
    on.exit({
        if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    }, add = TRUE)
    expect_no_error(Gviz::plotTracks(
        list(pointed, rectangular), chromosome = "chr1", from = 1, to = 60
    ))
    expect_no_error(Gviz::plotTracks(
        list(pointed, rectangular), chromosome = "chr1", from = 1, to = 60,
        reverseStrand = TRUE
    ))
    grDevices::dev.off()
    expect_gt(file.info(file)$size, 0)
})

test_that("DirectionalRangeTrack is exported", {
    expect_true(
        "DirectionalRangeTrack" %in% getNamespaceExports("GvizRegulatoryTracks")
    )
})
