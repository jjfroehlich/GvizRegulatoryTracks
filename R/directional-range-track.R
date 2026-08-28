.assign_directional_lanes <- function(hits, stacking) {
    if (identical(stacking, "dense")) return(rep.int(1L, length(hits)))
    ord <- order(
        IRanges::start(hits), IRanges::end(hits),
        as.character(GenomicRanges::strand(hits)), seq_along(hits)
    )
    lanes <- integer(length(hits))
    lane_ends <- integer()
    for (idx in ord) {
        available <- which(IRanges::start(hits)[idx] > lane_ends)
        lane <- if (length(available)) available[1L] else length(lane_ends) + 1L
        lanes[idx] <- lane
        if (lane > length(lane_ends)) {
            lane_ends <- c(lane_ends, IRanges::end(hits)[idx])
        } else {
            lane_ends[lane] <- IRanges::end(hits)[idx]
        }
    }
    lanes
}

.resolve_directional_fills <- function(features, fill) {
    features <- as.character(features)
    feature_levels <- unique(features)
    fill <- .validate_color_vector(fill, "fill")
    if (length(fill) == 1L) {
        return(rep.int(unname(fill), length(features)))
    }
    if (is.null(names(fill))) {
        if (length(fill) < length(feature_levels)) {
            stop("Unnamed `fill` needs at least one color per feature.", call. = FALSE)
        }
        fill <- fill[seq_along(feature_levels)]
        names(fill) <- feature_levels[seq_along(fill)]
    } else if (anyNA(names(fill)) || any(!nzchar(names(fill))) ||
               anyDuplicated(names(fill))) {
        stop("Named `fill` colors must have unique, non-empty names.",
             call. = FALSE)
    }
    missing_features <- setdiff(feature_levels, names(fill))
    if (length(missing_features)) {
        stop("`fill` lacks colors for features: ",
             paste(missing_features, collapse = ", "), call. = FALSE)
    }
    unname(fill[features])
}

.directional_tip_width <- function(rangeWidth, fixedWidth) {
    pmin(as.numeric(fixedWidth), pmax(0, as.numeric(rangeWidth) * 0.45))
}

.directional_polygon <- function(start, end, strand, lane, height, tipWidth) {
    start <- as.numeric(start) - 0.5
    end <- as.numeric(end) + 0.5
    midpoint <- lane
    lower <- midpoint - height / 2
    upper <- midpoint + height / 2
    strand <- as.character(strand)
    if (strand == "+") {
        body_end <- max(start, end - tipWidth)
        return(list(
            x = c(start, body_end, end, body_end, start),
            y = c(lower, lower, midpoint, upper, upper)
        ))
    }
    if (strand == "-") {
        body_start <- min(end, start + tipWidth)
        return(list(
            x = c(end, body_start, start, body_start, end),
            y = c(lower, lower, midpoint, upper, upper)
        ))
    }
    list(
        x = c(start, end, end, start),
        y = c(lower, lower, upper, upper)
    )
}

.directional_rectangle <- function(start, end, lane, height) {
    start <- as.numeric(start) - 0.5
    end <- as.numeric(end) + 0.5
    lower <- lane - height / 2
    upper <- lane + height / 2
    list(
        x = c(start, end, end, start),
        y = c(lower, lower, upper, upper)
    )
}

.directional_glyphs <- function(strand, reverseAxis = FALSE) {
    strand <- as.character(strand)
    plus <- if (isTRUE(reverseAxis)) "<" else ">"
    minus <- if (isTRUE(reverseAxis)) ">" else "<"
    ifelse(strand == "+", plus, ifelse(strand == "-", minus, ""))
}

.draw_directional_range_track <- function(object) {
    vars <- methods::slot(object, "variables")
    xscale <- grid::current.viewport()$xscale
    visible <- sort(as.numeric(xscale))
    keep <- IRanges::end(vars$hits) >= visible[1L] &
        IRanges::start(vars$hits) <= visible[2L]
    hits <- vars$hits[keep]
    fills <- vars$fills[keep]
    lanes <- S4Vectors::mcols(hits)$lane
    nlanes <- max(1L, if (length(lanes)) max(lanes) else 1L)
    grid::pushViewport(grid::viewport(
        xscale = xscale, yscale = c(0.4, nlanes + 0.6), clip = "on"
    ))
    on.exit(grid::popViewport(), add = TRUE)
    if (vars$showBoundary) {
        .draw_track_left_boundary(
            vars$boundaryColor, vars$boundaryWidth, vars$boundaryInset
        )
    }
    if (!length(hits)) return(invisible(NULL))
    fixed_width_native <- abs(grid::convertWidth(
        grid::unit(vars$triangleWidth, "mm"), "native", valueOnly = TRUE
    ))
    for (i in seq_along(hits)) {
        tip_width <- .directional_tip_width(
            IRanges::width(hits)[i], fixed_width_native
        )
        polygon <- if (identical(vars$geometry, "rectangle")) {
            .directional_rectangle(
                IRanges::start(hits)[i], IRanges::end(hits)[i],
                lanes[i], vars$rangeHeight
            )
        } else {
            .directional_polygon(
                IRanges::start(hits)[i], IRanges::end(hits)[i],
                GenomicRanges::strand(hits)[i], lanes[i], vars$rangeHeight,
                tip_width
            )
        }
        grid::grid.polygon(
            x = grid::unit(polygon$x, "native"),
            y = grid::unit(polygon$y, "native"),
            gp = grid::gpar(
                fill = grDevices::adjustcolor(fills[i], alpha.f = vars$alpha),
                col = grDevices::adjustcolor(vars$border, alpha.f = vars$alpha),
                lwd = vars$borderWidth
            )
        )
    }
    if (!identical(vars$directionIndicator, "none")) {
        strand_text <- as.character(GenomicRanges::strand(hits))
        glyph <- .directional_glyphs(
            strand_text, reverseAxis = xscale[1L] > xscale[2L]
        )
        indicator_width_native <- abs(grid::convertWidth(
            grid::unit(1.1, "mm"), "native", valueOnly = TRUE
        ))
        show_indicator <- nzchar(glyph) &
            IRanges::width(hits) >= indicator_width_native
        if (any(show_indicator)) {
            grid::grid.text(
                glyph[show_indicator],
                x = grid::unit(
                    (IRanges::start(hits)[show_indicator] +
                     IRanges::end(hits)[show_indicator]) / 2,
                    "native"
                ),
                y = grid::unit(lanes[show_indicator], "native"),
                gp = grid::gpar(
                    col = vars$directionColor,
                    fontsize = vars$directionFontsize,
                    fontface = "bold"
                )
            )
        }
    }
    if (vars$showLabels) {
        label_y <- if (identical(vars$labelPosition, "above")) {
            lanes + vars$rangeHeight / 2 + 0.08
        } else lanes
        midpoint <- (IRanges::start(hits) + IRanges::end(hits)) / 2
        label_layout <- .viewport_label_layout(midpoint, xscale)
        vertical_just <- if (identical(vars$labelPosition, "above")) {
            "bottom"
        } else "center"
        for (i in seq_along(hits)) {
            grid::grid.text(
                S4Vectors::mcols(hits)$label[i],
                x = grid::unit(label_layout$x[i], "native"),
                y = grid::unit(label_y[i], "native"),
                just = c(label_layout$just[i], vertical_just),
                gp = grid::gpar(
                    col = vars$labelColor, fontsize = vars$labelFontsize
                )
            )
        }
    }
    invisible(NULL)
}

#' Fixed-tip directional genomic ranges for Gviz
#'
#' Draw stranded genomic ranges as compact rectangles with a shallow triangular
#' point. The point has a fixed physical width and the same height as the body,
#' avoiding the oversized, span-dependent arrowheads used by standard Gviz
#' annotation arrows. Overlapping ranges are assigned deterministic lanes.
#'
#' @param hits A non-empty, single-chromosome `GRanges`.
#' @param name Gviz track title.
#' @param featureColumn Optional metadata column used to map `fill` colors.
#'   When `NULL`, `feature` is used if present; otherwise all ranges use one
#'   feature.
#' @param labelColumn Optional metadata column used for labels. When `NULL`,
#'   the feature values are used.
#' @param fill One color, a named feature-color vector, or an unnamed vector
#'   containing at least one color per feature.
#' @param border,boundaryColor,labelColor Colors for range outlines, the left
#'   track boundary, and optional labels.
#' @param triangleWidth Fixed triangle width in millimetres. Very short ranges
#'   cap the triangle at 45 percent of their displayed width.
#' @param geometry Range geometry. `"pointed"` retains the fixed-tip shape;
#'   `"rectangle"` draws a simple rectangular range.
#' @param directionIndicator Use `"arrow"` to add a subtle `>` or `<` strand
#'   glyph inside each sufficiently wide range, or use `"none"`.
#' @param directionColor,directionFontsize Appearance of the strand glyph.
#' @param rangeHeight Range height relative to one overlap lane.
#' @param stacking Either `"squish"` for deterministic non-overlapping lanes or
#'   `"dense"` for one lane.
#' @param showLabels Draw range labels.
#' @param labelPosition Draw labels `"inside"` or immediately `"above"` ranges.
#' @param labelFontsize Label size in points.
#' @param alpha Range alpha.
#' @param borderWidth Outline width.
#' @param showBoundary Draw a vertical boundary at the left edge of the panel.
#' @param boundaryWidth,boundaryInset Width and fractional end inset of the
#'   left boundary.
#' @param ... Additional Gviz display parameters.
#' @return A `Gviz::CustomTrack`.
#' @examples
#' hits <- GenomicRanges::GRanges(
#'     "chr1", IRanges::IRanges(c(10, 35), width = c(20, 15)),
#'     strand = c("+", "-"), feature = c("enhancer", "promoter")
#' )
#' DirectionalRangeTrack(
#'     hits,
#'     fill = c(enhancer = "#0072B2", promoter = "#D55E00")
#' )
#' @export
DirectionalRangeTrack <- function(
        hits, name = "Directional ranges", featureColumn = NULL,
        labelColumn = NULL, fill = "#4C78A8", border = "#333333",
        triangleWidth = 1.5, rangeHeight = 0.62,
        geometry = c("pointed", "rectangle"),
        directionIndicator = c("none", "arrow"),
        directionColor = "#F3F3F1", directionFontsize = 5,
        stacking = c("squish", "dense"), showLabels = FALSE,
        labelPosition = c("inside", "above"),
        labelFontsize = 6, labelColor = "#111111", alpha = 1,
        borderWidth = 0.6, showBoundary = TRUE,
        boundaryColor = "#222222", boundaryWidth = 0.7,
        boundaryInset = 0.06, ...) {
    if (!methods::is(hits, "GRanges") || !length(hits)) {
        stop("`hits` must be a non-empty GRanges object.", call. = FALSE)
    }
    .validate_one_based_granges(hits, "hits")
    chromosomes <- unique(as.character(GenomicRanges::seqnames(hits)))
    if (length(chromosomes) != 1L) {
        stop("`hits` must contain one chromosome.", call. = FALSE)
    }
    stacking <- match.arg(stacking)
    geometry <- match.arg(geometry)
    directionIndicator <- match.arg(directionIndicator)
    labelPosition <- match.arg(labelPosition)
    columns <- names(S4Vectors::mcols(hits))
    featureColumn <- .validate_optional_column_name(
        featureColumn, "featureColumn"
    )
    labelColumn <- .validate_optional_column_name(labelColumn, "labelColumn")
    if (is.null(featureColumn)) {
        featureColumn <- if ("feature" %in% columns) "feature" else NULL
    }
    if (!is.null(featureColumn) && !featureColumn %in% columns) {
        stop("`hits` lacks feature column `", featureColumn, "`.", call. = FALSE)
    }
    features <- if (is.null(featureColumn)) {
        rep.int("range", length(hits))
    } else as.character(S4Vectors::mcols(hits)[[featureColumn]])
    if (anyNA(features) || any(!nzchar(features))) {
        stop("Directional range features must be non-missing strings.", call. = FALSE)
    }
    if (!is.null(labelColumn) && !labelColumn %in% columns) {
        stop("`hits` lacks label column `", labelColumn, "`.", call. = FALSE)
    }
    labels <- if (is.null(labelColumn)) features else
        as.character(S4Vectors::mcols(hits)[[labelColumn]])
    if (anyNA(labels) || any(!nzchar(labels))) {
        stop("Directional range labels must be non-missing strings.",
             call. = FALSE)
    }
    triangleWidth <- .validate_scalar_number(
        triangleWidth, "triangleWidth", minimum = 0, minimumInclusive = FALSE
    )
    rangeHeight <- .validate_scalar_number(
        rangeHeight, "rangeHeight", minimum = 0, maximum = 1,
        minimumInclusive = FALSE
    )
    labelFontsize <- .validate_scalar_number(
        labelFontsize, "labelFontsize", minimum = 0, minimumInclusive = FALSE
    )
    directionFontsize <- .validate_scalar_number(
        directionFontsize, "directionFontsize", minimum = 0,
        minimumInclusive = FALSE
    )
    alpha <- .validate_scalar_number(alpha, "alpha", minimum = 0, maximum = 1)
    borderWidth <- .validate_scalar_number(
        borderWidth, "borderWidth", minimum = 0
    )
    border <- .validate_color(border, "border")
    labelColor <- .validate_color(labelColor, "labelColor")
    directionColor <- .validate_color(directionColor, "directionColor")
    showLabels <- .validate_flag(showLabels, "showLabels")
    .validate_track_boundary(
        showBoundary, boundaryColor, boundaryWidth, boundaryInset
    )
    S4Vectors::mcols(hits)$label <- labels
    hits <- sort(hits, ignore.strand = TRUE)
    S4Vectors::mcols(hits)$lane <- .assign_directional_lanes(hits, stacking)
    features <- if (is.null(featureColumn)) rep.int("range", length(hits)) else
        as.character(S4Vectors::mcols(hits)[[featureColumn]])
    fills <- .resolve_directional_fills(features, fill)
    vars <- list(
        hits = hits, fills = fills, triangleWidth = triangleWidth,
        rangeHeight = rangeHeight, geometry = geometry,
        directionIndicator = directionIndicator,
        directionColor = directionColor,
        directionFontsize = directionFontsize,
        showLabels = showLabels, labelPosition = labelPosition,
        labelFontsize = labelFontsize, labelColor = labelColor,
        alpha = alpha, border = border, borderWidth = borderWidth,
        showBoundary = showBoundary, boundaryColor = boundaryColor,
        boundaryWidth = boundaryWidth, boundaryInset = boundaryInset
    )
    Gviz::CustomTrack(
        plottingFunction = function(GdObject, prepare = FALSE, ...) {
            if (!prepare) .draw_directional_range_track(GdObject)
            invisible(GdObject)
        },
        variables = vars, name = name, ...
    )
}
