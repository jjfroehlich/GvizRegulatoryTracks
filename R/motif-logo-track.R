.as_ppm_matrix <- function(x, matrixType = "PPM") {
    if (inherits(x, "universalmotif")) {
        if (!requireNamespace("universalmotif", quietly = TRUE)) {
            stop("Package `universalmotif` is required for universalmotif inputs.", call. = FALSE)
        }
        x <- universalmotif::convert_type(x, "PPM")["motif"]
    }
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("Each motif must be a numeric matrix or universalmotif object.", call. = FALSE)
    }
    dna <- c("A", "C", "G", "T")
    if (!is.null(rownames(x)) && all(dna %in% rownames(x))) {
        x <- x[dna, , drop = FALSE]
    } else if (!is.null(colnames(x)) && all(dna %in% colnames(x))) {
        x <- t(x[, dna, drop = FALSE])
    } else if (nrow(x) == 4L) {
        rownames(x) <- dna
    } else {
        stop("Motif matrices need A/C/G/T rows (or columns).", call. = FALSE)
    }
    if (any(!is.finite(x)) || any(x < 0)) {
        stop("PPM/PCM entries must be finite and non-negative.", call. = FALSE)
    }
    matrixType <- toupper(matrixType)
    if (!matrixType %in% c("PPM", "PCM")) {
        stop("Plain matrices support `matrixType = 'PPM'` or `'PCM'`.", call. = FALSE)
    }
    totals <- colSums(x)
    if (any(totals <= 0)) stop("Every motif column must have positive mass.", call. = FALSE)
    x <- sweep(x, 2L, totals, "/")
    rownames(x) <- dna
    x
}

.reverse_complement_ppm <- function(ppm) {
    out <- ppm[c("T", "G", "C", "A"), ncol(ppm):1L, drop = FALSE]
    rownames(out) <- c("A", "C", "G", "T")
    out
}

.oriented_motif_ppm <- function(ppm, hitStrand = "+", reverseAxis = FALSE,
                                reverseComplementMinus = TRUE) {
    hit_minus <- identical(as.character(hitStrand), "-")
    if (xor(hit_minus && isTRUE(reverseComplementMinus), isTRUE(reverseAxis))) {
        ppm <- .reverse_complement_ppm(ppm)
    }
    ppm
}

.information_heights <- function(ppm) {
    entropy_terms <- ifelse(ppm > 0, ppm * log2(ppm), 0)
    information <- pmax(0, 2 + colSums(entropy_terms))
    sweep(ppm, 2L, information, "*")
}

.assign_motif_lanes <- function(hits) {
    ord <- order(IRanges::start(hits), IRanges::end(hits),
                 as.character(S4Vectors::mcols(hits)$motif_id))
    lanes <- integer(length(hits))
    lane_ends <- integer()
    for (idx in ord) {
        available <- which(IRanges::start(hits)[idx] > lane_ends)
        lane <- if (length(available)) available[1L] else length(lane_ends) + 1L
        lanes[idx] <- lane
        if (lane > length(lane_ends)) lane_ends <- c(lane_ends, IRanges::end(hits)[idx])
        else lane_ends[lane] <- IRanges::end(hits)[idx]
    }
    lanes
}

.normalize_motif_collection <- function(motifs) {
    if (is.character(motifs) && length(motifs) == 1L) {
        if (!requireNamespace("universalmotif", quietly = TRUE)) {
            stop("Package `universalmotif` is required to read MEME files.", call. = FALSE)
        }
        if (!file.exists(motifs)) stop("Motif file does not exist: ", motifs, call. = FALSE)
        motifs <- universalmotif::read_meme(motifs)
    }
    if (inherits(motifs, "universalmotif")) {
        motifs <- list(motifs)
    }
    if (!is.list(motifs) || !length(motifs)) {
        stop("`motifs` must be a motif object, named list, or MEME file path.",
             call. = FALSE)
    }
    motif_names <- names(motifs)
    if (is.null(motif_names)) motif_names <- rep("", length(motifs))
    for (i in seq_along(motifs)) {
        if (!nzchar(motif_names[i]) && inherits(motifs[[i]], "universalmotif")) {
            motif_names[i] <- as.character(methods::slot(motifs[[i]], "name"))
        }
    }
    if (any(!nzchar(motif_names))) {
        stop("Plain motif matrices must be supplied in a named list.", call. = FALSE)
    }
    if (anyDuplicated(motif_names)) stop("Motif IDs must be unique.", call. = FALSE)
    names(motifs) <- motif_names
    motifs
}

.coerce_motif_hits <- function(hits) {
    if (is.character(hits) && length(hits) == 1L) {
        if (!requireNamespace("rtracklayer", quietly = TRUE)) {
            stop("Package `rtracklayer` is required to read BED motif hits.", call. = FALSE)
        }
        if (!file.exists(hits)) stop("Motif-hit file does not exist: ", hits, call. = FALSE)
        hits <- rtracklayer::import(hits)
    }
    hits
}

.resolve_motif_id_column <- function(hits, motifIdColumn) {
    columns <- names(S4Vectors::mcols(hits))
    if (!is.null(motifIdColumn)) {
        if (!is.character(motifIdColumn) || length(motifIdColumn) != 1L ||
            is.na(motifIdColumn) || !nzchar(motifIdColumn)) {
            stop("`motifIdColumn` must be one non-empty column name.",
                 call. = FALSE)
        }
        return(motifIdColumn)
    }
    for (candidate in c("motif_id", "name", "pattern_name")) {
        if (candidate %in% columns) return(candidate)
    }
    if (!is.null(names(hits)) && all(nzchar(names(hits)))) {
        S4Vectors::mcols(hits)$motif_id <- names(hits)
        return(list(hits = hits, column = "motif_id"))
    }
    stop("Could not infer motif IDs; provide `motifIdColumn`.", call. = FALSE)
}

.prepare_motif_data <- function(hits, motifs, motifIdColumn, matrixType) {
    hits <- .coerce_motif_hits(hits)
    if (!methods::is(hits, "GRanges") || !length(hits)) {
        stop("`hits` must be a non-empty GRanges object.", call. = FALSE)
    }
    .validate_one_based_granges(hits, "hits")
    chromosomes <- unique(as.character(GenomicRanges::seqnames(hits)))
    if (length(chromosomes) != 1L) stop("`hits` must contain one chromosome.", call. = FALSE)
    resolved <- .resolve_motif_id_column(hits, motifIdColumn)
    if (is.list(resolved)) {
        hits <- resolved$hits
        motifIdColumn <- resolved$column
    } else {
        motifIdColumn <- resolved
    }
    if (!motifIdColumn %in% names(S4Vectors::mcols(hits))) {
        stop("`hits` lacks motif ID column `", motifIdColumn, "`.", call. = FALSE)
    }
    ids <- as.character(S4Vectors::mcols(hits)[[motifIdColumn]])
    motifs <- .normalize_motif_collection(motifs)
    absent <- setdiff(unique(ids), names(motifs))
    if (length(absent)) stop("Missing motif matrices: ", paste(absent, collapse = ", "), call. = FALSE)
    motifs <- motifs[unique(ids)]
    ppms <- lapply(motifs, .as_ppm_matrix, matrixType = matrixType)
    expected <- vapply(ids, function(id) ncol(ppms[[id]]), integer(1))
    bad <- IRanges::width(hits) != expected
    if (any(bad)) {
        stop("Motif-hit widths must exactly equal their motif matrix widths: ",
             paste(unique(ids[bad]), collapse = ", "), call. = FALSE)
    }
    S4Vectors::mcols(hits)$motif_id <- ids
    S4Vectors::mcols(hits)$lane <- .assign_motif_lanes(hits)
    list(hits = hits, ppms = ppms, chromosome = chromosomes)
}

.resolve_motif_colors <- function(ids, motifColorBy, motifPalette) {
    motifColorBy <- match.arg(motifColorBy, c("nucleotide", "motif_id"))
    ids <- sort(unique(as.character(ids)))
    if (motifColorBy == "nucleotide") return(structure(character(), names = character()))
    if (is.null(motifPalette)) {
        colors <- grDevices::hcl.colors(length(ids), palette = "Dark 3")
        names(colors) <- ids
        return(colors)
    }
    motifPalette <- .validate_color_vector(motifPalette, "motifPalette")
    if (!is.null(names(motifPalette)) && anyNA(names(motifPalette))) {
        stop("`motifPalette` names must not be missing.", call. = FALSE)
    }
    if (!is.null(names(motifPalette)) && any(nzchar(names(motifPalette)))) {
        missing_ids <- setdiff(ids, names(motifPalette))
        if (length(missing_ids)) {
            stop("`motifPalette` lacks colors for motif IDs: ",
                 paste(missing_ids, collapse = ", "), call. = FALSE)
        }
        colors <- unname(motifPalette[ids])
    } else {
        if (length(motifPalette) < length(ids)) {
            stop("An unnamed `motifPalette` needs at least one color per motif ID.",
                 call. = FALSE)
        }
        colors <- motifPalette[seq_along(ids)]
    }
    names(colors) <- ids
    colors
}

.mix_with_white <- function(colors, amount) {
    if (!length(colors)) return(character())
    amount <- rep_len(as.numeric(amount), length(colors))
    rgb <- grDevices::col2rgb(colors) / 255
    mixed <- sweep(rgb, 2L, 1 - amount, "*") +
        matrix(rep(amount, each = 3L), nrow = 3L)
    grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ])
}

.score_data_limits <- function(values, fallback = c(0, 1)) {
    values <- as.numeric(values)
    values <- values[is.finite(values)]
    if (!length(values)) return(fallback)
    limits <- range(values)
    if (limits[1] < limits[2]) return(limits)
    value <- limits[1]
    if (value >= 0 && value <= 1) {
        delta <- 0.05
        limits <- c(max(0, value - delta), min(1, value + delta))
        if (limits[1] == limits[2]) {
            limits <- if (value == 0) c(0, 0.1) else c(0.9, 1)
        }
    } else {
        delta <- max(abs(value) * 0.05, 0.5)
        limits <- value + c(-delta, delta)
    }
    limits
}

.resolve_score_limit_spec <- function(scoreLimits) {
    if (is.null(scoreLimits)) return(list(mode = "global", limits = NULL))
    if (is.character(scoreLimits) && length(scoreLimits) == 1L) {
        mode <- match.arg(tolower(scoreLimits), c("view", "global"))
        return(list(mode = mode, limits = NULL))
    }
    if (!is.numeric(scoreLimits) || length(scoreLimits) != 2L ||
        any(!is.finite(scoreLimits)) || scoreLimits[1] >= scoreLimits[2]) {
        stop("`scoreLimits` must contain two increasing finite numbers.",
             call. = FALSE)
    }
    list(mode = "fixed", limits = scoreLimits)
}

.motif_score_styles <- function(hits, baseFill, scoreColumn, scoreAesthetic,
                                scoreLimits, scorePalette, scoreColor,
                                scoreBrightness, scoreOpacity) {
    if (!is.character(scoreAesthetic) || anyNA(scoreAesthetic)) {
        stop("`scoreAesthetic` must be a character vector.", call. = FALSE)
    }
    aesthetics <- unique(tolower(scoreAesthetic))
    if (!length(aesthetics)) aesthetics <- "none"
    allowed <- c("none", "fill", "border", "opacity", "brightness")
    if (any(!aesthetics %in% allowed) ||
        ("none" %in% aesthetics && length(aesthetics) > 1L)) {
        stop("`scoreAesthetic` must be 'none' or any of 'fill', 'border', 'opacity', and 'brightness'.",
             call. = FALSE)
    }
    style <- data.frame(
        fill = baseFill,
        border = rep(NA_character_, length(hits)),
        alpha = rep(1, length(hits)),
        brightness = rep(0, length(hits)),
        score = rep(NA_real_, length(hits)),
        scaledScore = rep(NA_real_, length(hits)),
        stringsAsFactors = FALSE
    )
    if (!is.null(scorePalette)) {
        scorePalette <- .validate_color_vector(
            scorePalette, "scorePalette", minimumLength = 2L
        )
    }
    scoreColor <- .validate_color(scoreColor, "scoreColor")
    if (!is.numeric(scoreBrightness) || length(scoreBrightness) != 2L ||
        any(!is.finite(scoreBrightness)) || any(scoreBrightness < 0) ||
        any(scoreBrightness > 1)) {
        stop("`scoreBrightness` must contain two white-mixing fractions within 0--1.",
             call. = FALSE)
    }
    if (!is.numeric(scoreOpacity) || length(scoreOpacity) != 2L ||
        any(!is.finite(scoreOpacity)) || scoreOpacity[1] < 0 ||
        scoreOpacity[2] > 1 || scoreOpacity[1] > scoreOpacity[2]) {
        stop("`scoreOpacity` must be an increasing two-number range within 0--1.",
             call. = FALSE)
    }
    if (identical(aesthetics, "none")) {
        attr(style, "scoreLimits") <- NULL
        attr(style, "scoreAesthetic") <- aesthetics
        return(style)
    }
    if ("fill" %in% aesthetics && any(!is.na(baseFill))) {
        stop("Score-mapped `fill` cannot be combined with `motifColorBy = 'motif_id'`; use brightness, opacity, or border for the score.",
             call. = FALSE)
    }
    if (is.null(scoreColumn) || !scoreColumn %in% names(S4Vectors::mcols(hits))) {
        stop("A valid `scoreColumn` is required for motif score styling.", call. = FALSE)
    }
    values <- .score_values(
        S4Vectors::mcols(hits)[[scoreColumn]],
        paste0("motif score column `", scoreColumn, "`")
    )
    if (is.null(scoreLimits)) {
        scoreLimits <- if (all(values >= 0 & values <= 1)) c(0, 1) else range(values)
    }
    if (!is.numeric(scoreLimits) || length(scoreLimits) != 2L ||
        any(!is.finite(scoreLimits)) || scoreLimits[1] >= scoreLimits[2]) {
        stop("`scoreLimits` must contain two increasing finite numbers.", call. = FALSE)
    }
    scaled <- pmin(1, pmax(0, (values - scoreLimits[1]) / diff(scoreLimits)))
    brightness <- scoreBrightness[1] + scaled * diff(scoreBrightness)
    if (is.null(scorePalette)) {
        mapped <- .mix_with_white(rep(scoreColor, length(hits)), brightness)
    } else {
        palette <- grDevices::colorRampPalette(scorePalette)(256L)
        mapped <- palette[pmax(1L, pmin(256L, floor(scaled * 255) + 1L))]
    }
    if ("fill" %in% aesthetics) style$fill <- mapped
    if ("border" %in% aesthetics) style$border <- mapped
    if ("brightness" %in% aesthetics) style$brightness <- brightness
    if ("opacity" %in% aesthetics) {
        style$alpha <- scoreOpacity[1] + scaled * diff(scoreOpacity)
    }
    style$score <- values
    style$scaledScore <- scaled
    attr(style, "scoreLimits") <- scoreLimits
    attr(style, "scoreAesthetic") <- aesthetics
    style
}

.draw_motif_logo <- function(hit, ppm, baseline, motifHeight, reverse_axis,
                             colors, font, style, scoreBorderWidth,
                             reverseComplementMinus) {
    ppm <- .oriented_motif_ppm(
        ppm, GenomicRanges::strand(hit), reverse_axis, reverseComplementMinus
    )
    heights <- .information_heights(ppm)
    bit_scale <- motifHeight / 2
    for (j in seq_len(ncol(ppm))) {
        position <- IRanges::start(hit) + j - 1L
        h <- heights[, j] * bit_scale
        order_bases <- names(sort(h, decreasing = FALSE))
        bottom <- baseline
        for (base in order_bases) {
            top <- bottom + h[base]
            fill <- if (is.na(style$fill)) colors[base] else style$fill
            if (style$brightness > 0) fill <- .mix_with_white(fill, style$brightness)
            .draw_letter_glyph(base, position - 0.5, position + 0.5,
                               bottom, top, fill, alpha = style$alpha, xpad = 0.02,
                               font = font, borderColor = style$border,
                               borderWidth = scoreBorderWidth)
            bottom <- top
        }
    }
    invisible(NULL)
}

.motif_lane_layout <- function(nlanes, trackHeight, motifHeight, laneGap,
                               showLabels) {
    nlanes <- max(1L, as.integer(nlanes))
    label_height <- if (isTRUE(showLabels)) 3.2 else 0
    outer_padding <- 0.8
    gap_height <- laneGap * motifHeight / 2
    requested <- 2 * outer_padding +
        nlanes * (motifHeight + label_height) +
        (nlanes - 1L) * gap_height
    scale <- min(1, trackHeight / requested)
    actual_motif_height <- motifHeight * scale
    actual_label_height <- label_height * scale
    actual_gap_height <- gap_height * scale
    actual_outer_padding <- outer_padding * scale
    used <- 2 * actual_outer_padding +
        nlanes * (actual_motif_height + actual_label_height) +
        (nlanes - 1L) * actual_gap_height
    offset <- max(0, (trackHeight - used) / 2)
    baselines <- offset + actual_outer_padding +
        (seq_len(nlanes) - 1L) *
        (actual_motif_height + actual_label_height + actual_gap_height)
    list(
        baselines = baselines,
        motifHeight = actual_motif_height,
        labelHeight = actual_label_height,
        gapHeight = actual_gap_height,
        labelFontsize = 6 * scale,
        scale = scale,
        requiredHeight = requested
    )
}

.motif_hit_labels <- function(hits, showStrand) {
    labels <- if ("label" %in% names(S4Vectors::mcols(hits))) {
        as.character(S4Vectors::mcols(hits)$label)
    } else {
        as.character(S4Vectors::mcols(hits)$motif_id)
    }
    if (isTRUE(showStrand)) {
        labels <- paste0(labels, " (", as.character(GenomicRanges::strand(hits)), ")")
    }
    labels
}

.motif_label_layout <- function(hits, xscale, edgeFraction = 0.12) {
    midpoint <- (IRanges::start(hits) + IRanges::end(hits)) / 2
    .viewport_label_layout(midpoint, xscale, edgeFraction = edgeFraction)
}

.score_legend_styles <- function(vars, n = 64L) {
    scaled <- seq(0, 1, length.out = n)
    brightness <- vars$scoreBrightness[1] + scaled * diff(vars$scoreBrightness)
    mapped <- if (is.null(vars$scorePalette)) {
        .mix_with_white(rep(vars$scoreColor, n), brightness)
    } else {
        grDevices::colorRampPalette(vars$scorePalette)(n)
    }
    aesthetics <- vars$scoreAesthetic
    fill <- if (identical(aesthetics, "opacity")) {
        rep(vars$scoreLegendColor, n)
    } else mapped
    alpha <- if ("opacity" %in% aesthetics) {
        vars$scoreOpacity[1] + scaled * diff(vars$scoreOpacity)
    } else rep(1, n)
    data.frame(fill = fill, alpha = alpha, scaled = scaled)
}

.score_legend_width <- function(title, minimumWidth, fontsize = 5.5,
                                padding = 1.5) {
    title <- if (is.null(title) || !length(title)) "" else as.character(title)[1L]
    lines <- strsplit(title, "\n", fixed = TRUE)[[1L]]
    title_width <- max(vapply(lines, function(line) {
        grid::convertWidth(
            grid::grobWidth(grid::textGrob(line, gp = grid::gpar(
                fontsize = fontsize, fontface = "bold"
            ))),
            "mm", valueOnly = TRUE
        )
    }, numeric(1L)))
    max(minimumWidth, title_width + 2 * padding)
}

.score_legend_layout <- function(title) {
    title <- if (is.null(title) || !length(title)) "" else as.character(title)[1L]
    title_lines <- if (nzchar(title)) {
        length(strsplit(title, "\n", fixed = TRUE)[[1L]])
    } else 0L
    # Reserve a title band inside the legend border. In particular, anchoring
    # multi-line titles from their top edge prevents their first line from
    # crossing the rounded rectangle at short track heights.
    title_band <- if (title_lines) min(0.48, 0.12 + 0.16 * title_lines) else 0
    list(
        title_y = 0.95,
        title_just = c("center", "top"),
        gradient_min = 0.14,
        gradient_max = 0.91 - title_band
    )
}

.draw_score_legend <- function(vars, scoreLimits) {
    if (!isTRUE(vars$showScoreLegend) || is.null(scoreLimits)) {
        return(invisible(NULL))
    }
    styles <- .score_legend_styles(vars)
    layout <- .score_legend_layout(vars$scoreLegendTitle)
    n <- nrow(styles)
    track_height <- grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
    legend_height <- min(vars$scoreLegendHeight, max(3, track_height - 2))
    legend_width <- .score_legend_width(
        vars$scoreLegendTitle, vars$scoreLegendWidth
    )
    grid::pushViewport(grid::viewport(
        x = grid::unit(1, "npc") - grid::unit(1, "mm"), just = "right",
        y = grid::unit(0.5, "npc"),
        width = grid::unit(legend_width, "mm"),
        height = grid::unit(legend_height, "mm"),
        xscale = c(0, 1), yscale = c(0, 1), clip = "off"
    ))
    on.exit(grid::popViewport(), add = TRUE)
    grid::grid.roundrect(
        r = grid::unit(0.8, "mm"),
        gp = grid::gpar(fill = grDevices::adjustcolor("white", alpha.f = 0.96),
                        col = "#777777", lwd = 0.7)
    )
    ymin <- layout$gradient_min
    ymax <- layout$gradient_max
    edges <- seq(ymin, ymax, length.out = n + 1L)
    for (i in seq_len(n)) {
        grid::grid.rect(
            x = grid::unit(0.25, "native"),
            y = grid::unit(mean(edges[c(i, i + 1L)]), "native"),
            width = grid::unit(0.22, "native"),
            height = grid::unit(diff(edges)[i] * 1.02, "native"),
            gp = grid::gpar(
                fill = grDevices::adjustcolor(styles$fill[i], alpha.f = styles$alpha[i]),
                col = NA
            )
        )
    }
    grid::grid.rect(
        x = grid::unit(0.25, "native"), y = grid::unit(mean(c(ymin, ymax)), "native"),
        width = grid::unit(0.22, "native"), height = grid::unit(ymax - ymin, "native"),
        gp = grid::gpar(fill = NA, col = "#555555", lwd = 0.55)
    )
    ticks <- seq(scoreLimits[1], scoreLimits[2], length.out = 3L)
    tick_y <- seq(ymin, ymax, length.out = 3L)
    grid::grid.segments(
        x0 = grid::unit(0.37, "native"), x1 = grid::unit(0.44, "native"),
        y0 = grid::unit(tick_y, "native"), y1 = grid::unit(tick_y, "native"),
        gp = grid::gpar(col = "#444444", lwd = 0.5)
    )
    grid::grid.text(
        formatC(ticks, digits = vars$scoreLegendDigits, format = "fg", flag = "#"),
        x = grid::unit(0.48, "native"), y = grid::unit(tick_y, "native"),
        just = "left", gp = grid::gpar(col = "#333333", fontsize = 5.5)
    )
    grid::grid.text(
        vars$scoreLegendTitle, x = grid::unit(0.5, "native"),
        y = grid::unit(layout$title_y, "native"), just = layout$title_just,
        gp = grid::gpar(col = "#222222", fontsize = 5.5, fontface = "bold")
    )
    invisible(NULL)
}

.visible_motif_styles <- function(vars, keep) {
    hits <- vars$hits[keep]
    if (identical(vars$scoreAesthetic, "none")) {
        return(list(styles = vars$hitStyles[keep, , drop = FALSE],
                    limits = NULL))
    }
    limits <- switch(
        vars$scoreLimitMode,
        fixed = vars$scoreLimits,
        global = vars$globalScoreLimits,
        view = .score_data_limits(
            vars$hitStyles$score[keep], vars$globalScoreLimits
        )
    )
    styles <- .motif_score_styles(
        hits, vars$baseFill[keep], vars$scoreColumn, vars$scoreAesthetic,
        limits, vars$scorePalette, vars$scoreColor, vars$scoreBrightness,
        vars$scoreOpacity
    )
    list(styles = styles, limits = limits)
}

.motif_display_mode <- function(display, span, maxBases) {
    if (display != "auto") return(display)
    if (span <= maxBases) "logo" else "ranges"
}

.resolve_motif_range_styles <- function(styles, vars) {
    fill <- ifelse(is.na(styles$fill), vars$rangeFill, styles$fill)
    brighten <- is.finite(styles$brightness) & styles$brightness > 0
    if (any(brighten)) {
        fill[brighten] <- .mix_with_white(
            fill[brighten], styles$brightness[brighten]
        )
    }
    border <- ifelse(is.na(styles$border), vars$rangeBorder, styles$border)
    data.frame(
        fill = fill,
        border = border,
        alpha = vars$rangeAlpha * styles$alpha,
        stringsAsFactors = FALSE
    )
}

.draw_motif_ranges <- function(hits, styles, vars, xscale) {
    nlanes <- max(1L, if (length(hits)) max(S4Vectors::mcols(hits)$lane) else 1L)
    labels_visible <- vars$showLabels && length(hits) <= vars$maxRangeLabels
    upper_padding <- if (labels_visible) 0.85 else 0.5
    grid::pushViewport(grid::viewport(
        xscale = xscale, yscale = c(0.5, nlanes + upper_padding), clip = "on"
    ))
    on.exit(grid::popViewport(), add = TRUE)
    if (!length(hits)) return(invisible(NULL))
    lanes <- S4Vectors::mcols(hits)$lane
    range_styles <- .resolve_motif_range_styles(styles, vars)
    for (i in seq_along(hits)) {
        grid::grid.rect(
            x = grid::unit((IRanges::start(hits)[i] + IRanges::end(hits)[i]) / 2,
                           "native"),
            y = grid::unit(lanes[i], "native"),
            width = grid::unit(IRanges::width(hits)[i], "native"),
            height = grid::unit(0.64, "native"),
            gp = grid::gpar(
                fill = grDevices::adjustcolor(
                    range_styles$fill[i], alpha.f = range_styles$alpha[i]
                ),
                col = grDevices::adjustcolor(
                    range_styles$border[i], alpha.f = range_styles$alpha[i]
                ),
                lwd = 0.6
            )
        )
    }
    if (labels_visible) {
        ids <- .motif_hit_labels(hits, vars$showStrand)
        label_layout <- .motif_label_layout(hits, xscale)
        for (i in seq_along(ids)) {
            grid::grid.text(
                ids[i],
                x = grid::unit(label_layout$x[i], "native"),
                y = grid::unit(lanes[i] + 0.34, "native"),
                just = c(label_layout$just[i], "bottom"),
                gp = grid::gpar(col = "#333333", fontsize = 6)
            )
        }
    }
    invisible(NULL)
}

.draw_motif_track <- function(object) {
    vars <- methods::slot(object, "variables")
    xscale <- grid::current.viewport()$xscale
    visible <- sort(as.numeric(xscale))
    span <- floor(diff(visible)) + 1L
    display <- .motif_display_mode(vars$display, span, vars$maxBases)
    hits <- vars$hits
    keep <- IRanges::end(hits) >= visible[1] & IRanges::start(hits) <= visible[2]
    hits <- hits[keep]
    visible_styles <- .visible_motif_styles(vars, keep)
    styles <- visible_styles$styles
    if (display == "ranges") {
        .draw_motif_ranges(hits, styles, vars, xscale)
        if (vars$showBoundary) {
            .draw_track_left_boundary(
                vars$boundaryColor, vars$boundaryWidth, vars$boundaryInset
            )
        }
        .draw_score_legend(vars, visible_styles$limits)
        return(invisible(object))
    }
    nlanes <- max(1L, if (length(hits)) max(S4Vectors::mcols(hits)$lane) else 1L)
    track_height <- grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
    layout <- .motif_lane_layout(
        nlanes, track_height, vars$motifHeight, vars$laneGap, vars$showLabels
    )
    reverse_axis <- xscale[1] > xscale[2]
    grid::pushViewport(grid::viewport(
        xscale = xscale, yscale = c(0, track_height), clip = "on"
    ))
    on.exit(grid::popViewport(), add = TRUE)
    if (vars$showBoundary) {
        .draw_track_left_boundary(
            vars$boundaryColor, vars$boundaryWidth, vars$boundaryInset
        )
    }
    if (vars$showAxis) {
        axis_at <- layout$baselines[1L] +
            c(0, 0.5, 1) * layout$motifHeight
        grid::grid.yaxis(at = axis_at, label = c("0", "1", "2"), main = FALSE,
                         gp = grid::gpar(col = "#555555", fontsize = 7))
    }
    if (length(hits)) {
        for (i in seq_along(hits)) {
            id <- as.character(S4Vectors::mcols(hits)$motif_id[i])
            lane <- S4Vectors::mcols(hits)$lane[i]
            baseline <- layout$baselines[lane]
            .draw_motif_logo(hits[i], vars$ppms[[id]], baseline,
                             layout$motifHeight, reverse_axis,
                             vars$colors, vars$font, styles[i, ],
                             vars$scoreBorderWidth, vars$reverseComplementMinus)
            if (vars$showLabels) {
                y <- baseline + layout$motifHeight + layout$labelHeight / 2
                label_layout <- .motif_label_layout(hits[i], xscale)
                grid::grid.text(
                    .motif_hit_labels(hits[i], vars$showStrand),
                    x = grid::unit(label_layout$x, "native"),
                    y = grid::unit(y, "native"), just = label_layout$just,
                    gp = grid::gpar(
                        col = "#333333", fontsize = layout$labelFontsize
                    )
                )
            }
        }
    }
    .draw_score_legend(vars, visible_styles$limits)
    invisible(object)
}

#' Genomically positioned motif-logo track for Gviz
#'
#' Place information-content sequence logos over exact motif-hit ranges. Logos
#' are reverse-complemented for minus-strand hits, and overlapping hits are
#' placed into deterministic lanes.
#'
#' @param hits Single-chromosome `GRanges` of motif matches or a BED file path.
#' @param motifs Named list of PPM/PCM matrices, `universalmotif` object(s), or
#'   a MEME file path.
#' @param motifIdColumn Metadata column in `hits` matching `names(motifs)`.
#'   `NULL` automatically checks `motif_id`, BED `name`, and `pattern_name`.
#' @param labelColumn Optional metadata column used for displayed hit labels.
#'   `NULL` retains motif IDs, preserving the existing default.
#' @param matrixType Interpretation of plain matrices, `"PPM"` or `"PCM"`.
#' @param name Gviz track title.
#' @param colors Named nucleotide color vector.
#' @param motifColorBy Fill detailed logos and overview ranges by nucleotide or
#'   by resolved motif ID. Motif IDs may originate from `motif_id`, BED `name`,
#'   `pattern_name`, or `motifIdColumn`.
#' @param motifPalette Named motif-ID colors, or an unnamed vector with at least
#'   one color per motif. `NULL` uses a deterministic categorical HCL palette.
#' @param laneGap Gap between overlapping-hit lanes, expressed relative to the
#'   two-bit logo height. The default `0.30` becomes 15 percent of
#'   `motifHeight`.
#' @param motifHeight Maximum physical height in millimetres of a full two-bit
#'   motif stack. Logos retain this height when the Gviz track is taller and
#'   shrink uniformly only when the allocated track is too short for all lanes.
#'   Use `sizes` in `Gviz::plotTracks()` to allocate more room to tracks with
#'   several overlapping lanes.
#' @param showLabels Label each logo with its motif ID.
#' @param showStrand Append `(+)` or `(-)` to motif labels, making orientation
#'   explicit in dense multi-motif examples.
#' @param reverseComplementMinus Reverse-complement PPMs for minus-strand hits.
#'   The default displays each match in reference-genome orientation and also
#'   accounts for Gviz's `reverseStrand = TRUE` axis reversal.
#' @param showAxis Draw a 0--2 bit y axis.
#' @param showBoundary Draw a vertical boundary at the left edge of the data
#'   panel.
#' @param boundaryColor,boundaryWidth Color and line width of that boundary.
#' @param boundaryInset Fraction of track height inset at both ends of the
#'   boundary, leaving a visible gap between adjacent track segments.
#' @param maxBases Largest plotting viewport that will render logos.
#' @param display Rendering mode: automatic logo/range switching, forced logo,
#'   or forced ranges.
#' @param font A font name returned by `ggseqlogo::list_fonts(FALSE)`.
#' @param rangeFill,rangeBorder,rangeAlpha Appearance of overview ranges.
#' @param maxRangeLabels Maximum number of overview ranges that may be labeled.
#' @param scoreColumn Numeric hit metadata column used for confidence styling.
#'   Numeric factor and character columns are parsed strictly.
#' @param scoreAesthetic `"none"` or one or more of `"fill"`, `"border"`,
#'   `"opacity"`, and `"brightness"`. Styling is applied consistently to logos
#'   and overview ranges. Score fill and motif-ID fill are mutually exclusive;
#'   combine motif-ID fill with brightness, opacity, or border instead.
#' @param scoreLimits Numeric fixed score limits, `"view"` to use the minimum
#'   and maximum among currently visible hits, or `"global"` to use all hits in
#'   the track. View scaling is the default; use fixed limits when comparing
#'   score encodings across panels.
#' @param scorePalette Optional explicit sequential low-to-high colors for fill
#'   or border mapping. `NULL` derives shades from `scoreColor`.
#' @param scoreColor Base color used for score-mapped fill or border when
#'   `scorePalette = NULL`.
#' @param scoreBrightness Fractions of white mixed into a base color at the low
#'   and high score limits. This controls `"brightness"` and the default
#'   single-color fill/border gradient.
#' @param scoreOpacity Low-to-high alpha range for opacity mapping.
#' @param scoreBorderWidth Border width for logo glyphs.
#' @param showScoreLegend Draw a compact continuous score legend at the right
#'   of the data panel whenever score styling is active.
#' @param scoreLegendTitle Legend title. `NULL` uses `scoreColumn`.
#' @param scoreLegendWidth,scoreLegendHeight Minimum legend width and fixed
#'   legend height in mm. Width expands automatically to contain the title;
#'   height is capped when a track is shorter and does not grow with the track.
#' @param scoreLegendDigits Significant digits used for legend tick labels.
#' @param scoreLegendColor Neutral legend fill when score controls brightness,
#'   opacity, or border while motif identity controls the logo fill.
#' @param ... Additional Gviz display parameters.
#' @return A `Gviz::CustomTrack`.
#' @export
MotifLogoTrack <- function(hits, motifs, motifIdColumn = NULL,
                           labelColumn = NULL,
                           matrixType = c("PPM", "PCM"), name = "Motifs",
                           colors = .nucleotide_colors(), laneGap = 0.30,
                           motifHeight = 10,
                           motifColorBy = c("nucleotide", "motif_id"),
                           motifPalette = NULL,
                           showLabels = FALSE, showStrand = FALSE,
                           reverseComplementMinus = TRUE, showAxis = TRUE,
                           showBoundary = TRUE,
                           boundaryColor = "#222222",
                           boundaryWidth = 0.7,
                           boundaryInset = 0.06,
                           maxBases = 1000L,
                           display = c("auto", "logo", "ranges"),
                           font = "roboto_medium",
                           rangeFill = "#7A5E8E", rangeBorder = "#4A3A57",
                           rangeAlpha = 0.75, maxRangeLabels = 50L,
                           scoreColumn = NULL, scoreAesthetic = "none",
                           scoreLimits = "view",
                           scorePalette = NULL, scoreColor = "#111111",
                           scoreBrightness = c(0.85, 0),
                           scoreOpacity = c(0.25, 1),
                           scoreBorderWidth = 0.7,
                           showScoreLegend = TRUE, scoreLegendTitle = NULL,
                           scoreLegendWidth = 16, scoreLegendHeight = 23,
                           scoreLegendDigits = 3L,
                           scoreLegendColor = "#444444", ...) {
    matrixType <- match.arg(matrixType)
    display <- match.arg(display)
    motifColorBy <- match.arg(motifColorBy)
    colors <- .validate_colors(colors)
    if (!is.null(scoreColumn) &&
        (!is.character(scoreColumn) || length(scoreColumn) != 1L ||
         is.na(scoreColumn) || !nzchar(scoreColumn))) {
        stop("`scoreColumn` must be NULL or one non-empty column name.",
             call. = FALSE)
    }
    if (!is.null(labelColumn) &&
        (!is.character(labelColumn) || length(labelColumn) != 1L ||
         is.na(labelColumn) || !nzchar(labelColumn))) {
        stop("`labelColumn` must be NULL or one non-empty column name.",
             call. = FALSE)
    }
    showLabels <- .validate_flag(showLabels, "showLabels")
    showStrand <- .validate_flag(showStrand, "showStrand")
    reverseComplementMinus <- .validate_flag(
        reverseComplementMinus, "reverseComplementMinus"
    )
    showAxis <- .validate_flag(showAxis, "showAxis")
    showScoreLegend <- .validate_flag(showScoreLegend, "showScoreLegend")
    maxBases <- .validate_scalar_number(
        maxBases, "maxBases", minimum = 0, whole = TRUE,
        minimumInclusive = FALSE
    )
    maxRangeLabels <- .validate_scalar_number(
        maxRangeLabels, "maxRangeLabels", minimum = 0, whole = TRUE
    )
    rangeAlpha <- .validate_scalar_number(
        rangeAlpha, "rangeAlpha", minimum = 0, maximum = 1
    )
    scoreBorderWidth <- .validate_scalar_number(
        scoreBorderWidth, "scoreBorderWidth", minimum = 0
    )
    rangeFill <- .validate_color(rangeFill, "rangeFill")
    rangeBorder <- .validate_color(rangeBorder, "rangeBorder")
    if (!is.numeric(laneGap) || length(laneGap) != 1L ||
        !is.finite(laneGap) || laneGap < 0) {
        stop("`laneGap` must be one non-negative number.", call. = FALSE)
    }
    if (!is.numeric(motifHeight) || length(motifHeight) != 1L ||
        !is.finite(motifHeight) || motifHeight <= 0) {
        stop("`motifHeight` must be one positive number of mm.", call. = FALSE)
    }
    .validate_track_boundary(
        showBoundary, boundaryColor, boundaryWidth, boundaryInset
    )
    .load_glyph_font(font)
    prepared <- .prepare_motif_data(hits, motifs, motifIdColumn, matrixType)
    if (!is.null(labelColumn)) {
        if (!labelColumn %in% names(S4Vectors::mcols(prepared$hits))) {
            stop("`hits` lacks label column `", labelColumn, "`.", call. = FALSE)
        }
        labels <- as.character(S4Vectors::mcols(prepared$hits)[[labelColumn]])
        if (anyNA(labels) || any(!nzchar(labels))) {
            stop("Motif labels must be non-missing strings.", call. = FALSE)
        }
        S4Vectors::mcols(prepared$hits)$label <- labels
    }
    motif_colors <- .resolve_motif_colors(
        S4Vectors::mcols(prepared$hits)$motif_id, motifColorBy, motifPalette
    )
    base_fill <- if (motifColorBy == "motif_id") {
        unname(motif_colors[as.character(S4Vectors::mcols(prepared$hits)$motif_id)])
    } else {
        rep(NA_character_, length(prepared$hits))
    }
    score_limit_spec <- .resolve_score_limit_spec(scoreLimits)
    score_values <- if (!is.null(scoreColumn) &&
                        scoreColumn %in% names(S4Vectors::mcols(prepared$hits))) {
        .score_values(
            S4Vectors::mcols(prepared$hits)[[scoreColumn]],
            paste0("motif score column `", scoreColumn, "`")
        )
    } else numeric()
    global_score_limits <- if (length(score_values)) {
        .score_data_limits(score_values)
    } else NULL
    initial_score_limits <- if (score_limit_spec$mode == "fixed") {
        score_limit_spec$limits
    } else global_score_limits
    hit_styles <- .motif_score_styles(
        prepared$hits, base_fill, scoreColumn, scoreAesthetic, initial_score_limits,
        scorePalette, scoreColor, scoreBrightness, scoreOpacity
    )
    resolved_score_limits <- if (score_limit_spec$mode == "fixed") {
        attr(hit_styles, "scoreLimits")
    } else NULL
    resolved_score_aesthetic <- attr(hit_styles, "scoreAesthetic")
    score_legend_active <- showScoreLegend &&
        !identical(resolved_score_aesthetic, "none")
    if (is.null(scoreLegendTitle)) scoreLegendTitle <- scoreColumn
    if (!is.null(scoreLegendTitle) &&
        (!is.character(scoreLegendTitle) || length(scoreLegendTitle) != 1L ||
         is.na(scoreLegendTitle))) {
        stop("`scoreLegendTitle` must be NULL or one character string.",
             call. = FALSE)
    }
    if (!is.numeric(scoreLegendWidth) || length(scoreLegendWidth) != 1L ||
        !is.finite(scoreLegendWidth) || scoreLegendWidth <= 0) {
        stop("`scoreLegendWidth` must be one positive number of mm.", call. = FALSE)
    }
    if (!is.numeric(scoreLegendHeight) || length(scoreLegendHeight) != 1L ||
        !is.finite(scoreLegendHeight) || scoreLegendHeight <= 0) {
        stop("`scoreLegendHeight` must be one positive number of mm.", call. = FALSE)
    }
    scoreLegendDigits <- .validate_scalar_number(
        scoreLegendDigits, "scoreLegendDigits", minimum = 0, whole = TRUE,
        minimumInclusive = FALSE
    )
    scoreLegendColor <- .validate_color(scoreLegendColor, "scoreLegendColor")
    vars <- c(prepared, list(
        colors = colors, laneGap = laneGap, motifHeight = motifHeight,
        showLabels = showLabels,
        showStrand = showStrand,
        reverseComplementMinus = reverseComplementMinus,
        motifColorBy = motifColorBy, motifPalette = motifPalette,
        motifColors = motif_colors, baseFill = base_fill,
        showAxis = showAxis, maxBases = maxBases,
        showBoundary = isTRUE(showBoundary), boundaryColor = boundaryColor,
        boundaryWidth = boundaryWidth, boundaryInset = boundaryInset,
        display = display, font = font, rangeFill = rangeFill,
        rangeBorder = rangeBorder, rangeAlpha = rangeAlpha,
        maxRangeLabels = maxRangeLabels, hitStyles = hit_styles,
        scoreColumn = scoreColumn, scoreAesthetic = resolved_score_aesthetic,
        scoreLimits = resolved_score_limits, scorePalette = scorePalette,
        scoreLimitMode = score_limit_spec$mode,
        globalScoreLimits = global_score_limits,
        scoreColor = scoreColor, scoreBrightness = scoreBrightness,
        scoreOpacity = scoreOpacity, scoreBorderWidth = scoreBorderWidth,
        showScoreLegend = score_legend_active,
        scoreLegendTitle = scoreLegendTitle,
        scoreLegendWidth = scoreLegendWidth,
        scoreLegendHeight = scoreLegendHeight,
        scoreLegendDigits = scoreLegendDigits,
        scoreLegendColor = scoreLegendColor
    ))
    Gviz::CustomTrack(
        plottingFunction = function(GdObject, prepare = FALSE, ...) {
            if (!prepare) .draw_motif_track(GdObject)
            invisible(GdObject)
        },
        variables = vars, name = name, ...
    )
}
