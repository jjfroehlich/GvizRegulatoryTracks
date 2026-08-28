#' Symmetric y limits for contribution tracks
#'
#' Compute common symmetric limits for one or more numeric score vectors or
#' `GRanges` objects. Supplying the same result to several contribution tracks
#' makes their nucleotide heights directly comparable.
#'
#' @param ... Numeric vectors, `GRanges` objects, or a single list containing
#'   either. For `GRanges`, values are taken from `scoreColumn`; numeric factor
#'   and character columns are parsed strictly.
#' @param scoreColumn Metadata column containing scores in `GRanges` inputs.
#' @param padding Non-negative fractional padding added above the largest
#'   absolute score.
#' @return A numeric vector of length two.
#' @export
contributionYlim <- function(..., scoreColumn = "score", padding = 0.05) {
    if (!is.character(scoreColumn) || length(scoreColumn) != 1L ||
        is.na(scoreColumn) || !nzchar(scoreColumn)) {
        stop("`scoreColumn` must be one non-empty column name.", call. = FALSE)
    }
    padding <- .validate_scalar_number(padding, "padding", minimum = 0)
    inputs <- list(...)
    if (length(inputs) == 1L && is.list(inputs[[1L]]) &&
        !methods::is(inputs[[1L]], "GRanges")) {
        inputs <- inputs[[1L]]
    }
    values <- unlist(lapply(inputs, function(x) {
        if (methods::is(x, "GRanges")) {
            if (!scoreColumn %in% names(S4Vectors::mcols(x))) {
                stop("A GRanges input lacks score column `", scoreColumn, "`.", call. = FALSE)
            }
            return(.score_values(
                S4Vectors::mcols(x)[[scoreColumn]],
                paste0("GRanges score column `", scoreColumn, "`"),
                allowNumericNonFinite = TRUE
            ))
        }
        if (!is.numeric(x)) {
            stop("Inputs must be numeric vectors or GRanges objects.", call. = FALSE)
        }
        as.numeric(x)
    }), use.names = FALSE)
    values <- values[is.finite(values)]
    if (!length(values)) {
        return(c(-1, 1))
    }
    lim <- max(abs(values))
    if (lim == 0) lim <- 1
    lim <- lim * (1 + padding)
    c(-lim, lim)
}

#' @rdname contributionYlim
#' @export
scoreYlim <- contributionYlim

.validate_one_based_granges <- function(x, argument, requireOne = FALSE) {
    if (!methods::is(x, "GRanges")) {
        stop("`", argument, "` must be a GRanges object.", call. = FALSE)
    }
    if (requireOne && length(x) != 1L) {
        stop("`", argument, "` must be one GRanges range.", call. = FALSE)
    }
    if (length(x) && any(IRanges::start(x) < 1L)) {
        stop("`", argument, "` must use positive, one-based GRanges coordinates.",
             call. = FALSE)
    }
    if (length(x) && any(IRanges::width(x) < 1L)) {
        stop("`", argument, "` must contain positive-width ranges.",
             call. = FALSE)
    }
    invisible(x)
}

.one_based_coordinate <- function(x, argument) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 1 || x != floor(x) || x > .Machine$integer.max) {
        stop("`", argument, "` must be one positive, one-based integer coordinate.",
             call. = FALSE)
    }
    as.integer(x)
}

.sequence_letters <- function(sequence) {
    if (methods::is(sequence, "DNAString")) {
        sequence <- as.character(sequence)
    }
    if (!is.character(sequence) || length(sequence) != 1L || is.na(sequence)) {
        stop("`sequence` must be one DNAString or one character string.", call. = FALSE)
    }
    letters <- strsplit(toupper(sequence), "", fixed = TRUE)[[1L]]
    letters[!letters %in% c("A", "C", "G", "T")] <- "N"
    letters
}

.expand_granges_scores <- function(scores, scoreColumn, missing, targetRange = NULL) {
    if (!is.null(targetRange) &&
        (!methods::is(targetRange, "GRanges") || length(targetRange) != 1L)) {
        stop("`region` must be one GRanges range.", call. = FALSE)
    }
    .validate_one_based_granges(scores, "scores")
    if (!is.null(targetRange)) {
        .validate_one_based_granges(targetRange, "region", requireOne = TRUE)
    }
    if (!length(scores) && is.null(targetRange)) {
        stop("`scores` GRanges must not be empty.", call. = FALSE)
    }
    if (!scoreColumn %in% names(S4Vectors::mcols(scores))) {
        stop("`scores` lacks metadata column `", scoreColumn, "`.", call. = FALSE)
    }
    chromosomes <- unique(c(
        as.character(GenomicRanges::seqnames(scores)),
        if (!is.null(targetRange)) as.character(GenomicRanges::seqnames(targetRange))
    ))
    if (length(chromosomes) != 1L) {
        stop("`scores` must contain exactly one chromosome.", call. = FALSE)
    }
    if (!IRanges::isDisjoint(GenomicRanges::ranges(scores))) {
        stop("`scores` contains overlapping genomic ranges.", call. = FALSE)
    }
    score_values <- .score_values(
        S4Vectors::mcols(scores)[[scoreColumn]],
        paste0("scores column `", scoreColumn, "`")
    )
    lo <- if (is.null(targetRange)) min(IRanges::start(scores)) else IRanges::start(targetRange)
    hi <- if (is.null(targetRange)) max(IRanges::end(scores)) else IRanges::end(targetRange)
    positions <- seq.int(lo, hi)
    values <- rep(NA_real_, length(positions))
    for (i in seq_along(scores)) {
        left <- max(lo, IRanges::start(scores)[i])
        right <- min(hi, IRanges::end(scores)[i])
        if (left <= right) {
            idx <- seq.int(left, right) - lo + 1L
            values[idx] <- score_values[i]
        }
    }
    if (anyNA(values)) {
        if (identical(missing, "error")) {
            stop("`scores` does not cover every base in its genomic span.", call. = FALSE)
        }
        values[is.na(values)] <- 0
    }
    list(chromosome = chromosomes, start = lo, positions = positions, scores = values)
}

.prepare_contribution_data <- function(sequence, scores, chromosome, start,
                                       scoreColumn, missing, normalization,
                                       targetRange = NULL) {
    if (methods::is(scores, "GRanges")) {
        expanded <- .expand_granges_scores(scores, scoreColumn, missing, targetRange)
        if (!missing(chromosome) && !is.null(chromosome) &&
            !identical(as.character(chromosome), expanded$chromosome)) {
            stop("`chromosome` disagrees with the GRanges seqname.", call. = FALSE)
        }
        if (!missing(start) && !is.null(start) &&
            .one_based_coordinate(start, "start") != expanded$start) {
            stop("`start` disagrees with the GRanges start coordinate.", call. = FALSE)
        }
        chromosome <- expanded$chromosome
        start <- expanded$start
        values <- expanded$scores
        positions <- expanded$positions
    } else {
        if (!is.numeric(scores)) {
            stop("`scores` must be numeric or a GRanges object.", call. = FALSE)
        }
        if (!length(scores)) {
            stop("`scores` must contain at least one basewise value.", call. = FALSE)
        }
        if (missing(chromosome) || is.null(chromosome) || length(chromosome) != 1L) {
            stop("`chromosome` is required for numeric scores.", call. = FALSE)
        }
        if (missing(start) || is.null(start) || length(start) != 1L || !is.finite(start)) {
            stop("`start` is required for numeric scores.", call. = FALSE)
        }
        values <- .score_values(scores, "scores")
        start <- .one_based_coordinate(start, "start")
        if (!is.null(targetRange)) {
            .validate_one_based_granges(targetRange, "region", requireOne = TRUE)
            target_chromosome <- as.character(GenomicRanges::seqnames(targetRange))
            if (!identical(as.character(chromosome), target_chromosome) ||
                start != IRanges::start(targetRange)) {
                stop("Numeric score coordinates disagree with `region`.", call. = FALSE)
            }
            if (length(values) != IRanges::width(targetRange)) {
                stop("The number of numeric scores must equal the `region` width.",
                     call. = FALSE)
            }
        }
        positions <- seq.int(start, length.out = length(values))
    }
    has_sequence <- !is.null(sequence)
    if (has_sequence) {
        bases <- .sequence_letters(sequence)
        if (length(values) != length(bases)) {
            stop("Sequence length and number of basewise scores must match.", call. = FALSE)
        }
    } else {
        bases <- rep("N", length(values))
    }
    if (identical(normalization, "maxabs")) {
        denom <- max(abs(values))
        if (denom > 0) values <- values / denom
    }
    list(
        data = data.frame(position = positions, base = bases, score = values,
                          stringsAsFactors = FALSE),
        hasSequence = has_sequence,
        chromosome = as.character(chromosome),
        start = as.integer(start)
    )
}

.score_file_region <- function(sequence, region, chromosome, start, end) {
    if (!is.null(region)) {
        if (!methods::is(region, "GRanges") || length(region) != 1L) {
            stop("`region` must be one GRanges range.", call. = FALSE)
        }
        .validate_one_based_granges(region, "region", requireOne = TRUE)
        return(region)
    }
    if (is.null(chromosome) || is.null(start)) {
        stop("A score file requires `region`, or `chromosome` plus `start` and `end`.",
             call. = FALSE)
    }
    if (is.null(end)) {
        if (is.null(sequence)) {
            stop("`end` is required for a score file when `sequence` is NULL.", call. = FALSE)
        }
        start <- .one_based_coordinate(start, "start")
        end <- start + length(.sequence_letters(sequence)) - 1L
    } else {
        start <- .one_based_coordinate(start, "start")
        end <- .one_based_coordinate(end, "end")
    }
    if (end < start) {
        stop("`end` must be greater than or equal to `start`.", call. = FALSE)
    }
    GenomicRanges::GRanges(
        as.character(chromosome), IRanges::IRanges(as.integer(start), as.integer(end))
    )
}

.score_input_region <- function(scores, region, chromosome, start, end) {
    if (!is.null(region)) {
        if (!methods::is(region, "GRanges") || length(region) != 1L) {
            stop("`region` must be one GRanges range.", call. = FALSE)
        }
        .validate_one_based_granges(region, "region", requireOne = TRUE)
        return(region)
    }
    if (methods::is(scores, "GRanges") && length(scores)) {
        .validate_one_based_granges(scores, "scores")
        chromosomes <- unique(as.character(GenomicRanges::seqnames(scores)))
        if (length(chromosomes) != 1L) {
            stop("`scores` must contain one chromosome.", call. = FALSE)
        }
        return(GenomicRanges::GRanges(
            chromosomes,
            IRanges::IRanges(min(IRanges::start(scores)), max(IRanges::end(scores)))
        ))
    }
    if (is.numeric(scores) && !is.null(chromosome) && !is.null(start)) {
        start <- .one_based_coordinate(start, "start")
        return(GenomicRanges::GRanges(
            as.character(chromosome),
            IRanges::IRanges(start, start + length(scores) - 1L)
        ))
    }
    if (!is.null(chromosome) && !is.null(start) && !is.null(end)) {
        start <- .one_based_coordinate(start, "start")
        end <- .one_based_coordinate(end, "end")
        if (end < start) {
            stop("`end` must be greater than or equal to `start`.", call. = FALSE)
        }
        return(GenomicRanges::GRanges(
            as.character(chromosome), IRanges::IRanges(start, end)
        ))
    }
    stop("Automatic reference retrieval requires `region`, scored GRanges, or numeric scores with `chromosome` and `start`.",
         call. = FALSE)
}

.sequence_from_reference <- function(reference, region, genome = NULL) {
    if (!methods::is(region, "GRanges") || length(region) != 1L) {
        stop("Reference retrieval requires one genomic `region`.", call. = FALSE)
    }
    .validate_one_based_granges(region, "region", requireOne = TRUE)
    chromosome <- as.character(GenomicRanges::seqnames(region))
    if (methods::is(reference, "DNAString")) {
        sequence <- reference
    } else {
        if (methods::is(reference, "SequenceTrack")) {
            sequence_track <- reference
            sequence_track <- Gviz::`chromosome<-`(sequence_track, value = chromosome)
        } else {
            args <- list(sequence = reference, chromosome = chromosome)
            if (!is.null(genome)) args$genome <- genome
            sequence_track <- do.call(Gviz::SequenceTrack, args)
        }
        sequence <- Biostrings::subseq(
            sequence_track, start = IRanges::start(region), end = IRanges::end(region)
        )
    }
    if (!methods::is(sequence, "DNAString") ||
        length(sequence) != IRanges::width(region)) {
        stop("Reference sequence length does not match `region`.", call. = FALSE)
    }
    sequence
}

.import_score_file <- function(path, region, importFunction = NULL) {
    if (!is.null(importFunction)) {
        if (!is.function(importFunction)) stop("`importFunction` must be a function.", call. = FALSE)
        out <- importFunction(path, region)
        if (!methods::is(out, "GRanges")) {
            stop("`importFunction` must return a GRanges object.", call. = FALSE)
        }
        .validate_one_based_granges(out, "imported scores")
        return(out)
    }
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop("Package `rtracklayer` is required for BigWig/bedGraph score files.",
             call. = FALSE)
    }
    if (!file.exists(path)) stop("Score file does not exist: ", path, call. = FALSE)
    out <- rtracklayer::import(path, which = region)
    .validate_one_based_granges(out, "imported scores")
    out
}

.contribution_display_mode <- function(display, span, maxBases, hasSequence) {
    if (display != "auto") return(display)
    if (hasSequence && span <= maxBases) "letters" else "signal"
}

.make_signal_track <- function(data, chromosome, genome, ylim, showAxis,
                               baselineColor, signalType, signalAggregation,
                               signalWindow, signalWindowSize, signalColor,
                               signalFill) {
    args <- list(
        start = data$position, end = data$position,
        data = matrix(data$score, nrow = 1L), chromosome = chromosome,
        genome = if (is.null(genome)) "unknown" else as.character(genome),
        name = "", type = signalType, ylim = ylim, baseline = 0,
        showAxis = showAxis, col = signalColor, fill = signalFill,
        col.histogram = signalColor, fill.histogram = signalFill,
        col.baseline = baselineColor, aggregation = signalAggregation,
        # Gviz::collapseTrack() can drop the matrix dimension of a one-series
        # DataTrack when only some adjacent ranges collapse, which makes the
        # resulting object fail slot validation. Window aggregation already
        # provides the required overview downsampling, so range collapsing is
        # both redundant and unsafe here.
        window = signalWindow, collapse = FALSE
    )
    if (!is.null(signalWindowSize)) args$windowSize <- signalWindowSize
    do.call(Gviz::DataTrack, args)
}

.draw_gviz_signal_track <- function(object, prepare) {
    vars <- methods::slot(object, "variables")
    xscale <- as.numeric(grid::current.viewport()$xscale)
    signal_track <- vars$signalTrack
    Gviz::displayPars(signal_track) <- list(reverseStrand = xscale[1] > xscale[2])
    signal_track <- Gviz::drawGD(
        signal_track, minBase = min(xscale), maxBase = max(xscale), prepare = prepare
    )
    if (prepare) {
        needed <- Gviz::displayPars(signal_track)$neededVerticalSpace
        if (!is.null(needed)) {
            Gviz::displayPars(object) <- list(neededVerticalSpace = needed)
        }
    }
    invisible(object)
}

.draw_contribution_letters <- function(object) {
    vars <- methods::slot(object, "variables")
    xscale <- grid::current.viewport()$xscale
    visible <- sort(as.numeric(xscale))
    keep <- vars$data$position >= visible[1] & vars$data$position <= visible[2]
    dat <- vars$data[keep, , drop = FALSE]
    reverse_axis <- xscale[1] > xscale[2]
    if (reverse_axis) dat$base <- .complement_bases(dat$base)

    grid::pushViewport(grid::viewport(xscale = xscale, yscale = vars$ylim, clip = "on"))
    on.exit(grid::popViewport(), add = TRUE)
    if (vars$showBoundary) {
        .draw_track_left_boundary(
            vars$boundaryColor, vars$boundaryWidth, vars$boundaryInset
        )
    }
    grid::grid.lines(
        x = grid::unit(xscale, "native"), y = grid::unit(c(0, 0), "native"),
        gp = grid::gpar(col = vars$baselineColor, lwd = 0.7)
    )
    if (nrow(dat)) {
        for (i in seq_len(nrow(dat))) {
            b <- dat$base[i]
            col <- vars$colors[if (b %in% names(vars$colors)) b else "N"]
            .draw_letter_glyph(
                b, dat$position[i] - 0.5, dat$position[i] + 0.5,
                0, dat$score[i], col,
                alpha = if (dat$score[i] < 0) vars$negativeAlpha else 1,
                font = vars$font
            )
        }
    }
    invisible(object)
}

.draw_score_sequence_track <- function(object, prepare = FALSE) {
    vars <- methods::slot(object, "variables")
    xscale <- grid::current.viewport()$xscale
    span <- floor(diff(sort(as.numeric(xscale)))) + 1L
    display <- .contribution_display_mode(
        vars$display, span, vars$maxBases, vars$hasSequence
    )
    if (display == "signal") {
        .draw_gviz_signal_track(object, prepare)
        if (!prepare && vars$showBoundary) {
            .draw_track_left_boundary(
                vars$boundaryColor, vars$boundaryWidth, vars$boundaryInset
            )
        }
        return(invisible(object))
    }
    if (!prepare) .draw_contribution_letters(object)
    invisible(object)
}

#' Score-scaled sequence track for Gviz
#'
#' Construct a track whose nucleotide glyph heights encode signed,
#' base-resolution scores. The result inherits from both `Gviz::CustomTrack`
#' and `Gviz::NumericTrack`: the sequence glyphs use Gviz's custom Grid callback
#' while Gviz itself draws the quantitative y axis. At wider views the track
#' delegates rendering to a normal `Gviz::DataTrack`, including Gviz's window
#' aggregation. Supply `from`, `to`, and `chromosome` to `Gviz::plotTracks()`;
#' Gviz does not infer plotting limits from tracks with a custom callback.
#'
#' @param sequence One `Biostrings::DNAString` or character sequence. May be
#'   `NULL` for a signal-only track.
#' @importClassesFrom Biostrings DNAString
#' @param reference Optional Gviz-compatible reference source used when
#'   `sequence` is `NULL`: a `Gviz::SequenceTrack`, `BSgenome`, `DNAStringSet`,
#'   indexed FASTA, or 2bit file.
#' @param genome Optional genome identifier passed to `Gviz::SequenceTrack`
#'   and the overview `Gviz::DataTrack`.
#' @param scores Numeric basewise scores, a single-chromosome `GRanges`, or a
#'   BigWig/bedGraph file path. A range wider than one base represents the same
#'   score at every covered base.
#' @param chromosome Chromosome for numeric scores; inferred from `GRanges`.
#' @param start One-based genomic start for numeric scores; inferred from
#'   `GRanges`.
#' @param end One-based genomic end for score files. It can be inferred from a
#'   supplied sequence.
#' @param region Optional one-range `GRanges` selecting a score-file interval.
#' @param importFunction Optional function with arguments `(file, selection)`
#'   returning scored `GRanges`; useful for external large-BigWig readers.
#' @param scoreColumn Score metadata column for `GRanges` input. Numeric factor
#'   and character columns are parsed strictly.
#' @param name Gviz track title.
#' @param normalization Either `"none"` or within-track `"maxabs"` scaling.
#' @param ylim Numeric y limits. By default symmetric limits are derived from
#'   the displayed scores.
#' @param missing How uncovered positions inside a `GRanges` span are handled.
#' @param colors Named nucleotide color vector.
#' @param maxBases Largest plotting viewport that will render letters.
#' @param display Rendering mode: automatic letters/signal switching, forced
#'   letters, or forced signal.
#' @param font A font name returned by `ggseqlogo::list_fonts(FALSE)`.
#' @param showAxis Let Gviz draw its standard numeric y axis in the title panel.
#'   Standard Gviz display parameters such as `col.axis`, `cex.axis`, and
#'   `yTicksAt` can be supplied through `...` or changed with
#'   `Gviz::displayPars()`.
#' @param baselineColor Color of the zero baseline.
#' @param showBoundary Draw a vertical boundary at the left edge of the data
#'   panel, aligned with the plotting area of other Gviz tracks.
#' @param boundaryColor,boundaryWidth Color and line width of that boundary.
#' @param boundaryInset Fraction of track height inset at both ends of the
#'   boundary, leaving a visible gap between adjacent track segments.
#' @param negativeAlpha Alpha multiplier for negative glyphs.
#' @param signalType Gviz DataTrack type used for the overview.
#' @param signalAggregation Gviz aggregation function for overview windows.
#' @param signalWindow Gviz window setting, normally `"auto"`.
#' @param signalWindowSize Optional explicit Gviz window size.
#' @param signalColor,signalFill Overview signal border and fill colors.
#' @param ... Additional Gviz display parameters.
#' @return An object inheriting from both `Gviz::CustomTrack` and
#'   `Gviz::NumericTrack`.
#' @export
ScoreSequenceTrack <- function(sequence = NULL, scores, reference = NULL,
                               genome = NULL, chromosome = NULL,
                               start = NULL, end = NULL, region = NULL,
                               importFunction = NULL,
                               scoreColumn = "score",
                               name = "Score sequence",
                               normalization = c("none", "maxabs"),
                               ylim = NULL,
                               missing = c("error", "zero"),
                               colors = .nucleotide_colors(),
                               maxBases = 250L,
                               display = c("auto", "letters", "signal"),
                               font = "roboto_medium", showAxis = TRUE,
                               baselineColor = "#555555",
                               showBoundary = TRUE,
                               boundaryColor = "#222222",
                               boundaryWidth = 0.7,
                               boundaryInset = 0.06,
                               negativeAlpha = 0.82,
                               signalType = "histogram",
                               signalAggregation = "mean",
                               signalWindow = "auto",
                               signalWindowSize = NULL,
                               signalColor = "#59636B",
                               signalFill = "#7B858D", ...) {
    normalization <- match.arg(normalization)
    display <- match.arg(display)
    missing_was_default <- missing(missing)
    missing <- match.arg(missing)
    colors <- .validate_colors(colors)
    showAxis <- .validate_flag(showAxis, "showAxis")
    maxBases <- .validate_scalar_number(
        maxBases, "maxBases", minimum = 0, whole = TRUE,
        minimumInclusive = FALSE
    )
    negativeAlpha <- .validate_scalar_number(
        negativeAlpha, "negativeAlpha", minimum = 0, maximum = 1
    )
    baselineColor <- .validate_color(baselineColor, "baselineColor")
    signalColor <- .validate_color(signalColor, "signalColor")
    signalFill <- .validate_color(signalFill, "signalFill")
    if (!is.null(signalWindowSize)) {
        signalWindowSize <- .validate_scalar_number(
            signalWindowSize, "signalWindowSize", minimum = 0, whole = TRUE,
            minimumInclusive = FALSE
        )
    }
    if (!is.null(sequence) && !is.null(reference)) {
        stop("Supply only one of `sequence` and `reference`.", call. = FALSE)
    }
    if (!is.null(region)) {
        .validate_one_based_granges(region, "region", requireOne = TRUE)
        region_chromosome <- as.character(GenomicRanges::seqnames(region))
        if (!is.null(chromosome) &&
            !identical(as.character(chromosome), region_chromosome)) {
            stop("`chromosome` disagrees with `region`.", call. = FALSE)
        }
        if (!is.null(start) &&
            .one_based_coordinate(start, "start") != IRanges::start(region)) {
            stop("`start` disagrees with `region`.", call. = FALSE)
        }
        if (!is.null(end) &&
            .one_based_coordinate(end, "end") != IRanges::end(region)) {
            stop("`end` disagrees with `region`.", call. = FALSE)
        }
    }
    if (is.numeric(scores) && is.null(region) && !is.null(end)) {
        numeric_start <- .one_based_coordinate(start, "start")
        numeric_end <- .one_based_coordinate(end, "end")
        if (numeric_end != numeric_start + length(scores) - 1L) {
            stop("For numeric scores, `end` must equal `start + length(scores) - 1`.",
                 call. = FALSE)
        }
    }
    if (methods::is(scores, "GRanges") && !is.null(end) && length(scores) &&
        .one_based_coordinate(end, "end") != max(IRanges::end(scores))) {
        stop("`end` disagrees with the GRanges end coordinate.", call. = FALSE)
    }
    if (is.null(sequence) && !is.null(reference)) {
        reference_region <- .score_input_region(
            scores, region, chromosome, start, end
        )
        sequence <- .sequence_from_reference(reference, reference_region, genome)
        if (is.null(region)) region <- reference_region
    }
    target_range <- region
    score_source <- NULL
    if (is.character(scores) && length(scores) == 1L) {
        target_range <- .score_file_region(sequence, region, chromosome, start, end)
        score_source <- normalizePath(scores, winslash = "/", mustWork = TRUE)
        scores <- .import_score_file(score_source, target_range, importFunction)
        if (missing_was_default) missing <- "zero"
    } else if (!is.null(region)) {
        if (!methods::is(region, "GRanges") || length(region) != 1L) {
            stop("`region` must be one GRanges range.", call. = FALSE)
        }
        chromosome <- as.character(GenomicRanges::seqnames(region))
        start <- IRanges::start(region)
    }
    prepared <- .prepare_contribution_data(
        sequence, scores, chromosome, start, scoreColumn, missing, normalization,
        targetRange = target_range
    )
    data <- prepared$data
    if (is.null(ylim)) ylim <- contributionYlim(data$score, padding = 0.04)
    if (!is.numeric(ylim) || length(ylim) != 2L || any(!is.finite(ylim)) || ylim[1] >= ylim[2]) {
        stop("`ylim` must contain two increasing finite numbers.", call. = FALSE)
    }
    .validate_track_boundary(
        showBoundary, boundaryColor, boundaryWidth, boundaryInset
    )
    .load_glyph_font(font)
    signal_track <- .make_signal_track(
        data, prepared$chromosome, genome, as.numeric(ylim), showAxis,
        baselineColor, signalType, signalAggregation, signalWindow,
        signalWindowSize, signalColor, signalFill
    )

    vars <- list(
        data = data,
        chromosome = prepared$chromosome, source = score_source,
        ylim = as.numeric(ylim), colors = colors, maxBases = maxBases,
        display = display, hasSequence = prepared$hasSequence, font = font,
        showAxis = showAxis, baselineColor = baselineColor,
        showBoundary = isTRUE(showBoundary), boundaryColor = boundaryColor,
        boundaryWidth = boundaryWidth, boundaryInset = boundaryInset,
        negativeAlpha = negativeAlpha, signalTrack = signal_track,
        reference = reference, genome = genome
    )
    custom_track <- Gviz::CustomTrack(
        plottingFunction = function(GdObject, prepare = FALSE, ...) {
            .draw_score_sequence_track(GdObject, prepare)
        },
        variables = vars, name = name, ylim = as.numeric(ylim),
        showAxis = showAxis, ...
    )
    .as_gviz_numeric_custom_track(
        custom_track, data, prepared$chromosome, genome
    )
}

#' @rdname ScoreSequenceTrack
#' @export
DynSeqTrack <- function(...) ScoreSequenceTrack(...)
