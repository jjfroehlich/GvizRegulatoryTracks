.nucleotide_colors <- function() {
    c(A = "#2E8B57", C = "#2878B5", G = "#E69F00", T = "#D64B40", N = "#777777")
}

.glyph_cache <- new.env(parent = emptyenv())

.draw_track_left_boundary <- function(color = "#222222", width = 0.7,
                                      inset = 0.06) {
    grid::grid.lines(
        x = grid::unit(c(0, 0), "npc"),
        y = grid::unit(c(inset, 1 - inset), "npc"),
        gp = grid::gpar(col = color, lwd = width)
    )
    invisible(NULL)
}

.validate_track_boundary <- function(show, color, width, inset) {
    if (!is.logical(show) || length(show) != 1L || is.na(show)) {
        stop("`showBoundary` must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.character(color) || length(color) != 1L || is.na(color)) {
        stop("`boundaryColor` must be one valid color.", call. = FALSE)
    }
    tryCatch(grDevices::col2rgb(color), error = function(e) {
        stop("`boundaryColor` must be one valid color.", call. = FALSE)
    })
    if (!is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0) {
        stop("`boundaryWidth` must be one positive number.", call. = FALSE)
    }
    if (!is.numeric(inset) || length(inset) != 1L || !is.finite(inset) ||
        inset < 0 || inset >= 0.5) {
        stop("`boundaryInset` must be one number from zero up to, but not including, 0.5.",
             call. = FALSE)
    }
    invisible(NULL)
}

.load_glyph_font <- function(font = "roboto_medium") {
    available <- ggseqlogo::list_fonts(FALSE)
    if (!font %in% available) {
        stop("Unknown ggseqlogo font `", font, "`. Available fonts: ",
             paste(available, collapse = ", "), call. = FALSE)
    }
    if (!exists(font, envir = .glyph_cache, inherits = FALSE)) {
        path <- system.file("extdata", paste0(font, ".font"), package = "ggseqlogo")
        if (!nzchar(path)) stop("Could not locate ggseqlogo font `", font, "`.", call. = FALSE)
        dat <- readRDS(path)
        dat <- dat[dat$letter %in% c("A", "C", "G", "T", "N"), , drop = FALSE]
        glyphs <- split(dat[order(dat$letter, dat$order), c("x", "y")], dat$letter)
        assign(font, glyphs, envir = .glyph_cache)
    }
    get(font, envir = .glyph_cache, inherits = FALSE)
}

.validate_colors <- function(colors) {
    if (is.null(names(colors)) || !all(c("A", "C", "G", "T") %in% names(colors))) {
        stop("`colors` must be a named vector containing A, C, G, and T.", call. = FALSE)
    }
    if (!"N" %in% names(colors)) {
        colors <- c(colors, N = "#777777")
    }
    colors
}

.draw_letter_glyph <- function(base, xleft, xright, ybottom, ytop, color,
                               alpha = 1, xpad = 0.04, ypad = 0,
                               font = "roboto_medium", borderColor = NA,
                               borderWidth = 0.6) {
    if (!is.finite(ybottom) || !is.finite(ytop) || identical(ybottom, ytop)) {
        return(invisible(NULL))
    }
    base <- toupper(base)
    glyphs <- .load_glyph_font(font)
    if (!base %in% names(glyphs)) {
        grid::grid.rect(
            x = grid::unit((xleft + xright) / 2, "native"),
            y = grid::unit((ybottom + ytop) / 2, "native"),
            width = grid::unit(abs(xright - xleft) * (1 - 2 * xpad), "native"),
            height = grid::unit(abs(ytop - ybottom) * (1 - 2 * ypad), "native"),
            gp = grid::gpar(
                fill = grDevices::adjustcolor(color, alpha.f = alpha),
                col = if (is.na(borderColor)) NA else
                    grDevices::adjustcolor(borderColor, alpha.f = alpha),
                lwd = borderWidth
            )
        )
        return(invisible(NULL))
    }

    dx <- xright - xleft
    dy <- ytop - ybottom
    glyph <- glyphs[[base]]
    xx <- xleft + dx * (xpad + glyph$x * (1 - 2 * xpad))
    yy <- ybottom + dy * (ypad + glyph$y * (1 - 2 * ypad))
    grid::grid.path(
        x = grid::unit(xx, "native"), y = grid::unit(yy, "native"),
        rule = "evenodd",
        gp = grid::gpar(
            fill = grDevices::adjustcolor(color, alpha.f = alpha),
            col = if (is.na(borderColor)) NA else
                grDevices::adjustcolor(borderColor, alpha.f = alpha),
            lwd = borderWidth
        )
    )
    invisible(NULL)
}

.complement_bases <- function(x) {
    map <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
    out <- unname(map[toupper(x)])
    out[is.na(out)] <- "N"
    out
}
