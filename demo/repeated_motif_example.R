#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_path <- function() {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(arg)) return(normalizePath(sub("^--file=", "", arg[1]), winslash = "/"))
    normalizePath("demo/repeated_motif_example.R", winslash = "/", mustWork = FALSE)
}

package_root <- normalizePath(file.path(dirname(script_path()), ".."), winslash = "/")
if (file.exists(file.path(package_root, "DESCRIPTION")) &&
        requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(package_root, quiet = TRUE)
} else if (!requireNamespace("GvizRegulatoryTracks", quietly = TRUE)) {
    stop("Install devtools or GvizRegulatoryTracks.")
}

suppressPackageStartupMessages({
    library(Biostrings)
    library(GenomicRanges)
    library(Gviz)
    library(GvizRegulatoryTracks)
    library(grid)
    library(rtracklayer)
})

opened_null_device <- grDevices::dev.cur() == 1L
if (opened_null_device) grDevices::pdf(NULL)

fixture_dir <- file.path(package_root, "inst", "extdata", "example_repeated_motif")
if (!dir.exists(fixture_dir)) {
    fixture_dir <- system.file(
        "extdata", "example_repeated_motif", package = "GvizRegulatoryTracks"
    )
}
output_dir <- Sys.getenv("GVIZ_DEMO_OUTPUT_DIR", file.path(getwd(), "demo_output"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- file.path(fixture_dir, c(
    "locus_sequence.fa", "atac_seq_signal.bw",
    "chrombpnet_count_contribution.bw", "tfmodisco_seqlet.bed",
    "tfmodisco_pattern.meme", "fimo_hits.bed", "fimo_motifs.meme"
))
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing repeated-motif fixtures:\n", paste(missing, collapse = "\n"))

chromosome <- "chr11"
view_from <- 103102780L
view_to <- 103102890L
view_region <- GRanges(chromosome, IRanges(view_from, view_to))
nucleotide_colors <- c(
    A = "#2E8B57", C = "#2878B5", G = "#E69F00", T = "#D64B40", N = "#777777"
)
axis_color <- "#666666"
coordinate_color <- "#999999"
axis_cex <- 0.65
boundary_color <- "#222222"
boundary_inset <- 0.06

boundary_track <- CustomTrack(
    plottingFunction = function(GdObject, prepare = FALSE, ...) {
        if (!prepare) {
            grid::grid.lines(
                x = grid::unit(c(0, 0), "npc"),
                y = grid::unit(c(boundary_inset, 1 - boundary_inset), "npc"),
                gp = grid::gpar(col = boundary_color, lwd = 0.7)
            )
        }
        invisible(GdObject)
    },
    name = ""
)

atac_gr <- rtracklayer::import(
    file.path(fixture_dir, "atac_seq_signal.bw"), which = view_region
)
stopifnot(sum(width(atac_gr)) == width(view_region))
atac_track <- DataTrack(
    atac_gr, genome = "mm10", chromosome = chromosome,
    name = "ATAC-seq\n[fragment CPM]", type = "histogram", baseline = 0,
    ylim = c(0, max(atac_gr$score) * 1.05),
    col = "#5B6570", fill = "#737E89", col.baseline = "#777777",
    col.axis = axis_color, cex.axis = axis_cex, lwd = 0.45
)
atac_track <- OverlayTrack(
    trackList = list(atac_track, boundary_track),
    name = "ATAC-seq\n[fragment CPM]"
)

sequence <- as.character(Biostrings::readDNAStringSet(
    file.path(fixture_dir, "locus_sequence.fa")
)[[1]])
stopifnot(nchar(sequence) == width(view_region))

dna_track <- CustomTrack(
    plottingFunction = function(GdObject, prepare = FALSE, ...) {
        if (prepare) return(invisible(GdObject))
        vars <- methods::slot(GdObject, "variables")
        xscale <- grid::current.viewport()$xscale
        keep <- vars$position >= min(xscale) & vars$position <= max(xscale)
        grid::pushViewport(grid::viewport(xscale = xscale, yscale = c(0, 1), clip = "on"))
        grid::grid.lines(
            x = grid::unit(c(0, 0), "npc"),
            y = grid::unit(c(boundary_inset, 1 - boundary_inset), "npc"),
            gp = grid::gpar(col = boundary_color, lwd = 0.7)
        )
        grid::grid.text(
            vars$base[keep], x = grid::unit(vars$position[keep], "native"),
            y = grid::unit(0.5, "native"),
            gp = grid::gpar(col = vars$colors[vars$base[keep]], fontsize = 7)
        )
        grid::popViewport()
        invisible(GdObject)
    },
    variables = list(
        position = view_from:view_to,
        base = strsplit(sequence, "", fixed = TRUE)[[1]],
        colors = nucleotide_colors
    ),
    name = "DNA sequence"
)

contribution_track <- DynSeqTrack(
    sequence = sequence,
    scores = file.path(fixture_dir, "chrombpnet_count_contribution.bw"),
    region = view_region, genome = "mm10",
    name = "ChromBPNet\ncount contribution\n[contribution score]",
    colors = nucleotide_colors, showAxis = TRUE,
    yTicksAt = c(-0.04, 0, 0.04), col.axis = axis_color,
    cex.axis = axis_cex,
    maxBases = 250L
)

tfmodisco_hits <- rtracklayer::import(file.path(fixture_dir, "tfmodisco_seqlet.bed"))
tfmodisco_meme <- file.path(fixture_dir, "tfmodisco_pattern.meme")
hits <- rtracklayer::import(file.path(fixture_dir, "fimo_hits.bed"))
hits$relative_fimo_score <- hits$score / 1000
fimo_meme <- file.path(fixture_dir, "fimo_motifs.meme")
stopifnot(
    identical(as.character(tfmodisco_hits$name), "pattern_0"),
    identical(start(tfmodisco_hits), 103102793L),
    identical(end(tfmodisco_hits), 103102842L),
    identical(as.character(strand(tfmodisco_hits)), "-"),
    identical(as.character(hits$name), rep("CTCF.H14CORE.0.P.B", 2L)),
    identical(start(hits), c(103102812L, 103102851L)),
    identical(end(hits), c(103102831L, 103102870L)),
    identical(as.character(strand(hits)), c("+", "+")),
    identical(hits$relative_fimo_score, c(0.940, 0.884))
)

tfmodisco_track <- MotifLogoTrack(
    tfmodisco_hits, tfmodisco_meme,
    name = "TF-MoDISco pattern\n[nucleotide]",
    motifColorBy = "nucleotide", colors = nucleotide_colors,
    scoreAesthetic = "none", showLabels = TRUE, showStrand = TRUE,
    motifHeight = 8, showAxis = FALSE, maxBases = 150L
)

fimo_track <- MotifLogoTrack(
    hits, fimo_meme, name = "Repeated motif matches\n[relative FIMO score]",
    scoreColumn = "relative_fimo_score", scoreAesthetic = "fill",
    scoreLimits = "view", scoreColor = "#111111",
    scoreBrightness = c(0.82, 0), scoreLegendTitle = "Relative FIMO score",
    motifHeight = 8, showLabels = TRUE, showStrand = TRUE,
    showAxis = FALSE, maxBases = 150L
)

tracks <- list(
    GenomeAxisTrack(
        name = "", col = coordinate_color, fontcolor = coordinate_color,
        fontsize = 8, fontface = 1
    ),
    atac_track, dna_track, contribution_track,
    tfmodisco_track, fimo_track
)
sizes <- c(0.35, 0.85, 0.48, 0.85, 0.85, 0.85)

draw_demo <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
        x = 0.50, y = 0.50, width = 0.98, height = 0.94
    ))
    Gviz::plotTracks(
        tracks, chromosome = chromosome, from = view_from, to = view_to,
        sizes = sizes,
        main = "Example 2: repeated matches of one motif",
        background.title = "transparent", col.border.title = "transparent",
        fontcolor.title = "#222222", fontface.title = 1,
        fontsize.title = 9, cex.title = 0.72,
        fontsize = 8, cex.axis = axis_cex, col.axis = axis_color,
        title.width = 1.15, margin = c(2, 6, 2, 4), innerMargin = 1,
        cex.main = 0.78, add = TRUE
    )
    grid::popViewport()
}

pdf_file <- file.path(output_dir, "repeated_motif_example.pdf")
png_file <- file.path(output_dir, "repeated_motif_example.png")
grDevices::pdf(pdf_file, width = 10, height = 7.2, useDingbats = FALSE)
draw_demo()
grDevices::dev.off()
grDevices::png(png_file, width = 2000, height = 1440, res = 200,
               type = "cairo", bg = "white")
draw_demo()
grDevices::dev.off()

if (any(file.info(c(pdf_file, png_file))$size <= 0)) stop("Repeated-motif example failed.")
message("Wrote:\n  ", normalizePath(pdf_file, winslash = "/"),
        "\n  ", normalizePath(png_file, winslash = "/"))
if (opened_null_device) grDevices::dev.off()
