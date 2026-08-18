#' @noRd
#' @importClassesFrom Gviz CustomTrack NumericTrack
methods::setClass(
    "GvizRegulatoryNumericTrack",
    contains = c("CustomTrack", "NumericTrack")
)

.as_gviz_numeric_custom_track <- function(track, data, chromosome, genome) {
    track_range <- GenomicRanges::GRanges(
        chromosome,
        IRanges::IRanges(min(data$position), max(data$position))
    )
    out <- methods::new(
        "GvizRegulatoryNumericTrack",
        plottingFunction = methods::slot(track, "plottingFunction"),
        variables = methods::slot(track, "variables"),
        name = methods::slot(track, "name"),
        imageMap = methods::slot(track, "imageMap"),
        range = track_range,
        chromosome = as.character(chromosome),
        genome = if (is.null(genome)) "unknown" else as.character(genome)
    )
    # The inherited Gviz initializer treats a named `dp` argument as another
    # display parameter. Assigning the slot after initialization preserves the
    # exact CustomTrack DisplayPars object instead of nesting it under `$dp`.
    methods::slot(out, "dp") <- methods::slot(track, "dp")
    out
}
