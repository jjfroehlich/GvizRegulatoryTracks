.validate_flag <- function(x, argument) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
    }
    x
}

.validate_optional_column_name <- function(x, argument) {
    if (is.null(x)) return(NULL)
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
        stop("`", argument,
             "` must be NULL or one non-empty column name.", call. = FALSE)
    }
    x
}

.validate_scalar_number <- function(x, argument, minimum = -Inf,
                                    maximum = Inf, whole = FALSE,
                                    minimumInclusive = TRUE,
                                    maximumInclusive = TRUE) {
    if (isTRUE(whole)) maximum <- min(maximum, .Machine$integer.max)
    valid <- is.numeric(x) && length(x) == 1L && is.finite(x)
    if (valid && isTRUE(whole)) valid <- x == floor(x)
    if (valid) {
        valid <- if (minimumInclusive) x >= minimum else x > minimum
    }
    if (valid) {
        valid <- if (maximumInclusive) x <= maximum else x < maximum
    }
    if (!valid) {
        qualifier <- if (whole) " whole number" else " number"
        stop("`", argument, "` must be one finite", qualifier,
             " in the supported range.", call. = FALSE)
    }
    if (whole) as.integer(x) else as.numeric(x)
}

.validate_color <- function(x, argument) {
    if (!is.character(x) || length(x) != 1L || is.na(x)) {
        stop("`", argument, "` must be one valid color.", call. = FALSE)
    }
    tryCatch(grDevices::col2rgb(x), error = function(e) {
        stop("`", argument, "` must be one valid color.", call. = FALSE)
    })
    x
}

.validate_color_vector <- function(x, argument, minimumLength = 1L) {
    if (!is.character(x) || length(x) < minimumLength || anyNA(x)) {
        stop("`", argument, "` must contain valid colors.", call. = FALSE)
    }
    tryCatch(grDevices::col2rgb(x), error = function(e) {
        stop("`", argument, "` must contain valid colors.", call. = FALSE)
    })
    x
}

.score_values <- function(x, argument, allowNumericNonFinite = FALSE) {
    text_input <- is.factor(x) || is.character(x)
    if (is.factor(x)) x <- as.character(x)
    if (!is.numeric(x) && !is.character(x)) {
        stop("`", argument, "` must contain numeric values or numeric text.",
             call. = FALSE)
    }

    if (text_input) {
        missing_input <- is.na(x)
        values <- suppressWarnings(as.numeric(x))
        malformed <- !missing_input & (is.na(values) | is.infinite(values))
        if (any(malformed)) {
            stop("`", argument, "` contains non-numeric or non-finite text.",
                 call. = FALSE)
        }
    } else {
        values <- as.numeric(x)
    }

    if (!allowNumericNonFinite && any(!is.finite(values))) {
        stop("`", argument, "` must contain only finite, non-missing scores.",
             call. = FALSE)
    }
    values
}
