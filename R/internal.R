# Start of internal.R ##############################################################################

# get_indices --------------------------------------------------------------------------------------
#' Get indices of object elements.
#'
#' Get the indices of \code{x}, as constrained by the \code{requested} parameter. If using this
#' function without \code{requested} constraints, consider using the faster primitive R function
#' \code{seq_along} instead.
#'
#' @param x An object with elements that are accessible by an index.
#' @param requested A character vector of names, a logical vector of the same length as \code{x}, or
#' a numeric vector containing indices of \code{x}. If this parameter is not specified, all indices
#' are returned.
#' @param strict Option indicating that \code{requested}, if specified, must request indices that
#' are unique and in the same order as the corresponding elements of \code{x}.
#'
#' @return Integer vector containing indices of \code{x}.
#'
#' @keywords internal
#' @rdname get_indices
get_indices <- function(x, requested=NULL, strict=FALSE) {

    stopifnot( isTRUE(strict) || identical(FALSE, strict) )

    object.length <- length(x)

    if ( ! is_single_nonnegative_whole_number(object.length) ) {
        stop("cannot get object indices - object does not have length")
    }

    indices <- seq_along(x)

    if ( ! is.null(requested) ) {

        if ( is.numeric(requested) ) {

            nonintegers <- requested[ ! is_whole_number(requested) ]
            if ( length(nonintegers) > 0L ) {
                stop("requested indices are not integers - '", toString(nonintegers), "'")
            }

            exrange <- requested[ ! requested %in% indices ]
            if ( length(exrange) > 0L ) {
                stop("requested indices out of range - '", toString(exrange), "'")
            }

            indices <- indices[requested]

        } else if ( is.logical(requested) ) {

            if ( length(requested) != object.length ) {
                stop("cannot resolve indices by logical vector - length mismatch")
            }

            indices <- unname(which(requested))

        } else if ( is.character(requested) ) {

            object.names <- names(x)

            if ( anyNA(requested) ) {
                stop("cannot resolve indices by name - requested names are incomplete")
            }

            if ( is.null(object.names) ) {
                stop("cannot resolve indices by name - no object names found")
            }

            if ( anyDuplicated(object.names) ) {
                stop("cannot resolve indices by name - duplicate object names found")
            }

            unfound <- requested[ ! requested %in% object.names ]
            if ( length(unfound) > 0L ) {
                stop("requested names not found - '", toString(unfound), "'")
            }

            indices <- match(requested, object.names)

        } else {
            stop("requested indices must be specified by index, logical mask, or name")
        }

        if (strict) { # NB: also ensures no duplicates
            if ( is.unsorted(indices, strictly=TRUE) ) {
                stop("requested indices not specified in strictly increasing order")
            }
        }
    }

    return(indices)
}

# get_lod_col_index --------------------------------------------------------------------------------
#' Get LOD column index.
#'
#' @param x An \pkg{R/qtl} \code{scanone}, \code{scanoneperm}, or \code{summary.scanoneperm} object.
#' @param lodcolumn This parameter indicates for which LOD column an index should be returned. This
#' must be either a LOD column name or an index \emph{with respect to the set of LOD columns}. If no
#' LOD column is specified and one such column is found, the index of that column is returned by
#' default; otherwise a LOD column must be specified.
#'
#' @return LOD column index.
#'
#' @keywords internal
#' @rdname get_lod_col_index
get_lod_col_index <- function(x, lodcolumn=NULL) {

    lodcol.indices <- get_lod_col_indices(x, lodcolumns=lodcolumn)

    if ( length(lodcol.indices) > 1L ) {
        stop("object has multiple LOD columns - please choose one")
    } else if ( length(lodcol.indices) == 0L ) {
        stop("no LOD column found")
    }

    return(lodcol.indices[1L])
}

# get_lod_col_indices ------------------------------------------------------------------------------
#' Get LOD column indices.
#'
#' @param x An \pkg{R/qtl} \code{scanone}, \code{scanoneperm}, or \code{summary.scanoneperm} object.
#' @param lodcolumns This parameter indicates for which LOD columns indices should be returned. The
#' specified LOD columns must be a character vector of LOD column names, a logical vector with
#' length equal to the number of LOD columns, or a numeric vector of indices \emph{with respect
#' to the set of LOD columns}. If no LOD columns are specified, all LOD column indices are returned.
#' @param strict Option indicating that \code{lodcolumns}, if specified, must request LOD column
#' indices that are unique and in the same order as the LOD columns of \code{x}.
#'
#' @return Vector of LOD column indices.
#'
#' @keywords internal
#' @rdname get_lod_col_indices
get_lod_col_indices <- function(x, lodcolumns=NULL, strict=FALSE) {
    UseMethod('get_lod_col_indices', x)
}

# get_lod_col_indices.scanone ----------------------------------------------------------------------
#' @method get_lod_col_indices scanone
#' @rdname get_lod_col_indices
get_lod_col_indices.scanone <- function(x, lodcolumns=NULL, strict=FALSE) {
    stopifnot( identical(colnames(x)[1:2], c('chr', 'pos')) )
    available <- seq_len(ncol(x))[-(1:2)]
    names(available) <- colnames(x)[available]
    resolved <- get_indices(available, requested=lodcolumns, strict=strict)
    indices <- unname(available[resolved])
    return(indices)
}

# get_lod_col_indices.scanoneperm ------------------------------------------------------------------
#' @method get_lod_col_indices scanoneperm
#' @rdname get_lod_col_indices
get_lod_col_indices.scanoneperm <- function(x, lodcolumns=NULL, strict=FALSE) {
    stopifnot( ncol(x) >= 1L )
    available <- seq_len(ncol(x))
    names(available) <- colnames(x)
    resolved <- get_indices(available, requested=lodcolumns, strict=strict)
    indices <- unname(available[resolved])
    return(indices)
}

# get_lod_col_indices.summary.scanoneperm ----------------------------------------------------------
#' @export
#' @method get_lod_col_indices summary.scanoneperm
#' @rdname get_lod_col_indices
get_lod_col_indices.summary.scanoneperm <- function(x, lodcolumns=NULL, strict=FALSE) {
    stopifnot( ncol(x) >= 1L )
    available <- seq_len(ncol(x))
    names(available) <- colnames(x)
    resolved <- get_indices(available, requested=lodcolumns, strict=strict)
    indices <- unname(available[resolved])
    return(indices)
}

# get_run_index_list -------------------------------------------------------------------------------
#' Get index list of successive runs in a vector.
#'
#' @param x A vector.
#'
#' @return List of integer vectors, each containing indices for a run of repeated values in
#' \code{x}. Each list element takes its name from the corresponding repeated value. Returns
#' an empty list if \code{x} is of length zero.
#'
#' @keywords internal
#' @rdname get_run_index_list
get_run_index_list <- function(x, na.rm=FALSE) {

    stopifnot( is.vector(x) )
    stopifnot( isTRUE(na.rm) || identical(FALSE, na.rm) )

    if ( length(x) > 0L ) {

        # Get run-length encoding of vector.
        runs <- rle(x)

        # Set run names from RLE values.
        run.names <- runs$values

        # Get number of runs in RLE.
        num.runs <- unique(lengths(runs))

        # Get last index of each run.
        J <- cumsum(runs$lengths)

        # Get first index of each run.
        if ( num.runs > 1L ) {
            I <- c(1L, sapply(J[1:(length(J)-1)], function(j) j + 1))
        } else {
            I <- 1L
        }

        # Remove NA values, if specified.
        if (na.rm) {
            mask <- ! is.na(runs$values)
            run.names <- run.names[mask]
            I <- I[mask]
            J <- J[mask]
        }

        # Set index list from run index ranges.
        index.list <- mapply(function(i, j) i:j, I, J, SIMPLIFY=FALSE)

        # Set names of index list from run values.
        names(index.list) <- run.names

    } else {

        index.list <- list()
    }

    return(index.list)
}

# infer_map_precision ------------------------------------------------------------------------------
#' Infer precision of map positions.
#'
#' @param x An \pkg{R/qtl} \code{scanone} or \code{map} object.
#' @param tol Tolerance for map position equality.
#'
#' @return Inferred precision of map in \code{x}. Map precision is taken from the map step size, if
#' the majority of the map steps in \code{x} are the same. Otherwise, map precision is set from the
#' smallest distance between any two adjacent loci.
#'
#' @keywords internal
#' @rdname infer_map_precision
infer_map_precision <- function(x, tol=.Machine$double.eps^0.5) {

    stopifnot( any( c('map', 'scanone') %in% class(x) ) )
    stopifnot( is_single_nonnegative_number(tol) )

    if ( 'scanone' %in% class(x) ) {
        map.steps <- unlist(lapply(unique(x$chr), function(map.seq)
                            diff(x[x$chr == map.seq, 'pos'])))
    } else {  # 'map' %in% class(x)
        map.steps <- unlist(unname(lapply(x, function(map.pos) diff(as.numeric(map.pos)))))
    }

    # Get frequency table of map steps.
    step.freqs <- table(map.steps)

    # Get numeric step sizes.
    step.sizes <- as.numeric(names(step.freqs))

    # Get differences between step sizes.
    step.diffs <- diff(step.sizes)

    # Group step-size values that are very similar.
    size.groups <- vector('list')
    i <- 1L
    while ( i <= length(step.freqs) ) {

        j <- i

        while ( j < length(step.freqs) && step.diffs[j] < tol ) {
            j <- j + 1L
        }

        size.groups <- append(size.groups, list(i:j))

        i <- j + 1L
    }

    # Merge similar step values.
    merged.freqs <- integer(length=length(size.groups))
    for ( i in seq_along(size.groups) ) {
        g <- unlist(size.groups[i])
        merged.freqs[i] <- sum(step.freqs[g])
        names(merged.freqs)[i] <- names(sort(step.freqs[g], decreasing=TRUE))[1L]
    }

    # Sort steps by decreasing frequency.
    sorted.freqs <- sort(merged.freqs, decreasing=TRUE)

    # If the most frequent step is more frequent than all others combined, set as map precision..
    if ( length(sorted.freqs) == 1L || (length(sorted.freqs) > 1L &&
         sorted.freqs[1L] >= sum(sorted.freqs[2:length(sorted.freqs)])) ) {
        map.precision <- as.numeric(names(sorted.freqs)[1L])
    # ..otherwise set smallest step as map precision.
    } else {
        map.precision <- min(map.steps, na.rm=TRUE)
    }

    return(map.precision)
}

# is_single_char -----------------------------------------------------------------------------------
#' Test for a single character.
#'
#' @param x Test object.
#'
#' @return \code{TRUE} if \code{x} is a single character; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_char
is_single_char <- function(x) {
    return( is.character(x) && length(x) == 1L && ! is.na(x) && nchar(x) == 1L )
}

# is_single_nonnegative_number ---------------------------------------------------------------------
#' Test for a single non-negative number.
#'
#' @param n Test object.
#'
#' @return \code{TRUE} if \code{x} is a single non-negative number; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_nonnegative_number
is_single_nonnegative_number <- function(n) {
    return( is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 0.0 )
}

# is_single_nonnegative_whole_number ---------------------------------------------------------------
#' Test for a single non-negative whole number.
#'
#' @param n Test object.
#' @param tol Numeric tolerance.
#'
#' @return \code{TRUE} if \code{x} is a single non-negative whole number; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_nonnegative_whole_number
is_single_nonnegative_whole_number <- function(n, tol=.Machine$double.eps^0.5) {
    return( is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 0.0 &&
            abs(n - round(n)) < abs(tol) )
}

# is_single_probability ----------------------------------------------------------------------------
#' Test for a single valid probability.
#'
#' @param n Test object.
#'
#' @return \code{TRUE} if \code{x} is a single valid probability; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_probability
is_single_probability <- function(n) {
    return( is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 0.0 && n <= 1.0 )
}

# is_single_string ---------------------------------------------------------------------------------
#' Test for a single character string.
#'
#' @param x Test object.
#'
#' @return \code{TRUE} if \code{x} is a single string; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_string
is_single_string <- function(x) {
    return( is.character(x) && length(x) == 1L && ! is.na(x) )
}

# is_whole_number ----------------------------------------------------------------------------------
#' Test for whole numbers.
#'
#' @param n Test vector.
#' @param tol Numeric tolerance.
#'
#' @return Logical vector indicating which elements of \code{n} are whole numbers.
#'
#' @keywords internal
#' @rdname is_whole_number
is_whole_number <- function(n, tol=.Machine$double.eps^0.5) {
    return( is.numeric(n) & is.finite(n) & abs(n - round(n)) < abs(tol) )
}

# End of internal.R ################################################################################