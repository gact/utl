# Start of is_pseudomarker_id.R ####################################################################

# is_pseudomarker_id -------------------------------------------------------------------------------
#' Test for \pkg{R/qtl} pseudomarker IDs.
#'
#' @param x Character vector of locus IDs.
#'
#' @return Logical vector indicating which elements of \code{x} are \pkg{R/qtl} pseudomarker IDs.
#' Pseudomarker IDs are used by \pkg{R/qtl} for inter-marker loci. They indicate the reference
#' sequence and genetic map position of the locus (e.g. \code{'c4.loc33'} for a locus at position
#' 33cM on chromosome 4).
#'
#' @export
#' @include internal.R
#' @rdname is_pseudomarker_id
is_pseudomarker_id <- function(x) {
    return( is.character(x) & nzchar(x) &
            grepl('^c([[:alnum:]]+)[.]loc(-?[[:digit:]]+(?:[.][[:digit:]]+)?)$', x) )
}

# End of is_pseudomarker_id.R ######################################################################