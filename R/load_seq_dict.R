# Start of load_seq_dict.R

# load_seq_dict
#' Load sequence dictionary from file.
#'
#' Load a sequence dictionary for a given reference genome, such as may
#' be output by the Picard Tools command \code{CreateSequenceDictionary}.
#'
#' @param file A sequence dictionary file path.
#'
#' @return A \pkg{GenomeInfoDb} \code{Seqinfo} object containing
#' information from the input sequence dictionary.
#'
#' @seealso \href{https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html}{GenomeInfoDb package}
#' @seealso \href{https://broadinstitute.github.io/picard/}{Picard Tools documentation}
#'
#' @export
#' @rdname load_seq_dict
load_seq_dict <- function(file) {
    temp_path <- tempfile(pattern='seqdict_', fileext='.bam')
    on.exit( if ( file.exists(temp_path) ) { file.remove(temp_path) } )
    temp_stem <- substr(temp_path, 0, nchar(temp_path) - nchar('.bam'))
    bam_path <- Rsamtools::asBam(file, destination=temp_stem)
    bam <- Rsamtools::BamFile(temp_path)
    seqinfo <- Rsamtools::seqinfo(bam)
    return(seqinfo)
}

# End of load_seq_dict.R