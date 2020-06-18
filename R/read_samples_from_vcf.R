# Start of read_samples_from_vcf.R

# read_samples_from_vcf
#' Read sample IDs from a VCF file.
#'
#' @param file A VCF file path.
#'
#' @return Vector of the sample IDs in the given VCF file.
#'
#' @export
#' @rdname read_samples_from_vcf
read_samples_from_vcf <- function(file) {
    stopifnot( is_single_string(file) )
    stopifnot( file.exists(file) )
    header <- VariantAnnotation::scanVcfHeader(file)
    return( VariantAnnotation::samples(header) )
}

# End of read_samples_from_vcf.R