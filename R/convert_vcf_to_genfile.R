# Start of convert_vcf_to_genfile.R

# convert_vcf_to_genfile
#' Convert VCF data to an \pkg{R/qtl} genotype file.
#'
#' Read SNP genotype data from one or more VCF \code{infiles}, and output
#' these to a \code{genfile} that may be accepted as input by \pkg{R/qtl}.
#'
#' @section What does this function do?:
#'
#' Given one or more VCF \code{infiles}, along with a set of identifiers for \code{samples} and
#' \code{founders} from a two-founder cross, this function first reads SNP allele and genotype
#' calls. Then, for the set of SNPs that are common to samples and founders, where the base
#' call of a sample allele matches that of a founder, the sample is assigned the corresponding
#' founder allele. The resulting sample genotype is made up of the combination of its founder
#' alleles. A SNP is retained only if each founder has a distinct homozygous genotype, and there
#' are at least two distinct genotypes among the samples for that SNP. Finally, the set of SNP
#' markers satisfying these conditions are output to \code{genfile}.
#'
#' @section What does this function \emph{not} do?:
#'
#' This function does not convert genotypes for a multi-founder cross.
#'
#' Variants are taken at face value, and there is no consideration of variant or genotype quality.
#' To exclude low-quality variants and genotypes, you may wish to perform a hard filter beforehand
#' using a suitable VCF or BCF toolkit.
#'
#' No attempt is made to collapse redundant markers. Depending on SNP density and linkage block
#' size, the output \code{genfile} may contain groups of linked SNPs that could be collapsed to
#' a single marker. If appropriate, this can be done with the data in the \code{genfile}.
#'
#' While care has been taken to ensure that the output \code{genfile} accurately reflects
#' the genotype data in the input VCF file(s), the genotype data is not validated with
#' respect to a specific cross type or experimental design. In general, you should review
#' the output file to ensure that the results are appropriate for your dataset.
#'
#' @section SNP marker IDs:
#' Each SNP marker is assigned an identifier of the form \code{'chr11:002160000'}, where
#' \code{'chr11'} is the chromosome identifier taken from the VCF \code{CHROM} field, and
#' \code{'002160000'} is a zero-padded number giving the position of the given SNP in the
#' reference genome, taken from the VCF \code{POS} field.
#'
#' The width to which the genomic position is padded is determined by the maximum sequence
#' length. Ideally reference sequence lengths are already specified in contig fields in the
#' header of the VCF. In that case the maximum sequence length is taken from the input VCF
#' file, and there is no need to specify a \code{max.seqlength} parameter.
#'
#' If the VCF header does not contain the necessary reference sequence length information,
#' then the \code{max.seqlength} parameter should be set to the appropriate value for the
#' given reference genome to ensure consistent formatting of SNP marker IDs across datasets.
#'
#' Sequence length information may be obtained by using one of the \code{ChromInfo} functions
#' provided by the \pkg{GenomeInfoDb} package. Alternatively, it may be obtained from a Picard
#' Tools sequence dictionary file using the function \code{\link{load_seq_dict}}, or directly
#' from a reference genome FASTA file using the \pkg{Biostrings} function \code{fasta.seqlengths}.
#' For more, see the \strong{Examples} section below.
#'
#' @param infiles Input VCF file paths.
#' @param genfile Output \pkg{R/qtl} genotype data file.
#' @param samples Cross sample IDs.
#' @param founders Founder sample IDs.
#' @param alleles Optional mapping of founder IDs to allele symbols (e.g. \code{utl::mapping(
#' c(DBVPG6044 = 'W', Y12 = 'S') )}). If this parameter is not specified, allele symbols are
#' taken from the letters of the alphabet (i.e. \code{'A'}, \code{'B'} etc.).
#' @param max.seqlength Optional parameter to indicate the maximum reference sequence length, which
#' is used to determine the zero-padded width of genomic positions in SNP marker IDs. Without this
#' information, SNP marker IDs may be formatted inconsistently in different datasets. For more
#' details, see the \strong{SNP marker IDs} section.
#' @param na.string String to replace \code{NA} values.
#'
#' @examples
#' \dontrun{
#'
#' # To get maximum sequence length from a sequence dictionary.
#' seqinfo <- utl::load_seq_dict('reference_genome.dict')
#' max.seqlength <- max(GenomeInfoDb::seqlengths(seqinfo))
#'
#' # To get maximum sequence length from reference sequence FASTA file.
#' max.seqlength <- max(Biostrings::fasta.seqlengths('reference_genome.fasta'))
#'
#' # Convert VCF to R/qtl genfile.
#' infiles <- c('samples.vcf', 'founders.vcf')
#' genfile <- 'geno.csv'
#' sample.ids <- utl::read_samples_from_vcf('samples.vcf')
#' founder.ids <- utl::read_samples_from_vcf('founders.vcf')
#' alleles <- utl::mapping(c(FOUNDER1='A', FOUNDER2='B'))
#' utl::convert_vcf_to_genfile(infiles, genfile, sample.ids, founder.ids,
#'                             alleles=alleles, max.seqlength=max.seqlength)
#'
#' }
#'
#' @seealso \href{https://bioconductor.org/packages/release/bioc/html/Biostrings.html}{Biostrings package}
#' @seealso \href{https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html}{GenomeInfoDb package}
#' @seealso \href{https://broadinstitute.github.io/picard/}{Picard Tools documentation}
#' @seealso \href{http://www.rqtl.org}{R/qtl website}
#'
#' @export
#' @rdname convert_vcf_to_genfile
convert_vcf_to_genfile <- function(infiles, genfile, samples, founders, alleles=NULL,
                                   max.seqlength=NULL, na.string=c('-', 'NA')) {

    stopifnot( is_single_string(genfile) )
    na.string <- match.arg(na.string)

    clashing.ids <- intersect(samples, founders)
    if ( length(clashing.ids) > 0 ) {
        stop("sample/founder ID clash - '", toString(clashing.ids), "'")
    }

    sample.data <- read_snps_from_vcf(infiles, samples=samples, max.seqlength=max.seqlength,
                                      require.any=TRUE, require.polymorphic=TRUE)

    founder.data <- read_snps_from_vcf(infiles, samples=founders, max.seqlength=max.seqlength,
                                       require.all=TRUE, require.polymorphic=TRUE)

    geno.mat <- make_geno_matrix(sample.data, founder.data, alleles=alleles)
    snp.loc <- parse_snp_marker_ids(colnames(geno.mat))

    id.col <- c('id', '', rownames(geno.mat))
    geno.mat <- rbind(colnames(geno.mat), snp.loc$chr, geno.mat)
    geno.mat <- cbind(id.col, geno.mat)

    utils::write.table(geno.mat, file=genfile, na=na.string, sep=',',
                       quote=FALSE, row.names=FALSE, col.names=FALSE)

    return( invisible() )
}

# End of convert_vcf_to_genfile.R