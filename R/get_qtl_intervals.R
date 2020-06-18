# Start of get_qtl_intervals.R

# get_qtl_intervals
#' Get list of QTL intervals.
#'
#' Given a specified LOD \code{threshold} and associated significance level \code{alpha}, this
#' function gets a list of approximate QTL intervals for an \pkg{R/qtl} \code{scanone} object.
#'
#' @param x An \pkg{R/qtl} \code{scanone} object.
#' @param chr Vector indicating which sequences to consider. If no sequences are specified, the QTL
#' interval list is created with respect to every sequence in \code{x}.
#' @param lodcolumn This parameter indicates which LOD column to consider. This must be either a LOD
#' column name or an index \emph{with respect to the set of LOD columns}. If no LOD column is
#' specified and one such column is found, that column is used by default; otherwise a LOD column
#' must be specified.
#' @param threshold A single \code{numeric} LOD significance threshold, or a
#' \code{summary.scanoneperm} object containing one such threshold and its
#' associated significance level.
#' @param ci.function Option to indicate which \pkg{R/qtl} function should be used for estimating
#' approximate confidence intervals for QTL location. Set to \code{'lodint'} for LOD support
#' intervals (adjusting stringency with the \code{drop} parameter), or to \code{'bayesint'} for
#' Bayesian credible intervals (adjusting stringency with the \code{prob} parameter). For more
#' information on the QTL interval methods used, see functions \code{lodint} and \code{bayesint}
#' in the \pkg{R/qtl} manual, as well as Section 4.5 of Broman and Sen (2009).
#' @param drop LOD units that the LOD profile must drop to form the interval.
#' This is used only if \code{ci.function} is set to \code{'lodint'}.
#' @param prob Desired probability coverage for the Bayesian credible interval.
#' This is used only if \code{ci.function} is set to \code{'bayesint'}.
#' @param expandtomarkers Expand the LOD interval to the nearest flanking
#' markers, or to the respective terminal loci.
#'
#' @return A list of \code{data.frame} objects, each containing three rows of information about the
#' lower interval limit, peak, and upper interval limit (respectively) of a QTL. Returns an empty
#' list if there are no significant QTLs.
#'
#' @examples
#' \dontrun{
#' # Load R/qtl hyper dataset.
#' data(hyper, package='qtl')
#'
#' # Estimate genetic map of hyper data.
#' gmap <- qtl::est.map(hyper, offset=0.0)
#'
#' # Set newly estimated genetic map.
#' hyper <- qtl::replace.map(hyper, gmap)
#'
#' # Calculate genotype probabilities.
#' hyper <- qtl::calc.genoprob(hyper, step=1.0)
#'
#' # Do single-QTL analysis.
#' scanone.result <- qtl::scanone(hyper, pheno='bp')
#'
#' # Do single-QTL permutation analysis.
#' scanone.perms <- qtl::scanone(hyper, pheno='bp', n.perm=1000L)
#'
#' # Get LOD threshold values from single-QTL permutation results.
#' threshold.vals <- summary(scanone.perms, alpha=0.05)
#'
#' # Get LOD threshold value for first LOD column.
#' threshold.val <- qtl:::subset.scanoneperm(threshold.vals, lodcolumn=1L)
#'
#' # Get QTL intervals for first LOD column.
#' qtl.intervals <- get_qtl_intervals(scanone.result, lodcolumn=1L, threshold=threshold.val,
#'                                    ci.function='bayesint')
#' }
#'
#' @author Thomas A. Walsh
#' @author Yue Hu
#' @references Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping
#' in experimental crosses. \emph{Bioinformatics} \bold{19}:889-890.
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/12724300}{PubMed})
#' @references Broman KW, Sen S (2009) \emph{A guide to QTL mapping with R/qtl.}
#' New York: Springer. (\href{http://www.rqtl.org/book/}{Website})
#' @seealso \href{http://www.rqtl.org}{R/qtl website}
#'
#' @export
#' @include internal.R
#' @rdname get_qtl_intervals
get_qtl_intervals <- function(x, chr=NULL, lodcolumn=NULL, threshold=NULL, ci.function=c('lodint',
    'bayesint'), drop=1.5, prob=0.95, expandtomarkers=FALSE) {

    # validate arguments

    # Validate scanone object.
    stopifnot( 'scanone' %in% class(x) )
    stopifnot( nrow(x) > 0L )
    if ( is.factor(x$chr) ) {
        x$chr <- as.character(x$chr)
    } else {
        stopifnot( is.character(x$chr) )
    }

    # Resolve which sequences to consider.
    if ( ! is.null(chr) ) {
        if ( is.factor(chr) ) {
            chr <- as.character(chr)
        } else {
            stopifnot( is.character(chr) )
        }
        stopifnot( length(chr) > 0L )
        stopifnot( ! anyNA(chr) )
        stopifnot( ! anyDuplicated(chr) )
        stopifnot( all( chr %in% x$chr ) )
    } else {
        chr <- unique(x$chr)
    }

    # Resolve which LOD column to consider.
    lodcol.index <- get_lod_col_index(x, lodcolumn=lodcolumn)

    # Validate LOD threshold.
    stopifnot( ! is.null(threshold) )
    if ( 'summary.scanoneperm' %in% class(threshold) ) {
        stopifnot( nrow(threshold) == 1L )
        stopifnot( ncol(threshold) == 1L )
        stopifnot( ! anyNA(threshold) )
        tinfo <- list(threshold = threshold[1L, 1L], alpha = 0.01 *
                      as.numeric(sub('%', '', rownames(threshold)[1L])))
    } else if ( 'numeric' %in% class(threshold) ) {
        stopifnot( is_single_nonnegative_number(threshold) )
        tinfo <- list(threshold=unname(threshold), alpha=NULL)
    } else {
        stop("threshold must have class 'summary.scanoneperm' or 'numeric', not '",
             toString(class(threshold)), "'")
    }

    # Resolve which interval function to use.
    ci.function <- match.arg(ci.function)

    # If using LOD support intervals, validate LOD support interval drop..
    if ( ci.function == 'lodint' ) {
        stopifnot( is_single_nonnegative_number(drop) )
    } else {  # ..otherwise validate Bayesian credible interval probability.
        stopifnot( is_single_probability(prob) )
    }

    # Validate marker expansion option.
    stopifnot( isTRUE(expandtomarkers) || identical(FALSE, expandtomarkers) )

    # prepare data

    # Get LOD profile for the given LOD column.
    poscol.indices <- which( colnames(x) %in% c('chr', 'pos') )
    x <- x[, c(poscol.indices, lodcol.index)]
    colnames(x) <- c('chr', 'pos', 'lod')

    # Get map sequences.
    map.seqs <- unique(x$chr)

    # Init QTL interval list.
    intervals <- list()
    attr(intervals, 'threshold') <- tinfo[['threshold']]
    attr(intervals, 'alpha') <- tinfo[['alpha']]
    if ( ci.function == 'lodint' ){
        attr(intervals, 'drop') <- drop
    } else {  # ci.function == 'bayesint'
        attr(intervals, 'prob') <- prob
    }

    # Return empty QTL interval list if no relevant LOD data.
    if( nrow(x) == 0L || all( is.na(x[, 'lod']) ) ) {
        return(intervals)
    }

    # find QTL peaks

    # Create LOD character mask. Regions with significant LOD values are marked
    # with the sequence name, while all other loci are marked with NA.
    # NB: this is used instead of a simple logical mask so as to prevent
    # significant regions 'spilling over' across sequence boundaries.
    lod.mask <- vector('character', nrow(x))
    for ( chr.seq in chr ) {
        indices <- which( x$chr == chr.seq )
        seq.lod <- x[indices, 'lod']
        bool.mask <- ! is.na(seq.lod) & seq.lod >= tinfo[['threshold']]
        char.mask <- sapply(bool.mask, function(sig.lod) if (sig.lod) {chr.seq} else { NA })
        lod.mask[indices] <- char.mask
    }

    # Get list of vectors, each containing the row indices of the loci in a significant region.
    significant.rows <- get_run_index_list(lod.mask, na.rm=TRUE)

    # Return empty QTL interval list if no significant regions found.
    if ( length(significant.rows) == 0L ) {
        return(intervals)
    }

    # Get row index of peak for each significant region, taking first peak in case of ties.
    peak.indices <- unlist(lapply(significant.rows, function(rows) rows[which.max(x[rows, 'lod'])]))

    # get QTL intervals

    # Get list of sequence ranges.
    index.ranges <- lapply(chr, function (chr.seq) range(which( x$chr == chr.seq )))
    names(index.ranges) <- chr

    # Set standard QTL interval indices.
    ci <- c(low=1L, peak=2L, high=3L)

    # Init matrix of raw intervals, each row of which will contain the row
    # indices of the lower limit, peak, and upper limit of a significant region.
    raw.intervals <- matrix(nrow=length(peak.indices), ncol=3L, dimnames=list(NULL, names(ci)))

    # Estimate LOD support intervals, if specified..
    if ( ci.function == 'lodint' ) {

        for ( i in seq_along(peak.indices) ) {

            # Init interval indices to peak index.
            l <- u <- peak.index <- peak.indices[i]

            # Get sequence endpoint indices.
            peak.seq <- x$chr[peak.index]
            seq.start <- index.ranges[[peak.seq]][1L]
            seq.end <- index.ranges[[peak.seq]][2L]

            # Set LOD cutoff.
            cutoff <- max(x$lod[peak.index] - drop, 0.0)

            # Decrement lower limit within this sequence while LOD score meets cutoff.
            while ( x$lod[l] >= cutoff && l > seq.start && ! is.na(x$lod[l-1L]) ) {
                l <- l - 1L
            }

            # Increment upper limit within this sequence while LOD score meets cutoff.
            while ( x$lod[u] >= cutoff && u < seq.end && ! is.na(x$lod[u+1L]) ) {
                u <- u + 1L
            }

            raw.intervals[i, ] <- c(l, peak.index, u)
        }

    # ..otherwise estimate Bayesian credible intervals.
    } else {  # ci.function == 'bayesint'

        for ( i in seq_along(peak.indices) ) {

            # Init interval indices to peak index.
            l <- u <- peak.index <- peak.indices[i]

            # Get sequence endpoint indices.
            peak.seq <- x$chr[peak.index]
            seq.start <- index.ranges[[peak.seq]][1L]
            seq.end <- index.ranges[[peak.seq]][2L]

            # Call R/qtl function to do Bayesian credible interval.
            bayes.interval <- qtl::bayesint(results=x, chr=peak.seq, qtl.index=1L, prob=prob,
                                            expandtomarkers=FALSE)

            # Decrement lower limit within this sequence until position matches
            # Bayesian interval lower endpoint (within numeric tolerance).
            while ( ! isTRUE(all.equal(x$lod[l], bayes.interval$pos[1L])) &&
                    l > seq.start && ! is.na(x$lod[l-1L]) ) {
                l <- l - 1L
            }

            # Increment upper limit within this sequence until position matches
            # Bayesian interval upper endpoint (within numeric tolerance).
            while ( ! isTRUE(all.equal(x$lod[u], bayes.interval$pos[3L])) &&
                    u < seq.end && ! is.na(x$lod[u+1L]) ) {
                u <- u + 1L
            }

            raw.intervals[i, ] <- c(l, peak.index, u)
        }
    }

    # Create matrix of merged intervals, formed by combining overlapping raw intervals.
    merged.intervals <- matrix(nrow=0L, ncol=3L, dimnames=list(NULL, names(ci)))
    i <- 1L
    while ( i <= nrow(raw.intervals) ) {

        # Merged range will include this interval.
        j <- i

        # Add any intervals whose limits overlap the merged range.
        while ( j < nrow(raw.intervals) && raw.intervals[j+1L, 'low'] <= raw.intervals[j, 'high'] ) {
            j <- j + 1L
        }

        # Set the lower limit to that of the first interval.
        l <- raw.intervals[i, 'low']

        # Set the merged interval peak to the highest peak across the
        # merged intervals, taking the first peak in case of ties.
        peak.indices <- raw.intervals[i:j, 'peak']
        peak.index <- peak.indices[which.max(x$lod[peak.indices])]

        # Set the upper limit to that of the last interval.
        u <- raw.intervals[j, 'high']

        # Add to merged intervals.
        merged.intervals <- rbind(merged.intervals, c(l, peak.index, u))

        i <- j + 1L
    }

    # If expanding to markers, expand each interval limit
    # to next flanking marker or sequence endpoint.
    if (expandtomarkers) {

        # Create mask indicating which loci are markers.
        x.rownames <- attr(x, 'row.names')
        if ( ! is.null(x.rownames) && ! identical(x.rownames, seq_len(nrow(x))) ) {
            marker.mask <- ! is_pseudomarker_id(rownames(x))
        } else {
            marker.mask <- rep_len(FALSE, nrow(x))
        }

        # Expand each interval to next flanking marker or sequence endpoint.
        for ( r in seq_len(nrow(merged.intervals)) ) {

            # Get interval limit indices.
            l <- merged.intervals[r, 'low']
            u <- merged.intervals[r, 'high']

            # Get sequence endpoint indices.
            s <- as.character(x$chr[l])  # NB: sequence identical across interval.
            seq.start <- index.ranges[[s]][1L]
            seq.end <- index.ranges[[s]][2L]

            # Decrement lower limit within this sequence while locus is not a marker.
            while ( ! marker.mask[l] && l > seq.start && ! is.na(x$lod[l-1L]) ) {
                l <- l - 1L
            }

            # Increment upper limit within this sequence while locus is not a marker.
            while ( ! marker.mask[u] && u < seq.end && ! is.na(x$lod[u+1L]) ) {
                u <- u + 1L
            }

            # Set expanded interval limit indices.
            merged.intervals[r, 'low'] <- l
            merged.intervals[r, 'high'] <- u
        }
    }

    # Rename interval endpoints if they coincide with peak.
    for ( r in seq_len(nrow(merged.intervals)) ) {

        merged.interval <- merged.intervals[r, ]

        interval <- as.data.frame(x[merged.interval, ])
        colnames(interval)[2L] <- 'pos'

        # If peak coincides with start of interval,
        # rename start of interval to 'ci.low'.
        if ( merged.interval['low'] == merged.interval['peak'] ) {
            temp.id <- rownames(interval)[ ci['low'] ]
            rownames(interval)[ ci['low'] ] <- 'ci.low'
            rownames(interval)[ ci['peak'] ] <- temp.id
        }

        # If peak coincides with end of interval,
        # rename end of interval to 'ci.high'.
        if ( merged.interval['peak'] == merged.interval['high'] ) {
            rownames(interval)[ ci['high'] ] <- 'ci.high'
        }

        intervals[[r]] <- interval
    }

    # Get interval peaks.
    interval.peaks <- do.call(rbind, lapply(intervals, function(interval) interval[ci['peak'], ]))

    # Get locus IDs at interval peaks.
    peak.ids <- rownames(interval.peaks)

    # Identify unnamed interval peaks.
    unnamed <- is_pseudomarker_id(peak.ids)

    # Ensure every QTL interval has a name.
    if ( any(unnamed) ) {
        map.precision <- infer_map_precision(x)
        digits <- max(1L, -floor(log10(map.precision)))  # from R/qtl::makeqtl
        unnamed.pos <- sprintf(paste0('%.', digits, 'f'), interval.peaks[unnamed, 'pos'])
        peak.ids[unnamed] <- paste(interval.peaks[unnamed, 'chr'], unnamed.pos, sep='@')
    }

    # Set QTL interval names.
    names(intervals) <- peak.ids

    return(intervals)
}

# End of get_qtl_intervals.R