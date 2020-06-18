# Start of plot_qtl_scanone.R

# plot_qtl_scanone
#' Draw a plot of an \pkg{R/qtl} \code{scanone} result.
#'
#' Plotting function to plot either a LOD curve plot or Manhattan plot, using methods based on
#' \pkg{R/qtl} \code{plot.scanone} (Broman \emph{et al.} 2003) and \pkg{qqman} \code{manhattan}
#' (Turner 2018), respectively. If you find either useful, please give credit where it is due
#' to the original author(s).
#'
#' @param x An \pkg{R/qtl} \code{scanone} object.
#' @param chr Vector indicating which sequences to plot. If none are specified, all are plotted.
#' @param lodcolumn This parameter indicates which LOD column of \code{x} to consider. This must
#' be either a LOD column name or an index \emph{with respect to the set of LOD columns}. If no
#' LOD column is specified and one such column is found, that column is used by default; otherwise
#' a LOD column must be specified.
#' @param threshold A single \code{numeric} LOD significance threshold, or a
#' \code{summary.scanoneperm} object containing one such threshold and its
#' associated significance level.
#' @param qtl.intervals A \code{list} of \code{data.frame} objects created from the same data as
#' \code{x}, such that each \code{data.frame} contains three rows of information about the lower
#' interval limit, peak, and upper interval limit (respectively) of a given QTL. These intervals
#' are only included in the plot if they would be visually distinguishable.
#' @param col Analogous to the standard \code{'col'} plotting parameter. This is recycled to match
#' the number of sequences being plotted. As in the package \pkg{qqman}, this defaults to two
#' alternating monochrome shades.
#' @param gap Gap (in centiMorgans) between sequences in multi-sequence plots.
#' @param phenotype Name of the phenotype, to be shown in plot information.
#' @param type Type of plot. Set to \code{'l'} to output a standard LOD curve, or to \code{'p'} for
#' a Manhattan plot of individual LOD scores. If no plot type is specified, this is automatically
#' set based on the number of markers being plotted.
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
#'
#' # Plot complete LOD profile, including any QTL intervals.
#' plot_qtl_scanone(scanone.result, lodcolumn=1L, threshold=threshold.val,
#'                  qtl.intervals=qtl.intervals, phenotype='Blood pressure')
#'
#' # Plot LOD profile of chromosomes with significant QTL intervals.
#' plot_qtl_scanone(scanone.result, chr=c('1', '4'), lodcolumn=1L, threshold=threshold.val,
#'                  qtl.intervals=qtl.intervals, phenotype='Blood pressure')
#' }
#'
#' @references Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping in experimental
#' crosses. \emph{Bioinformatics} \bold{19}:889-890.
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/12724300}{PubMed})
#' @references Turner SD (2018) qqman: an R package for visualizing GWAS results using Q-Q and
#' manhattan plots. (\href{https://doi.org/10.21105/joss.00731}{Journal of Open Source Software})
#' @seealso \href{http://www.rqtl.org}{R/qtl website}
#' @seealso \href{https://cran.r-project.org/web/packages/qqman/}{qqman package}
#'
#' @export
#' @importFrom grDevices graphics.off
#' @include internal.R
#' @rdname plot_qtl_scanone
plot_qtl_scanone <- function(x, chr=NULL, lodcolumn=NULL, threshold=NULL, qtl.intervals=NULL,
    col=c('gray10', 'gray60'), gap=25L, phenotype=NULL, type=NULL) {

    opar <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(opar))

    # validate arguments

    # Validate scanone object.
    stopifnot( 'scanone' %in% class(x) )
    stopifnot( nrow(x) > 0L )
    if ( is.factor(x$chr) ) {
        x$chr <- as.character(x$chr)
    } else {
        stopifnot( is.character(x$chr) )
    }

    # Resolve which sequences to plot.
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

    # Resolve which LOD column to plot.
    lodcol.index <- get_lod_col_index(x, lodcolumn=lodcolumn)

    # Validate LOD threshold, if specified.
    if ( ! is.null(threshold) ) {
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
    } else {
        tinfo <- NULL
    }

    # Validate QTL intervals, if specified.
    if ( ! is.null(qtl.intervals) ) {
        stopifnot( is.list(qtl.intervals) )
        stopifnot( length(qtl.intervals) > 0L )
        for ( i in seq_along(qtl.intervals) ) {
            stopifnot( is.data.frame(qtl.intervals[[i]]) )
            stopifnot( nrow(qtl.intervals[[i]]) == 3L )
            stopifnot( 'chr' %in% colnames(qtl.intervals[[i]]) )
            if ( is.factor(qtl.intervals[[i]]$chr) ) {
                qtl.intervals[[i]]$chr <- as.character(qtl.intervals[[i]]$chr)
            } else {
                stopifnot( is.character(qtl.intervals[[i]]$chr) )
            }
            stopifnot( ! anyNA(qtl.intervals[[i]]$chr) )
            stopifnot( length(unique(qtl.intervals[[i]]$chr)) == 1L )
            stopifnot( unique(qtl.intervals[[i]]$chr) %in% x$chr )
            stopifnot( 'pos' %in% colnames(qtl.intervals[[i]]) )
            stopifnot( is.numeric(qtl.intervals[[i]]$pos) )
            stopifnot( ! anyNA(qtl.intervals[[i]]$pos) )
            stopifnot( ! is.unsorted(qtl.intervals[[i]]$pos) )
            stopifnot( min(qtl.intervals[[i]]$pos) >= 0.0 )
            stopifnot( max(qtl.intervals[[i]]$pos) <=
                       max(x[x$chr == unique(qtl.intervals[[i]]$chr), ]$pos) )
        }
    }

    # Validate plot 'col' parameter.
    stopifnot( is.vector(col) )
    stopifnot( length(col) > 0L )
    stopifnot( ! anyNA(col) )

    # Validate inter-chromosomal gap.
    stopifnot( is_single_nonnegative_number(gap) )

    # Resolve phenotype display name, if possible.
    if ( ! is.null(phenotype) ) {
        stopifnot( is_single_string(phenotype) )
    } else if ( colnames(x)[lodcol.index] != 'lod' ) {
        phenotype <- colnames(x)[lodcol.index]
    }

    # Validate plot type, if specified.
    if ( ! is.null(type) ) {
        stopifnot( is_single_char(type) )
        stopifnot( type %in% c('l', 'p') )
    }

    # prepare data

    # Subset scanone result by specified sequences.
    x <- x[x$chr %in% chr, ]

    # Replace any NA values with zero.
    x[is.na(x[, lodcol.index]), lodcol.index] <- 0.0

    # Get list containing scanone row indices for each sequence.
    seq.indices <- lapply(chr, function(s) which( x$chr == s ))
    names(seq.indices) <- chr

    # Set maximum y-value, ensure greater than or equal to one.
    max.lod <- max(x[, lodcol.index], tinfo[['threshold']], 1.0, na.rm=TRUE)

    # assemble sequence plotting info

    # Set width of gaps between sequences.
    gap.width <- ifelse(length(chr) > 1L, gap, 0.0)

    # Init cumulative plot width.
    cum.plot.width <- 0.0

    # Init sequence plotting info.
    seq.par <- matrix(NA_real_, nrow=length(chr), ncol=3L, dimnames=
                      list(chr, c('offset', 'midpoint', 'length')))

    for ( i in seq_along(chr) ) {

        if ( length(chr) > 1L ) {
            seq.start <- x[utils::head(seq.indices[[i]], 1L), 'pos']
        } else {
            seq.start <- 0.0
        }

        seq.end <- x[utils::tail(seq.indices[[i]], 1L), 'pos']

        seq.par[i, 'offset'] <- cum.plot.width

        seq.par[i, 'length'] <- seq.end - seq.start

        seq.par[i, 'midpoint'] <- seq.par[i, 'offset'] + (0.5 * seq.par[i, 'length'])

        cum.plot.width <- cum.plot.width + seq.par[i, 'length']

        if ( i < length(chr) ) {
            cum.plot.width <- cum.plot.width + gap.width
        }
    }

    # assemble general plot info

    plot.info <- list()

    # Get phenotype from argument or LOD column, if possible.
    if ( ! is.null(phenotype) ) {
        stopifnot( is_single_string(phenotype) )
        plot.info['Phenotype'] <- phenotype
    } else {
        profile.phename <- colnames(x)[lodcol.index]
        if ( profile.phename != 'lod' ) {
            plot.info['Phenotype'] <- profile.phename
        }
    }

    # For single-sequence plots, add sequence label to plot info.
    if ( length(chr) == 1L ) {
        plot.info['Chromosome'] <- chr
    }

    # draw plot

    # Set top margin line counts, accounting for any plot info lines.
    adj.top.mar <- length(plot.info) + 4.0 # NB: minimum 4 margin lines
    graphics::par(mar=(c(5.0, 4.0, adj.top.mar, 2.0) + 0.1))

    # Set plot info margin indices.
    plot.info.indices <- rev(seq_along(plot.info))
    plot.title.index <- adj.top.mar - 2.0

    # Set x-axis plotting parameters; values taken from R/qtl.
    if ( length(chr) > 1L ) {
        xlim <- c(-(0.5 * gap.width), cum.plot.width + (0.5 * gap.width))
        xlab <- 'Chromosome'
        xaxt <- 'n'
    } else {
        xlim <- c(0.0, cum.plot.width)
        xlab <- 'Map position (cM)'
        xaxt <- 's'
    }

    # Set y-axis plotting parameters.
    headspace <- 0.2
    ylim <- c(0.0, max.lod / (1.0 - headspace))
    ylab <- 'LOD'

    # Set palette for sequences.
    col <- rep_len(col, length(chr))

    # Start new plot.
    graphics::plot(NA, family='sans', xaxs='i', xpd=FALSE, yaxs='i', bg='white', cex=1.0, las=1L,
                   lty='solid', xaxt=xaxt, xlab=xlab, xlim=xlim, ylab=ylab, ylim=ylim)

    # If plot type not specified, choose based on number of marker loci.
    if ( is.null(type) ) {
        num.markers <- length( ! is_pseudomarker_id(rownames(x)) )
        if ( num.markers > 10000L ) {  # marker count threshold
            type <- 'p'
        } else {
            type <- 'l'
        }
    }

    # Plot scanone result for each sequence.
    for ( i in seq_along(chr) ) {

        seq.x <- x[seq.indices[[i]], 'pos'] + seq.par[i, 'offset']
        seq.y <- x[seq.indices[[i]], lodcol.index]

        if (type == 'l') {
            graphics::lines(seq.x, seq.y, col=col[i], lwd=2.0, lty='solid')
        } else if (type == 'p') {
            graphics::points(seq.x, seq.y, col=col[i], lwd=0.5, pch=20L)
        }
    }

    # If QTL intervals available, add to plot if they will be visually distinguishable.
    if ( ! is.null(qtl.intervals) ) {

        # Subset QTL intervals by specified sequences.
        qtl.intervals <- qtl.intervals[ sapply(qtl.intervals, function(qtl.interval)
                                        unique(qtl.interval$chr)) %in% chr ]

        if ( length(qtl.intervals) > 0L ) {

            # Get size of QTL interval in centiMorgans.
            interval.widths <- sapply(qtl.intervals, function(qtl.interval)
                                      diff(qtl.interval[c(1L, 3L), 'pos']))

            # Try to get size of bullet symbol (pch=20) in user-space..
            bullet.width <- tryCatch({
                result <- suppressWarnings(graphics::strwidth('\u2022',
                                           cex=args$cex, family ='sans'))
            }, error=function(e) { # ..otherwise get default bullet width.
                bullet.inches <- 0.06944444 # default bullet width in inches (family='sans',cex=1.0)
                plot.inches <- graphics::par('pin')[1L] # plot width in inches
                plot.space <- diff(xlim) # plot width in user-space
                result <- (bullet.inches / plot.inches) * plot.space
            })

            # Plot QTL intervals if plotting LOD curves and typical
            # QTL intervals would be visually distinguishable.
            if ( type == 'l' && stats::median(interval.widths) >= bullet.width ) {

                # Get sequences corresponding to QTL intervals.
                interval.seqs <- sapply(qtl.intervals, function(obj) unique(obj[, 'chr']))

                # Set vertical offset of 1.5-LOD interval from
                # LOD peak in terms of the overall plot height.
                y.offset <- 0.1 * diff(ylim)

                # Display each QTL interval, remember the plot regions it occupies.
                for ( i in seq_along(qtl.intervals) ) {

                    qtl.interval <- qtl.intervals[[i]]

                    # Get the horizontal offset for this sequence (i.e. cumulative
                    # length of previous sequences and gaps between.).
                    interval.offset <- seq.par[interval.seqs[[i]], 'offset']

                    # Get horizontal positions of QTL interval parts.
                    interval.xpos <- qtl.interval[, 'pos'] + interval.offset

                    # Get maximum LOD value for this QTL interval.
                    interval.start <- qtl.interval[1L, 'pos']
                    interval.end <- qtl.interval[3L, 'pos']
                    interval.lods <- x[x$chr == interval.seqs[[i]] & x$pos >= interval.start &
                                       x$pos <= interval.end, lodcol.index]
                    interval.max.lod <- max(interval.lods)

                    # Get interval line vertical position from interval max LOD and vertical offset.
                    iline <- interval.max.lod + y.offset

                    # Get vertical positions of QTL interval parts.
                    interval.ypos <- c(
                        iline - 0.25 * y.offset,
                        iline,
                        iline + 0.25 * y.offset
                    )

                    # Get QTL interval line segment endpoints.
                    x0 <- interval.xpos[c(1L, 1L, 3L)]
                    x1 <- interval.xpos[c(1L, 3L, 3L)]
                    y0 <- interval.ypos[c(1L, 2L, 1L)]
                    y1 <- interval.ypos[c(3L, 2L, 3L)]

                    # Draw QTL interval line segments.
                    graphics::segments(x0, y0, x1, y1, col='black', lwd=0.5, lty='solid')

                    # Draw point at position of LOD peak.
                    graphics::points(interval.xpos[2L], iline, col='black', lwd=0.5, pch=20L)
                }
            }
        }
    }

    # Plot LOD threshold if available.
    if ( ! is.null(tinfo[['threshold']]) ) {

        # Draw horizontal red dashed line to indicate LOD threshold.
        graphics::abline(tinfo[['threshold']], 0.0, col='red', lwd=1.0, lty='dotted')

        # If significance level available, add to threshold line.
        if ( ! is.null(tinfo[['alpha']]) ) {

            # Set threshold label text.
            thresh.label.text <- bquote(bold(alpha ~ '=' ~ .(tinfo[['alpha']])))

            # Set threshold label width.
            thresh.label.width <- graphics::strwidth(thresh.label.text)

            # Set threshold label position.
            thresh.label.pos <- cum.plot.width - thresh.label.width - graphics::xinch(0.025)

            # Draw LOD significance threshold label.
            graphics::text(thresh.label.pos, tinfo[['threshold']] + graphics::yinch(0.025),
                           thresh.label.text, col='red', adj=c(0.0, 0.0), cex=0.8)
        }
    }

    # If multiple sequences, add sequence ticks and labels.
    if ( length(chr) > 1L ) {

        # Adjust size of sequence labels to ensure that they don't overlap.
        seq.label.size <- max(graphics::strwidth(chr))
        seq.plot.size <- min(seq.par[, 'length'])
        if ( seq.plot.size < seq.label.size ) {
            cex.axis <- seq.plot.size / seq.label.size
        } else {
            cex.axis <- 1.0
        }

        # Plot x-axis with given ticks and labels.
        for ( i in seq_along(chr) ) {
            graphics::axis(side=1L, at=seq.par[i, 'midpoint'], labels=chr[i], cex.axis=cex.axis)
        }
    }

    # Plot box around graph.
    graphics::box(lwd=3.0)

    # Write plot title.
    graphics::title(main='Scanone', line=plot.title.index, cex=1.2, col='black', family='sans')

    # If any plot info, add to top margin.
    if ( length(plot.info) > 0L ) {

        plot.info.lines <- sapply(names(plot.info), function(k) paste0(k, ': ', plot.info[[k]]))

        for ( i in seq_along(plot.info.lines) ) {
            graphics::mtext(plot.info.lines[i], line=plot.info.indices[i], side=3L, adj=0.0,
                            cex=1.0, col='black', family='sans', font=1L)
        }
    }

    return(invisible())
}

# End of plot_qtl_scanone.R