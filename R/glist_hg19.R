#' Gene positions - hg19
#'
#' Gene range lists downloaded from PLINK 1.9 resources, generated from UCSC Table Browser RefSeq track
#'
#' @format A dataframe with 26,292 rows and 4 variables:
#' \describe{
#'  \item{chr}{Chromosome code}
#'  \item{start}{Start of gene (base-pair units, 1-based)}
#'  \item{end}{End of gene (this position is included in the interval)}
#'  \item{name}{Gene ID}
#' }
#' @source \url{https://www.cog-genomics.org/static/bin/plink/glist-hg19}
"glist_hg19"