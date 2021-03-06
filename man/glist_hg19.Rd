% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{glist_hg19}
\alias{glist_hg19}
\title{Gene positions - hg19}
\format{
A dataframe with 26,292 rows and 4 variables:
\describe{
 \item{chr}{Chromosome code}
 \item{start}{Start of gene (base-pair units, 1-based)}
 \item{end}{End of gene (this position is included in the interval)}
 \item{name}{Gene ID}
}
}
\source{
\url{https://www.cog-genomics.org/static/bin/plink/glist-hg19}
}
\usage{
glist_hg19
}
\description{
Gene range lists downloaded from PLINK 1.9 resources, generated from
UCSC Table Browser RefSeq track
}
\keyword{datasets}
