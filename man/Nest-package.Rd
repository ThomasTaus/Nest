\name{Nest-package}
\alias{Nest-package}
\alias{Nest}
\docType{package}
\title{
Study Random Genetic Drift
}
\description{
The \code{Nest} package allows you to simulate allele frequency trajectories under random genetic drift using the Wright-Fisher model and estimate the effective population size (Ne) applying different moment-based methods. Pool-seq data can be simulated by re-sampling allele count data. Additionally, allele frequency data can be read from file and Ne can be readily estimated for windows along the genome, or for whole chromosomes instead.
}
\details{

}
\author{
\packageAuthor{Nest}

Maintainer: \packageMaintainer{Nest}
}
\references{
Waples R. S.: A generalized approach for estimating effective population size from temporal changes in allele frequency, \emph{Genetics} \bold{1989}, 121, 379–391.

Jorde P. E. and Ryman N.: Unbiased estimator for genetic drift and effective population size, \emph{Genetics} \bold{2007}, 177 927–935.

Frick K., Munk, A. and Sieling, H.: Multiscale Change-Point Inference, \emph{Journal of the Royal Statistical Society: Series B} \bold{2014}, 76, 495-580.

Futschik A., Hotz T., Munk A. and Sieling H.: Multiresolution DNA partitioning: statistical evidence for segments, \emph{Bioinformatics} \bold{2014}, 30, 2255-2262.

Jónás A., Taus T., Kosiol C., Schlötterer C. & Futschik A.: Estimating effective population size from temporal allele frequency changes in experimental evolution, manuscript in preparation.
}
\keyword{ package }
\seealso{
\code{\link{checkSNP}}, \code{\link{estimateNe}}, \code{\link{estimateWndNe}}, \code{\link{read.sync}}, \code{\link{sample.alleles}} and \code{\link{wf.traj}}.
}
\examples{

}
