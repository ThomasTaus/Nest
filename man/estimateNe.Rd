\name{estimateNe}
\alias{estimateNe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Effetive Population Size
}
\description{
\code{estimateNe} estimates the effective population size (Ne) based on tempral allele frequency changes.
}
\usage{
estimateNe(p0, pt, cov0, covt, t, ploidy = 2, truncAF = NA, method = "P_AF2_planI", poolSize = NA, Ncensus = NA, asList = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p0}{
numeric vector of starting allele frequencies.
}
  \item{pt}{
numeric vector of allele frequencies after \code{t} generations.
}
  \item{cov0}{
numeric vector of sequencing coverage at the first time point.
}
  \item{covt}{
numeric vector of sequencing coverage at the second time point.
}
  \item{t}{
integer specifying the number of generations between the two allele frequency measurements.
}
  \item{ploidy}{
integer specifying the ploidy of the organisms in the population. Diploids (\code{ploidy = 2} are the default.)
}
  \item{truncAF}{
numeric indicating whether SNPs with extreme allele frequencies should be excluded from Ne estimation, see \code{\link{checkSNP}}.
}
  \item{method}{ 
character vector defining the method(s) that should be used to estimate Ne. Possible values are \code{"P_AF"}, \code{"P_AF2_planI"} (the default), \code{"P_AF2_planII"}, \code{"P_AJ_planI"}, \code{"P_AJ_planII"}, \code{"JR_planI"}, \code{"JR_planII"}, \code{"W_planI"} and \code{"W_planII"}, see 'Details'.
}
  \item{poolSize}{
numeric vector indicating the number of individuals pooled at the first and second time point. If the entire population has been subjected to PoolSeq (default) then \code{poolSize = NA}.
}
  \item{Ncensus}{
integer specifying the census population size, which is required if \code{method} contains \code{"JR_planI"}.
}
  \item{asList}{
logical asking whether the estimates should be returned as a list. By default (\code{asList = FALSE}) a named vector is returned.
}
}
\details{
MISSING
}
\value{
\code{estimateNe} returns either a named vector (default) or a list (if \code{asList = TRUE}) containing estimates of Ne applying all the different method(s) specified in \code{method}.
}
\references{
Jonas, A., Taus, T., Kosiol, C., Schlötterer, C. & Futschik, A.: Estimating effective population size from temporal allele frequency changes in experimental evolution, manuscript in preparation
}
\author{
Thomas Taus & Agnes Jonas
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wf.traj}}, \code{\link{sample.alleles}}, \code{\link{checkSNP}} and \code{\link{estimateWndNe}}
}
\examples{
# simulate allele frequency trajectories of individual loci using the Wright-Fisher model
Ne <- 1000
afData <- wf.traj(p0=runif(1000, 0, 1), Ne=2*Ne, t=c(0, 10, 20))
head(afData)

# introduce sampling error to mimic Pool-seq data
cov <- 200
beg <- sample.alleles(afData[,"F0"], size=cov, mode="coverage")
end <- sample.alleles(afData[,"F20"], size=cov, mode="coverage")


# estimate Ne after 20 generations using different methods
estNe <- estimateNe(p0=beg$p.smpld, pt=end$p.smpld, cov0=beg$size, covt=end$size, t=20, method=c("P_AF", "JR_planI", "JR_planII"), Ncensus=Ne)
estNe
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
