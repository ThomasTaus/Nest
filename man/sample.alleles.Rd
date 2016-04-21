\name{sample.alleles}
\alias{sample.alleles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample alleles
}
\description{
Introduce sampling error in allele frequency data to mimic the Pool-seq approach. On one hand, subjecting only a subset of individuals in a population to Pool-seq is modeled with hypergeometric sampling (\code{mode = "individuals"}). On the other hand, sampling error introduced by sequencing only a fraction of all DNA fragments is modeled with binomial sampling (\code{mode = "coverage"}).
}
\usage{
sample.alleles(p, size, mode = c("coverage", "individuals"), Ncensus = NA, ploidy = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
numeric vector defining allele frequencies, which are used as success probabilities in the sampling process.
  }
  \item{size}{
numeric indicating the sample size to be used for binomial (\code{mode = "coverage"}) or Poisson sampling (\code{mode = "individuals"}). If \code{mode = "coverage"} and \code{length(size) == 1} then for each allele frequency an individual sequence coverage value will be drawn from a Poisson distribution with \code{lambda = size}. Otherwise (\code{length(size) > 1}) the values in \code{size} will be used directly and recycled if necessary. If \code{mode = "individuals"} then \code{size} has to be an integer specifying the number of individuals that were sampled from the population.
  }
  \item{mode}{
character string specifying the sampling mode. Possible values are \code{"coverage"} and \code{"individuals"}, see 'Details'.
  }
  \item{Ncensus}{
numeric specifying the census size of the entire population (before sampling).
  }
  \item{ploidy}{
numeric, the ploidy of the individuals.
  }
}
\details{
MISSING
}
\value{
A numeric vector of allele frequencies or (if \code{mode = "coverage"} and \code{length(size) == 1}) a list containing the following items 
  \item{p.smpld}{allele frequencies after sampling}
  \item{size}{sequence coverage for each position, drawn from a Poisson distribution with \code{lambda = size}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
