\name{checkSNP}
\alias{checkSNP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtaine Appropriate SNPs for Ne Estimation
}
\description{
\code{checkSNP} helps to filter SNPs for those that are best suited for estimating Ne.
}
\usage{
checkSNP(p0, pt, cov0, covt, truncAF = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p0}{
numeric vector of starting allele frequencies.
}
  \item{pt}{
numeric vector of end allele frequencies.
}
  \item{cov0}{
numeric vector of sequence coverages at the earlier time point.
}
  \item{covt}{
numeric vector of sequence coverages at the later time point.
}
  \item{truncAF}{
numeric indicating whether SNPs with extreme allele frequencies should be excluded from Ne estimation. If \code{truncAF = NA} (default) even extreme allele frequencies will be considered. If \code{truncAF = F} then only SNPs with an allele frequency \code{>= F} and \code{<= 1-F} will be considered for Ne estimation.   
}
}
\details{
MISSING
}
\value{
A logical vector indicating which SNPs are suited for Ne estimation.
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
