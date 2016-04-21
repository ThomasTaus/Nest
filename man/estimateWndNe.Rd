\name{estimateWndNe}
\alias{estimateWndNe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Effetive Population Size for Genome-Wide Data
}
\description{
\code{estimateNe} splits the genome-wide data into non-overlapping windows of size \code{wndSize} and estimates for each window the effective population size (Ne) based on tempral allele frequency changes.
}
\usage{
estimateWndNe(chr, pos, wndSize, p0, pt, cov0, covt, t, unit = c("bp", "SNP"), ploidy = 2, truncAF = NA, method = "P_AF2_planI", poolSize = NA, Ncensus = NA)
}
\arguments{
  \item{chr}{
vector of chromosome identifiers.
}
  \item{pos}{
numeric vector of sequence positions.
}
  \item{wndSize}{
integer indicating the window size to be used for Ne estimation. If \code{NA}, Ne will be estimated for entire chromosomes. The interpretation of this parameter dependes on how \code{unit} is specified.
}
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
  \item{unit}{
character string indicating the unit of \code{wndSize}. Possible values are \code{"bp"} (number of base pairs) and \code{"SNP"} (number of single-nucleotide polymorphisms).
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
}
\details{
MISSING
}
\value{
\code{estimateWndNe} returns an object of type \code{data.table} that contains Ne estimates for all windows of each chromosome using the method(s) specified in \code{method}. The resulting \code{data.table} contains the following columns:
  \item{Nx}{Ne estimate for the current window using method x, see 'Details'}
  \item{chr}{chromosome name}
  \item{wndStart}{start position of the current window}
  \item{wndStop}{stop position}
  \item{SNPs}{number of SNPs within the window}
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
\code{\link{estimateNe}} and \code{\link{read.sync}}
}
\examples{
#MISSING: use example dataset (returned after reading small sync file) and plot some results (barplot & genome-wide Ne)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
