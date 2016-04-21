\name{read.sync}
\alias{read.sync}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Sync-File Format
}
\description{
Reads a synchronized file generated e.g. by the \href{http://sourceforge.net/projects/popoolation2/}{Popoolation2} software and extracts biallelic count data.
}
\usage{
read.sync(file, gen, repl, rising = FALSE)
}
\arguments{
  \item{file}{
the name of the file, which the data is to be read from.
}
  \item{gen}{
numeric vector specifying the number of generations for each population, see 'Details'.
}
  \item{repl}{
numeric vector defining how populations should be grouped by replicate, see 'Details'.
}
  \item{rising}{
a logical value. If \code{FALSE} (default) allele counts are polarized for the minor allele. If \code{TRUE} allele counts are polarized for the rising allele.
}
}
\details{
MISSING
}
\value{
A list containing the following items:
  \item{chr}{named vector containing ID and name of each chromosome in the sync-file.}
  \item{alleles}{\code{data.table} that contains major and minor allele for each genome position.}
  \item{cnts}{list of numeric matrices with allele frequency count data, one matrix for each replicate (as specified in \code{repl}). Each matrix contains the chromosome ID (\code{chrID}), sequence position (\code{pos}), biallelic minor allele count (\code{Fx_cnt}) and sequence coverage (\code{Fx_cov}) for each SNP (and generation x). The format of this object is similar to the one returned by \code{\link{wf.traj}()}.}
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
#MISSING: include example dataset (returned after reading small sync file)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
