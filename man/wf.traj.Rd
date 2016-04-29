\name{wf.traj}
\alias{wf.traj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Random Genetic Drift
}
\description{
\code{wf.traj} simulates allele frequency trajectories based on the Wright-Fisher model.
}
\usage{
wf.traj(p0, Ne, t, s = 0, h = 0.5)
}
\arguments{
  \item{p0}{
numeric vector with starting allele frequencies.
}
  \item{Ne}{
integer specifying the desired effective population size of haploid individuals.
}
  \item{t}{
numeric vector defining for which generations allele frequencies should be returned.
}
  \item{s}{
  numeric specifying the selection coefficient used. By default (\code{s = 0}) no selection is applied.
  }
  \item{h}{
  numeric indicating the dominance coefficient.
  }
}
\details{
MISSING

If a custom function is specified in \code{startAF} it might be necessary to export some objects or functions to the cluster for parallelization. \code{exportVars} allows to the specify object names that need to be exported. See \code{\link{clusterExport}} in the \code{\link{parallel}} package for further details.
}
\value{
\code{wf.traj} returns a list of matrices (one for each replicate) that contain the simulated allele frequency data. In each matrix rows correspond to individual SNPs and the columns provide the following information for each SNP:
  \item{Fx_cnt}{allele counts after x generations of random genetic drift.}
  \item{Fx_cov}{target sequence coverage. Only available if re-sampling is performed to model Pool-seq data (\code{seqCov != NA}).}
  \item{Fx_sCnt}{re-sampled allele counts after x generations of drift. Only available if \code{seqCov != NA}.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Thomas Taus
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# simulate allele frequency trajectories of individual loci using the Wright-Fisher model
Ne <- 200
gen <- seq(0, 100, by=10)
afData <- wf.traj(p0=rep(0.5, times=500), Ne=Ne, t=gen)
# look at a subset of the generated data
head(afData)

# plot allele frequency trajectories
plot(1, type="n", xlim=c(0, max(gen)), ylim=c(0, 1), main="Allele frequency", xlab="Generation", ylab="Allele frequency change (\%)")
for(r in 1:nrow(afData)) {
  lines(gen, afData[r,])
}

# simulate allele frequency trajectories including selection
afData <- wf.traj(p0=rep(0.05, times=500), Ne=Ne, t=gen, s=0.1, h=0.5)

# plot results
plot(1, type="n", xlim=c(0, max(gen)), ylim=c(0, 1), main="Allele frequency", xlab="Generation", ylab="Allele frequency change (\%)")
for(r in 1:nrow(afData)) {
  lines(gen, afData[r,])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
