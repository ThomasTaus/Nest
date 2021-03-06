\name{sync}
\alias{sync-class}
\alias{is.sync}



\title{Sync Object - Creation and Usage}

\description{\code{sync} is an S4 class that stores all information of a synchronized file generated e.g. by the \href{http://sourceforge.net/projects/popoolation2/}{Popoolation2} software. }

\usage{
is.sync(x)
}

\arguments{
  \item{x}{
any \R object.
}
}

\details{
\code{is.sync} checks if an object is a \code{sync} object.
}

\value{
\code{is.sync} returns \code{TRUE} if its argument is a \code{sync} object and \code{FALSE} otherwise.
}

\author{
Thomas Taus
}

\seealso{
\code{\link{read.sync}}, \code{\link{alleles}}, \code{\link{af.traj}}, \code{\link{af}} and \code{\link{coverage}}
}
