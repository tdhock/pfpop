\name{pfpop}
\alias{pfpop}

\title{
  Piecewise Functional Pruning Optimal Partitioning
}
\description{
  Compute best change-points using a piecewise loss function, for a data
  sequence, and a penalty.
}
\usage{
pfpop(degrees_vec, penalty, weight_vec = rep(1, length(degrees_vec)))
}
\arguments{
  \item{degrees_vec}{
    Numeric vector of data to segment, degrees in [0,360).
}
  \item{penalty}{
    non-negative penalty, possibly Inf.
}
  \item{weight_vec}{
    Numeric vector of positive weights, used to determine the contribution of
  each data point toward the total loss. Default: weight of one for each
  data point.
}
}
\value{
Named list of two elements: segments is a data frame with one row per
  segment, iterations is a data frame with one row per iteration of
  dynamic programming.
}
\references{
  TODO
}
\author{
  Toby Dylan Hocking
}
\examples{

pfpop::pfpop(c(1.1, 5, 350), Inf)

pfpop::pfpop(c(10:12, 40:42, 80:82), 10)

}
