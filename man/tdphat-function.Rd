\name{TDP point estimate}

\alias{tdphat}

\title{TDP point estimate}

\description{ Provides a point estimate for the number of active features in the specified feature-set
}

\usage{
tdphat(hommel, ix)
}

\arguments{
  \item{hommel}{Hommel object. The hommel object for the complete dataset. For more details refer to hommel package.}
  \item{ix}{The selection of features defining the feature-set. Any selector that works on the original vector of p-values (index, negative index, logical, names) is allowed. If ix is missing, the total proportion  of active features in the dataset is estimated.}
}

\value{Returns the median point estimate of TDP.}

\references{
  Goeman, J., Meijer, R., Krebs, T., Solari, A. (2017). Simultaneous Control of All False Discovery Proportions in Large-Scale Multiple Hypothesis Testing. arXiv preprint  arXiv:1611.06739v2.
}

\author{Mitra Ebrahimpoor.}

\examples{
  #Generate a vector of pvalues
  set.seed(102)
  m <- 10
  pvalues <- c(runif(0.5*m,0,0.02), runif(0.5*m,0,1))

  # First step: create a hommel object.
  hom <- hommel(pvalues, simes = TRUE)

  # Get the median point estimate for the feature-set of interest:
  set <- 1:5
  tdphat(hom, set)
}
