\name{Top feature-sets}

\alias{sea.top}

\title{Top feature-sets from SEA_chart}

\description{ Provides a sorted subset of SEA_chart with selected properties.
}

\usage{
sea.top(SEA_chart, pathid=NULL, thresh=0.5, num=10, coverage=NULL)

}

\arguments{
  \item{SEA_chart}{The SEA chart output from sea.test function.}
  \item{pathid}{A vector of feature-set names.
  This is used to subset the SEA_chart to feature-sets of interest.}
  \item{thresh}{Any value which is [0,1]. This threshold is applied on estimated proportion of active features in each set. Feature-sets with smaller estimated TDP are removed.}
  \item{num}{Integer. Number of raws for the output.}
  \item{coverage}{Any value which is [0,1]. A threshold for the coverage of the feature-sets with your data. feature-sets with a coverage lower than the specified value are removed.}


}

\value{Returns a subset of SEA_chart sorted by the proportion of active features.}

\references{
  Goeman, J., Meijer, R., Krebs, T., Solari, A. (2017). Simultaneous Control of All False Discovery Proportions in Large-Scale Multiple Hypothesis Testing. arXiv preprint  arXiv:1611.06739v2.
}

\author{Mitra Ebrahimpoor}

\examples{
  #Generate a vector of pvalues
  set.seed(102)
  m <- 10
  pvalues <- c(runif(0.5*m,0,0.02), runif(0.5*m,0,1))

  # First create a hommel object:
  hom <- hommel(pvalues, simes = TRUE)

  # Estimate te proportion of active features in the set by:
  set <- 1:5
  tdphat(hom, set)
}
