\name{Testing the Unified null hypothesis for a feature-set}

\alias{set.test}

\title{Feature-set testing}

\description{ Based on a given hommel object, returns an FWER-adjusted p-value for any feature-set of interest for the null hypothesis: "Is the proportion of active features less than c?" for a selected threshold c.
}

\usage{
set.test(hommel, ix, testype, testvalue)}

\arguments{
  \item{hommel}{The hommel object for the complete multiple testing problem.}
  \item{ix}{The selection of features defining the feature-set. Any selector that works on the original vector of p-values (index, negative index, logical, names) is allowed. If ix is missing, the total proportion  of active features in the dataset is estimated.}
  \item{testype}{Type of the test, namely "selfcontained" or "competitive". Choosing the self-contained option will set the threshold to zero. Choosing the competitive option without a testvalue, will set the threshold to the overall estimated proportion of true hypotheses.}
  \item{testvalue}{Setting this value to c along with testype="competitive" will lead to testing the null hypothesis against a threshold c. Note: this value needs to be a proportion.}

}

\value{Returns the adjusted p-value of the unified null hypothesis.}

\references{
Goeman, J., Meijer, R., Krebs, T., Solari, A. (2017). Simultaneous Control of All False Discovery Proportions in Large-Scale Multiple Hypothesis Testing. arXiv preprint  arXiv:1611.06739v2.}

\author{Mitra Ebrahimpoor.}

\examples{
  #Generate a vector of pvalues
  set.seed(102)
  m <- 10
  pvalues <- c(runif(0.5*m,0,0.02), runif(0.5*m,0,1))

  # First step: create a hommel object.
  hom<-hommel(pvalues )
  # Use Unified null to test the proportion of active features against c=0
  set.test (hom, 1:3, "selfcontained")

}

