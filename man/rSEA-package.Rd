\name{rSEA-package}
\alias{rSEA-package}
\docType{package}
\title{Simultaneous Enrichment Analysis (SEA) of all possible feature-sets using the unified null hypothesis
}
\description{
The unified null hypothesis is tested using closed testing procedure and all-resolutions inference. In short, using p-values of the individual features as input, the package can provide an FWER-adjusted p-value along with a lower bound and a point estimate for the proportion of true discoveries per feature-setas output. The flexibility in revising the choice of feature-sets without inflating type-I error is the most important property of SEA.
}
\details{
The starting point of the package is the hommel function. Then the created object can be used to test any feature-set of interest. Pre-defined subsets from GO, Reactome, and Wikipathways are accessible through package.
}
\author{
Mitra Ebrahimpoor.

Maintainer: Mitra Ebrahimpoor<m.ebrahimpoor@lumc.nl>
}
\references{
  Goeman, J., Meijer, R., Krebs, T., Solari, A. (2017). Simultaneous Control of All False Discovery Proportions in Large-Scale Multiple Hypothesis Testing. arXiv preprint  arXiv:1611.06739v2.
}

\keyword{ SEA, package }




