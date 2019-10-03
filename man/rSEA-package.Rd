\name{rSEA-package}
\alias{rSEA-package}
\docType{package}
\title{Simultaneous Enrichment Analysis (SEA) of all possible feature-sets using the unified null hypothesis
}
\description{
The unified null hypothesis is tested using closed testing procedure and all-resolutions inference. In short, using p-values of the individual features as input, the package can provide an FWER-adjusted p-value along with a lower bound and a point estimate for the proportion of true discoveries per feature-set as output. The flexibility in revising the choice of feature-sets without inflating type-I error is the most important property of SEA. For more details see: Ebrahimpoor, M (2019) <doi:10.1093/bib/bbz074>
}

\details{
The package uses raw p-values of features as input and evaluates any given set of hypotheses. For each selected set the adjusted p-value and TDP lower-bound are calculated. The type of test can be defined by arguments and can be refined as necessary. The p-values are corrected for every possible set of features, making the method flexible in ch0ice of set and test. A list of sets can also be passed to get simultaneous analysis of sets. Any genomic database can be adopted.
}
\author{
Mitra Ebrahimpoor.

Maintainer: Mitra Ebrahimpoor<m.ebrahimpoor@lumc.nl>
}
\references{
Goeman, J., Meijer, R., Krebs, T., Solari, A. (2017). Simultaneous Control of All False Discovery Proportions in Large-Scale Multiple Hypothesis Testing. arXiv preprint  arXiv:1611.06739v2.
}

\keyword{rSEA}
