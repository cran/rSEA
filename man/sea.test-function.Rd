\name{Simultaneous Testing of Many Feature-sets}

\alias{sea.test}

\title{Simultaneous Testing of Many Feature-sets}

\description{ Based on a given hommel object, returns an FWER-adjusted p-value for the unified null hypothesis corresponding to a list of feature-sets from a chosen database.
}

\usage{
sea.test(hommel, geneid, database, pathid, source=NULL,

tdp=TRUE, tdphat=TRUE, selfcontained=TRUE, competitive=TRUE, thresh=NULL, no.name=TRUE)}

\arguments{
  \item{hommel}{Hommel object. The hommel object for the complete dataset. For more details refer to hommel package.}
  \item{geneid}{A vector. List of all features that are included in the analysis.}
  \item{database}{A list of pathway lists. Current options are GO, Reactome, or Wikipathways. Any other list can be provided, check the databases in the package for correct structure.}
  \item{pathid}{A vector. List of pathways of interest from the database of choice. If missing, all pathways of the database will be included.}
  \item{source}{One of the following options are possible:  molecular_function, biological_process or cellular_component. In case GO is selected as the database, this argument will subset the GO data base by your choice. In case no source is not specified, all GO terms are included. This argument is ignored for other databases.}
  \item{tdp}{Logical. If true, the lower_bound of the confidence interval for the proportion of true discoveries for each pathway will be estimated.}
  \item{tdphat}{Logical. If true, the point estimate for the proportion of true discoveries for each pathway  will be calculated.}
  \item{selfcontained}{Logical. If true, the self-contained null hypothesis will be tested for each pathway  and the corresponding adj. p-value is returned.}
  \item{competitive}{Logical. If true, competitive null hypothesis will be tested for each pathway  with default threshold and the corresponding adj. p-value is returned.}
  \item{thresh}{Logical. If specified, the unified null hypothesis will be tested for each pathway with the specified threshold and the corresponding adj. p-value is returned.}
  \item{no.name}{Logical. If TRUE, the names of the feature-sets are not published in the SEA chart.}
}


\value{Returns the SEA-chart containing a list of pathways with corresponding TDP bound estimate, TDP point estimate and adjusted p-values for the specified tests.}

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
  hom<-hommel(pvalues)
  # Use Unified null to test the proportion of true discoveries against 0
  set.test (hom, 1:3, "selfcontained")

}
