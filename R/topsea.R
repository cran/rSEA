#' @title topSEA
#'
#' @description returns a permutation  of SEA-chart which rearranges
#' the feature-sets according to the selected argument into ascending or
#' descending order.
#'
#' @param object A SEA-chart object which is the output of \code{SEA} function
#'
#' @param by Variable name by which the ordering should happen. It should be a column of SEA-chart.
#' The default is TDP_bound.
#'
#' @param thresh A real number between 0 and 1. If specified the values of the variable defined in \code{by}
#' will be thresholded accordingly.
#'
#' @param n Integer. Number of raws of the output chart
#'
#' @param descending Logical. If \code{TRUE} The output chart is organized in a descending order
#'
#' @param cover An optional threshold for coverage, whcih muct be a real number between 0 and 1.
#' If specified, feature-sets with a coverage lower than or equal to this value are removed.
#'
#' @param digits Optional integer value. Number of decimal places to be shown in the output for numeric columns.
#' Default is the value defined by the local R options
#'
#' @param scientific Logical. Optional argument, if TRUE, the values are returned in scientific notofication as defined in
#' \code{\link{format}}  function. Default is \code{FALSE}
#'
#'
#' @return Returns a subset of SEA_chart sorted and refined according to the arguments
#'
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso \code{\link{SEA}}
#'
#' @references
#' Mitra Ebrahimpoor, Pietro Spitali, Kristina Hettne, Roula Tsonaka, Jelle Goeman,
#' Simultaneous Enrichment Analysis of all Possible Gene-sets: Unifying Self-Contained
#' and Competitive Methods, Briefings in Bioinformatics,bbz074
#'
#' @examples
#'
#' \dontrun{
#' ##Generate a vector of pvalues
#' set.seed(102)
#'
#' m<- 100
#' pvalues <- c(runif(0.5 * m,0,0.02), runif(0.5*m,0,1))
#' # Use Unified null to test the proportion of active features against c=0
#' set.test (hom, 1:3, "selfcontained")
#' }
#' @export
#'
#' @importFrom hommel hommel tdp localtest
#'
#'
topSEA=topsea<- function(object, by, thresh=NULL, descending=TRUE, n=20, cover,
                         digits= getOption("digits"), scientific=FALSE){

  #save the call to function
  cl <- match.call()

  #evaluate by argument
  if(missing(by)) byName <- "TDP.estimate"
  else byName <- deparse(cl$by)

  #check the arguments
    if(sum(colnames(object) %in% c("ID","Name","Size","Coverage","TDP.bound","TDP.estimate",
                                   "SC.adjP","Comp.adjP"))< 5)

       stop('Maybe the SEA-chart object is not specified correctly.')

    if(! byName %in% colnames(object))
          stop('The argument by should match the colnames of SEA-chart!')

    if( byName %in% c("ID","Name"))
          stop('The chart can not be reordered by ID or Name!')

  #remove low cover
      if(!missing(cover)){
        object<-with(object, object[round(coverage,4) >= cover,])
          if(nrow(object)==0)
             stop('Nothing selected, modify Cover value!')}

  #remove thoes below threshold
      if(!missing(thresh)){
        object <- with(object, object[object[,byName]>=thresh,])
          if(nrow(object)==0)
             stop('Nothing selected, modify threshold!')}

      if(descending==TRUE)
          SEA_top<-object[order(object[,byName], decreasing = TRUE),]

      if(descending==FALSE)
          SEA_top<-object[order(object[,byName], decreasing = FALSE),]

  #select the top
      if(nrow(object)>n)
        SEA_top<-SEA_top[1:n,]

      SEA_top<-round(SEA_top,digits)

      if(scientific==TRUE) SEA_top<-format(SEA_top, scientific=TRUE)

      SEA_top<-apply(SEA_top, 2, function(a) as.numeric(a))


  return(SEA_top)
}
