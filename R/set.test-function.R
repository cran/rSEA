set.test <- function(hommel, ix, testype, testvalue)
{

 m <- length(hommel@p)
  if (missing(ix)) {
    p = hommel@p
    n = m
  }
  else {
    p <- hommel@p[ix]
    n <- length(p)
  }

total_tdp<- ceiling(hommel::tdp(hommel)*n)/n

  if (any(is.na(p)))
    stop('NAs produced by selecting with ix')

  if (n == 0) {
    warning('empty selection')
    return(p=1, adjusted=1)
  }

  if(testype!="selfcontained" & testype!="competitive")
    stop("Test type is not correctly specifed")

  if (testype=="selfcontained")
    thr<-0
  if (testype=="competitive" & !missing(testvalue))
    thr<-testvalue*n
  if (testype=="competitive" & missing(testvalue))
    {
    thr=total_tdp*n
    }

  adjustedp <- hommel::localtest(hommel, ix, tdp=thr)

  return(adjustedp)
 }
