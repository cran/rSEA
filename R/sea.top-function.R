sea.top<- function(SEA_chart, pathid=NULL, thresh=0.5, num=10, coverage=NULL){
  
  if (!missing(pathid))
  SEA_chart<-subset(SEA_chart,SEA_chart$Set %in% pathid)
  
  if (!missing(coverage))
  SEA_chart<-subset(SEA_chart,SEA_chart$Coverage >=coverage)
  
  if(missing(SEA_chart))
    stop('SEA_chart not specified correctly!')
  
  SEA_chart<-subset(SEA_chart,SEA_chart$TDP_bound >=thresh)
  
  SEA_top<-SEA_chart[order(-SEA_chart$TDP_bound), , drop = FALSE]
  
  #SEA_top<-SEA_top[1:num,]
  
  return(SEA_top)
}