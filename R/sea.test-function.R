sea.test<- function(hommel, geneid, database, pathid, source=NULL, tdp=TRUE, tdphat=TRUE, selfcontained=TRUE,
                    competitive=TRUE, thresh=NULL, no.name=TRUE ){

if (missing(database))
  stop('Database is not specified!')

if(database=="GO"){
  ifelse(!missing(source),pathlist<-GO[which(names(GO) %in% GO_names$ID[which(GO_names$Source==source)])],pathlist<-GO)
  ifelse(!missing(pathid),pathlist<-pathlist[which(names(pathlist) %in% pathid )],pathlist<-GO) }

if(database=="Reactome"){
  ifelse(!missing(pathid),pathlist<-Reactome[which(names(Reactome) %in% pathid )] ,pathlist<-Reactome)}

if(database=="Wiki"){
  ifelse(!missing(pathid),pathlist<-Wiki[which(names(Wiki) %in% pathid )] ,pathlist<-Wiki)}

#Some functions
tdpfunc<-function(x){ifelse(length(which(geneid  %in% x ))>0, tdp(hommel, geneid  %in% x),NA)}

tdp50func<-function(x){ifelse(length(which(geneid  %in% x ))>0, tdphat(hommel, geneid  %in% x),NA)}

SC_adjPfunc<-function(x){ifelse(length(which(geneid  %in% x ))>0,
                                set.test (hommel, which(geneid  %in% x ), testype="selfcontained"),NA)}

C_adjPfunc<-function(x){ifelse(length(which(geneid  %in% x ))>0,
                               set.test (hommel,testype="competitive", which(geneid  %in% x )),NA)}

thr_adjPfunc<-function(x){ifelse(length(which(geneid  %in% x ))>0,
                                 set.test(hommel,testype="competitive",testvalue=thresh,
                                               which(geneid  %in% x )),NA)}

#Filling in the SEA_chart
SEA_chart <-data.frame(matrix(vector(), length(names(pathlist)), 9, dimnames=list(c(),
                       c("ID","Name","Size","Coverage","TDP_bound","TDP_estimate","SC_adjP","C_adjP","custom_adjP"))))

SEA_chart$ID<-names(pathlist)

if(database=="GO" && no.name==FALSE) SEA_chart$Name<-GO_names$Name[which(GO_names$ID %in% SEA_chart$ID )]
if(database=="Reactome" && no.name==FALSE) SEA_chart$Name<-Reactome_names$Name[which(Reactome_names$ID %in% SEA_chart$ID )]
if(database=="Wiki" && no.name==FALSE) SEA_chart$Name<-Wiki_names$Name[which(Wiki_names$ID %in% SEA_chart$ID )]

SEA_chart$Size<-sapply(pathlist,FUN=function(x){length(x)},simplify = TRUE)

SEA_chart$Coverage<-sapply(pathlist,FUN=function(x){round(length(which(geneid  %in% x ))/length(x),2)},
                           simplify = TRUE)
if(tdp==TRUE)
  SEA_chart$TDP_bound<-sapply(pathlist,tdpfunc,simplify = TRUE)
else
  SEA_chart$TDP_bound<-NA

if(tdphat==TRUE)
  SEA_chart$TDP_estimate<-sapply(pathlist,tdp50func,simplify = TRUE)
else
  SEA_chart$TDP_estimate<-NA

if(selfcontained==TRUE)
  SEA_chart$SC_adjP<-sapply(pathlist,SC_adjPfunc,simplify = TRUE)

else
  SEA_chart$SC_adjP<-NA

if(competitive==TRUE)
  SEA_chart$C_adjP<-sapply(pathlist,C_adjPfunc,simplify = TRUE)
else
  SEA_chart$C_adjP<-NA

if(!missing(thresh))
  SEA_chart$custom_adjP<-sapply(pathlist,thr_adjPfunc,simplify = TRUE)
else
  SEA_chart$custom_adjP<-NA

SEA_chart <- Filter(function(x) !(all(x=="")), SEA_chart)

return(SEA_chart)

}
