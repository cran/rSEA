## ----echo = FALSE, message = FALSE--------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(rSEA)

## ----data_sim-----------------------------------------------------------------

set.seed(753) #to get the same result as you repeat
names<-sapply(1:300, function(i) paste0( combn(letters, 2)[,i], collapse = "")) # a vector of two-letter names
p<-c(runif(100, 0, 1)^6,runif(200, 0, 1)) # a vector of random p-values, the first 100 are smaller
simdat<-data.frame(pvals=p, ids=names, stringsAsFactors=FALSE) #the dataframe
head(simdat) #take a look at the dataset

pathsize<-floor(runif(50, 10, 60)) #generate a vector of random pathway sizes between 10 and 60

mocklist<-lapply(pathsize,function(x) sample(names,size=x)) # create a pseudo list of pathways


## ----setTDP_chunk1------------------------------------------------------------
require(rSEA) #load rSEA
#setTDP(pvals, ids, data=simdat) 

setTDP(simdat$pvals, simdat$ids) 


## ----setTDP_chunk2------------------------------------------------------------
require(rSEA) #load rSEA

setTDP(simdat$pvals, simdat$ids, set=mocklist[[3]]) 


## ----setTest_chunk1-----------------------------------------------------------
require(rSEA) 

#selfcontained test of all features
setTest(pvals, ids, data=simdat, testype = "selfcontained")  

#default comp. test
setTest(pvals, ids, data=simdat, set=mocklist[[3]], testype = "competitive")  

#custom comp. test
setTest(pvals, ids, data=simdat, set=mocklist[[3]], testype = "competitive", testvalue = 0.5) 

## ----SEA_chunk1---------------------------------------------------------------
require(rSEA) #load rSEA
testchart1<-SEA(simdat$pvals, simdat$ids, pathlist = mocklist)
head(testchart1)

## ----SEA_chunk2---------------------------------------------------------------
require(rSEA) #load rSEA
testchart2<-SEA(pvals, ids, data=simdat, pathlist = mocklist, select =11:30)
testchart2


## ----setTest_chunk2-----------------------------------------------------------
require(rSEA) #load rSEA
lap<-union(mocklist[[8]], mocklist[[9]]) # the overlapping feature-set
length(lap)

setTDP(pvals, ids, data=simdat, set =lap) 
setTest(pvals, ids, data=simdat, set =lap , testype = "selfcontained")
setTest(pvals, ids, data=simdat, set =lap , testype = "competitive", testvalue = 0.1) 


## ----SEA_chunk3---------------------------------------------------------------
require(rSEA) #load rSEA
testchart3<-topSEA(testchart1) #sorted by large TDP.estimates
head(testchart3)

testchart4<-topSEA(testchart1, by=Comp.adjP, descending=FALSE) #sorted by smallest competitive adj.p-values 
head(testchart3)

sigchart<-topSEA(testchart1, by=Comp.adjP, thresh = 0.05) #keep only significant self-contained p-values
sigchart2<-topSEA(sigchart, by=Size, descending=TRUE) #sorted by pathway size
head(sigchart2)

## ----pathlist_chunk1, eval = FALSE--------------------------------------------
#  
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("hgu133a.db")
#  
#  ###if you are running an older version of R
#  #source("http://bioconductor.org/biocLite.R")
#  #biocLite(c("limma","edgeR"))

## ----pathlist_chunk2, eval = FALSE--------------------------------------------
#  library(hgu133a.db)
#  ls("package:hgu133a.db")
#  columns(hgu133a.db)
#  
#  gobimap<-toTable(hgu133aGO)
#  gobimap<-gobimap[gobimap$Ontology=="CC",]
#  
#  head(gobimap)
#  

## ----pathlist_chunk3, eval = FALSE--------------------------------------------
#  
#  GOIDs<-unique(gobimap$go_id)
#  GOList<-lapply(GOIDs,
#                 function(id)
#                     gobimap$probe_id[gobimap$go_id==id])
#  
#  names(GOList)<-GOIDs
#  #head(GOList)
#  
#  #make sure the Ids are unique and the paths are non-empty
#  GOList<-lapply(GOList, unique)
#  GOList <- lapply(GOList, function(path) if (all(is.na(path))) character(0) else path)
#  
#  #save(GOList, file="GOList.RData")

## ----pathlist_chunk4, eval = FALSE--------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("org.Mm.eg.db")
#  
#  library(org.Mm.eg.db)
#  ls("package:org.Mm.eg.db")
#  columns(org.Hs.eg.db)

## ----pathlist_chunk5, eval = FALSE--------------------------------------------
#  
#  
#  uniKeys <- keys(org.Hs.eg.db, keytype="ENTREZID")
#  cols <- c("SYMBOL", "PATH")
#  
#  keggbimap<-select(org.Hs.eg.db, keys=uniKeys,
#                    columns=cols, keytype="ENTREZID")
#  
#  head(keggbimap)
#  
#  keggIDs<-unique(keggbimap$PATH)
#  KEGGList<-lapply(keggIDs,
#                   function(id)
#                      keggbimap$ENTREZID[keggbimap$PATH==id])
#  names(KEGGList)<-keggIDs
#  head(KEGGList)
#  
#  #make sure the Ids are unique and the paths are non-empty
#  KEGGList<-lapply(KEGGList, unique)
#  KEGGList<-lapply(KEGGList, function(path)
#                    path[!is.na(path)])
#  KEGGList <- lapply(KEGGList,
#                     function(path)
#                         if (all(is.na(path))) character(0) else path)
#  
#  #save(KEGGList, file="KEGGList.RData")
#  

## ----pathlist_chunk6, eval = FALSE--------------------------------------------
#  
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("rWikiPathways")
#  

## ----pathlist_chunk7, eval = FALSE--------------------------------------------
#  
#  #Matching the IDs to create list of wikipathways for metabs
#  library(rWikiPathways)
#  
#  homoPathways<-listPathwayIds(organism="Homo sapiens")
#  wikiList<-lapply(homoPathways,
#                  function(x)
#                      getXrefList(pathway = x, systemCode="En"))
#  
#  names(wikiList)<-homoPathways
#  
#  #make sure the Ids are unique and the paths are non-empty
#  wikiList<-lapply(wikiList, unique)
#  wikiList <- lapply(wikiList, function(path) if (all(is.na(path))) character(0) else path)
#  
#  #save(wikiList, file="wikiList.RData")
#  

