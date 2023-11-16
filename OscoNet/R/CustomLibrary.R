# Author: Luisa Cutillo
community.significance.test <- function(graph, Nodeids, ...) {
  if (is.directed(graph)) stop("The graph is not undirected")
  subgraph <- induced.subgraph(graph, Nodeids)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(graph, Nodeids) - in.degrees
  res=wilcox.test(in.degrees, out.degrees, alternative = "greater",paired=FALSE,exact=FALSE)
  return(res)
}


zerofilter<-function(counts,alpha){
  #this function will keep the genes that have at most alpha% zeros
  #
  #input:
  #counts, matrix of counts
  #output:Matrix with genes that have at least (1-alpha)% non zero values
  zerosperraw=rowSums(counts ==0)
  minZero=min(zerosperraw)
  maxZero=max(zerosperraw)
  Nc=dim(counts)[2]
  
  #alpha=0.3
  ###Keeping genes that have at most alpha % of zeros
  p=alpha*Nc
  DATA=counts[zerosperraw<p,]
  return(DATA) 
}

quantilefilter<-function(counts,alpha){
  #this function will remove all the genes (raws) for which the tails quantiles are the same: q(a)=q(1-a).The value a=0.95 is suggested.	
  #alpha= tails quantiles for the data cut-off
  DataNorm=as.matrix(counts)
  Q5 <- apply(DataNorm, 1, function(i) quantile(i, (1-alpha)))
  Q95 <- apply(DataNorm, 1, function(i) quantile(i, alpha))
  Rg <- Q95 - Q5
  id0=which(Rg==0)
  DataCut=DataNorm[-id0,]
  return(DataCut)
}

MVfilter<-function(counts,meancutlow,alpha){
  library(Oscope)
  MV<-CalcMV(Data= counts, Sizes= NULL, NormData=TRUE,MeanCutLow= meancutlow)
  #This is consistent with what we did in oscoNet 
  DataSubset<- counts[MV$GeneToUse,]####THIS SHOULD BE DATACUT!!!!
  return(DataSubset)
}


Allfilter<-function(counts,meancutlow, alpha){
  #apply all the filtering + [-1,1] range transformation
  library(Oscope)
  #countsq=quantilefilter(counts,alphacut)
  countsq=counts
  countsqz=zerofilter(countsq,alpha)
  countsqzmv=MVfilter(countsqz, meancutlow)
  countsnorm <- NormForSine(countsqzmv)
  return(countsnorm)
}


Cluster_filtering<-function(label,metadata,ALLDATA,alpha,meancutlow){
  id=which(metadata$x==label)
  genenames=ALLDATA[,1]
  counts=as.matrix(ALLDATA[,c(2:dim(ALLDATA)[2])])
  rownames(counts)= genenames
  allcases=colnames(counts)
  cellnames=metadata$X
  posid=match(cellnames[id], allcases)
  D=counts[, posid]
  Dnorm=Allfilter(D,meancutlow, alpha)
  return(Dnorm)
}