####Required packages in the following. Please uncomment if you need to install them.
#install.packages('xml2')
#install.packages('neat')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("biomaRt")

#####Description######
#This script runs a diagnostic and a clustering on the graphs produced as a result of OscoNet.
#You can specify which level of 'alpha' and which 'cases' you ran OscoNet on and are now ready for disgnostic. 
#Defaults are cases=c('D1','D2','D4','D5','D7','D8') and alpha=0.1
#For each case Dx:
# we provide the number of co-oscillating genes and the number of total genes in input
# we extract the communities using the Walktrap algorithm from package igraph. A plot of the Network with overlaied communities is provided in the folder Plots. 
# For each identified community we run a linearity test and flag the community accordingly with a 1 (default value for communitied not affected by the linearity is 0). only groups having within-group phase differences are further considered in order recovery. Specifically, for any pair of genes gi,gj within a group, we define the phase-shift residual as υgigj = min((π – mm_gigj,mm_gigj), mm_gigj), in which mm_gigj = (ψ_gigj mod π). Our default takes groups whose 80th quantile of υgigj's is greater than π/4 for further order recovery. Note that this is the same as the linearity filtering used in Oscope but we are estimating our ψ_gigj differently and also originally in Oscope the 90th quantile is used.
#Output: The output summary is stored in the folder ./Results in the file:
#'/Report',experiment,case,'_',alpha,'.csv' where again case is can be 'D1','D2','D4','D5','D7','D8' and alpha equals to the chosen value (0.1 or 0.2 in our exploration).
#This output is a data frame showing for each community, the gene names, number of nodes, number of edges, Significance (p-value of Wilcox test to compare the distribution of the in-degree against the out degree, details in OsconNet), Density (Number of edges divided by the total number of edges) ,RelativeDensity (number of edges divided by the maximum number of possible edges), Psi values, average psi values, variance of the psi, number of different psi (num_diff_psi) and the binary flag value for linearity (LinFlag).



rm(list=ls())
library('igraph')
library('neat')
customlib='/OscoNet/R/CustomLibrary.R'
source(customlib)
##############INPUT SPECIFICATIONS
alpha=0.2
experiment="YourDataset" # Change this to your datset name/experiment name/subpopulation name
cases=c('D0','D1','D2','D3','D4','D5','D6', 'D7', 'D8') #Change this to refelct the number of 'D' clusters you have in Input Directory

outdir='/OscoNet/YourDataset/Results'
plotroot='/OscoNet/YourDataset/Plots/Community'

if(!dir.exists('/OscoNet/YourDataset/Plots/')) {
  dir.create('/OscoNet/YourDataset/Plots/')}
if(!dir.exists(plotroot)) {
  dir.create(plotroot)}

if(!dir.exists(outdir)) {
  dir.create(outdir)}

########################
for(case in cases){
	print(case)
  plotdir=file.path(plotroot,paste0('Hong_',experiment,case))
  if(!dir.exists(plotdir)) {
    dir.create(plotdir)}
#datafile= file containing the preprocessed data, output of the filtering step, before oscoNet
datafile=paste0('/OscoNet/YourDataset/OsconetInput/filter',experiment,case,'_',alpha,'.csv')
#filename= graph (network) estimated by OscoNet -Only significant co-oscillators
filename=paste0('/OscoNet/YourDataset/OsconetOutput/filter',experiment,case,'_',alpha,'.out.csv')
#filepsi=Matrix of all the estimated psi for the genes in datafile. 
filepsi=paste0('/OscoNet/YourDataset/OsconetOutput/filter',experiment,case,'_',alpha,'.psi_ng.csv')

print(filename)
print(filepsi)

dataorig=read.csv(file=datafile, header= TRUE)
datafilt=read.csv(file=filename, header= TRUE)
psi=read.csv(file= filepsi, header= FALSE)
genes=dataorig[,1]
rownames(psi)=genes
colnames(psi)=genes
mpsi=as.matrix(psi)
mpsi1=mpsi
mpsi1[mpsi1==0]=0.0001
psiG=graph_from_adjacency_matrix(mpsi1,mode="undirected", diag=FALSE,weighted=TRUE)
#> E(psiG)$weight
E(psiG)$psi= E(psiG)$weight
psiG=delete_edge_attr(psiG,'weight')
Gf=graph_from_data_frame(datafilt,directed=FALSE)
#E(psiG)[E(G)]$weight
#G=set.edge.attribute(Gf,'psi',index=E(Gf),value=E(psiG)[E(Gf)]$weight)
G=intersection(psiG,Gf,keep.all.vertices = FALSE)
#ww=get.edge.attribute(G,'cost')
###########
#MISSING GENEs!!!!
fD=datafile#paste0('../OsconetInput/filter',case,'_',alpha,'.csv')
ID=read.csv(file=fD, header= TRUE)
gn=get.vertex.attribute(G,'name')
og=length(ID$X)
rg=length(gn)
mg=length(setdiff(ID$X,gn))
message=paste0('co-oscillating genes: ', rg, '; missing ',mg,' genes,', ' out of ', og)
print(message)
#imiss=which(ID$X==setdiff(ID$X,gn))
#ID[imiss,]



###############WALKTRAP COMMUNITY 
wt=walktrap.community(G)
modularity(wt)
mwt=membership(wt)
sz=sizes(wt)
#sz[which(sz>=10)]
Mwt=as.matrix(membership(wt))
fn=paste('Full_Plot',experiment,case,'_',alpha,'.pdf',sep='')
pdf(file=file.path(plotdir,fn))
plot(wt,G,vertex.label.dist=1,vertex.size=3,vertex.label.cex=0.5)
dev.off()
#genes=read.delim(file.path('../casestudy/Whitfield/ImageIdConversions.txt'),sep="\t",header=TRUE,quote="")
#details=genes[,c(1:3)]
#colnames(details)=c('CloneID','Ref','GeneSymbol')
colnames(Mwt)='CommunityID'
dg=degree(G,rownames(Mwt))
Mwtdg=merge(Mwt,dg,by.x=0,by.y=0)
colnames(Mwtdg)=c('GeneName','CommunityID','Degree')
outfile=paste0(outdir,'/','MGH_',experiment,case,'_',alpha,'Comm.csv')
write.csv(Mwtdg,file=outfile)
##########################################
IdGE10=which(sz>1)
n10=length(IdGE10)
NumberOfEdges=matrix('NA',n10,1)
CommDensity=matrix('NA',n10,1)
RelativeCommDensity=matrix('NA',n10,1)
PsigCom=matrix('NA',n10,1)
actpsi=matrix('NA',n10,1)
av_actpsi=matrix('NA',n10,1)
var_actpsi=matrix('NA',n10,1)
num_diff_psi=matrix('NA',n10,1)
actgenes=matrix('NA',n10,1)
TotalEdges=ecount(G)
V=rep('NA',n10)
qt=0.9
thre=pi/4
for (i in c(1:n10)){
  actCom=IdGE10[i]
  actNodes=which(Mwt==actCom)
  actNumVertex=length(actNodes)
  subG=induced_subgraph(G,actNodes)
  fn=paste(case,'_',alpha,'_Comm', actCom,'.pdf',sep='')
  pdf(file=file.path(plotdir,fn))
  plot(subG,vertex.label.dist=2,vertex.size=6,vertex.label.cex=1)
  title(paste0(case,' ',alpha,' Comm', actCom,sep=''))
  dev.off()
  actEdge=E(subG)
  actNumEdges=ecount(subG)
  NumberOfEdges[i]=actNumEdges
  CommDensity[i]=actNumEdges/TotalEdges
  RelativeCommDensity[i]=2*actNumEdges/(actNumVertex^2-actNumVertex)
  temppsi=E(subG)$psi
  actpsi[i]=paste(temppsi,collapse=',')
  actgenes[i]=paste(rownames(Mwt)[actNodes],collapse=',')
  num_diff_psi[i]=dim(table(temppsi))
  av_actpsi[i]=mean(temppsi)
  var_actpsi[i]=var(temppsi)
  res=community.significance.test(G, actNodes)
  PsigCom[i]=res$p.value
  #%%%%%%%Modulus operation in each cluster%%%%%
  mm= temppsi%%pi
  M0=abs(mm-0)
  M90=abs(mm-pi)
  out=ifelse(M0<M90,M0,M90)
  #%%%%%%%Flag the clusters for linearity%%%%%%%
  # print(paste0('Checking Community=', actCom))
  tmpNotNa= out[which(!is.na(out))]
  #print(summary(tmpNotNa))
  V[i]=quantile(tmpNotNa,qt)
  
}
####Diagnostic on Linearity
Vmark=which(V<thre)
for (kk in (1:length(Vmark))){
  message("flagged Cluster", Vmark[kk])}
Flag=rep(0,n10)
Flag[which(V<thre)]=1

sizes10=sz[IdGE10]

Report10=data.frame(Community=names(sizes10),Genes= actgenes,NumberOfNodes= as.vector(sizes10), NumberOfEdges=NumberOfEdges,Significance=PsigCom,Density=CommDensity,RelativeDensity=RelativeCommDensity, Psi= actpsi,av_psi=av_actpsi, var_psi= var_actpsi, num_diff_psi= num_diff_psi,LinFlag=Flag)
#Report10sig=Report10[which(as.numeric(as.vector(Report10$Significance))<0.01),]
###########
###basicInfo contains the basic info for the Significant communities:
#ID, Significance, Number of nodes, Density and relative density
#we have 8 communities (25  38  61  70  74  92  124 139) significant
#basicInfo=Report10sig[,c(1:5)]#data.frame(sizes10,NumberOfNodes=NumberOfNodes,NumberOfEdges=NumberOfEdges,Significance=PsigCom)
#valCC=merge(outCC,basicInfo,by.y='Community',by.x=1)
#AAA=merge(basicInfo,outCC,by.x='Community',by.y=1,all.x=TRUE)

#Table1=data.frame(Size=as.vector(AAA$NumberOfNodes),CC=as.vector(AAA$'size CC Overlapping'),Density=as.vector(AAA$Density),Realative_density=as.vector(AAA$RelativeDensity),Neat=as.vector(AAA$adjusted_p),Hypergeometic=as.vector(AAA$'CC sig BH corr'), Significance=AAA$Significance)

Report10 = Report10[order(Report10$Significance,decreasing = FALSE),]

#Report10
print('Communities identified:')
print(length(Report10$Community))

ToCutComm=Report10$Community[Report10$LinFlag==1]
print('Communities flagged as linear:')
print(length(ToCutComm))
print('ids:')
print(ToCutComm)
#communities to retain, removing the linear ones:
ToPlotComm=Report10$Community[Report10$LinFlag==0]
print('Communities to retain:')
print(length(ToPlotComm))
print('ids:')
print(ToPlotComm)

write.csv(Report10,file=paste(outdir,'/Report',experiment,case,'_',alpha,'.csv',sep=''))


pdf(file=paste(plotdir,'/Cut_Plot',experiment,case,'_',alpha,'.pdf',sep=''))

g3<- delete_vertices(G, V(G)[wt$membership %in% ToCutComm])
#wtcut=wt
#wtcut$membership=wt$membership[ V(G)[wt$membership %in% ToPlotComm]]
#wtcut$names=wt$names[ V(G)[wt$membership %in% ToPlotComm]]
#wtcut$count=length(which(wt$membership %in% ToPlotComm))

cutmembership=wt$membership[ V(G)[wt$membership %in% ToPlotComm]]
cutnames=wt$names[ V(G)[wt$membership %in% ToPlotComm]]
cutcount=length(which(wt$membership %in% ToPlotComm))
wtcut=make_clusters(g3, membership = cutmembership, algorithm = 'walktrap cut',
                    merges = NULL, modularity = FALSE)
#plot(wt,G,vertex.label.dist=1,vertex.size=3,vertex.label.cex=0.5)
plot(wtcut,g3,vertex.label.dist=1,vertex.size=3,vertex.label.cex=0.5)
dev.off()


}

## capture all the output to a file.


