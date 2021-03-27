require(XML)
require(stringi)
require(factoextra)
require(NbClust)
require(Rtsne)
require(ggseqlogo)
require(png)
require(imager)
require(Biostrings)
require(ggplot2)
mim7lindb=xmlToList(xmlParse(file="D:\\Documents\\BG\\FNI2016\\work\\PavlinaPap\\mimodb\\AdvanceSearch_By_monoclonal antibody(TargetType)+7(LibrarySeqLength)+Linear(LibraryTopology).xml"), simplify = T)
m7lpep=lapply(mim7lindb[[1]], function(x){
  y=unlist(stri_extract_all_regex(x$Peptides,"[A-Z]\\w+"))
  y=y[sapply(y,nchar)==7]
  return(unlist(y))
})
m7lpep=m7lpep[lengths(m7lpep)>0]
names(m7lpep)=NULL
iis=sapply(seq_along(m7lpep), function(i){rep(i,length(m7lpep[[i]]))})
ls=lengths(m7lpep)
j=ls!=1
ls=lengths(m7lpep[j])
m7pep=data.frame(id=unlist(iis[j]),pep=unlist(m7lpep[j]), stringsAsFactors = F)
write.table(m7pep[,2], file="m7pep.txt", quote = F,row.names = F)
# Clean at http://immunet.cn/bdb/index.php/site/tools?type=TUPScan for parasitic sequences
clTUP=unlist(read.table(file="cleaned_1565435748.txt", header = F, stringsAsFactors = F))
m7pep=m7pep[m7pep[,2] %in% c(clTUP),]
j=table(m7pep[,1])
m7pep=m7pep[j>1,]
m7pep=unique(m7pep)

m7blsq=array(0,dim = c(8000,50))
pid=array(0, dim=c(8000,2))
co=1
for (i in 1:nrow(m7pep)) {
  for (j in 1:nrow(m7pep)){
    if ((i!=j) & ((m7pep[i,1]==m7pep[j,1])|(sample(50,1)==50))) {
      m7blsq[co,2:50]=c(blosq(m7pep[i,2],m7pep[j,2]))
      if (m7pep[i,1]==m7pep[j,1]) m7blsq[co,1]=1 else m7blsq[co,1]=0
      pid[co,]=c(i,j)
      co=co+1
    }
  }
  print(i)
}

print(which(rowSums(m7blsq)==0))
m7blsq=m7blsq[1:6321,]
m7blsq_in=m7blsq[m7blsq[,1]==1,2:50]
pid=pid[1:6321,]
pidin=pid[m7blsq[,1]==1,]

PCAm7blsq=prcomp(m7blsq[,2:50])
fviz_pca_biplot(PCAm7blsq,axes = c(5,6),geom = "point",col.ind = m7blsq[,1]+1)
PCAm7blsq_pos=prcomp(m7blsq_in)
fviz_pca_biplot(PCAm7blsq_pos,axes = c(3,4),geom = "point")
fviz_eig(PCAm7blsq_pos,ncp=45)

hcm7blsq=hclust(dist(m7blsq_in),method="ward.D2")
NbClm7blsq=NbClust(m7blsq_in, max.nc = 400, method = "ward.D2", index = "ccc")
plot(diff(NbClm7blsq$All.index))
print(which(diff(NbClm7blsq$All.index)>2))
clstrs=c(8,19,134,188)
bestcl=sapply(clstrs, function(n){cutree(hcm7blsq,k=n)})
cl8sums=aggregate(m7blsq_in, by=list(bestcl[,1]), FUN="mean")
cl8sums=cl8sums[,-1]
blsq8=lapply(1:nrow(cl8sums),function(i){
  s=matrix(unlist(cl8sums[i,]),7)
  image(s,main=i)
  return(s)
})

#peps10_2=pidin[bestcl[,1]==2,]

tsne_m7blsq_in=Rtsne(m7blsq_in, initial_dim=49)
for (i in 1:10){
  plot(tsne_m7blsq_in$Y,col=1+1*(bestcl[,1]==i), main=i)
}
