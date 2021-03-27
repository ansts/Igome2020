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
require(fields)
require(reshape2)
require(igraph)
require(ggseqlogo)
require(ggplot2)
require(stringdist)

cpl=colorRampPalette(c("#0000A000","#0000FF00","#00FF0000","#FFFF0000","#FFFFFF00","#FF000000"))
# Load search results from BDB
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
dj=which(duplicated(m7pep$pep))
m7pep=m7pep[-dj,]
j=table(m7pep[,1])
jn=as.double(names(j))
jn=jn[j==1]
m7pep=m7pep[!(m7pep$id %in% jn),]

###### Recalculate everything from here#######


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
nc=NbClm7blsq$All.index
plot(diff(nc))
clp=diff(nc)
print(clp[clp>2])
clstrs=c(10,21,134,188,303)
bestcl=sapply(clstrs, function(n){cutree(hcm7blsq,k=n)})

cl10sums=aggregate(m7blsq_in, by=list(bestcl[,1]), FUN="mean")
cl10sums=cl10sums[,-1]
br=c(seq(from=min(cl10sums),to=0,by=(0-min(cl10sums))/50),seq(from=0.01,to=max(cl10sums),by=max(cl10sums)/100))
blsq10=lapply(1:nrow(cl10sums),function(i){
  s=matrix(unlist(cl10sums[i,]),7)
  
  fnm=paste("Pattern",i,".pdf", sep="", collapse="")
  pdf(file=fnm,width=8, height=8)
  opar=par(mai=c(2,2,2,2))
  
  image.plot(s,main=i,col=cpl(150), bigplot = c(0.2,0.8,0.18,0.7), breaks =br, nlevel=150, xaxt="n", yaxt="n")
  axis(1,at=c(0,0.1667,0.33,0.5,0.66667,0.83333,1), labels=c("P1.1","P1.2","P1.3","P1.4","P1.5","P1.6","P1.7"))
  axis(2,at=c(0,0.1667,0.33,0.5,0.66667,0.83333,1), labels=c("P2.1","P2.2","P2.3","P2.4","P2.5","P2.6","P2.7"))
  
  par(opar)
  dev.off()
  return(s)
})

m7pepcl=lapply(seq_along(unique(bestcl[,1])),function(i){
  ij=pidin[bestcl[,1]==i,]
  x=list(unique(m7pep[ij[,1],2]),unique(m7pep[ij[,2],2]))
  return(x)
})

dir.create("Cl10_logos")
for (i in seq_along(m7pepcl)){                                            
  png(file=paste("Cl10_logos\\logo_",i,"1.png",sep=""))
  seq=m7pepcl[[i]][1]
  a=ggseqlogo(seq)+ggtitle(paste("Cluster ",i,"1",sep=""))
  print(a)
  dev.off()
  png(file=paste("Cl10_logos\\logo_",i,"2.png",sep=""))
  seq=m7pepcl[[i]][2]
  a=ggseqlogo(seq)+ggtitle(paste("Cluster ",i,"2",sep=""))
  print(a)
  dev.off()
}

cl21sums=aggregate(m7blsq_in, by=list(bestcl[,2]), FUN="mean")
cl21sums=cl21sums[,-1]
br=c(seq(from=min(cl21sums),to=0,by=(0-min(cl21sums))/50),seq(from=0.01,to=max(cl21sums),by=max(cl21sums)/100))
blsq21=lapply(1:nrow(cl21sums),function(i){
  s=matrix(unlist(cl21sums[i,]),7)
  image.plot(s,main=i,col=cpl(150), bigplot = c(0.2,0.8,0.18,0.7), breaks =br, nlevel=150, xaxt="n", yaxt="n")
  axis(1,at=c(0,0.1667,0.33,0.5,0.66667,0.83333,1), labels=c("P1.1","P1.2","P1.3","P1.4","P1.5","P1.6","P1.7"))
  axis(2,at=c(0,0.1667,0.33,0.5,0.66667,0.83333,1), labels=c("P2.1","P2.2","P2.3","P2.4","P2.5","P2.6","P2.7"))
  
  return(s)
})


treee=rbind(table(bestcl[,1],bestcl[,2]),table(bestcl[,3],bestcl[,2]))
rn1=rownames(treee)[1:10]
rn2=colnames(treee)
rn3=rownames(treee)[11:144]
rownames(treee)[1:10]=paste("r1",rn1, sep=".")
rownames(treee)[11:144]=paste("r3",rn3, sep=".")
colnames(treee)=paste("r2",rn2, sep=".")
tredges=melt(treee)
tredges=tredges[tredges[,3]>0,]
s1=table(bestcl[,1])
s2=table(bestcl[,2])
s3=table(bestcl[,3])
names(s1)=paste("r1",names(s1), sep=".")
names(s2)=paste("r2",names(s2), sep=".")
names(s3)=paste("r3",names(s3), sep=".")
r3grsizes=c(s1,s2,s3)
r3graph=graph_from_edgelist(as.matrix(tredges[,1:2]), directed=F)
x=V(r3graph)$name
xs=r3grsizes[x]
V(r3graph)$size=xs*0.01

tkplot(r3graph, vertex.color=0,arrow.mode="-")



#peps10_2=pidin[bestcl[,1]==2,]

tsne_m7blsq_in=Rtsne(m7blsq_in, initial_dim=49)
for (i in 1:10){
  plot(tsne_m7blsq_in$Y,col=1+1*(bestcl[,1]==i), main=i)
}

mabsin=sapply(1:max(bestcl[,1]), function(cl){
  x=pidin[bestcl[,1]==cl,1]
  y=cbind(rep(cl,length(x)),m7pep[x,1])
  return(y)
})

for (i in 2:max(bestcl[,1])){
  mabsin[[1]]=rbind(mabsin[[1]],mabsin[[i]])
}
mabsin=mabsin[[1]]
tamabsin=table(mabsin[,1], mabsin[,2])
cpl1=colorRampPalette(c("#80808000","#0009FF00","#00FFA000","#A0FF0000","#FFFF0000","#FF000000"))
br=c(0,0.5,0.75,1,2,4,8,16,32,64,128,256,305)
image.plot(tamabsin,col=cpl1(12), bigplot = c(0.2,0.4,0.18,0.9), breaks =br, xaxt="n", yaxt="n", legend.shrink = 0.5, legend.line=3.5, legend.lab = "Number of MAb Pairs")
axis(1,at=(0:(max(bestcl[,1])-1))/(max(bestcl[,1])-1), las=2, labels=paste("Pattern",1:max(bestcl[,1])))
axis(2,at=(0:55)/55, las=2, cex.axis=0.75, labels=colnames(tamabsin))
title(ylab="Monoclonal Antibodies")

MAb=unique(m7pep[unique(c(pidin)),1])
mabsimx=lapply(MAb, function(m){
  i=which(m7pep[,1]==m)
  x=which(pidin[,1] %in% i & pidin[,2] %in% i)
  y=bestcl[x,1]
  xi=pidin[x,1]-min(pidin[x,1])+1
  xj=pidin[x,2]-min(pidin[x,2])+1
  z=cbind(xi,xj,y)
  z0=matrix(0,max(xi),max(xj)) 
  for (i in 1:nrow(z)) z0[z[i,1],z[i,2]]=z[i,3]
  
  image(z0, col=0:8)
  return(z0)
})
names(mabsimx)=MAb

# Preparing a scrambled db (by column)

m7pepmx=t(sapply(m7pep[,2],function(p){unlist(strsplit(p, split=""))}))
m7scrmb=sapply(1:100,function(s){
  m=apply(m7pepmx,2,sample)
  m=t(apply(m,1,function(l){paste(l,sep="",collapse="")}))
  return(m)
})

bsqrmb=apply(m7scrmb,2,function(x){
  clblsq=patex(x,s=30)[[2]]
  print(which(m7scrmb==x[1], arr.ind = T)[2])
  x=lapply(clblsq,function(i){
    y=mean(i[rank(i)<11])
    return(y)
  })
  clrt=clblsq[which.max(x)]
  return(clrt)
})
x=unlist(bsqrmb, recursive = F)
bsqrmba=array(unlist(x), dim=c(7,7,100))
x=apply(bsqrmba,c(1,2),max)

stdm7p=stringdistmatrix(m7pep$pep,method = "hamming")
stdm7p=as.matrix(stdm7p)
din=c()
dout=c()
for (i in 1:(length(m7pep$id)-1)){
  for (j in (i+1):length(m7pep$id)){
    if (m7pep$id[i]==m7pep$id[j]) din=c(din,stdm7p[i,j]) else dout=c(dout,stdm7p[i,j])
  }
}
