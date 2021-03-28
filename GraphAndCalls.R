require(rgl)
require(car)
require(reshape2)
require(umap)
require(igraph)
require(parallel)
require(qualV)
require(PST)
require(stringdist)

pth="E:\\Documents\\BG\\FNI2016\\work\\NewSeqs"
d=dir(pth)
d=d[grep("bigAL_",d)]
load(paste(pth,d[[1]],sep="\\"))
L0=L
for (f in d[2:7]){
  load(paste(pth,f[[1]],sep="\\"))
  L0=c(L0,L)
}
rm(L)     

load(paste(pth,"glo_gbm_ordered", sep="\\"))
load(paste(pth,"fl_ordered", sep="\\"))

proct=proc.time()
x=melt(L0)
Linv=aggregate(x[,2], by=list(x[,1]), "list")
print(proc.time()-proct)
save(L0, file ="L0" )
rm(L0)
x=Linv[,2]
names(x)=Linv[,1]
Linv=x
pp=names(Linv)
fl1_3=fl[glo_gbm %in% pp]
gg1_3=glo_gbm[glo_gbm %in% pp]
names(fl1_3)=gg1_3
Linv=Linv[gg1_3]
save(Linv, "Linv")
x=unlist(Linv)
xl=lengths(Linv)

Linvf=fl1_3[x]

proct=proc.time()
j=1
co=1
y=array(0, dim = c(length(xl),3))
colnames(y)=1:3
for (i in xl) {
  tf=table(Linvf[j:(j+i-1)])
  y[co,names(tf)]=tf
  j=j+i
  co=co+1
}
print(proc.time()-proct)

rownames(y)=names(Linv)

t0Linv=y
rm(y)

tLinv=sweep(t0Linv, 2, colMeans(tLinv), "/")

LEinv=melt(Linv)
Gall1third=graph_from_edgelist(as.matrix(LEinv), directed = F) #partial graph of the mimotope library
Gall1third=simplify(Gall1third)
V=names(V(Gall1third))
Gall1third=set_vertex_attr(Gall1third, "Group", value=fl[glo_gbm %in% V])  # add a vertex attribute indicating the GBm, healthy or common category
save(Gall1third, file="Gall1third")
cGall1=components(Gall1third)
fcGall1=which(cGall1$csize>5 & cGall1$csize<1e6)
Gall1_frags=lapply(fcGall1,function(i){
  v=names(cGall1$membership)[cGall1$membership==i]
})
pG=table(fl)[2]/length(fl)
Gall1_frags_pG=sapply(seq_along(Gall1_frags), function(i){   # isolate the small components of the graph 
  l=Gall1_frags[[i]]                                         # and test if the number of GBM igome sequences 
  f=fl[glo_gbm %in% l]                                       # is significantly different from the overall distribution
  n=sum(f==2)
  print(c(n, length(l)))
  p=binom.test(n,length(l),pG)
  return(p$p.value)
})
Gall1_frags_pG=p.adjust(Gall1_frags_pG)                # no small fragment has significantly different number of GBM mimotopes

Gall1_3main=induced.subgraph(Gall1third, names(cGall1$membership[cGall1$membership==1])) # take the giant component
write.graph(Gall1_3main, file = "Gall1_3main.graphml", format = "graphml")             
save(Gall1_3main, file = "Gall1_3main")

clqG1m=max_cliques(Gall1_3main,min=4)       # find all maximal cliques in the giant component
save(clqG1m, file="clqG1m")

length(clqG1m)
table(lengths(clqG1m))
length(unique(unlist(clqG1m)))

lclqG1m=lengths(clqG1m)
flclq=unlist(sapply(seq_along(lclqG1m), function(i){rep(i,lclqG1m[i])}))
clqG1mn=split(names(unlist(clqG1m)), flclq)
names(flclq)=unlist(clqG1mn)
names(fl)=glo_gbm
fflclq=fl[names(flclq)]
tff=t(table(fflclq,flclq))

clGBM=tff[,2]>((tff[,1]+tff[,3])*3)
clGBM=clqG1mn[clGBM]
clGBM=clGBM[order(lengths(clGBM),decreasing = T)]

clovrlp=tff[,1]==tff[,2]
clovrlp=names(clovrlp[clovrlp])
clovrlp=clqG1mn[clovrlp]
clovrlpup=clovrlp[lengths(clovrlp)>6]

clnonGBM=clqG1mn[(tff[,2]+tff[,3])==0]
clnonGBMup=clnonGBM[lengths(clnonGBM)>4]
GBMcalls=bestSeq(clGBM)
GBMclls_d=stayAway(GBMcalls,d=5)

x0=unique(unlist(clnonGBMup))
x=Linv[x0]
x=lapply(x, function(l){
  l[!l %in% x0]
})
x=x[lengths(x)>0]
allx=unique(unlist(x))
flax=fl[allx]
t0=t(sapply(x,function(l){ # table of neighbors outside of the cliques
  x=rep(0,3)
  names(x)=1:3
  i=table(flax[l])
  x[names(i)]=i
  return(x)
}))
x=names(x) # x0 which have neighbors outside of the clique
x1=t(sapply(clnonGBMup, function(ii) {  # summarize distribution of neighbors 
  j=ii[ii %in% x]                       # outside of the cliques for each clique
  if (length(j)>1) colSums(t0[j,]) else t0[j,]
}))
critnG=x1[,1]/(x1[,2]+x1[,3]+0.5)
clnonGBMups=clnonGBMup[order(critnG, decreasing = T)]
clnonGBMups=clnonGBMups[1:1000]
clnonGBMups=clnonGBMups[order(lengths(clnonGBMups), decreasing = T)]
nonGBMcalls=bestSeq(clnonGBMups)
l=length(GBMclls_d)
GBMnonGBMclls_d=stayAway(c(GBMclls_d,nonGBMcalls),d=5,startat=l+1)
length(GBMnonGBMclls_d)

l=length(GBMnonGBMclls_d)
x0=unique(unlist(clovrlpup))
x=Linv[x0]
x=lapply(x, function(l){
  l[!l %in% x0]
})
x=x[lengths(x)>0]
allx=unique(unlist(x))
flax=fl[allx]
t0=t(sapply(x,function(l){
  x=rep(0,3)
  names(x)=1:3
  i=table(flax[l])
  x[names(i)]=i
  return(x)
}))
x=names(x)
x1=t(sapply(clovrlpup, function(ii) {
  j=ii[ii %in% x]
  colSums(t0[j,])
}))
corf=colSums(tff)/sum(tff)
critclo=abs(log((x1[,3]/corf[3]+x1[,2]/corf[2])/(x1[,3]/corf[3]+x1[,1]/corf[1]+0.5)))
clovrlpups=clovrlpup[order(critclo)]
clovrlpups=clovrlpups
clovrlpups=clovrlpups[order(lengths(clovrlpups), decreasing = T)]
clovrlpcalls=bestSeq(clovrlpups)

clovrlGBMnonGBMclls_d=stayAway(c(GBMnonGBMclls_d,clovrlpcalls),d=5,startat=l+1)
length(clovrlGBMnonGBMclls_d)

x=rep(0, length(clovrlGBMnonGBMclls_d))
x[clovrlGBMnonGBMclls_d %in% GBMclls_d]="GBM"
x[clovrlGBMnonGBMclls_d %in% nonGBMcalls]="nonGBM"
x[clovrlGBMnonGBMclls_d %in% clovrlpcalls]="Common"

seqtotest=cbind(clovrlGBMnonGBMclls_d,x)
seqtoorder=rbind(seqtotest, seqtotest)

i=seq_along(seqtoorder[,1])
i=sample(i)
seqtoorder=seqtoorder[i,]
colnames(seqtoorder)=c("Sequence","Group")
write.csv(seqtoorder, file="seqtoorder.csv")

#Visualization

x=tLinv[!(rownames(tLinv) %in% seqtotest[,1]),]
x=x[sample(nrow(x),5e5),]
x=rbind(x,tLinv[clovrlGBMnonGBMclls_d,])
noise=array(rnorm(length(x),0,1), dim = dim(x))
noise=sweep(noise, 2, (c(0.3051517,1.054103, 0.6368911 ))*0.2, "*")

uconf=umap.defaults
uconf$n_components=3
uconf$n_neighbors=12
uconf$verbose=T
proct=proc.time()
utLinvf=umap(x+noise, config = uconf)
print(proc.time()-proct)

cols=apply(tLinv[rownames(utLinvf$layout),],1,function(i){rgb(sum(i[2:3])/sum(i),sum(i[c(1,3)])/sum(i),0,.7)})         #sum(i)/47
cols1=cols
x=rep(0, length(clovrlGBMnonGBMclls_d))
x[clovrlGBMnonGBMclls_d %in% GBMclls_d]=1
x[clovrlGBMnonGBMclls_d %in% nonGBMcalls]=2
x[clovrlGBMnonGBMclls_d %in% clovrlpcalls]=3

cols1[x==2]=rgb(0.1,0.5,0,1)
cols1[x==1]=rgb(0.75,0,0,1)
cols1[x==3]=rgb(1,1,0.3,1)

plot3d(utLinvf$layout, col=cols)
plot3d(utLinvf$layout[clovrlGBMnonGBMclls_d,], add=T, col=cols1[clovrlGBMnonGBMclls_d], type="s", size=0.5)
browseURL(paste("file://", writeWebGL(dir=file.path(getwd(), "webGL"), width=1500), sep=""))

uconf$n_components=2
proct=proc.time()
utLinv2df=umap(x, config = uconf)
print(proc.time()-proct)


sz=rowSums(st0Linv)/max(rowSums(st0Linv))
grffl=as.factor(groupingfinal2)
grffl=as.double(grffl)
grffl[grffl==2]=rgb(0.1,0.5,0,1)
grffl[grffl==1]=rgb(0.5,0,0.2,1)
grffl[grffl==3]=rgb(0.8,0.8,0,1)
pdf("utlinv2d_finaloverlay.pdf", width=10, height=10)
xr=range(utLinv2df$layout[,1])
yr=range(utLinv2df$layout[,2])
plot(utLinv2df$layout, cex=0.1, pch=16, col=cols, xlim=xr, ylim=yr)
par(new=T)
plot(utLinv2df$layout[allsel,], cex=1,  col=grffl, xlim=xr, ylim=yr)
dev.off()

