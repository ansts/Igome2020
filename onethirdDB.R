
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

uconf=umap.defaults
uconf$n_components=3
uconf$n_neighbors=15
uconf$verbose=T
utLinv=umap(tLinv)


#noise=array(rnorm(length(tLinv),0,1), dim = dim(tLinv))
cols=apply(tLinv,1,function(i){rgb(sum(i[c(1,3)])/sum(i),sum(i[2:3])/sum(i),0,0.3)})

proct=proc.time()
utLinv=umap(tLinv[1:3.4e5,]+noise[1:3.4e5,], config = uconf)
print(proc.time()-proct)

uconf$n_components=2
uconf$n_neighbors=15
proct=proc.time()
utLinv=umap(tLinv[1:1e6,]+noise[1:1e6,], config = uconf)
print(proc.time()-proct)

pdf("utlinv12.pdf", width=10, height=10)
plot(utLinv$layout, cex=0.3, pch=16, col=cols[1:300000])
dev.off()
plot3d(utLinv$layout, col=cols, cex=.3)
browseURL(paste("file://", writeWebGL(dir=file.path(getwd(), "webGL"), width=1500), sep=""))

tLinv=round(tLinv,2)
stLinv=apply(tLinv,1,paste, collapse="_")
ttLinv=aggregate(1:nrow(tLinv), by=list(tLinv[,1],tLinv[,2],tLinv[,3]), "length")

chsqt0=chisq.test(t0Linv, simulate.p.value = T)
pt0=p.adjust(2*pnorm(chsqt0$stdres, lower.tail = F))
pt0=array(pt0, dim=dim(chsqt0$stdres))
calls=apply(pt0,1,function(l){any(l<0.05)})
sum(calls)
pcalls=rownames(t0Linv[calls,])
ALcalls=Linv[pcalls]
ELcalls=as.matrix(melt(ALcalls))
Gcalls=graph_from_edgelist(ELcalls, directed = F)
Gcalls=simplify(Gcalls)
write.graph(Gcalls, file="Gcalls.graphml",format = "graphml")
write(paste(">",seq_along(pcalls),"\n",pcalls,"\n",collapse=""), "pcalls.txt")

call3_0=tLinv[(abs(chsqt0$stdres[,2])>2 | abs(chsqt0$stdres[,1])>2) & chsqt0$stdres[,3]<(-2),]
dim(call3_0)
calls3=rownames(call3_0)
cols3=rownames(utLinv$layout) %in% calls3
plot3d(utLinv$layout, col=cols3*1, cex=.3)
pdf("utlinv3_0.pdf", width=10, height=10)
plot(utLinv$layout, cex=0.3, pch=16, col=cols3*1+1)
dev.off()
plot3d(call3_0,  cex=.3)

calls3i=list(calls3[1:1370], calls3[1371:2742])
fcalls3=fl[glo_gbm %in% calls3]
names(fcalls3)=glo_gbm[glo_gbm %in% calls3]
fcalls3=fcalls3[calls3]


cl=makeCluster(2)
clusterExport(cl,list("glo_gbm","vstrdsum","calls3i"))
clusterExport(cl,list("fl"), envir = environment())
proct=proc.time()
calls3AL=parSapply(cl,calls3i, function(p){
  x=t(vstrdsum(p,glo_gbm, f=fl))
  return(list(x))
})
stopCluster(cl)
print(proc.time()-proct)

dummy=rbind(calls3AL[[1]], calls3AL[[2]])
calls3AL=dummy[,2]
calls3t3=t(array(unlist(dummy[,1]), dim=c(3,length(dummy[,1])), dimnames = list(1:3, names(dummy[,1]))))
plot(calls3t3[,1]/(calls3t3[,2]+0.5),t0Linv[calls3,1]/(t0Linv[calls3,1]+0.5))


calls3tAL=sapply(calls3AL,function(l){
  x0=rep(0,3)
  names(x0)=1:3
  x=table(fl[glo_gbm %in% l])
  x0[as.double(names(x))]=x
  return(x0)
})
calls3tAL=t(calls3tAL)
chsqcalls3t=chisq.test(calls3tAL)
pt1=p.adjust(2*pnorm(chsqcalls3t$stdres, lower.tail = F))
pt1=array(pt1, dim=dim(chsqcalls3t$stdres))
clls=apply(pt1,1,function(l){any(l<0.05)})
sum(clls)
clls=calls3tAL[clls,]
ranker=(chsqcalls3t$stdres[,2]-chsqcalls3t$stdres[,1])/abs(chsqcalls3t$stdres[,3])

LEinv=melt(Linv)
Gall1third=graph_from_edgelist(as.matrix(LEinv), directed = F)
Gall1third=simplify(Gall1third)
V=names(V(Gall1third))
Gall1third=set_vertex_attr(Gall1third, "Group", value=fl[glo_gbm %in% V])
save(Gall1third, file="Gall1third")
cGall1=components(Gall1third)
fcGall1=which(cGall1$csize>5 & cGall1$csize<1e6)
Gall1_frags=lapply(fcGall1,function(i){
  v=names(cGall1$membership)[cGall1$membership==i]
})
pG=table(fl)[2]/length(fl)
Gall1_frags_pG=sapply(seq_along(Gall1_frags), function(i){
  l=Gall1_frags[[i]]
  f=fl[glo_gbm %in% l]
  n=sum(f==2)
  print(c(n, length(l)))
  p=binom.test(n,length(l),pG)
  return(p$p.value)
})
Gall1_frags_pG=p.adjust(Gall1_frags_pG)

Gall1_3main=induced.subgraph(Gall1third, names(cGall1$membership[cGall1$membership==1]))
write.graph(Gall1_3main, file = "Gall1_3main.graphml", format = "graphml")
save(Gall1_3main, file = "Gall1_3main")

trGall1m=transitivity(Gall1_3main,type="local")
badTr=trGall1m==0 | is.na(trGall1m)
dGall1m=degree(Gall1_3main, mode="all")

pdf(file="TrxDeg_Gall1m.pdf", width = 10, height = 10)
plot((dGall1m[!badTr]),(trGall1m[!badTr]), cex=0.3)
dev.off()




#clqG1m=max_cliques(Gall1_3main,min=4,file = "G1_3main_cliques.txt")
#clqG1m=scan(file = "G1_3main_cliques.txt", what="", sep="\n")
clqG1m=max_cliques(Gall1_3main,min=4)
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

# goes to postFinal_3xchips




x=tLinv[!(rownames(tLinv) %in% final),]
x=x[sample(nrow(x),5e5),]
x=rbind(x,tLinv[final,])
noise=array(rnorm(length(x),0,1), dim = dim(x))
noise=sweep(noise, 2, (c(0.3051517,1.054103, 0.6368911 ))*0.2, "*")

x=x+noise[1:nrow(x),]

uconf=umap.defaults
uconf$n_components=3
uconf$n_neighbors=12
uconf$verbose=T
proct=proc.time()
utLinv=umap(x, config = uconf)
print(proc.time()-proct)

cols=apply(tLinv[rownames(utLinv$layout),],1,function(i){rgb(sum(i[2:3])/sum(i),sum(i[c(1,3)])/sum(i),0,0.3)})
pdf("utlinv12_finaloverlay.pdf", width=10, height=10)
xr=range(utLinv$layout[,1])
yr=range(utLinv$layout[,2])
plot(utLinv$layout, cex=0.1, pch=16, col=cols, xlim=xr, ylim=yr)
par(new=T)
plot(utLinv$layout[final,], cex=1,  col=1, xlim=xr, ylim=yr)
dev.off()


labs=apply(t0Linv[rownames(utLinv$layout),],1,function(l){paste(l,collapse="_")})
pdf("utlinv12_text.pdf", width=10, height=10)
xr=range(utLinv$layout[,1])
yr=range(utLinv$layout[,2])
ij=sample(nrow(utLinv$layout),10000)
plot(utLinv$layout[ij,], cex=0, col=cols, xlim=xr, ylim=yr)
text(utLinv$layout[ij,],labels = labs[ij], cex=0.05, col=cols[ij])
dev.off()

cols1=cols
cols1[500001:500600]=rgb(0,0,0,0.3)
plot3d(utLinv$layout, type="p", col=cols1, radius=0.1)   #c(rep(0.1,500000),rep(1,600))

uconf$n_components=2
proct=proc.time()
utLinv2d=umap(x, config = uconf)
print(proc.time()-proct)



ppinclq=unique(names(flclq))
jppin=rownames(utLinv2d$layout) %in% ppinclq
pdf("utlinv2d_incliques.pdf", width=10, height=10)
xr=range(utLinv2d$layout[,1])
yr=range(utLinv2d$layout[,2])
plot(utLinv2d$layout, cex=sz, pch=16, col=jppin*1+1, xlim=xr, ylim=yr)
dev.off()




