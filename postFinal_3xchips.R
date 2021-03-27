# from onethirdDB

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

