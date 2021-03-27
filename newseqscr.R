rm(list=ls())
ls()

require(stringi)
require(reshape2)

cpl=colorRampPalette(c("#00004000","#00FFFF00","#00FF0000","#FFFF0000","#FFFFFF00","#FF000000"))
cpl1=colorRampPalette(c("#0000AF00","#FF000000"))

folders=dir("E:\\Documents\\BG\\FNI2016\\work\\OsloSeq2020\\Fastq_Alpha-GBM-batch2")
foli=sapply(stri_locate_all(folders, fixed="_UNI"), function(i){!is.na(i[1])})
folders=folders[foli]

GBMseq=list()
for (d in folders) {
  qf=read.table(paste("E:\\Documents\\BG\\FNI2016\\work\\OsloSeq2020\\Fastq_Alpha-GBM-batch2",d,"uniqueP_QF.txt", sep = "\\", collapse=""), stringsAsFactors=F, header=F)
  qr=read.table(paste("E:\\Documents\\BG\\FNI2016\\work\\OsloSeq2020\\Fastq_Alpha-GBM-batch2",d,"uniqueP_QR.txt", sep = "\\", collapse=""), stringsAsFactors=F, header=F)
  GBMseq=c(GBMseq, qf=list(qf),qr=list(qr))
}
n=paste(names(GBMseq),rep(1:16, each=2), sep="")
names(GBMseq)=n

folders=dir("E:\\Documents\\BG\\FNI2016\\work\\OsloSeq2020\\Fastq_Alphal-vlgM-batch2\\Fastq_Alphal-vlgM-batch2")
foli=sapply(stri_locate_all(folders, fixed="_UNI"), function(i){!is.na(i[1])})
folders=folders[foli]
IVIgMseq=list()
for (d in folders) {
  qf=read.table(paste("E:\\Documents\\BG\\FNI2016\\work\\OsloSeq2020\\Fastq_Alphal-vlgM-batch2\\Fastq_Alphal-vlgM-batch2",d,"uniqueP_QF.txt", sep = "\\", collapse=""), stringsAsFactors=F, header=F)
  qr=read.table(paste("E:\\Documents\\BG\\FNI2016\\work\\OsloSeq2020\\Fastq_Alphal-vlgM-batch2\\Fastq_Alphal-vlgM-batch2",d,"uniqueP_QR.txt", sep = "\\", collapse=""), stringsAsFactors=F, header=F)
  IVIgMseq=c(IVIgMseq, qf=list(qf),qr=list(qr))
}
n=paste(names(IVIgMseq),rep(1:16, each=2), sep="")
names(IVIgMseq)=n
comx=sapply(GBMseq, function(g){
        sapply(IVIgMseq, function(i){
          length(intersect(g[,1],i[,1]))
        })
})
image(comx)

GBMseq=lapply(1:16, function(i){
  x=rbind(GBMseq[[i*2-1]],GBMseq[[i*2]])
  x=aggregate(x[,2], by=list(x[,1]), "sum")
  return(x)
})

IVIgMseq=lapply(1:16, function(i){
  x=rbind(IVIgMseq[[i*2-1]],IVIgMseq[[i*2]])
  x=aggregate(x[,2], by=list(x[,1]), "sum")
  return(x)
})

X=c(GBMseq,IVIgMseq)
comx=sapply(X, function(g){
  sapply(X, function(i){
    length(intersect(g[,1],i[,1]))
  })
})
comxi=array(0, dim=dim(comx))
for (i in 1:32) {
  for (j in 1:32) {
    comxi[i,j]=comx[i,j]/(comx[i,i]+comx[j,j]-comx[i,j])
  }
}

image(comxi, col = cpl(2000))

save(IVIgMseq,file="IVIgMseq")
save(GBMseq,file="GBMseq")
x=c()
for (i in 1:16){
  x=rbind(x,GBMseq[[i]])
}
GBMseq=aggregate(x[,2], by=list(x[,1]), "sum")
colnames(GBMseq)=c("Seq","N")
GBMseq=GBMseq[order(GBMseq[,2], decreasing=T),]

x=c()
for (i in 1:16){
  x=rbind(x,IVIgMseq[[i]])
}
IVIgMseq=aggregate(x[,2],by=list(x[,1]),FUN="sum")
colnames(IVIgMseq)=c("Seq","N")
IVIgMseq=IVIgMseq[order(IVIgMseq[,2], decreasing=T),]

length(GBMseq[GBMseq[,2]>2&GBMseq[,2]<11,2])
length(IVIgMseq[IVIgMseq[,2]>2&IVIgMseq[,2]<11,2])
tg=table(GBMseq[,2])
tiv=table(IVIgMseq[,2])
plot(log10(as.double(names(tg))), log10(tg), xlim=c(0,4), ylim = c(0,7))
par(new=T)
plot(log10(as.double(names(tiv))), log10(tiv), xlim=c(0,4), ylim = c(0,7), col=2)

# First experiment reselected
P1_QF=read.table("uniqueP1_QF.txt")
P1_QR=read.table("uniqueP1_QR.txt")
P2_QF=read.table("uniqueP2_QF.txt")
P2_QR=read.table("uniqueP2_QР.txt")
P3_QF=read.table("uniqueP3_QF.txt")
P3_QR=read.table("uniqueP3_QR.txt")
oldmims=rbind(P1_QF,P1_QR,P2_QR,P2_QF)
oldmims=aggregate(oldmims[,2],by=list(oldmims[,1]),FUN="sum")
colnames(oldmims)=c("Seq","N")
oldmims=oldmims[order(oldmims[,2], decreasing=T),]
to=table(oldmims[,2])

plot(log10(as.double(names(tg))), log10(tg), xlim=c(0,4), ylim = c(0,7))
par(new=T)
plot(log10(as.double(names(tiv))), log10(tiv), xlim=c(0,4), ylim = c(0,7), col=2)
par(new=T)
plot(log10(as.double(names(to))), log10(to), xlim=c(0,4), ylim = c(0,7), col=3)

allIgM=rbind(oldmims,IVIgMseq)
allIgM=aggregate(allIgM[,2],by=list(allIgM[,1]),FUN="sum")
colnames(allIgM)=c("Seq","N")

rm(P1_QF,P1_QR,P2_QF,P2_QR,P3_QF,P3_QR,qf,qr)

save(oldmims,file="oldmims")
save(IVIgMseq, file="IVIgMseqbeforeTUPclean")
save(GBMseq, file="GBMseqbeforeTUPclean")
write(oldmims$Seq, file = "oldmims.txt")
write(allIgM[,1],"allIgM.txt")
write(GBMseq[,1], "GBMseq.txt")
allIgMTUPs=read.table("exp1.txt",sep = '\t', header = T)
GBMseqTUPs=read.table("exp2.txt",sep = '\t', header = T)
oldmimsTUPs=read.table("oldmimsTUPs.txt",sep = '\t', header = T)
typbad=unique(c(allIgMTUPs$Brief.Description,GBMseqTUPs$Brief.Description,oldmimsTUPs$Brief.Description))
aIbad=allIgMTUPs[allIgMTUPs$Brief.Description %in% typbad[c(1,3,4,5,7)], 1]
gbmbad=GBMseqTUPs[GBMseqTUPs$Brief.Description %in% typbad[c(1,3,4,5,7)], 1]
oldbad=oldmimsTUPs[oldmimsTUPs$Brief.Description %in% typbad[c(1,3,4,5,7)], 1]

aIgM=IVIgMseq[!(IVIgMseq$Seq %in% aIbad),]
Gbm=GBMseq[!(GBMseq$Seq %in% gbmbad),]
oldm=oldmims[!(oldmims$Seq %in% oldbad),]
save(list = c("aIgM","Gbm","oldm"), file="aIgMGbmOldm")


allseq=rbind(aIgM,Gbm,oldm)
allseq=aggregate(allseq$N, by=list(allseq$Seq), "sum")
colnames(allseq)=c("Seq","N")
allseq=allseq[order(allseq$N, decreasing=T),]
tall=table(allseq[,2])
plot(log10(as.double(names(tall))), log10(tall), xlim=c(0,4), ylim = c(0,7))
seqToRetain=allseq$Seq[allseq$N>2]
length(seqToRetain)
aIgM=aIgM$Seq[aIgM$Seq %in% seqToRetain]
Gbm=Gbm$Seq[Gbm$Seq %in% seqToRetain]
oldm=oldm$Seq[oldm$Seq %in% seqToRetain]
rm(allseq, tall)

allIgM=union(oldm,aIgM)

glo_gbm=union(Gbm,allIgM)
fl=(glo_gbm %in% allIgM)*1+(glo_gbm %in% Gbm)*2
rm(allIgM,Gbm,X,GBMseq,IVIgMseq,oldmims,aIgM,oldm,seqToRetain,xx,xxx)
save(glo_gbm,file="glo_gbm_ordered")
save(fl, file = "fl_ordered")
i=sample(length(fl))
glo_gbm=glo_gbm[i]
fl=fl[i]
