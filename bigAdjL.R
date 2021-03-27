require(stringdist)
require(parallel)
require(pbapply)

proct=proc.time()
cl=makeCluster(4, outputfile="")
clusterExport(cl,list("glo_gbm","stringdist", "fl"))
dummy=parSapply(cl,glo_gbm, function(p){
  x=stringdist(p,glo_gbm, method = "lcs")
  x=table(fl[x<5])
  write(paste(c(p,x), collapse =" "), file="bigAL.txt", append = T)
  return("")
})
stopCluster(cl)
print(proc.time()-proct)

vstrdsum=Vectorize(strdsum, "p")

N=length(glo_gbm) %/% 20
N=sapply(1:20, function(i){
  (i-1)*N+1:N
})
n=nrow(N)

########################
x=apply(N,2,function(i){
  j=list(i[1:(n %/% 3)],i[((n %/% 3)+1):(2*(n %/% 3))],i[(2*(n %/% 3)+1):n])
  cc=sum(lengths(j))
  cl=makeCluster(3)
  clusterExport(cl,list("glo_gbm","vstrdsum"))
  clusterExport(cl,list("fl","cc"), envir = environment())
  proct=proc.time()
  dummy=parSapply(cl,j, function(ii){
    l=glo_gbm[ii]
    x=t(vstrdsum(l,glo_gbm, f=fl))
    return(list(x))
  })
  stopCluster(cl)
  print(proc.time()-proct)
  dummy=rbind(dummy[[1]], dummy[[2]], dummy[[3]])
  L=dummy[,2]
  dummy=t(array(unlist(dummy[,1]), dim=c(3,cc), dimnames = list(1:3, names(dummy[,1]))))
  fnm=paste("bigAL",i[1],i[n],sep = "_")
  save(L, file=fnm) 
  fnm=paste("nn3",i[1],i[n],sep = "_")
  save(dummy, file=fnm)
  
  gc()
  return(i[n])
})
