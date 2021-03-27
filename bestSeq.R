# Given a list of lists of similar strings
# builds probabilistic suffix trees 
# for the separate lists. 
# Based on the PST model, the highest probability 
# sequence is selected from each list as
# a prototype representative.

bestSeq=function(L){
  require(PST)
  require(parallel)
  
  cl=makeCluster(4)
  clusterEvalQ(cl, library(PST))
  clusterExport(cl,list("L"))  
  mxpL=parLapply(cl,L,function(l){
    l=t(as.data.frame(strsplit(l,split=""),colnames=seq_along(l), stringsAsFactors=F))
    sts=seqdef(l)
    tree=pstree(sts,weighted=F)
    x=which.max(predict(tree,tree@data, decomp=F))
    return(x)
  })
  stopCluster(cl)
  
  calls=sapply(seq_along(L), function(i){
    cl=L[[i]]
    j=mxpL[[i]]
    return(cl[j])
  })

  return(unique(calls))
}