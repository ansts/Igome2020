# Selects from a list of strings
# a subset of strings no two of which 
# are closer to each other than d using 
# the LCS metric.

stayAway=function(l, d=3, startat=1){
  require(stringdist)
  xres=l[1:(startat-1)]
  l=l[startat:length(l)]
  for (p in xres){
    x1=stringdist(p,l, method = "lcs")
    print(p)
    l=l[x1>d]
    print(length(l))
  }
  x=l[1]
  l=l[-1]
  repeat{
    x1=stringdist(x,l, method = "lcs")
    l=l[x1>d]
    x=l[1]
    l=l[-1]
    xres=c(xres,x)
    print(length(l))
    if (length(l)==0) break
  }
  return(xres)
}