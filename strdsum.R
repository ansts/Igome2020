strdsum=function(p,m,f=fl,d="lcs"){
  require(stringdist)
  #m=m[m!=p]
  x=stringdist(p,m, method = d)
  y=table(f[x<5])
  return(list(y,m[x<3]))
}