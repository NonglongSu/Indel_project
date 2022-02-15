base_freq = function(dat,cod){
  #codon and nuc freq
  f1     = rowSums(dat) + colSums(dat)  
  f1.id  = sapply(seq(64), function(x){match(cod[x, ], DNA_BASES)})
  base.f = c()
  for (i in seq(4)) {
    base.id = sapply(seq(64), function(x){length(which(f1.id[, x] == i))})
    base.f  = c(base.f, sum(base.id*f1))
  }
  f  = base.f/sum(base.f)
  if(sum(f)!=1){print("Warning:sum of nuc freq is not 1!")}
  
  cf = sapply(seq(64), function(x){prod(f[match(cod[x, ], DNA_BASES)])})
  if(round(sum(cf),3)!=1){print("Warning:sum of codon freq is not 1!")}
  
  res=list(f,cf)
  return(res)
}