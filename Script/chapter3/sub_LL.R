
#4*4 GTR matrix
GTR = function(Si, pai){
  r1 = matrix(0, 4, 4)
  r1[lower.tri(r1)] = Si
  r  = r1 + t(r1)
  r  = t(r*pai)
  diag(r) = -rowSums(r)
  return(r)
}

#Construct a 64*64 MG94 matrix
MG94 = function(g, oga, cod, syn, codonstrs){
  R  = matrix(0, 64, 64)
  for (i in 1:64) {
    for (j in 1:64) {
      if(i == j){
        R[i, j] = 0
      }else if(sum(cod[i, ] != cod[j, ]) > 1){#more than 1 nucleotide changes
        R[i, j] = 0
      }else{
        if(codonstrs[j] %in% syn[[codonstrs[i]]]){#syn_subs
          w = 1
        }else{#nonsyn_subs
          w = oga
        }
        pos = which(cod[i, ] != cod[j, ])
        x   = which(DNA_BASES == cod[i, pos])
        y   = which(DNA_BASES == cod[j, pos])
        R[i, j] = w * g[x, y]
      }
    }
  }
  diag(R) = -rowSums(R)
  
  return(R)
}


sub_weight = function(dat, codons, syn, codonstrs, p0){
  #codon and nuc freq
  f1     = rowSums(dat) + colSums(dat)  
  f1.id  = sapply(seq(64), function(x){match(codons[x, ], DNA_BASES)})
  base.f = c()
  for (i in seq(4)) {
    base.id = sapply(seq(64), function(x){length(which(f1.id[, x] == i))})
    base.f  = c(base.f, sum(base.id*f1))
  }
  f  = base.f/sum(base.f)
  if(sum(f)!=1){print("Warning:sum of nuc freq is not 1!")}
  
  cf = sapply(seq(64), function(x){prod(f[match(codons[x, ], DNA_BASES)])})
  if(round(sum(cf),3)!=1){print("Warning:sum of codon freq is not 1!")}
  
  #gtr + mg94 
  r   = GTR(p0[1:6],f)
  R   = MG94(r,p0[7],codons,syn,codonstrs)
  ll  = sum(log(expm(R)*cf)*dat)
  
  return(ll)
}