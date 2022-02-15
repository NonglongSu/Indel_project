
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

sub_cost = function(dat, cod, syn, codonstrs, omega){
  #nucletide freq (Yang (1994) Estimating the pattern of nucleotide substitution)
  nuc_freqs = c(0.308, 0.185, 0.199, 0.308)
  nuc_q     = c(-0.818, 0.132,  0.586,  0.1,            
                0.221, -1.349,  0.231,  0.897,
                0.909,  0.215, -1.322,  0.198,
                0.1,    0.537,  0.128, -0.765)
  nuc_q     = matrix(nuc_q,4,4, byrow = T)
  cod_freqs = sapply(seq(64), function(x){prod(nuc_freqs[match(cod[x, ], DNA_BASES)])})
  if(sum(cod_freqs)!=1){print("Warning:sum of codon freq is not 1!")}
  #mg94 
  R   = MG94(nuc_q,omega,cod,syn,codonstrs)
  ll  = sum(log(expm(R)*cod_freqs)*dat)
  return(ll)
}