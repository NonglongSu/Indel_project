#Simulate the Prob of typeN and typeS in phase1 and phase 2 via 64 * 4 codon space.
#Assume symmetry between insertion and deletion, hence I only simulate insertion here.

library(Biostrings)
library(stringr)
library(seqinr)

#setwd("~/Dropbox (ASU)/Indel_project/Script")

sim_ins = function(s1, s2, pos){
  
  c1  = str_split(s1, "")[[1]]
  c2  = str_split(s2, "")[[1]]
  
  if(pos == 1){
    sub1   = paste0(c1[1], c2[1], c2[2], collapse = "")
    sub2   = paste0(c2[3], c1[2], c1[3], collapse = "")
    sec    = syn[[s1]]
    if(sub1 %in% sec || sub2 %in% sec){ 
      ty = "s"  
    }else{
      ty = "n"   
    }
  }else{
    sub1   = paste0(c1[1], c1[2], c2[1], collapse = "")
    sub2   = paste0(c2[2], c2[3], c1[3], collapse = "")
    sec    = syn[[s1]]
    if(sub1 %in% sec || sub2 %in% sec){ 
      ty = "s"  
    }else{
      ty = "n"   
    }
  }
  return(ty)
}


#########################################################MAIN
# Assume nucleotide frequencies: A, C, G, T
nuc_freqs = c(0.308, 0.185, 0.199, 0.308)    


# construct codons and its degeneracy
codons = cbind(rep(DNA_BASES, each = 16),
               rep(DNA_BASES, times = 4, each = 4),
               rep(DNA_BASES, 16))

codonstrs  = apply(codons, 1, stringr::str_c, collapse = "")
syn        = syncodons(codonstrs)
names(syn) = toupper(names(syn))
syn        <<- lapply(syn, toupper)


df           = data.frame(matrix(ncol = 64, nrow = 4))
colnames(df) = codonstrs
rownames(df) = c("typeS1", "typeN1", "typeS2", "typeN2")

for (i in 1:64) {
  refC   = codonstrs[i]
  phase1 = c()
  phase2 = c()
  j      = 1
  repeat{
    insC     = sample(DNA_BASES, 3, prob = nuc_freqs, replace = TRUE)
    insC     = paste0(insC, collapse = "")
    res1     = sim_ins(refC, insC, 1)
    res2     = sim_ins(refC, insC, 2)
    phase1   = c(phase1, res1)
    phase2   = c(phase2, res2)
    j        = j + 1
    if(j == 1000){
      break
    }
  }
  
  typeS_1  = round(length(which(phase1 == "s")) / 1000, 2)
  typeS_2  = round(length(which(phase2 == "s")) / 1000, 2)
  df[1, i] = typeS_1
  df[2, i] = 1 - typeS_1
  df[3, i] = typeS_2
  df[4, i] = 1 - typeS_2
}
  
  
  

print.data.frame(df)

