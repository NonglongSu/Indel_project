library(Biostrings)
library(stringr)
library(seqinr)

#setwd("~/Dropbox (ASU)/Indel_project/Script")

# Determine number of non-syn/syn sites
cal_sites = function(x, tag){
  
  unit.x  = substr(x, tag, tag + 2)
  toMatch = c('-', '\\+', '=', 'N')
  s.site  = 0
  
  if(grepl(paste0(toMatch, collapse = "|"), unit.x)){
    return(c(0, 0))
  }else{
    section = codon[[which(sapply(codon, function(X){unit.x %in% X}))]]
    for (k in 1:3) {
      Base     = DNA_BASES[substr(unit.x, k, k) != DNA_BASES]
      unit.sub = unit.x
      for (w in 1:length(Base)) {
        substr(unit.sub, k, k) = Base[w]
        if(unit.sub %in% section){
          s.site = s.site + 1/3
        }
      }
    }
  }
  n.site = 3 - s.site
  res = c(n.site, s.site)
  return(res)
}

# Determine number of non-syn/syn substitutions
cal_subs = function(x, y, tag){
  
  unit.x  = substr(x, tag, tag + 2)
  unit.y  = substr(y, tag, tag + 2)
  sub.s   = 0   
  sub.n   = 0  
  toMatch = c('-', '\\+', '=', 'N')
  
  if(grepl(paste0(toMatch, collapse = "|"), unit.x) || grepl(paste0(toMatch, collapse = "|"), unit.y)){
    return(c(sub.s, sub.n))
  }else{
    section = codon[[which(sapply(codon, function(X){unit.x %in% X}))]]
    mut.num = mapply(function(X, Y) sum(X != Y), strsplit(unit.x, ""), strsplit(unit.y, ""))
    mut.pos = mapply(function(X, Y) which(X!=Y), strsplit(unit.x, ""), strsplit(unit.y, ""))
    if(mut.num == 1){#ATG/TTG
      if(unit.y %in% section){
        sub.s = sub.s + 1
      }else{
        sub.n = sub.n + 1
      }
    }else if(mut.num == 2){#ATG/TCG
      unit.list = collect_allSubs(unit.x, unit.y, mut.pos)
      sub.rate  = cal_subsub(unit.list, sub.s, sub.n, section)
      sub.s = sub.rate[1] / 2
      sub.n = sub.rate[2] / 2
    }else if(mut.num == 3){#ATG/CCT
      unit.list = collect_allSubs(unit.x, unit.y, mut.pos)
      sub.rate  = cal_subsub(unit.list, sub.s, sub.n, section)
      sub.s = sub.rate[1] / 6
      sub.n = sub.rate[2] / 6
    }else{#no substitution ATG/ATG
      return(c(sub.s, sub.n))
    }
  }
  return(c(sub.s, sub.n))
}

# Find all possible substitution pathways.        
# unit.x = "TTT"
# unit.y = "TAC"
# unit.y = "GAC"
# mut.pos = 2:3
# mut.pos = 1:3
collect_allSubs = function(unit.x, unit.y, mut.pos){
  unit.set  = c()
  for (i in 1:length(mut.pos)) {
    unit.sub = unit.x
    substr(unit.sub, mut.pos[i], mut.pos[i]) = substr(unit.y, mut.pos[i], mut.pos[i])
    unit.set = c(unit.set, unit.sub)
    if (length(mut.pos) == 3) {
      for (j in mut.pos[mut.pos != mut.pos[i]]) {
        unit.sb               = unit.sub
        substr(unit.sb, j, j) = substr(unit.y, j, j)
        unit.set              = c(unit.set, unit.sb)
      }
    }
  }
  unit.list = list()
  if(length(unit.set) == 2){
    for (k in 1:length(unit.set)) {
      unit.list = c(unit.list, list(c(unit.x, unit.set[k], unit.y)))
    }
  }else{
    for (k in 1:length(unit.set)) {
      if (k %% 3 == 1) {
        for (w in (k + 1):(k + 2)) {
          unit.list = c(unit.list, list(c(unit.x, unit.set[k], unit.set[w], unit.y)))
        }
      }
    }
  }
  return(unit.list)
}

# Sum the number of substitutions of all possible pathways
# section = c("TTT","TTC")
# 12 , 6
cal_subsub = function(unit.list, sub.s, sub.n, sec){
  for (m in 1:length(unit.list)) {
    sec.var = sec
    for (n in 2:(length(unit.list[[m]])) ) {
      if(unit.list[[m]][n] %in% sec.var){
        sub.s = sub.s + 1
      }else{
        sub.n = sub.n + 1
      }
      sec.var = codon[[which(sapply(codon, function(X){unit.list[[m]][n] %in% X}))]]
    }
  }
  return(c(sub.s, sub.n))
}



# inDir  = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_cds"
# ouFile = "../test_human_mouse_rat/Data_6.2/Results/dNdS_mul.txt"
# name   = c("mouse", "rat")

main = function(inDir, ouFile, name){
  
  #Create a codon table
  codon <<- list (c("TTT","TTC"),
                  c("TTA","TTG","CTT","CTC","CTA","CTG"),
                  c("ATT","ATC","ATA"),
                  c("ATG"),
                  c("GTT","GTC","GTA","GTG"),
                  c("TCT","TCC","TCA","TCG","AGT","AGC"),
                  c("CCT","CCC","CCA","CCG"),
                  c("ACT","ACC","ACA","ACG"),
                  c("GCT","GCC","GCA","GCG"),
                  c("TAT","TAC"),
                  c("CAT","CAC"),
                  c("CAA","CAG"),
                  c("AAT","AAC"),
                  c("AAA","AAG"),
                  c("GAT","GAC"),
                  c("GAA","GAG"),
                  c("TGT","TGC"),
                  c("TGG"),
                  c("CGT","CGC","CGA","CGG","AGA","AGG"),
                  c("GGT","GGC","GGA","GGG"),
                  c("TAA","TGA","TAG")
  )
  
  Files = list.files(inDir, full.names = TRUE)
  
  n.m  = 0 # number of non-synonymous sites (relative to species)
  s.m  = 0 # number of synonymous sites (relative to species)
  Nd.m = 0 # number of non-synonymous mutations
  Sd.m = 0 # number of synonymouse mutations
  
  n.r  = 0
  s.r  = 0
  Nd.r = 0
  Sd.r = 0
  
  for (i in 1:length(Files)) {
    dna   = readBStringSet(Files[i] ,format = "fasta")
    Len   = length(dna)
    len   = width(dna)[Len]
    
    M = toString(dna[Len - 1]) 
    R = toString(dna[Len]) 
    
    # Cal. the dN/dS of each one 
    j = 1                            #reset to 1 each run
    while (j < len) {
      #Mouse as ref
      nSite_M  = cal_sites(M, j)
      Nd_M     = cal_subs(M, R, j)
      n.m      = n.m + nSite_M[1]
      s.m      = s.m + nSite_M[2]
      Sd.m     = Sd.m + Nd_M[1]
      Nd.m     = Nd.m + Nd_M[2]
      #Rat as ref
      nSite_R  = cal_sites(R, j) 
      Nd_R     = cal_subs(R, M, j)
      n.r      = n.r + nSite_R[1]
      s.r      = s.r + nSite_R[2]
      Sd.r     = Sd.r + Nd_R[1]
      Nd.r     = Nd.r + Nd_R[2]
      
      j = j + 3                     #use codon length as increment unit
    }
  }
  
  # Comparison of observed vs expected. (KaKs = n/s vs N/S) 
  # KaKs_M_obs = Nd.m / Sd.m 
  # KaKs_M_exp = n.m / s.m
  # 
  # KaKs_R_obs = Nd.r / Sd.r 
  # KaKs_R_exp = n.r / s.r
  
  # Cal. the W (dN/dS)
  P.n_m = Nd.m / n.m
  P.s_m = Sd.m / s.m
  
  P.n_r = Nd.r / n.r
  P.s_r = Sd.r / s.r
  
  # Jukes Cantor formula (1969)
  dNdS_M2R = log(1 - 4 * P.n_m / 3, exp(1)) / log(1 - 4 * P.s_m / 3, exp(1))
  dNdS_R2M = log(1 - 4 * P.n_r / 3, exp(1)) / log(1 - 4 * P.s_r / 3, exp(1))
  
  # Generate a database
  df = data.frame(round(dNdS_M2R, 3), round(dNdS_R2M, 3))
  colnames(df) = name
  write.table(df, file = ouFile,
              sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE )
}





args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], eval(parse(text = args[3])))