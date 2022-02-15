# Nei Gojobori (1986) method

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



# inDir  = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_anc"
# ouFile = "../test_human_mouse_rat/Data_6.2/Results/dNdS.txt"

main = function(inDir, ouFile){
  
  # Create a codon table
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
  # Filter pre-stop codons
  preStop = c("TAA","TGA","TAG")
  
  Files = list.files(inDir, full.names = TRUE)
  
  n.site = 0 # number of non-synonymous sites 
  s.site = 0 # number of synonymous sites 
  Nd.m   = 0 # number of non-synonymous substitutions
  Sd.m   = 0 # number of synonymouse substitutions
  Nd.r   = 0
  Sd.r   = 0
  
  
  for (i in 1:length(Files)) {
    dna   = readDNAStringSet(Files[i] ,format = "fasta")
    len   = width(dna)[1]
    
    A = toString(dna[[1]]) 
    M = toString(dna[[2]]) 
    R = toString(dna[[3]]) 
    
    # Cal. the dN/dS of each one 
    j = 1                              #reset to 1 each run
    while (j < len) {
      pre.aa = substr(A, j, j + 2)
      if(pre.aa %in% preStop){
        j = j + 3
        next()
      }else{
        #Expected sites
        nSite    = cal_sites(A, j)
        n.site   = n.site + nSite[1]
        s.site   = s.site + nSite[2]
        #Mouse
        Nd_M     = cal_subs(A, M, j)
        Sd.m     = Sd.m + Nd_M[1]
        Nd.m     = Nd.m + Nd_M[2]
        #Rat
        Nd_R     = cal_subs(A, R, j)
        Sd.r     = Sd.r + Nd_R[1]
        Nd.r     = Nd.r + Nd_R[2]
        j = j + 3                     #use codon length as increment unit
      }
    }
  }
  
  # Comparison of observed vs expected. (KaKs = n/s vs N/S) 
  KaKs_A_exp = n.site / s.site
  KaKs_M_obs = Nd.m / Sd.m
  KaKs_R_obs = Nd.r / Sd.r

  # Cal. the W (dN/dS)
  P.n_m = Nd.m / n.site
  P.s_m = Sd.m / s.site
  
  P.n_r = Nd.r / n.site
  P.s_r = Sd.r / s.site

  # Jukes Cantor formula (1969) 
  dNdS_M = log(1 - 4 * P.n_m / 3, exp(1)) / log(1 - 4 * P.s_m / 3, exp(1))
  dNdS_R = log(1 - 4 * P.n_r / 3, exp(1)) / log(1 - 4 * P.s_r / 3, exp(1))
  
  
  # Generate a database
  df  = data.frame("nonsyn_subs"  = c(Nd.m, Nd.r),             "syn_subs"  = c(Sd.m, Sd.r),
                   "nonsyn_sites" = c(n.site, n.site),         "syn_sites" = c(s.site, s.site), 
                   "obs_ratio"    = c(KaKs_M_obs, KaKs_R_obs), "exp_ratio" = c(KaKs_A_exp, KaKs_A_exp), 
                   "omega"        = c(dNdS_M, dNdS_R))
  df1 = as.data.frame(lapply(df, function(x){sprintf("%.3f", x)}))
  rownames(df1) = c("mouse", "rat")
  write.table(df1, file = ouFile,
              sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE )
}


args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])
