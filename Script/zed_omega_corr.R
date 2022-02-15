# Cal the znzs value based on separate bin groups
# Make a correlation plot of dnds vs znzs

library(Biostrings)
library(stringr)
library(seqinr)
library(dplyr)
library(readr)

#setwd("~/Dropbox (ASU)/Indel_project/Script")

#####################################################PARTI

# Number of non-syn/syn indel sites
cal_zed_sites = function(anc.seq, wid, len){
  D_mer    = len - wid + 1
  anc.char = str_split(anc.seq, "")[[1]] 
  
  s.d = 0
  n.d = 0
  i   = 1
  while (i <= D_mer) {
    if(i %% 3 == 1){#Phase-0
      unit.1 = substr(anc.seq, i, i + 2)
      unit.2 = substr(anc.seq, i + wid, i + wid + 2)
      if(grepl('N', unit.1) || grepl('N', unit.2)){
        i = i + 3
      }else{
        s.d = s.d + 1
        i   = i + 1
      }
    }else{
      if(i %% 3 == 0){#phase-2
        sub = paste0(c(anc.char[i - 2], anc.char[i - 1], anc.char[i + wid]), collapse = "")
        sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
        if(unit.1 %in% sec || unit.2 %in% sec){
          s.d = s.d + 1
        }else{
          n.d = n.d + 1  
        }
      }else{#phase-1
        sub = paste0(c(anc.char[i - 1], anc.char[i + wid], anc.char[i + wid + 1]), collapse = "")
        sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
        if(unit.1 %in% sec || unit.2 %in% sec){
          s.d = s.d  + 1
        }else{
          n.d = n.d + 1
        }
      }
      i = i + 1
    }
  }
  res = c(s.d, n.d) 
  return(res)
}

# Number of observed mutational indel events
cal_zed_subs = function(x, y, Gj){
  
  pos.sta = start(Gj)
  pos.end = end(Gj)
  
  y.char  = str_split(y, "")[[1]]    # y has gaps
  
  n.indel = s.indel = 0
  
  if(pos.sta %% 3 == 0){#phase-2
    sub = paste0(c(y.char[pos.sta - 2], y.char[pos.sta - 1], y.char[pos.end + 1]), collapse = "")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = substr(x, pos.sta - 2, pos.sta)
    unit.2 = substr(x, pos.end - 1, pos.end + 1)
    if(unit.1 %in% sec || unit.2 %in% sec){
      s.indel = 1
    }else{
      n.indel = 1
    }
  }else if(pos.sta %% 3 == 1){#phase-0
    s.indel = 1
  }else{#phase-1
    sub = paste0(c(y.char[pos.sta - 1], y.char[pos.end + 1], y.char[pos.end + 2]), collapse = "")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = substr(x, pos.sta - 1, pos.sta + 1)
    unit.2 = substr(x, pos.end, pos.end + 2)
    if(unit.1 %in% sec || unit.2 %in% sec){
      s.indel = 1
    }else{
      n.indel = 1
    }
  }
  
  res = c(s.indel, n.indel)
  return(res)
}

#include '+/='
get_znzs_par = function(dna, par){
  
  sm = par[1] 
  nm = par[2] 
  sr = par[3] 
  nr = par[4] 
  
  Sm = par[5] 
  Nm = par[6] 
  Sr = par[7] 
  Nr = par[8] 
  
  a1 = 1 
  a2 = 2
  M        = toString(dna[a1]) 
  R        = toString(dna[a2]) 
  toRemove = paste0(c('-', '\\+', '='), collapse = "|")
  M1       = str_remove_all(M, toRemove)
  R1       = str_remove_all(R, toRemove)
  len.M1   = nchar(M1)
  len.R1   = nchar(R1)
  
  dna.str  = str_split(as.character(dna[a1:a2]), '')
  g        = IRangesList(lapply(dna.str, function(x){IRanges(x == '-')}))
  gap      = find_gap(dna.str, "-")
  
  nSite.M  = cal_zed_sites(M1, 3, len.M1)    
  nSite.R  = cal_zed_sites(R1, 3, len.R1)
  s.Site   = (nSite.M[1] + nSite.R[1]) / 2
  n.Site   = (nSite.M[2] + nSite.R[2]) / 2
  
  if(gap == 0){
    res = c(s.Site, n.Site, 0, 0)
  }else{
    if(length(g[[1]]) > 0){# indels in mouse
      for (j in 1:length(g[[1]])) { 
        subs = cal_zed_subs(R, M, g[[1]][j])
        Sm   = Sm + subs[1]
        Nm   = Nm + subs[2]
      } 
    }
    if(length(g[[2]]) > 0){# indels in rat
      for (j in 1:length(g[[2]])) {
        subs = cal_zed_subs(M, R, g[[2]][j])
        Sr   = Sr + subs[1]
        Nr   = Nr + subs[2]
      }
    }
    s.indel = Sm + Sr 
    n.indel = Nm + Nr
    res     = c(s.Site, n.Site, s.indel, n.indel)
  }
  
  return(res)
}

#Find 3-mer gaps of -/=/+
find_gap = function(seq.c, sign){
  g        = lapply(seq.c, function(x){IRanges(x == sign)})
  g        = IRangesList(g)
  gaps     = length(which(unlist(width(g)) == 3))
  return(gaps)
}

######################################################PARTII

# Determine number of non-syn/syn substitution sites
cal_omega_sites = function(x, tag){
  
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

# Determine number of non-syn/syn substitution events
cal_omega_subs = function(x, y, tag){
  
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
    mut.pos = mapply(function(X, Y) which(X != Y), strsplit(unit.x, ""), strsplit(unit.y, ""))
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

# Get the dnds value for every single gene
get_dnds_par = function(dna, par){
  
  sm = par[1] 
  nm = par[2] 
  sr = par[3] 
  nr = par[4] 
  
  s.subs = par[5] 
  n.subs = par[6] 
 
  len = width(dna)[1]
  M   = toString(dna[1]) 
  R   = toString(dna[2]) 
  
  j = 1
  while (j < len) {
    Subs   = cal_omega_subs(M, R, j)
    s.subs = s.subs + Subs[1]
    n.subs = n.subs + Subs[2]
    
    #Mouse as ref
    nSite_M  = cal_omega_sites(M, j)
    nm       = nm + nSite_M[1]
    sm       = sm + nSite_M[2]
    
    #Rat as ref
    nSite_R  = cal_omega_sites(R, j) 
    nr       = nr + nSite_R[1]
    sr       = sr + nSite_R[2]
    
    j = j + 3    #use codon length as increment unit
  }
  s.Site = (sm + sr) / 2
  n.Site = (nm + nr) / 2
  res    = c(s.Site, n.Site, s.subs, n.subs)
  return(res)
}




##################################################################
inDir  = "../test_human_mouse_rat/Data_6.2/Mafft/mixed_cds/"
file1  = "../test_human_mouse_rat/Data_6.2/Results/omega_ladder.txt"
file2  = "../test_human_mouse_rat/Data_6.2/Results/omega_bin.txt"
ouFile = "../test_human_mouse_rat/Data_6.2/Results/omega_zed.txt"
ouFig  = "../test_human_mouse_rat/Data_6.2/Figure/omega_zed.pdf"

main = function(inDir, file1, file2, ouFile, ouFig){
  
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
  
  DF.rank = read_delim(file1, "\t", col_names  = TRUE)
  geneId  = DF.rank[[1]]
  dnds    = DF.rank[[2]]
  index   = read_delim(file2, "\t", col_names  = FALSE)[[1]]
  
  pars.1   = c(0, 0, 0, 0, 0, 0, 0, 0)
  pars.2   = c(0, 0, 0, 0, 0, 0)
  znzs.set = c()
  dnds.set = c()
  k        = 1
  
  for (i in index) {
    S.site.1 = N.site.1 = S.indels = N.indels = 0
    S.site.2 = N.site.2 = S.subs = N.subs = 0
    
    while (k <= i) {
      dna       = readBStringSet(paste0(inDir, geneId[k], ".fa"))
      
      zed.value = get_znzs_par(dna, pars.1)
      S.site.1  = S.site.1 + zed.value[1]
      N.site.1  = N.site.1 + zed.value[2]
      S.indels  = S.indels + zed.value[3]
      N.indels  = N.indels + zed.value[4]
      
      omega.value = get_dnds_par(dna, pars.2)
      S.site.2  = S.site.2 + omega.value[1]
      N.site.2  = N.site.2 + omega.value[2]
      S.subs    = S.subs + omega.value[3]
      N.subs    = N.subs + omega.value[4]
      k         = k + 1
    }
    znzs     = (N.indels / N.site.1) / (S.indels / S.site.1)
    dnds     = (N.subs / N.site.2) / (S.subs / S.site.2)
    znzs.set = c(znzs.set, znzs)
    dnds.set = c(dnds.set, dnds)
  }
  
  W.zed = data.frame("znzs" = znzs.set, "dnds" = dnds.set)
  write.table(W.zed, file = ouFile,
              sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE )
  # #Plot
  pdf(ouFig)
  plot(W.zed$dnds, W.zed$znzs,
       main = "The correlation between znzs and dnds ratio across the genome",
       xlab = "dnds ratio", ylab = "znzs ratio", pch = 20, col = "#88CCEE")
  dev.off()
  
}



args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4], args[5])
