# Cal the znzs value based on separate bin groups
# Make a correlation plot of dnds vs znzs

library(Biostrings)
library(stringr)
library(seqinr)
library(dplyr)
library(readr)

#setwd("~/Dropbox (ASU)/Indel_project/Script")

#####################################################
# Number of non-syn/syn indel sites
cal_zed_sites = function(seq, wid, len){
  D_mer = len - wid + 1
  seq.c = str_split(seq, "")[[1]] 
  
  s.d = 0
  n.d = 0
  i   = 1
  while (i <= D_mer) {
    if(i %% 3 == 1){#Phase-0
      unit.1 = substr(seq, i, i + 2)
      unit.2 = substr(seq, i + wid, i + wid + 2)
      if(grepl('N', unit.1) || grepl('N', unit.2)){
        i = i + 3
      }else{
        s.d = s.d + 1
        i   = i + 1
      }
    }else{
      if(i %% 3 == 0){#phase-2
        sub = paste0(c(seq.c[i - 2], seq.c[i - 1], seq.c[i + wid]), collapse = "")
        sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
        if(unit.1 %in% sec || unit.2 %in% sec){
          s.d = s.d + 1
        }else{
          n.d = n.d + 1  
        }
      }else{#phase-1
        sub = paste0(c(seq.c[i - 1], seq.c[i + wid], seq.c[i + wid + 1]), collapse = "")
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
    sub    = paste0(c(y.char[pos.sta - 2], y.char[pos.sta - 1], y.char[pos.end + 1]), collapse = "")
    sec    = codon[[which(sapply(codon, function(X){sub %in% X}))]]
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
    sub    = paste0(c(y.char[pos.sta - 1], y.char[pos.end + 1], y.char[pos.end + 2]), collapse = "")
    sec    = codon[[which(sapply(codon, function(X){sub %in% X}))]]
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

# Get the znzs value for every single gene
get_znzs_par = function(dna, par){
  
  sm = par[1] 
  nm = par[2] 
  sr = par[3] 
  nr = par[4] 
  
  Sm = par[5] 
  Nm = par[6] 
  Sr = par[7] 
  Nr = par[8] 
  
  a1 = 3 
  a2 = 4
  M        = toString(dna[a1]) 
  R        = toString(dna[a2]) 
  toRemove = paste0(c('-', '\\+', '='), collapse = "|")
  M1       = str_remove_all(M, toRemove)
  R1       = str_remove_all(R, toRemove)
  len.M1   = nchar(M1)
  len.R1   = nchar(R1)
  
  dna.str  = str_split(as.character(dna[a1:a2]), '')
  g        = lapply(dna.str, function(x){IRanges(x == '-')})
  g        = IRangesList(g)
  gaps     = length(which(unlist(width(g)) == 3))
  
  nSite.M = cal_zed_sites(M1, 3, len.M1)    
  nSite.R = cal_zed_sites(R1, 3, len.R1)
  s.Site  = (nSite.M[1] + nSite.R[1]) / 2
  n.Site  = (nSite.M[2] + nSite.R[2]) / 2
  
  if(gaps == 0){
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
    
    res = c(s.Site, n.Site, s.indel, n.indel)
  }
  return(res)
}


##################################################################
# inDir  = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_cds/"
# ouFile = "../test_human_mouse_rat/Data_6.2/Results/znzs_gene.txt"

main = function(inDir, ouFile){
  
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
  
 
  pars     = c(0, 0, 0, 0, 0, 0, 0, 0)
  S.site   = N.site = S.indels = N.indels = c()
  
  Files = list.files(inDir, full.names = TRUE)
  
  for (i in 1:length(Files)) {
      dna       = readBStringSet(Files[i])
      
      zed.value = get_znzs_par(dna, pars)
      S.site    = c(S.site, zed.value[1])
      N.site    = c(N.site, zed.value[2])
      S.indels  = c(S.indels, zed.value[3])
      N.indels  = c(N.indels, zed.value[4])
  }
  
  Df = data.frame("N-site"  = N.site,   "S-site"  = S.site,
                  "N-indel" = N.indels, "S-indel" = S.indels )
  
  N.ind = sum(N.indels) 
  S.ind = sum(S.indels) 
  N.s   = sum(N.site)
  S.s   = sum(S.site)
  
  # Comparison of observed vs expected. (KaKs = n/s vs N/S) 
  zed_exp = N.s / S.s
  zed_obs = N.ind / S.ind
  znzs    = zed_obs / zed_exp
  
  write.table(Df, file = "../test_human_mouse_rat/Data_6.2/Results/Tmp/zed_gene0.txt",
              sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE )
 
  
}



args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])
