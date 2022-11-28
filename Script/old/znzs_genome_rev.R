# Treat indels as an integrate part.
# Insertion in ancestor are considered as deletion in focal species

# insertion in A (ancestor) is not the same as deletion in B (derived) in expected number of sites (hard to fix)
# phase_0 indels will cause the expected ratio of n.site/s.site < 1 (which is the opposite of dnds scenario)

library(Biostrings)
library(stringr)
library(seqinr)

# setwd("~/Dropbox (ASU)/Indel_project/Script")

# Number of non-syn/syn sites
cal_sites = function(anc.seq, wid, len){
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

 
#Number of non-syn/syn indels
cal_subs = function(x, y, Gj){
  
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



#inDir  = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_anc"
#ouFile = "../test_human_mouse_rat/Data_6.2/Results/zNzS.txt"
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
  
  Files = list.files(inDir, full.names=TRUE)
  
  ssite.m = 0 #number of synonymous deletion sites in mouse.
  nsite.m = 0 #number of non-synonymous deletion sites in mouse.
  ssite.r = 0 #number of synonymous insertion sites in rat.
  nsite.r = 0 #number of non-synonymous insertion sites in rat. 
  
  S.del_m = 0 #number of synonymous deletions in mouse.
  N.del_m = 0 #number of non-synonymous deletions in mouse. 
  S.ins_m = 0 #number of synonymous insertions in mouse.
  N.ins_m = 0 #number of non-synonymous insertions in mouse.
  
  S.ins_r = 0 #number of synonymous insertions in rat.
  N.ins_r = 0 #number of non-synonymous insertions in rat. 
  S.del_r = 0 #number of synonymous deletions in rat. 
  N.del_r = 0 #number of non-synonymous deletions in rat. 
  
  for (i in 1:length(Files)) { 
    dna = readDNAStringSet(Files[i], format = "fasta")
    len = width(dna)[1]
    
    A = toString(dna[[1]])
    M = toString(dna[[2]])
    R = toString(dna[[3]])
  
    #Remove all gaps from ancestor sequence
    A1 = str_remove_all(A, '-')
    M1 = str_remove_all(M, '-')
    R1 = str_remove_all(R, '-')
    
    len.A1  = nchar(A1)
    len.M1  = nchar(M1)
    len.R1  = nchar(R1)
    
    dna.str = str_split(as.character(dna), '')
    g       = lapply(dna.str, function(x){IRanges(x == '-')})
    g       = IRangesList(g)
    
    m.width = width(g[[2]])
    r.width = width(g[[3]])
    
    #Count the number of expect indel sites
    print(paste0(basename(Files[i]), '-', i))
    for (k in c(3, 6, 9, 12)) {
      ssite.m = ssite.m + (cal_sites(A1, k, len.A1)[1] + cal_sites(M1, k, len.M1)[1])/8
      nsite.m = nsite.m + (cal_sites(A1, k, len.A1)[2] + cal_sites(M1, k, len.M1)[2])/8
      ssite.r = ssite.r + (cal_sites(A1, k, len.A1)[1] + cal_sites(R1, k, len.R1)[1])/8
      nsite.r = nsite.r + (cal_sites(A1, k, len.A1)[2] + cal_sites(R1, k, len.R1)[2])/8
    }
    
    
    if(length(g[[2]]) > 0){
      toMatch = lapply(g, function(x){g[[2]] %in% x})
      pat.1   = c(TRUE,  TRUE, FALSE)
      pat.2   = c(FALSE, TRUE, FALSE)
      
      for (j in 1:length(g[[2]])) {
        pat.obs = unlist(lapply(toMatch, `[[`, j))
        if(all(pat.obs == pat.1)){                        #insertion in rat
          Ni       = cal_subs(R, A, g[[2]][j])
          S.ins_r  = S.ins_r + Ni[1]
          N.ins_r  = N.ins_r + Ni[2]
        } 
        if(all(pat.obs == pat.2)){                        #deletion in mouse
            Nd       = cal_subs(A, M, g[[2]][j])
            S.del_m  = S.del_m + Nd[1]
            N.del_m  = N.del_m + Nd[2]
          }
        }
      }
    
    if(length(g[[3]]) > 0){
      toMatch = lapply(g, function(x){g[[3]] %in% x})
      pat.1 = c(TRUE,  FALSE, TRUE)
      pat.2 = c(FALSE, FALSE, TRUE)
     
      for (j in 1:length(g[[3]])) {
        pat.obs = unlist(lapply(toMatch, `[[`, j))
        if(all(pat.obs == pat.1)){                          #insertion in mouse
          Ni       = cal_subs(M, A, g[[3]][j])
          S.ins_m  = S.ins_m + Ni[1]
          N.ins_m  = N.ins_m + Ni[2]
        }
        if(all(pat.obs == pat.2)){                          #deletion in rat
          Nd       = cal_subs(A, R, g[[3]][j])
          S.del_r  = S.del_r + Nd[1]
          N.del_r  = N.del_r + Nd[2]
        }
      }
    }
  }
  
  #Comparison of observed vs expected. (KaKs = n/s vs N/S) 
  KaKs_M_obs = (N.ins_m + N.del_m) / (S.ins_m + S.del_m)
  KaKs_M_exp = (nsite.m) / (ssite.m)
  
  KaKs_R_obs = (N.ins_r + N.del_r) / (S.ins_r + S.del_r) 
  KaKs_R_exp = (nsite.r) / (ssite.r)
  
  # Cal. the zNzS without any model. 
  zNzS_M_indel = KaKs_M_obs / KaKs_M_exp
  zNzS_R_indel = KaKs_R_obs / KaKs_R_exp
   
  # Generate a database
  df  = data.frame("nonsyn_indels"  = c(N.ins_m + N.del_m, N.ins_r + N.del_r), "syn_indels" = c(S.ins_m + S.del_m, S.ins_r + S.del_r),
                    "nonsyn_sites"  = c(nsite.m, nsite.r),                     "syn_sites"  = c(ssite.m, ssite.r), 
                    "obs_ratio"     = c(KaKs_M_obs, KaKs_R_obs),               "exp_ratio"  = c(KaKs_M_exp, KaKs_R_exp), 
                    "omega"         = c(zNzS_M_indel, zNzS_R_indel))
  df1 = as.data.frame(lapply(df, function(x){sprintf("%.3f", x)}))
  rownames(df1) = c("mouse", "rat")
   
   
   write.table(df1, file = ouFile,
               sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names=TRUE )
}




args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])

  
  
  
 
  
  
  

  
