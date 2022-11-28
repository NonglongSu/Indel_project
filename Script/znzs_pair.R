#Treat indels as an integrate part.
#Insertion in ancestor are considered as deletion in focal species

#insertion in A (ancestor) is not the same as deletion in B (derived) in expected number of sites (hard to fix)
#phase_0 indels will cause the expected ratio of n.site/s.site < 1 (which is the opposite of dnds scenario)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")

#Number of non-syn/syn sites
cal_sites = function(anc.seq,wid,len){
  D_mer    = len-wid+1
  anc.char = str_split(anc.seq, "")[[1]] 
  
  s.d = 0
  n.d = 0
  i   = 1
  while (i <= D_mer) {
    if(i %% 3 == 1){#Phase-0
      unit.1 = substr(anc.seq, i, i+2)
      unit.2 = substr(anc.seq, i+wid, i+wid+2)
      s.d    = s.d + 1
    }else if(i %% 3 == 0){#phase-2
        sub = paste0(c(anc.char[i-2], anc.char[i-1], anc.char[i+wid]), collapse = "")
        sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
        if(unit.1%in%sec || unit.2%in%sec){
          s.d = s.d + 1
        }else{
          n.d = n.d + 1  
        }
      }else{#phase-1
        sub = paste0(c(anc.char[i-1], anc.char[i+wid], anc.char[i+wid+1]), collapse = "")
        sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
        if(unit.1%in%sec || unit.2%in%sec){
          s.d = s.d + 1
        }else{
          n.d = n.d + 1
        }
      }
      i = i + 1
    }
  
  res = c(s.d,n.d) 
  return(res)
}


#Number of non-syn/syn indels
cal_subs = function(x,Gj){
  pos.sta = start(Gj)
  pos.end = end(Gj)
  
  x.char   = str_split(x,'')[[1]]
  #y.char  = str_split(y,"")[[1]]    # y has gaps
  n.indel = s.indel = 0
  
  if(pos.sta %%3 == 0){#phase-2
    sub = paste0(c(x.char[pos.sta-2], x.char[pos.sta-1], x.char[pos.end+1]),collapse="")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = substr(x, pos.sta-2, pos.sta)
    unit.2 = substr(x, pos.end-1, pos.end + 1)
    if(unit.1 %in% sec || unit.2 %in% sec){
      s.indel=1
    }else{
      n.indel=1
    }
  }else if(pos.sta %%3 == 1){#phase-0
    s.indel = 1
  }else{#phase-1
    sub = paste0(c(x.char[pos.sta-1], x.char[pos.end+1], x.char[pos.end+2]),collapse="")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = substr(x, pos.sta-1, pos.sta+1)
    unit.2 = substr(x, pos.end, pos.end+2)
    if(unit.1 %in% sec || unit.2 %in% sec){
      s.indel=1
    }else{
      n.indel=1
    }
  }
  
  res = c(s.indel,n.indel)
  return(res)
}



#inD = "Data_6/Mafft/mapped_cds"
#ouF = "Data_6/Results/znzs.txt"
#k   = '3'
main = function(inD, ouF, k){
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
  
  Files = list.files(inD, full.names=TRUE)
  k     = as.numeric(k)
  
  ssite = 0 #number of synonymous deletion sites 
  nsite = 0 #number of non-synonymous deletion sites 
  S.id  = 0 #number of synonymous deletions.
  N.id  = 0 #number of non-synonymous deletions. 
 
  #pattern
  pat = c('-', '\\+', '=')
  pat = paste0(pat,collapse='|')
  
  for (i in 1:length(Files)) { 
    dna = readBStringSet(Files[i])
    M   = toString(dna[[1]])
    R   = toString(dna[[2]])
    
    #Remove all weird gaps 
    M1  = str_remove_all(M,pat)
    R1  = str_remove_all(R,pat)
    
    len.M1  = nchar(M1)
    len.R1  = nchar(R1)
    
    dna.str = str_split(dna, '')
    g       = lapply(dna.str, function(x){IRanges(x == '-')})
    g       = IRangesList(g)
    
    #number of expected indel sites
    site1 = cal_sites(M1,k,len.M1)
    site2 = cal_sites(R1,k,len.R1)
    ssite = ssite + (site1[1]+site2[1])/2
    nsite = nsite + (site1[2]+site2[2])/2
    
    if((sum(site1)+2 != len.M1) || (sum(site2)+2 != len.R1)){
      print(paste0(basename(Files[i]),':',i))
      break
    } 
    
    #number of observed indel events
    if(length(g[[1]])>0){#mouse
      for (j in 1:length(g[[1]])) {
          Subs  = cal_subs(R,g[[1]][j])
          S.id  = S.id + Subs[1]
          N.id  = N.id + Subs[2]
      }
    }
    
    if(length(g[[2]])>0){#rat
      for (j in 1:length(g[[2]])) {
          Subs  = cal_subs(M,g[[2]][j])
          S.id  = S.id + Subs[1]
          N.id  = N.id + Subs[2]
      }
    }
    print(paste0(basename(Files[i]),':',i))
  }
  
  #observed vs expected. (KaKs = n/s vs N/S) 
  KaKs_obs = N.id/S.id
  KaKs_exp = nsite/ssite
  cat(sprintf("obs.ratio:%.3f exp.ratio:%.3f\n",KaKs_obs,KaKs_exp))
  
  #cal.the zNzS
  zNzS = KaKs_obs/KaKs_exp
  write(zNzS,file=ouF)
  
}




args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2], args[3])










