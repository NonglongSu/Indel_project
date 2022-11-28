#Focus on indels of three
#Have not write in makefile yet.


library(Biostrings)
library(stringr)
library(seqinr)

#setwd("~/Dropbox (ASU)/Indel_project/Script")

#Cal. the freq of ancester seq
countN = function(A1){
  anc = str_split(A1,'')[[1]] 
  nuc.freq = c()
  for (i in DNA_BASES) {
    nuc.freq =c(nuc.freq,length(which(anc==i)))
  }
  nfreq = nuc.freq/sum(nuc.freq)
  return(nfreq)
}

#Number of non-syn/syn insertion sites
cal_isite = function(A1,k,Pi){
  ins = sample(DNA_BASES,k,prob=Pi,replace=TRUE)  
  anc = str_split(A1,'')[[1]] 
  
  s.d=0
  n.d=0
  i  =1
  while (i<=nchar(A1)) {
    if(i%%3 == 1){#Phase-0
      unit = substr(A1,i,i+2)
      if(grepl('N', unit)){
        i = i+3
      }else{
        s.d = s.d+1
        i   = i+1
      }
    }else{
      if(i%%3 == 0){#phase-2
        ref    = paste0(anc[i-2],  anc[i-1], anc[i], collapse="")
        unit.1 = paste0(anc[i-2],  anc[i-1], ins[1], collapse="")
        unit.2 = paste0(ins[k-1],  ins[k],   anc[i], collapse="")
        sec    = codon[[which(sapply(codon, function(X){ref %in% X}))]]
        if(unit.1%in%sec || unit.2%in%sec){
          s.d = s.d+1
        }else{
          n.d = n.d+1  
        }
      }else{#phase-1
        ref    = paste0(anc[i-1], anc[i], anc[i+1], collapse="")
        unit.1 = paste0(anc[i-1], ins[1], ins[2],   collapse="")
        unit.2 = paste0(ins[k],   anc[i], anc[i+1], collapse="")
        sec    = codon[[which(sapply(codon, function(X){ref %in% X}))]]
        if(unit.1%in%sec || unit.2%in%sec){
          s.d = s.d+1
        }else{
          n.d = n.d+1
        }
      }
      i = i+1
    }
  }
  res = c(s.d,n.d) 
  return(res)
}

#Number of non-syn/syn deletion sites
cal_dsite = function(A1,k,len){
  D_mer = len-k+1
  anc   = str_split(A1,'')[[1]] 
  
  s.d = 0
  n.d = 0
  i   = 1
  while (i <= D_mer) {
    if(i%%3 == 1){#Phase-0
      unit.1 = substr(A1, i, i+2)
      unit.2 = substr(A1, i+k, i+k+2)
      if(grepl('N',unit.1) || grepl('N',unit.2)){
        i = i+3
      }else{
        s.d = s.d+1
        i   = i+1
      }
    }else{
      if(i%%3 == 0){#phase-2
        sub = paste0(c(anc[i-2], anc[i-1], anc[i+k]), collapse='')
        sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
        if(unit.1%in%sec || unit.2%in%sec){
          s.d = s.d+1
        }else{
          n.d = n.d+1  
        }
      }else{#phase-1
        sub = paste0(c(anc[i-1], anc[i+k], anc[i+k+1]), collapse='')
        sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
        if(unit.1%in%sec || unit.2%in%sec){
          s.d = s.d+1
        }else{
          n.d = n.d+1
        }
      }
      i = i+1
    }
  }
  res = c(s.d,n.d) 
  return(res)
}


#Number of non-syn/syn insertions
cal_ins = function(x,y,Gj){
  pos.sta = start(Gj)
  pos.end = end(Gj)
  
  insc    = str_split(substr(y,pos.sta,pos.end), '')[[1]]
  leni    = length(insc)
  
  xchar   = str_split(x,'')[[1]] 
  n.indel = s.indel = 0
  
  if(pos.sta%%3 == 0){#phase-2
    sub = paste0(c(xchar[pos.sta-2], xchar[pos.sta-1], xchar[pos.end+1]), collapse="")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = paste0(c(xchar[pos.sta-2], xchar[pos.sta-1], insc[1]), collapse="")
    unit.2 = paste0(c(insc[leni-1], insc[leni], xchar[pos.end+1]),  collapse="")
    if(unit.1%in%sec || unit.2%in%sec){
      s.indel = 1
    }else{
      n.indel = 1
    }
  }else if(pos.sta%%3 == 1){#phase-0
    s.indel = 1
  }else{#phase-1
    sub = paste0(c(xchar[pos.sta-1], xchar[pos.end+1], xchar[pos.end+2]), collapse="")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = paste0(c(xchar[pos.sta-1], insc[1], insc[2]))
    unit.2 = paste0(c(insc[leni], xchar[pos.end+1], xchar[pos.end+2]))
    if(unit.1%in%sec || unit.2%in%sec){
      s.indel = 1
    }else{
      n.indel = 1
    }
  }
  
  res = c(s.indel,n.indel)
  return(res)
}

#Number of non-syn/syn deletions
cal_del = function(x,Gj){
  pos.sta = start(Gj)
  pos.end = end(Gj)
  
  xchar   = str_split(x,'')[[1]] 
  n.indel = s.indel = 0
  
  if(pos.sta%%3 == 0){#phase-2
    sub = paste0(c(xchar[pos.sta-2], xchar[pos.sta-1], xchar[pos.end+1]), collapse="")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = substr(x, pos.sta-2, pos.sta)
    unit.2 = substr(x, pos.end-1, pos.end+1)
    if(unit.1%in%sec || unit.2%in%sec){
      s.indel = 1
    }else{
      n.indel = 1
    }
  }else if(pos.sta%%3 == 1){#phase-0
    s.indel = 1
  }else{#phase-1
    sub = paste0(c(xchar[pos.sta-1], xchar[pos.end+1], xchar[pos.end+2]), collapse="")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = substr(x, pos.sta-1, pos.sta+1)
    unit.2 = substr(x, pos.end, pos.end+2)
    if(unit.1%in%sec || unit.2%in%sec){
      s.indel = 1
    }else{
      n.indel = 1
    }
  }
  
  res = c(s.indel,n.indel)
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
  
  #reading files
  Files = list.files(inDir, full.names=TRUE)
  
  dsite = c(0,0) #number of syn/non deletion sites 
  isite = c(0,0) #number of syn/non insertion sites 
  
  Del   = c(0,0) #number of syn/non deletions.
  Ins   = c(0,0) #number of syn/non insertions.
  
  
  for (i in 1:length(Files)) { 
    dna = readDNAStringSet(Files[i])
    len = width(dna)[1]
    
    A = toString(dna[[1]])
    M = toString(dna[[2]])
    R = toString(dna[[3]])
    
    #Remove all gaps from ancestor sequence
    A1      = str_remove_all(A,'-')
    len.A1  = nchar(A1)
    
    dna.str = str_split(as.character(dna), '')
    g       = lapply(dna.str, function(x){IRanges(x=='-')})
    g       = IRangesList(g)
    
    m.wid = width(g[[2]])
    r.wid = width(g[[3]])
    
    
    cat(sprintf(":%d  %s\n",i,paste0(basename(Files[i]))))
    
    #Count the # of expected indel sites
    dsite = dsite + cal_dsite(A1,3,len.A1) 
    isite = isite + cal_isite(A1,3,countN(A1))
      
    #Count the number of observed indels
    
    if(length(g[[2]])>0){
      toMatch = lapply(g, function(x){g[[2]] %in% x})
      pat.1   = c(TRUE,TRUE,FALSE)
      pat.2   = c(FALSE,TRUE,FALSE)
      for (j in 1:length(g[[2]])) {
        if(width(g[[2]][j])!=3){
          next
        }else{
          pat.obs = unlist(lapply(toMatch, `[[`, j))
          if(all(pat.obs==pat.1)){#ins in R
            Ins = Ins + cal_ins(A, R, g[[2]][j])
          } 
          if(all(pat.obs == pat.2)){#del
            Del = Del + cal_del(A, g[[2]][j])
          }
        }
      }
    }
    
    if(length(g[[3]])>0){
      toMatch = lapply(g, function(x){g[[3]] %in% x})
      pat.1 = c(TRUE,FALSE,TRUE)
      pat.2 = c(FALSE,FALSE,TRUE)
      for (j in 1:length(g[[3]])) {
        if(width(g[[3]][j])!=3){
          next
        }else{
          pat.obs = unlist(lapply(toMatch, `[[`, j))
          if(all(pat.obs == pat.1)){#ins in M
            Ins = Ins + cal_ins(A, M, g[[3]][j])
          }
          if(all(pat.obs == pat.2)){#deletion 
            Del = Del + cal_del(A, g[[3]][j])
          }
        }
      }
    }
  }
  
  #Comparison of observed vs expected. (KaKs = n/s vs N/S) 
  KaKs_obsI = Ins[2]/Ins[1]
  KaKs_expI = isite[2]/isite[1]
  
  KaKs_obsD = Del[2]/Del[1]
  KaKs_expD = dsite[2]/dsite[1]
  
  #Cal. zNzS 
  zNzS_I = KaKs_obsI / KaKs_expI
  zNzS_D = KaKs_obsD / KaKs_expD
  
  zNzS.sum = ((Ins[2]+Del[2])/(Ins[1]+Del[1])) / ((isite[2]+dsite[2])/(isite[1]+dsite[1]))
  
  cat(sprintf("znzs_I:%.6f  znzs_D:%.6f\n ", zNzS_I, zNzS_D))
  
}




args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])
