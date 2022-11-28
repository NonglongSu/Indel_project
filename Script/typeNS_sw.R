# Calculate the typeN/typeS

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

#setwd("~/Dropbox (ASU)/Indel_project")

eff_phase = function(seq1,pos,wid,len){#seq1:ref
  seq1.char = str_split(seq1,"")[[1]]
  eff.sub = 0
  non.sub = 0
  for(j in 1: len){
    if(pos[j]%%3 == 1){#phase-0
      non.sub = non.sub+1
    }else if(pos[j]%%3 == 0){#phase-2
      pos.ori = pos[j]-2 
      unit.1  = substr(seq1, pos.ori, pos.ori+2)
      unit.2  = substr(seq1, pos.ori+wid[j], pos.ori+wid[j]+2)
      sub     = paste0(c(seq1.char[pos[j]-2], seq1.char[pos[j]-1], seq1.char[pos[j]+wid[j]]), collapse="")
      sec     = codon[[which(sapply(codon, function(X){sub %in% X}))]]
      if(unit.1%in%sec || unit.2%in%sec){
        non.sub = non.sub+1
      }else{
        eff.sub = eff.sub+1
      }
    }else{#phase-1
      pos.ori = pos-1 
      unit.1  = substr(seq1, pos.ori, pos.ori+2)
      unit.2  = substr(seq1, pos.ori + wid[j], pos.ori+wid[j]+2)
      sub     = paste0(c(seq1.char[pos[j]-1], seq1.char[pos[j]+wid[j]], seq1.char[pos[j]+wid[j]+1]), collapse="")
      sec     = codon[[which(sapply(codon, function(X){sub %in% X}))]]
      if(unit.1%in%sec || unit.2%in%sec){
        non.sub = non.sub+1
      }else{
        eff.sub = eff.sub+1
      }
    }
  }  
  res = c(eff.sub,non.sub) 
  return(res)
}

#inD = "test_90_species/Data/sw_cds/07_yeast_aligned_cds"
#ouF = "test_90_species/Results/07_yeast_aligned_cds.sw.txt"

main = function(inD,ouF){
  
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
  
  eff.1 = 0
  non.1 = 0
  eff.2 = 0
  non.2 = 0
  
  File.total = list.files(inD,full.names=TRUE)
  for(i in 1:length(File.total)){
    dna = readBStringSet(File.total[i])
    len = length(dna)
    M   = toString(dna[[1]])
    R   = toString(dna[[2]])
    
    dna.1 = str_split(as.character(dna),'')
    g     = lapply(dna.1, function(x) {IRanges(x=='-')})
    g     = IRangesList(g)
    g1    = g[[1]][which(width(g[[1]])<=12)] 
    g2    = g[[2]][which(width(g[[2]])<=12)] 
    
    wid.m = width(g1)
    wid.r = width(g2)
    pos.m = start(g1)
    pos.r = start(g2)
    l.m   = length(wid.m)
    l.r   = length(wid.r)
    
    if(l.m > 0){#mouse
      eff_M = eff_phase(R,pos.m,wid.m,l.m)    
      eff.1 = eff.1 + eff_M[1]
      non.1 = non.1 + eff_M[2]
    }
    if(l.r > 0){#rat
      eff_R = eff_phase(M,pos.r,wid.r,l.r)
      eff.2 = eff.2 + eff_R[1] 
      non.2 = non.2 + eff_R[2] 
    }
    print(i)
  }
  
  typeN  = eff.1 + eff.2
  typeS  = non.1 + non.2
  typeNS = typeN/typeS
  
  #output
  write(typeNS,file=ouF)
}

#######################################
args = commandArgs(trailingOnly=TRUE)
main(args[1],args[2])
