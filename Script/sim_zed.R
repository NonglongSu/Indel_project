# The parameters from this simulation comes from the paper
#"An Empirical Codon Model for Protein Sequence Evolution (2007)"

library(Biostrings)
library(stringr)
library(seqinr)
library(dplyr)
library(readr)
#setwd("~/Dropbox (ASU)/Indel_project/Script")

sim_indel = function(seq){
  dna   = paste0(seq, collapse = "")
  dnas  = str_split(dna, "")[[1]]
  pos = sample(1:3, 1, replace = TRUE, prob = c(1/3, 1/3, 1/3)) 
  
  if(pos == 1){
    phase = 0 
     ty   = "s"  
  }else if(pos == 2){
    phase = 1
    sub   = paste0(dnas[1], dnas[5], dnas[6], collapse = "") 
    sec   = codon[[which(sapply(codon, function(X){sub %in% X}))]] 
    if(seq[[1]] %in% sec || seq[[2]] %in% sec){ 
      ty = "s"  
    }else{
      ty = "n"   
    }
  }else{
    phase = 2
    sub   = paste0(dnas[1], dnas[2], dnas[6], collapse = "") 
    sec   = codon[[which(sapply(codon, function(X){sub %in% X}))]] 
    if(seq[[1]] %in% sec || seq[[2]] %in% sec){ 
      ty = "s"  
    }else{
      ty = "n"   
    }
  }
  
  res = c(phase, ty)
  return(res)
}


#file = "../test_human_mouse_rat/Data_6.2/Simulation/codon_freq.txt"
main = function(file){
  
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
  
  # Input
  codon.df = read_delim(file, "\t", col_names  = FALSE)
  Codon    = codon.df[[1]]
  freq     = as.numeric(codon.df[[2]])
  
  Phase = c()
  Type  = c()
  i     = 0
  repeat{
    six_mer  = sample(Codon, 2, prob = freq, replace = TRUE)
    count    = sim_indel(six_mer)
    Phase    = c(Phase, count[1])
    Type     = c(Type, count[2])
    i        = i + 1
    if(i == 1000000){
      break
    }
  }
  
  DF   = data.frame("Phase" = Phase, "Type" = Type)
  P0.s  = length(which(DF$Phase == 0))
  
  P.s = DF %>% group_by(Phase) %>% filter(Type == "s")
  P.n = DF %>% group_by(Phase) %>% filter(Type == "n")
  
  P1.s = length(which(P.s$Phase == 1))
  P1.n = length(which(P.n$Phase == 1))
  P2.s = length(which(P.s$Phase == 2))
  P2.n = length(which(P.n$Phase == 2))
  
  Df = data.frame("phase" = 0:2, "S" = c(P0.s, P1.s, P2.s), "N" = c(0, P1.n, P2.n)) 
  
  print.data.frame(Df)
  
}

args = commandArgs(trailingOnly = TRUE)
main(args[1])
