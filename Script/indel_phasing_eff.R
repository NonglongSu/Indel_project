# Calculate the proportion of effective phase 1 & phase 2 INDELs. 
# Define effective as "trigger  a substitution"
# mututal anacestors for each other

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(stringr)
library(dplyr)

# setwd("~/Dropbox (ASU)/Indel_project/Script")

eff_phase = function(seq1, seq2, pos, wid, len){ #seq1---reference
  seq2.char = str_split(seq2, "")[[1]]
  eff.sub.1 = 0
  non.sub.1 = 0
  eff.sub.2 = 0
  non.sub.2 = 0
  for(j in 1: len){
    if(pos[j] %% 3 == 1){ #phase-0
      next
    }else if(pos[j] %% 3 == 0){ #phase-2
      pos.ori = pos[j] - 2 
      unit.1 = substr(seq1, pos.ori, pos.ori + 2)
      unit.2 = substr(seq1, pos.ori + wid[j], pos.ori + wid[j] + 2)
      sub = paste0(c(seq2.char[pos[j] - 2], seq2.char[pos[j] - 1], seq2.char[pos[j] + wid[j]]), collapse = "")
      sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
      if(unit.1 %in% sec || unit.2 %in% sec){
        non.sub.2 = non.sub.2 + 1
      }else{
        eff.sub.2 = eff.sub.2 + 1
      }
    }else{ #phase-1
      pos.ori = pos - 1 
      unit.1 = substr(seq1, pos.ori, pos.ori + 2)
      unit.2 = substr(seq1, pos.ori + wid[j], pos.ori + wid[j] + 2)
      sub = paste0(c(seq2.char[pos[j] - 1], seq2.char[pos[j] + wid[j]], seq2.char[pos[j] + wid[j] + 1]), collapse = "")
      sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
      if(unit.1 %in% sec || unit.2 %in% sec){
        non.sub.1 = non.sub.1 + 1
      }else{
        eff.sub.1 = eff.sub.1 + 1
      }
    }
  }  
 res = c(eff.sub.1, non.sub.1, eff.sub.2, non.sub.2) 
 return(res)
}

# dir    = "../test_human_mouse_rat/Data_6/Mafft/mapped_cds"
# ouFig  = "../test_human_mouse_rat/Data_6/Figure/phase.mafft.eff.pdf"

main = function(dir, ouFig){
  
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
  
  M.eff.1 = 0 ; R.eff.1 = 0
  M.non.1 = 0 ; R.non.1 = 0
  M.eff.2 = 0 ; R.eff.2 = 0
  M.non.2 = 0 ; R.non.2 = 0
  
  File.total = list.files(dir, full.names = TRUE)
  for(i in 1:length(File.total)){
    dna = readBStringSet(File.total[i], format = "fasta")
    len = length(dna)
    M   = toString(dna[[len - 1]])
    R   = toString(dna[[len]])
    
    dna.1 = str_split(as.character(dna), '')
    g     = lapply(dna.1, function(x) {IRanges(x == '-')})
    g     = IRangesList(g)
    
    m     = g[[len -1 ]]
    r     = g[[len]]
    wid.m = width(m)
    wid.r = width(r)
    pos.m = start(m)
    pos.r = start(r)
    l.m   = length(wid.m)
    l.r   = length(wid.r)
    
    if(l.m > 0){
      eff_M   = eff_phase(R, M, pos.m, wid.m, l.m)    
      M.eff.1 = M.eff.1 + eff_M[1]
      M.non.1 = M.non.1 + eff_M[2]
      M.eff.2 = M.eff.2 + eff_M[3] 
      M.non.2 = M.non.2 + eff_M[4] 
    }
    if(l.r > 0){
      eff_R   = eff_phase(M, R, pos.r, wid.r, l.r)
      R.eff.1 = R.eff.1 + eff_R[1]
      R.non.1 = R.non.1 + eff_R[2]
      R.eff.2 = R.eff.2 + eff_R[3] 
      R.non.2 = R.non.2 + eff_R[4] 
    }
  }
  
  # Generate an effective phase table
  Eff1     = M.eff.1 + R.eff.1 
  N.eff1   = M.non.1 + R.non.1
  Eff2     = M.eff.2 + R.eff.2
  N.eff2   = M.non.2 + R.non.2
  Phase.df = data.frame("Phase_1" = c(Eff1, N.eff1), "Phase_2" = c(Eff2, N.eff2))
  rownames(Phase.df) = c("Effective", "Non-effective")
  
  
  # Plot
  pdf(ouFig)
  safe_color  = c("#88CCEE", "#CC6677")
  upper.bound = max(colSums(Phase.df)) + 100
  barplot(as.matrix(Phase.df), main = "Effective vs Non-effective proportion of phases across focal species",
          names.arg =  colnames(Phase.df), col = safe_color, cex.main = 0.8, ylim = c(0, upper.bound), cex.main = 1)
  legend("topleft", legend = c("Effective", "Non-effective"), bty = "n", bg ="transparent", fill = safe_color)
  dev.off()
}

#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])
