#Calculate the proportion of effective phase 1 & phase 2 indels. 
#mututal anacestors for each other.

suppressWarnings(suppressMessages(library(Biostrings)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressPackageStartupMessages(library(stringr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Data_6")

#check if non-synonymous
eff_phase = function(seq1, pos, wid, len){#seq1---reference
  seq1.char = str_split(seq1, "")[[1]]
  eff.sub.1 = 0
  non.sub.1 = 0
  eff.sub.2 = 0
  non.sub.2 = 0
  for(j in 1: len){
    if(pos[j] %% 3 == 1){#phase-0
      next
    }else if(pos[j] %% 3 == 0){#phase-2
      pos.ori = pos[j]-2 
      unit.1  = substr(seq1, pos.ori, pos.ori+2)
      unit.2  = substr(seq1, pos.ori+wid[j], pos.ori+wid[j]+2)
      sub     = paste0(c(seq1.char[pos[j]-2], seq1.char[pos[j]-1], seq1.char[pos[j]+wid[j]]), collapse = "")
      sec     = codon[[which(sapply(codon, function(X){sub %in% X}))]]
      if(unit.1 %in% sec || unit.2 %in% sec){
        non.sub.2 = non.sub.2 + 1
      }else{
        eff.sub.2 = eff.sub.2 + 1
      }
    }else{ #phase-1
      pos.ori = pos-1 
      unit.1  = substr(seq1, pos.ori, pos.ori+2)
      unit.2  = substr(seq1, pos.ori + wid[j], pos.ori+wid[j]+2)
      sub     = paste0(c(seq1.char[pos[j]-1], seq1.char[pos[j]+wid[j]], seq1.char[pos[j]+wid[j]+1]), collapse = "")
      sec     = codon[[which(sapply(codon, function(X){sub %in% X}))]]
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

#generate a effective table
eff_gen = function(dir){
  eff.1 = 0
  non.1 = 0
  eff.2 = 0
  non.2 = 0
  
  File.total = list.files(dir,full.names=TRUE)
  for(i in 1:length(File.total)){
    dna = readBStringSet(File.total[i])
    len = length(dna)
    M   = toString(dna[[len-1]])
    R   = toString(dna[[len]])
    
    dna.1 = str_split(as.character(dna),'')
    g     = lapply(dna.1, function(x) {IRanges(x=='-')})
    g     = IRangesList(g)
    
    m     = g[[len-1]]
    r     = g[[len]]
    wid.m = width(m)
    wid.r = width(r)
    pos.m = start(m)
    pos.r = start(r)
    l.m   = length(wid.m)
    l.r   = length(wid.r)
    
    if(l.m > 0){#mouse
      eff_M = eff_phase(R,pos.m,wid.m,l.m)    
      eff.1 = eff.1 + eff_M[1]
      non.1 = non.1 + eff_M[2]
      eff.2 = eff.2 + eff_M[3] 
      non.2 = non.2 + eff_M[4] 
    }
    if(l.r > 0){#rat
      eff_R = eff_phase(M,pos.r,wid.r,l.r)
      eff.1 = eff.1 + eff_R[1]
      non.1 = non.1 + eff_R[2]
      eff.2 = eff.2 + eff_R[3] 
      non.2 = non.2 + eff_R[4] 
    }
    print(i)
  }
  
 ##generate an effective phase table
 data.frame("Phase_1" = c(eff.1,non.1), "Phase_2" = c(eff.2,non.2))
}

#ggplot of the eff phase
gg_eff = function(dat,win){
  safe_color  = c("#E69F00", "#56B4E9")
  upper.bound = max(colSums(dat)) + 100
  barplot(as.matrix(dat), main=NULL, xlab=paste0("window size of ",win),
          names.arg=colnames(dat), col=safe_color, cex.main=1, ylim=c(0,upper.bound), cex.main=1)
  legend("topleft", legend=c("typeN","typeS"), bty="n", bg="transparent", fill=safe_color)
}


###############################################
#ouFig = "Data_6/Figure/phase.mafft.eff.pdf"

main = function(ouFig){
  
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
  
  dir1  = "../Data_3/Mafft/mapped_cds"
  dir2  = "Mafft/mapped_cds"
  dir3  = "../Data_9/Mafft/mapped_cds"
  dir4  = "../Data_12/Mafft/mapped_cds"
  
  dat1 = eff_gen(dir1)
  dat2 = eff_gen(dir2)
  dat3 = eff_gen(dir3)
  dat4 = eff_gen(dir4)
  
  #ggplot
  pdf(ouFig,onefile=T)
  par(mfrow=c(2,2))
  gg_eff(dat1,3)
  gg_eff(dat2,6)
  gg_eff(dat3,9)
  gg_eff(dat4,12)
  
  dev.off()   
  
}

#######################################
args = commandArgs(trailingOnly=TRUE)
main(args[1])
