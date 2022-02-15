# Calculate the proportion of phase 0 / phase 1 / phase 2 INDELs. 

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(stringr)

# setwd("~/Dropbox (ASU)/Indel_project/Script")

phasing = function(file){
  dna = readBStringSet(file, format = "fasta")
  l   = length(dna)
 
  dna.1 = str_split(as.character(dna), '')
  g = lapply(dna.1, function(x) { IRanges(x == '-')})
  g = IRangesList(g)
  
  m = g[[l - 1]]
  r = g[[l]]
  
  wid.m = width(m)
  wid.r = width(r)
  
  l.m = length(wid.m)
  l.r = length(wid.r)
  
  pos.m = start(m)
  pos.r = start(r)
  
  
  df.m = data.frame("pos" = pos.m, "wid" = wid.m)
  df.r = data.frame("pos" = pos.r, "wid" = wid.r)
  
  if(l.m > 0 & l.r == 0){
    return(df.m)
  }else if(l.m == 0 & l.r > 0){
    return(df.r)
  }else if(l.m > 0 & l.r > 0){
    df = merge(df.m, df.r, all = TRUE)
    return(df)
  }else{
    return(NULL)
  }
  
}


Record = function(PST){
  if(is.null(PST)){
    return(PST)
    }
  rem   = PST %% 3
  phase.0 = length(which(rem == 1))
  phase.1 = length(which(rem == 2)) 
  phase.2 = length(which(rem == 0)) 
  
  df.phase = data.frame(Phase_0 = phase.0, Phase_1 = phase.1, Phase_2 = phase.2)
  return(df.phase)
}  


# file1  = "../test_human_mouse_rat/Data_6/Mafft/mapped_cds/ENSG00000000460.fa"
# ouFile = "../test_human_mouse_rat/Data_6/Results/Phase.mafft.txt"
# ouFig  = "../test_human_mouse_rat/Data_6/Figure/Phase.mafft.pdf"
# dir    = "../test_human_mouse_rat/Data_6/Mafft/mapped_cds"
# Window = 6

main = function(dir, ouFile, ouFig, num){
  
  File.total = list.files(dir, full.names = TRUE)
  Window     = as.numeric(num)
  
  
  PST = c()
  for(i in 1:length(File.total)){
    phase_score = phasing(File.total[i])
    PST         = rbind(PST, phase_score)
  }
  
  pos.3  = PST[PST$wid == 3, ]$pos
  pos.6  = PST[PST$wid == 6, ]$pos
  pos.9  = PST[PST$wid == 9, ]$pos
  pos.12 = PST[PST$wid == 12,]$pos
  
  phase.3  = Record(pos.3)
  phase.6  = Record(pos.6)
  phase.9  = Record(pos.9)
  phase.12 = Record(pos.12)
  
  Phase.df   = rbind(phase.3, phase.6, phase.9, phase.12)
  Phase.df   = rbind(Phase.df, colSums(Phase.df))
  gap_length = c('3', '6', '9', '12', 'T')
  Phase.df   = cbind(gap_length, Phase.df)
  
  #Output a tale
  write.table(Phase.df, file = ouFile, 
              sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  # Draw a figure
  pdf(ouFig)    
  Pdf      = as.table(noquote(t(Phase.df[-5, ])))
  col.name = as.matrix(Pdf[1, ])
  Pdf      = Pdf[-1, ]
  colnames(Pdf) = col.name
  safe_colorblind_palette = c("#88CCEE", "#CC6677", "#DDCC77")
  
  barplot(Pdf, main = "The proportion of indel phases across length of gaps", 
          xlab = paste0("Gap length (window size = ", Window, ")"), ylab = "Count", col = safe_colorblind_palette) 
  legend("topright", bg ="transparent", bty = "n", rownames(Pdf), fill = safe_colorblind_palette, cex = 0.75)
  dev.off()
}

#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4])
