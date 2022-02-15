# Calculate the Displacement distribution of INDELs of different lengths. 

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(stringr)
library(ggplot2)  

# setwd("~/Dropbox (ASU)/Indel_project/Script")

phasing = function(file){
  
  dna = readBStringSet(file, format = "fasta")
  l   = length(dna)
  
  dna.1 = str_split(as.character(dna), '')
  g1 = lapply(dna.1, function(x) { IRanges(x %in% c("-", "+"))})
  g2 = lapply(dna.1, function(x) { IRanges(x == "+")})
  g1 = IRangesList(g1)
  g2 = IRangesList(g2)
  
  mr     = g1[(l-1):l]
  wid.mr = unlist(width(mr))
  pos.mr = unlist(start(mr))
  
  # set up flag as "+++"
  flag  = rep(0, length(wid.mr))
  flag1 = unlist(start(g2[(l-1):l]))
  flag[which(pos.mr %in% flag1)] = rep(1, length(flag1))
  
  df.mr = data.frame("pos" = pos.mr, "wid" = wid.mr, "flag" = flag)
  return(df.mr)
}

# dir1   = "../test_human_mouse_rat/Data_6/Mafft/updated_cds"
# dir2   = "../test_human_mouse_rat/Data_6/Mafft/mapped_cds"
# ouFile = "../test_human_mouse_rat/Data_6/Results/dis.mafft.txt"
# ouFig  = "../test_human_mouse_rat/Data_6/Figure/dis.mafft.pdf"

main = function(dir1, dir2, ouFile, ouFig, num){
  
  Files = list.files(dir2, full.names = FALSE)
  
  PST.1 = c()
  PST.2 = c()
  for(i in 1:length(Files)){
    file1 = paste0(dir1, "/", Files[i])
    file2 = paste0(dir2, "/", Files[i])
    
    pScore.1 = phasing(file1)
    pScore.2 = phasing(file2)
    PST.1    = rbind(PST.1, pScore.1)
    PST.2    = rbind(PST.2, pScore.2)
  }
  
  # Filter out the "+++"
    PST.2f = PST.2[PST.2$flag == 0, ]
    PST.1f = PST.1[which(PST.2$flag == 0), ]        
    
  # updated_cds
  pos1.3  = PST.1f[PST.1f$wid == 3, ]$pos
  pos1.6  = PST.1f[PST.1f$wid == 6, ]$pos
  pos1.9  = PST.1f[PST.1f$wid == 9, ]$pos
  pos1.12 = PST.1f[PST.1f$wid == 12,]$pos
  
  # mapped_cds
  pos2.3  = PST.2f[PST.2f$wid == 3, ]$pos
  pos2.6  = PST.2f[PST.2f$wid == 6, ]$pos
  pos2.9  = PST.2f[PST.2f$wid == 9, ]$pos
  pos2.12 = PST.2f[PST.2f$wid == 12,]$pos
  
  diff.3  = pos2.3  - pos1.3
  diff.6  = pos2.6  - pos1.6
  diff.9  = pos2.9  - pos1.9
  diff.12 = pos2.12 - pos1.12
  
  #Create diplacement column
  Window = as.numeric(num)
  
  Dis  = seq(-Window, Window, 1)
  G.3  = c()
  G.6  = c()
  G.9  = c()
  G.12 = c()
  
  for(i in 1:length(Dis)){
    G.3[i]  = length(diff.3[diff.3 == Dis[i]])
    G.6[i]  = length(diff.6[diff.6 == Dis[i]])
    G.9[i]  = length(diff.9[diff.9 == Dis[i]])
    G.12[i] = length(diff.12[diff.12 == Dis[i]])
  }
  GD.df1 = rbind(G.3, G.6, G.9, G.12) 
  GD.df2 = rbind(G.3/sum(G.3), G.6/sum(G.6), G.9/sum(G.9), G.12/sum(G.12)) * 100
  GD.df2 = round(GD.df2, 2)
  colnames(GD.df2) = Dis
  rownames(GD.df2) = c(3, 6, 9, 12)
  
  # Output
  write.table(GD.df2, file = ouFile, sep = "\t", append = FALSE, quote = FALSE, row.names = TRUE, col.names = NA)
  
  
  # Stacked area chart
  distance = rep(Dis, each = 4)
  value    = as.vector(GD.df1)
  gap.len  = rep(c("A-3", "B-6", "C-9", "D-12"), time = length(Dis))
  data     = data.frame(distance, value, gap.len)
  
  dis.plot = ggplot(data, aes(x = distance, y = value, fill = gap.len)) + 
    geom_area()
  print(dis.plot + ggtitle("Displacement distribution across lengths of gaps") + theme(plot.title = element_text(face = "bold", hjust = 0.5)))
  
  # Output the plot
  ggsave(ouFig)
  
}

#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4], args[5])
