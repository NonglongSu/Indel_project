# Calculate the proportion of phase 0 / phase 1 / phase 2 INDELs. 

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(stringr)
library(ggplot2)  

# setwd("~/Dropbox (ASU)/Indel_project/Script/90")

phasing = function(file){
  dna = readBStringSet(file, format = "fasta")
  l   = length(dna)
 
  dna.1 = str_split(as.character(dna), '')
  g     = lapply(dna.1, function(x) { IRanges(x == '-')})
  g     = IRangesList(g)
  
  m = g[[l - 1]]
  r = g[[l]]
  
  wid.m = width(m)
  wid.r = width(r)
  l.m   = length(wid.m)
  l.r   = length(wid.r)
  pos.m = start(m)
  pos.r = start(r)
  
  df.m  = data.frame("pos" = pos.m, "wid" = wid.m)
  df.r  = data.frame("pos" = pos.r, "wid" = wid.r)
  
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


##############################################################

# dir    = "../../test_90_species/Data/sw_cds"
# ouFile = "../../test_90_species/Results/Phase.txt"
# ouFig  = "../../test_90_species/Figure/Phase.pdf"
# num    = 6

main = function(dir, ouFile, ouFig, num){
  
  # Read input 
  Files  = list.files(dir, full.names = TRUE)
  Window = as.numeric(num)
  
  
  PST = list()
  for(i in 1:length(Files)){
    Filei = list.files(Files[i], full.names = TRUE)
    pst   = c()
    for(j in Filei){
      phase_score = phasing(j)
      pst         = rbind(pst, phase_score)
    }
    PST[[i]] = pst
  }
  

  Gap = c()
    for (j in 1:length(PST)) {
      pos     = PST[[j]]$pos
      Gap     =  rbind(Gap,  Record(pos))
  }
  
  species.name = unlist(lapply(Files, function(x){basename(x)}))
  Phase.df     = cbind(species.name, Gap)
  
  #Output a tale
  capture.output(print(Phase.df, print.gap = 3), file = ouFile)
  # write.table(Phase.df, file = ouFile, 
  #             sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  # Draw a figure
  pdf(ouFig)    
  
  Species    = rep(1:90, each = 3)
  Count      = as.vector(t(Gap)) 
  Group      = rep(c("Phase 0", "Phase 1", "Phase 2"), times = 90)
  data       = data.frame(Species, Count, Group)
  phase.plot = ggplot(data, aes(x = Species, y = Count, fill = Group)) + 
    geom_area()
  print(phase.plot + ggtitle("The proportion of three indel phases across 90 species") + 
        theme(plot.title = element_text(face = "bold", hjust = 0.5)))
  ggsave(ouFig)
  dev.off()
}

#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4])
