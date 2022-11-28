# Calculate the Positional difference distribution of indels. 

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))  
suppressPackageStartupMessages(library(ggpubr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Data_6")

#return pos | wid | flag
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

#generate the distance dataframe
dis_gen = function(dir1,dir2){
  Files  = list.files(dir2, full.names=F)
  PST.1 = c()
  PST.2 = c() 
  for(i in 1:length(Files)){
    file1 = paste0(dir1, "/", Files[i])
    file2 = paste0(dir2, "/", Files[i])
    
    pScore.1 = phasing(file1)
    pScore.2 = phasing(file2)
    PST.1    = rbind(PST.1, pScore.1)
    PST.2    = rbind(PST.2, pScore.2)
    print(i)
  }
  
  ##Filter out the "+++"
  PST.2f = PST.2[PST.2$flag == 0,]
  PST.1f = PST.1[which(PST.2$flag == 0),]        
  
  ##updated_cds
  pos1.3  = PST.1f[PST.1f$wid == 3, ]$pos
  pos1.6  = PST.1f[PST.1f$wid == 6, ]$pos
  pos1.9  = PST.1f[PST.1f$wid == 9, ]$pos
  pos1.12 = PST.1f[PST.1f$wid == 12,]$pos
  
  ##mapped_cds
  pos2.3  = PST.2f[PST.2f$wid == 3, ]$pos
  pos2.6  = PST.2f[PST.2f$wid == 6, ]$pos
  pos2.9  = PST.2f[PST.2f$wid == 9, ]$pos
  pos2.12 = PST.2f[PST.2f$wid == 12,]$pos
  
  diff.3  = pos2.3  - pos1.3
  diff.6  = pos2.6  - pos1.6
  diff.9  = pos2.9  - pos1.9
  diff.12 = pos2.12 - pos1.12
  
  dat = data.frame(
    type     = c(rep(3,length(diff.3)),rep(6,length(diff.6)),rep(9,length(diff.9)),rep(12,length(diff.12))),
    distance = c(diff.3,diff.6,diff.9,diff.12)
  )
  dat$type = as.factor(dat$type)
  
  dat
}

#ggplot of each window size
gg_win = function(dat,win){
  g = dat %>% 
    ggplot(aes(x=distance,fill=type)) + 
    geom_density(alpha=0.4) + scale_fill_manual(values=c("#F0E442","#CC79A7","#009E73","#0072B2")) + 
    theme_bw() + labs(fill="gap-length",x="Distance",subtitle=paste0("Window size of ",win)) 
  g
}


####################################################################
#ouFig = "Figure/dis.mafft.pdf"
main = function(ouFig){
  dir1a   = "../Data_3/Mafft/updated_cds"
  dir1b   = "../Data_3/Mafft/mapped_cds"
  dir2a   = "Mafft/updated_cds"
  dir2b   = "Mafft/mapped_cds"
  dir3a   = "../Data_9/Mafft/updated_cds"
  dir3b   = "../Data_9/Mafft/mapped_cds"
  dir4a   = "../Data_12/Mafft/updated_cds"
  dir4b   = "../Data_12/Mafft/mapped_cds"
  
  dat1 = dis_gen(dir1a,dir1b)
  dat2 = dis_gen(dir2a,dir2b)
  dat3 = dis_gen(dir3a,dir3b)
  dat4 = dis_gen(dir4a,dir4b)
  
  gg1 = gg_win(dat1,3)
  gg2 = gg_win(dat2,6)
  gg3 = gg_win(dat3,9)
  gg4 = gg_win(dat4,12)
  
  #ggplot
  pdf(ouFig,onefile=T)
  gg    = ggarrange(gg1, gg2, gg3, gg4, labels=c("A","B","C","D"), ncol=2, nrow=2)
  print(gg)
  dev.off()
}

#######################################
args = commandArgs(trailingOnly=TRUE)
main(args[1])
