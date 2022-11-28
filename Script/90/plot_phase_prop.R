
suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(readr))  
suppressPackageStartupMessages(library(ggplot2))  
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggtern))

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#stack barplot
plt_bar = function(x){
  Species      = rep(1:90,each=3)
  Phases       = rep(c("phase-0","phase-1","phase-2"),90)
  Proportion   =  c(t(x[-1]))
  df           = data.frame(Species,Phases,Proportion)
  ggplot(df, aes(fill=Phases, y=Proportion, x=Species)) + geom_bar(position="stack", stat="identity")
}

#############################################
# inD   = "Results/phases"
# ouFig = "Figure/ZnZs/phase_prop.pdf"

main = function(inD,ouFig){
  Files  = list.files(inD,full.names=T)
  
  df.max = read.table(Files[1],header=T,sep='\t',row.names=NULL)
  df.sam = read.table(Files[2],header=T,sep='\t',row.names=NULL)
  df.sw  = read.table(Files[3],header=T,sep='\t',row.names=NULL)
  
  pdf(ouFig,onefile=T)
  
  gg.sw  = plt_bar(df.sw)
  gg.max = plt_bar(df.max)
  gg.sam = plt_bar(df.sam)
  
  gg    = ggarrange(gg.sw, gg.max, gg.sam, labels=c("A","B","C"), ncol=1,nrow=3)
  print(gg)
  #dev.off()
  
  ##ternary plot
  spec = rep(1:90,3)
  meth = rep(c('mafft+sw','coati-max','coati-sampling'),each=90)    
  df   = data.frame(phase0=c(df.sw$phase.0,df.max$phase.0,df.sam$phase.0), 
                    phase1=c(df.sw$phase.1,df.max$phase.1,df.sam$phase.1),
                    phase2=c(df.sw$phase.2,df.max$phase.2,df.sam$phase.2),
                    method = meth,
                    species=spec)
  
  #ggtern(data=df.sw, aes(x=phase.0,y=phase.1,z=phase.2)) + geom_point() 
  
  gg1 = ggtern(data=df, aes(x=phase0,y=phase1,z=phase2,color=method)) + geom_point(alpha=0.5) +
        scale_colour_manual(values = c("#CC79A7", "#009E73", "#56B4E9" ))
  print(gg1)
  
  dev.off()
  
  
  #########################find the abnormal species. 
  ph_comp = function(x) which(x$phase.1 > x$phase.2)
  
  ab_sw = ph_comp(df.sw)
  ab_max= ph_comp(df.max)
  ab_sam= ph_comp(df.sam)
  
  intersect(intersect(ab_sam,ab_max),ab_sw)
  
}


###############################
args=commandArgs(trailingOnly=T)
main(args[1],args[2])