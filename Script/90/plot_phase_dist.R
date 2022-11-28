suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(readr))  
suppressPackageStartupMessages(library(ggplot2))  
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra))

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#merge table
mertab = function(x,y,z){
  x0 = x[-1]
  y0 = y[-1]
  z0 = z[-1]
  c(unlist(x0),unlist(y0),unlist(z0))
}

#generate a mean/var table
mean_var_tb = function(x){
  vars= x %>% summarise_if(is.numeric,var) 
  avgs= colMeans(x) 
  tb  = rbind(avgs,vars) %>% round(4)
  rownames(tb) = c('mean','var')
  tb
}


#grouped boxplot
plt_gbox = function(x,y,z){
  Methods    = rep(c("mafft+sw","coati-max","coati-sampling"),each=270)
  Phases     = rep(c("phase-0","phase-1","phase-2"),each=90)
  Proportion = mertab(x,y,z)
  df         = data.frame(Methods,Phases,Proportion)
  df$Methods = factor(Methods,levels=unique(Methods))
  g          = ggplot(df, aes(x=Methods,y=Proportion,color=Phases)) + geom_boxplot() + xlab("") +
               geom_point(size=1,position=position_jitterdodge(jitter.width=0.1)) + 
               theme(axis.text.x = element_text(colour=c("red","green","blue")))
  
  #create a table
  tb         = cbind(mean_var_tb(x[-1]),mean_var_tb(y[-1]),mean_var_tb(z[-1]))
  cols       = matrix("black", nrow(tb), ncol(tb))
  col.head   = rep(c("red","green","blue"),each=3)
  mytab      = tableGrob(tb,theme=ttheme_default(base_size=7,core=list(fg_params = list(col=cols),
                                                                       bg_params = list(col=NA)),
                                                 rowhead=list(bg_params = list(col=NA)),
                                                 colhead=list(bg_params = list(col=col.head))))     
  
  g + annotation_custom(mytab, xmin=0.25,xmax=Inf,ymin=0.55,ymax=0.68)
}

#############################################
#inD   = "Results/phases"
#ouFig = "Figure/ZnZs/phase_dist.pdf"

main = function(inD,ouFig){
  Files  = list.files(inD,full.names=T)
  
  df.max = read.table(Files[1],header=T,sep='\t',row.names=NULL)
  df.sam = read.table(Files[2],header=T,sep='\t',row.names=NULL)
  df.sw  = read.table(Files[3],header=T,sep='\t',row.names=NULL)
  
 
  
  
  
  
  
  pdf(ouFig,onefile=T)
  gg  = plt_gbox(df.sw,df.max,df.sam)
  print(gg)
  dev.off()
  
}


###############################
args=commandArgs(trailingOnly=T)
main(args[1],args[2])