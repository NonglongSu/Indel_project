suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(readr))  
suppressPackageStartupMessages(library(ggplot2))  
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra))


#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#est dnds
dnds_est = function(x){
  Nd = sum(x$Nd)
  Sd = sum(x$Sd)
  N  = sum(x$N)
  S  = sum(x$S)
  
  PN = (Nd/N) 
  PS = (Sd/S) 
  
  #jukes cantor formula (1986)
  dNdS = log(1-4*PN/3)/log(1-4*PS/3)
  return(dNdS)
}

#est ZnZs
ZnZs_est = function(x){
  tyN = sum(x$Zn)
  tyS = sum(x$Zs)
  c(tyN,tyS)
}

#generate dnds,ZnZs matrix
ZnZs_gen = function(x){
  ZnZs_90 = matrix(0,90,3)
  for(i in 1:length(x)){
    spec.sum    = read_tsv(x[i],col_names=T,show_col_types=F)
    ZnZs_90[i,1:2] = ZnZs_est(spec.sum) 
    ZnZs_90[i,3]   = dnds_est(spec.sum) 
  }
  colnames(ZnZs_90) = c('Zn','Zs', 'dNdS')
  return(ZnZs_90)
}


######################################
#stack barplot
plt_bar = function(x){
  x            = x[,1:2]
  Species      = rep(1:90,each=2)
  ZnZs         = rep(c("Zn","Zs"),90)
  Proportion   =  c(t(x))
  df           = data.frame(Species,ZnZs,Proportion)
  ggplot(df, aes(fill=ZnZs, y=Proportion, x=Species)) + geom_bar(position="stack", stat="identity") + ylim(0,30000)
}

#grouped boxplot
plt_gbox = function(x,y,z){
  Methods    = rep(c("mafft+sw","coati-max","coati-sampling"),each=180)
  Stats      = rep(c("dN/dS","Zn/(Zn+Zs)"),each=90)
  Ratio      = c(c(x),c(y),c(z))
  df         = data.frame(Methods,Stats,Ratio)
  df$Methods = factor(Methods,levels=unique(Methods))
  g          = ggplot(df, aes(x=Methods,y=Ratio,color=Stats)) + geom_boxplot() + xlab("") +
    geom_point(size=1,position=position_jitterdodge(jitter.width=0.1)) + 
    theme(axis.text.x = element_text(colour=c("red","green","blue")))
  
  #create a table
  tb         = cbind(mean_var_tb(x),mean_var_tb(y),mean_var_tb(z))
  cols       = matrix("black", nrow(tb), ncol(tb))
  col.head   = rep(c("red","green","blue"),each=2)
  mytab      = tableGrob(tb,theme=ttheme_default(base_size=8,core=list(fg_params = list(col=cols),
                                                                        bg_params = list(col=NA)),
                                                 rowhead=list(bg_params = list(col=NA)),
                                                 colhead=list(bg_params = list(col=col.head))))     
  
  g + annotation_custom(mytab, xmin=-0.3,xmax=Inf,ymin=0.45,ymax=0.55)
}

#generate a mean/var table
#mean_var_vec = function(x){c(round(mean(x),4),round(var(x),4))}
mean_var_tb = function(x){
  xdf = data.frame(x)
  vars= xdf %>% summarise_if(is.numeric,var) 
  avgs= colMeans(x) 
  tb  = rbind(avgs,vars) %>% round(4)
  rownames(tb) = c('mean','var')
  colnames(tb) = c('dN/dS','Zn/(Zn+Zs)') 
  tb
}

#wilcox test 
plot_p_stats = function(before,after,m1,m2){
  Methods = rep(c(m1,m2),each=90)
  Ratio   = c(before,after)
  df      = data.frame(Methods,Ratio)
  df$Methods = factor(Methods,levels=unique(Methods))
  
  g = ggpaired(df,x='Methods',y='Ratio',color='Methods',line.color='gray',line.size=0.4,palette=NULL) + 
      xlab("") + ylab("Zn/(Zn+Zs)") + ylim(0,0.4)
  g + stat_compare_means(method='wilcox.test',paired=T,label = "p.format",vjust=2,hjust=-1) + scale_color_discrete(name="") + theme(legend.position="right") 
}


#######################################################
#inD  = "Results/ZD_sum"
#ouFig= "Figure/ZnZs/ZnZs_prop.pdf"
main = function(inD,ouFig){
  F.sw  = list.files(inD,full.names=T,pattern='sw')
  F.max = list.files(inD,full.names=T,pattern='max')
  F.sam = list.files(inD,full.names=T,pattern='sample')
  
  dmat.sw  = ZnZs_gen(F.sw)
  dmat.max = ZnZs_gen(F.max)
  dmat.sam = ZnZs_gen(F.sam)
  
  
  #############################################
  #proportion-histogram plot 
  gg.sw  = plt_bar(dmat.sw)
  gg.max = plt_bar(dmat.max)
  gg.sam = plt_bar(dmat.sam)
  gg     = ggarrange(gg.sw, gg.max, gg.sam, labels=c("A","B","C"), ncol=1,nrow=3)
  
  #proportion-distribution plot
  x0 = dmat.sw[,1]/rowSums(dmat.sw[,-3])
  y0 = dmat.max[,1]/rowSums(dmat.max[,-3])
  z0 = dmat.sam[,1]/rowSums(dmat.sam[,-3])
  
  dMat.sw  = cbind(dmat.sw[,3],x0)
  dMat.max = cbind(dmat.max[,3],y0)
  dMat.sam = cbind(dmat.sam[,3],z0)
  gg2      = plt_gbox(dMat.sw,dMat.max,dMat.sam)
  
  #Wilcoxon signed-rank test (all rejects!)
  # wilcox.test(x0,y0,paired=T)
  # wilcox.test(x0,z0,paired=T)
  # wilcox.test(y0,z0,paired=F)
  gg_1 = plot_p_stats(x0,y0,"mafft+sw","coati-max")
  gg_2 = plot_p_stats(x0,z0,"mafft+sw","coati-sampling")
  gg_3 = plot_p_stats(y0,z0,"coati-max","coati-sampling")
  gg3  = ggarrange(gg_1, gg_2, gg_3, ncol=1,nrow=3,labels=c("A","B","C")) + theme(plot.margin = margin(t=5,r=50,b=5,l=40))  
  #gg3  = ggarrange(gg_1, gg_2, gg_3, ncol=1,nrow=3,labels=c("A","B","C")) + theme(plot.margin = margin(t=1,r=1.25,b=1,l=1.25,unit='in'))  
  
  
  pdf(ouFig,onefile=T,paper='a4',width=8, height=15)
  print(gg)
  print(gg2)
  print(gg3)                                                                         
  dev.off()
}


###############################
args=commandArgs(trailingOnly=T)
main(args[1],args[2])



