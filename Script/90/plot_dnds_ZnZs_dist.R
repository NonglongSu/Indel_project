suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(utils))
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
znzs_est = function(x){
  tyN = sum(x$Zn)
  tyS = sum(x$Zs)
  return(tyN/tyS)
}

#generate dnds,ZnZs matrix
mat_gen = function(x){
  dnds_90 = c()
  ZnZs_90 = c()
  for(i in 1:length(x)){
    spec.sum   = read_tsv(x[i],col_names=T,show_col_types=F)
    dnds_90[i] = dnds_est(spec.sum)
    ZnZs_90[i] = znzs_est(spec.sum)
  }
  dmat = matrix(c(dnds_90,ZnZs_90),90,2)
  colnames(dmat) = c('dnds','ZnZs')
  cat(sprintf("dnds-mean:%.3f ZnZs-mean:%.3f \n dnds-var:%.3f ZnZs-var:%.3f",
              mean(dnds_90),mean(ZnZs_90),var(dnds_90),var(ZnZs_90)),"\n")
  return(dmat)
}



##################################
#generate a mean/var table
mean_var_tb = function(x){
  xdf = data.frame(x)
  vars= xdf %>% summarise_if(is.numeric,var) 
  avgs= colMeans(x) 
  tb  = rbind(avgs,vars) %>% round(4)
  rownames(tb) = c('mean','var')
  tb
}

#grouped boxplot
plt_gbox = function(x,y,z){
  Methods    = rep(c("mafft+sw","coati-max","coati-sampling"),each=180)
  Stats      = rep(c("dnds","ZnZs"),each=90)
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
  mytab      = tableGrob(tb,theme=ttheme_default(base_size=10,core=list(fg_params = list(col=cols),
                                                                        bg_params = list(col=NA)),
                                                                        rowhead=list(bg_params = list(col=NA)),
                                                                        colhead=list(bg_params = list(col=col.head))))     
  
  g + annotation_custom(mytab, xmin=-0.6,xmax=Inf,ymin=0.55,ymax=0.65)
}


#ggplot
# pltt = function(x){
#   tb  = mean_var_tb(x)
#   Pdf = as.data.frame(as.table(x))
#   Pdf = Pdf %>% rename(Ratio=Var2)
#   g = ggplot(Pdf, aes(x=Ratio,y=Freq,color=Ratio)) + geom_boxplot() + xlab("") + ylab("Ratio") + ylim(0,0.6)
#   g + geom_jitter(size=1,alpha=0.25,width=0.2) + stat_summary(fun=mean,geom='point',size=2) +
#       annotation_custom(tableGrob(tb,theme=ttheme_default(base_size=6)), xmin=-1.2,xmax=Inf,ymin=0.4,ymax=0.6)
# }

#######################################################
#inD  = "Results/ZD_sum"
#ouFig= "Figure/ZnZs/dnds_ZnZs.pdf"
main = function(inD,ouFig){
  F.sw  = list.files(inD,full.names=T,pattern='sw')
  F.max = list.files(inD,full.names=T,pattern='max')
  F.sam = list.files(inD,full.names=T,pattern='sample')
  
  dmat.sw  = mat_gen(F.sw)
  dmat.max = mat_gen(F.max)
  dmat.sam = mat_gen(F.sam)
  

  pdf(ouFig,onefile=T)
  gg = plt_gbox(dmat.sw,dmat.max,dmat.sam)
  print(gg)
  dev.off()
  
  
  #More test
  dnds = c(dmat.sw[,1],dmat.max[,1],dmat.sam[,1])
  ZnZs = c(dmat.sw[,2],dmat.max[,2],dmat.sam[,2])
  met  = rep(c("mafft+sw","coati-max","coati-sampling"),each=90)
  df   = data.frame('Methods'=met,'dnds'=dnds,'ZnZs'=ZnZs)
  
  
  ##one-way ANOVA of dnds
  aov.dnds = aov(dnds~Methods,data=df)
  summary(aov.dnds)
  
  ##normality test (violated)
  par(mfrow=c(1,2))
  res_dnds   = aov.dnds$residuals
  res_dnds_z = (res_dnds-mean(res_dnds))/sd(res_dnds)
  hist(res_dnds)
  qqnorm(res_dnds_z,xlim=c(-3,6),ylim=c(-3,6))
  abline(0,1,col='red')
  
  shapiro.test(res_dnds)
  
  ##try nonparameteric test instead of ANOVA (rejected)
  kruskal.test(dnds~Methods,data=df)
  
  ##post-hoc test
  library(FSA)
  dunnTest(dnds~Methods,data=df,method='holm')
  
  
  
  ##one-way ANOVA of ZnZs
  aov.ZnZs = aov(ZnZs~Methods,data=df)
  summary(aov.ZnZs)
  
  ##normality test (violated)
  par(mfrow=c(1,2))
  res_ZnZs   = aov.ZnZs$residuals
  res_ZnZs_z = (res_ZnZs-mean(res_ZnZs))/sd(res_ZnZs)
  hist(res_ZnZs)
  qqnorm(res_ZnZs_z,xlim=c(-3,6),ylim=c(-3,6))
  abline(0,1,col='red')
  
  shapiro.test(res_ZnZs)
  
  ##try nonparameteric test instead of ANOVA (rejected)
  kruskal.test(ZnZs~Methods,data=df)
  
  ##post-hoc test (max-sampling not rejected)
  dunnTest(ZnZs~Methods,data=df,method='holm')
  
  
  
}




##################################
args=commandArgs(trailingOnly=T)
main(args[1],args[2])
