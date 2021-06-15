#Data analysis of simulation(100) and real data(90)
#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))



##>>>est90
inD = "../../test_90_species/Results/est90"

main = function(inD){
  Files = list.files(inD, full.names=TRUE)
  
  dmat = matrix(0,length(Files),8)
  for(i in 1:length(Files)){
    Jtmp     = fromJSON(Files[i])
    dmat[i,] = c(Jtmp$par, Jtmp$p.intermed[nrow(Jtmp$p.intermed),10])
  }
  ssv  = dmat[,1:6]/dmat[,8]
  dMat = cbind(ssv,dmat[,7], dmat[,8]) 
  colnames(dMat)=c('s1','s2','s3','s4','s5','s6','omega','tau')
  
  #>boxplot
  boxplot(dMat, ylab='value')
  
  #>sorted boxplot
  dF = as.data.frame(as.table(dMat))
  colnames(dF) = c('haha','paras','value')
  ggplot(dF, aes(x=paras, y=value)) + geom_boxplot()
  dF_sorted = dF %>% mutate(paras = fct_reorder(paras,-value))
  ggplot(dF_sorted, aes(x=paras, y=value)) + geom_boxplot()
  
  #>flip he boxplot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_boxplot() + coord_flip() 
 
  #>violin plot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_violin() + coord_flip()
  
  #>lineplot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_line(size=1) + coord_flip()
   
  #>Dot strip plot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_point(size=2, alpha=0.6) + coord_flip()
  
  
  #create basic g
  g = ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + coord_flip()
  
  #>boxplot(median)+dotplot
  ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_boxplot() + coord_flip()
  g + geom_boxplot(outlier.alpha=0) + geom_point(size=2, alpha=0.6) 
  
  #>hightlight the mean
  set.seed(2021)
  g + geom_jitter(size=2,alpha=0.25,width=0.2) + stat_summary(fun=mean,geom='point',size=5)
  
  #>>>try w vs tau
  plot(dMat[,7]~dMat[,8], xlab='tau',ylab='omega')
  y = log(dMat[,7])
  x = dMat[,8]
  plot(y~x, xlab='tau',ylab='omega')
  abline(lm(y~x),col='red',lwd=2)
  plot(lowess(x,y),col='blue')
  
  cor(x,y)
  
}

dMat[order(dMat[,8]),]













