#Data analysis of real species(90)

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")


#inD  = "Results/est90"
#ouF1 = "Figure/EM/est90.pdf"
main = function(inD, ouF1, ouF2){
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
  #boxplot(dMat, ylab='value')
  
  #>sorted boxplot
  dF = as.data.frame(as.table(dMat))
  colnames(dF) = c('haha','paras','value')
  ggplot(dF, aes(x=paras, y=value)) + geom_boxplot()
  #dF_sorted = dF %>% mutate(paras = fct_reorder(paras,-value))
  #ggplot(dF_sorted, aes(x=paras, y=value)) + geom_boxplot()
  
  #>flip he boxplot
  # ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_boxplot() + coord_flip() 
  # 
  # #>violin plot
  # ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_violin() + coord_flip()
  # 
  # #>lineplot
  # ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_line(size=1) + coord_flip()
  #  
  # #>Dot strip plot
  # ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_point(size=2, alpha=0.6) + coord_flip()
  
  #>boxplot(median)+dotplot
  #ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_boxplot() + coord_flip()
  
  para = c(bquote(sigma*1),bquote(sigma*2),bquote(sigma*3),bquote(sigma*4),bquote(sigma*5),bquote(sigma*6),bquote(omega),bquote(tau))
  pdf(ouF1)
  gs    = ggplot(dF, aes(x=paras,y=value,color=paras)) + coord_flip() 
  plot1 = gs + geom_boxplot(outlier.alpha=0) + geom_point(size=2, alpha=0.6) + xlab("parameters") 
                                                                                                          
  plot2 = plot1 + scale_x_discrete(labels=para) + scale_color_hue(labels=para) 
  print(plot2)
  dev.off()
  
  #>hightlight the mean
  set.seed(2021)
  gs + geom_jitter(size=2,alpha=0.25,width=0.2) + stat_summary(fun=mean,geom='point',size=5)
  
  
  #>>>>>>>>>>>w vs tau
# df2  = data.frame(omega=dMat[,7], tau=dMat[,8])
# 
# pdf(ouF2)
# plot2 = ggplot(df2, aes(x=tau, y=omega)) + geom_point(alpha=0.6, col="#56B4E9") + 
#                                             geom_text_repel(label=rownames(df2), size=2.5)
# print(plot2)
# dev.off()
  
  ##
  # plot(dMat[,7]~dMat[,8], xlab='tau',ylab='omega')
  # y = log(dMat[,7])
  # x = dMat[,8]
  # plot(y~x, xlab='tau',ylab='omega')
  # abline(lm(y~x),col='red',lwd=2)
  # plot(lowess(x,y),col='blue')
  # cor(x,y)
  
}

#dMat[order(dMat[,8]),]


################################
args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2], args[3])










