#Data analysis of simulation(100) and real data(90)
#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))


inD    = "Results/est100"
inD    = "Results/nmkb100.1"
inF    = "Results/truePar_100.txt"
pat    = "5e"
ouFig1 = "Figure/EM/est100.5e.pdf"
ouFig2 = "Figure/EM/qqplot.5e.pdf"
ouF    = "Results/rmse.5e.txt"
##########################################################
main = function(inF, inD, ouFig1, ouFig2, ouF, pat){
  Files     = list.files(inD, pattern=pat, full.names=T)             
  pre.order = as.numeric(gsub('.*[\\/]([^.]+)[.].*','\\1',Files))
  Files     = Files[order(pre.order)]
  
  #true par
  trueP = read.table(inF, header=T, sep="")
  sst   = trueP[,5:10]/trueP[,11] 
  tMat  = as.matrix(cbind(sst,trueP[,12:11]))
  
  #est 
  dmat = matrix(0,length(Files),8)
  for(i in 1:length(Files)){
    Jtmp     = fromJSON(Files[i])
    if(!is.null(Jtmp$p.intermed)){
      dmat[i,] = c(Jtmp$par, Jtmp$p.intermed[nrow(Jtmp$p.intermed),10])
    }else{
      dmat[i,] = c(Jtmp$par[1:8])
    }
  }
  ssv  = dmat[,1:6]/dmat[,8]
  dMat = cbind(ssv, dmat[,7], dmat[,8]) 
  
  #>sorted boxplot
  # dF = as.data.frame(as.table(dMat))
  # dF = dF[,-1]
  # colnames(dF) = c('paras','value')
  # 
  # ggplot(dF, aes(x=paras, y=value)) + geom_boxplot()
  # dF_sorted = dF %>% mutate(paras = fct_reorder(paras,-value))
  # ggplot(dF_sorted, aes(x=paras, y=value)) + geom_boxplot()
  
  #>flip he boxplot
  #ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_boxplot() + coord_flip() 
 
  #>violin plot
  #ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_violin() + coord_flip()
  
  #>lineplot
  #ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_line(size=1) + coord_flip()
   
  #>Dot strip plot
  #ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_point(size=2, alpha=0.6) + coord_flip()
  

  #>boxplot(median)+dotplot
  #ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + geom_boxplot() + coord_flip()
  
  #>hightlight the mean
  #set.seed(2021)
  #g + geom_jitter(size=2,alpha=0.25,width=0.2) + stat_summary(fun=mean,geom='point',size=5)
  
  # pdf(ouFig1)
  # g     = ggplot(dF_sorted, aes(x=paras,y=value,color=paras)) + coord_flip()
  # plot1 = g + geom_boxplot(outlier.alpha=0) + geom_point(size=2, alpha=0.6) 
  # print(plot1)
  # dev.off()
  

  #>qqplot (true vs est)
  para = c(bquote(sigma*1),bquote(sigma*2),bquote(sigma*3),bquote(sigma*4),bquote(sigma*5),bquote(sigma*6),bquote(omega),bquote(tau))
  pdf(ouFig2,paper='a4', width=8, height=15, pointsize=18)
  par(mfrow=c(3,3),mai=c(1, 0.3, 0.1, 0.25), mgp=c(3,0.6,0))
  colorblind_palette1 = c("#CC6677")
  
  for (i in seq(8)) {
    qqplot(tMat[,i],dMat[,i],xlab=para[i],ylab="")
    abline(0,1,col=colorblind_palette1,lwd=2)
  }
  dev.off()
  
  # for (i in seq(6)) {
  #   qqplot(tMat[,i],dMat[,i],xlab=colnames(tMat)[i],ylab=colnames(dMat)[i])
  #   abline(0,1,col=colorblind_palette1,lwd=2)
  # }
  # par(mfrow=c(1,2),mar=c(15.1, 4.1, 1, 2.1))
  # for (i in 7:8) {
  #   qqplot(tMat[,i],dMat[,i],xlab=colnames(tMat)[i],ylab=colnames(dMat)[i])
  #   abline(0,1,col=colorblind_palette1,lwd=2)
  # }
  # dev.off()
  
  
  
  #>rmse output
  #column-wise
  rmse1 = rep(0,8)
  for (i in 1:8) {
    rmse1[i] = sqrt(crossprod(tMat[,i]-dMat[,i])/100)
  }
  #print(rmse1)
  
  #row-wise
  rmse2 = rep(0,100)
  for (i in 1:100) {
    rmse2[i] = sqrt(crossprod(tMat[i,]-dMat[i,])/8)
  }
  #print(rmse2)
  write(rmse2,ouF,sep='\t',ncolumns=1)
  plot(density(rmse2))
  rm2.avg = mean(rmse2)
  rm2.sd  = sd(rmse2)
  #95% CI
  rm2.lower  = rm2.avg-qnorm(.95)*rm2.sd/sqrt(100)
  rm2.higher = rm2.avg+qnorm(.95)*rm2.sd/sqrt(100)
  
  #print(c(rm2.lower,rm2.higher))
}


########################
args = commandArgs(trailingOnly=TRUE)
main(args[1],args[2],args[3],args[4], args[5],args[6])


#>>Paired t-test of theta: sample size-100. 
# se   = rep(0,8)
# pval = rep(0,8) 
# for (i in 1:8) {
#   ttest   = t.test(tMat[,i], dMat[,i], paired=TRUE)
#   se[i]   = ttest$stderr
#   pval[i] = ttest$p.value
# }
#print(pval)



#>>>>figure out how to draw LL contour 


