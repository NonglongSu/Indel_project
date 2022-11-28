#compare the omega (phylo-EM) with dnds (Nei Gojo (1986))

suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggpubr))

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")


# inD1 = "../chapter3/90/Results/PISE"
# inD2 = "Results/dnds"
# ouFig= "Figure/ZnZs/omega_dnds.pdf"
main = function(inD1,inD2,ouFig){
  Files1   = list.files(inD1,full.names=T,'est')
  Files2.1 = list.files(inD2,full.names=T,pattern='sw')
  Files2.2 = list.files(inD2,full.names=T,pattern='max')
  Files2.3 = list.files(inD2,full.names=T,pattern='sample')
  n        = length(Files1)
  
  omega  = rep(0,n)
  dnds1  = dnds2  = dnds3  = rep(0,n)
  
  for(i in 1:n){
    j.tab     = fromJSON(Files1[i])
    omega[i]  = j.tab$omega
    
    d1.tab    = read.table(Files2.1[i],header=F)
    d2.tab    = read.table(Files2.2[i],header=F)
    d3.tab    = read.table(Files2.3[i],header=F)
    
    dnds1[i]  = d1.tab[[1]]
    dnds2[i]  = d2.tab[[1]]
    dnds3[i]  = d3.tab[[1]]
  }
  
  
  # omat = cbind(omega,dnds1,dnds2,dnds3)
  # colnames(omat) = c('omega','dnds.sw','dnds.max','dnds.sampling')
  # write.csv(omat,"omega_dnds.csv",row.names=F,quote=F)
  
  #single-norm
  par(mfrow=c(3,1))
  
  plot(dnds1~omega,xlim=c(0,0.6),ylim=c(0,0.6),main='sw',ylab='dnds')
  lm.mod1 = lm(dnds1~omega)
  ypred1  = predict(lm.mod1,list(omega))
  lines(omega,ypred1,col="#661100",lwd=2)
  
  plot(dnds2~omega,xlim=c(0,0.6),ylim=c(0,0.6),main='max',ylab='dnds')
  lm.mod2 = lm(dnds2~omega)
  ypred2  = predict(lm.mod2,list(omega))
  lines(omega,ypred2,col="#661100",lwd=2)
  
  
  plot(dnds3~omega,xlim=c(0,0.6),ylim=c(0,0.6),main='sample',ylab='dnds')
  lm.mod3 = lm(dnds3~omega)
  ypred3  = predict(lm.mod3,list(omega))
  lines(omega,ypred3,col="#661100",lwd=2)
  dev.off()
  
  # sm1 = summary(lm.mod1)
  # sm2 = summary(lm.mod2)
  # sm3 = summary(lm.mod3)
  # sm1$r.squared
  # sm2$r.squared
  # sm3$r.squared
  
  slp1 = lm.mod1$coefficients[2]
  slp2 = lm.mod2$coefficients[2]
  slp3 = lm.mod3$coefficients[2]
  cat("the linear slope of sw,max,sample is: \n",slp1,slp2,slp3)
  
  pdf(ouFig)
  ##try with ggplot2
  df  = data.frame(omega=omega,sw=dnds1,max=dnds2,sample=dnds3)
  ddf = df %>% pivot_longer(cols=sw:max:sample, names_to="Methods", values_to="dnds")
  ddf = ddf %>% filter(is.finite(dnds))
  
  gg = ggplot(ddf,aes(x=omega,y=dnds,shape=Methods)) + geom_point(aes(colour=Methods)) + theme_classic() + xlim(0,0.6) + ylim(0,0.6)
  GG = gg + geom_smooth(method=lm,se=F,fullrange=T,aes(color=Methods)) + 
       stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~"),color=Methods),label.x=0.1,method='pearson')
  
  GG = GG + geom_abline(slope=1,show.legend=T,linetype='dashed') 
  
  print(GG)
  dev.off()
}





####################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3])
