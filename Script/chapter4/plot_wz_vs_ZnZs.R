#Usage: Rscript --vanilla ../Script/chapter4/plot_wz_vs_ZnZs.R 1 k12


suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))

######################################################
#setwd("~/Dropbox (ASU)/Indel_project/chapter4")

#num=1
#targ="k12"
main = function(num, targ){
  
  n = as.numeric(num)
  
  #read true input
  inF  = "trueP.100.txt"
  tP   = read.table(inF,header=T)
  tMat = as.matrix(tP)
  nvar <<- ncol(tMat)
  
  omegaV  = tMat[,11]
  omegazV = tMat[,18]
  
  #read simulated input
  inD1 = paste0(targ,"/Results/Gse")
  inD2 = paste0(targ,"/Results/Pise")
  #ouD1 = paste0(targ,"/Figure/500_sim")
  #ouD2 = paste0(targ,"/Figure/500_im")
  
  if(n==1){
    
  }else{
    
  }
  
  plot(omegaV,omegazV)
}