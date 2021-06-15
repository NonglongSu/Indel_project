#Find the linearship between rms and sample size.{5e..7e}

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(jsonlite))

#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

#rmsd = fromJSON("../test_90_species/Results/rms/1.json")

#######################################
#inD = "../../test_90_species/Results/rms"
main = function(inD, ouF){
  Files = list.files(inD, pattern='.txt', full.names=TRUE)
  dmat  = matrix(0,100,2)
  for(i in 1:length(Files)){
    dtmp      = read.delim(Files[i],header=F,sep='\t')
    dmat[i,1] = log(i*10^5)
    dmat[i,2] = dtmp[[2]]
  }
  
  #linearship
  colnames(dmat) = c("size","rms")
  dMat           = as.data.frame(dmat)
  lmmod = lm(rms~size, data=dMat)
  lmmod$coefficients
  summary(lmmod)
  plot(rms~size,data=dMat )
  
  glmod = glm(rms~size,data=dMat,family='gaussian')
  glmod$coefficients
  summary(glmod)
}



