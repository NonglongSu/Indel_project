#make supplementary table for chapter 3
#setwd("~/Dropbox (ASU)/Indel_project/chapter3")

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))


datgen = function(inD,pat){
  Files = list.files(inD,pattern=pat,full.names=T)
  index = as.numeric(str_extract(basename(Files),'[^.]+')) #rectify the order
  Files = Files[order(index)]
  
  dMat = matrix(0,100,nvar)
  colnames(dMat)=c('s1','s2','s3','s4','s5','s6','omega','tau','ext.I','ext.D','r0*t', 'r1*t','r2*t')
  for (i in 1:100) {
    Jmp      = fromJSON(Files[i])
    go       = Jmp$gap.openning
    go.comb  = c(go[1]+go[4],go[2]+go[5],go[3]+go[6])
    rt       = -log(1-go.comb)
    ge       = Jmp$gap.extension
    dMat[i,] = c(Jmp$sigmas, Jmp$omega, Jmp$branch.length, ge, rt) 
  }
  dMat = dMat %>% round(3)
  dMat = apply(dMat, 2, function(x){sprintf("%.3f",x)})
  dMat
}


############################################
inF    = "trueP.100.txt"
inD1   = "Results/Gse"
inD2   = "Results/Pise"

#true par
tP   = read.table(inF,header=T)
tMat = as.matrix(tP)
tmat = tMat[,-1:-4] %>% round(3)
tMat = apply(tmat, 2, function(x){sprintf("%.3f",x)})
nvar <<- ncol(tmat)

#est par
dMat1 = datgen(inD1,'est')
dMat2 = datgen(inD2,'est')


###output files
ouF      = "supp/true_par.csv"
ouF1     = "supp/sim_par.csv"
ouF2     = "supp/im_par.csv"


write.csv(tMat,ouF,row.names=F,quote=F)
write.csv(dMat1,ouF1,row.names=F,quote=F)
write.csv(dMat2,ouF2,row.names=F,quote=F)







