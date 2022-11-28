#make supplementary table for chapter 1
#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

datgen = function(inD,pat){
  Files     = list.files(inD, pattern=pat, full.names=T)        
  if(!pat==""){
    pre.order = as.numeric(gsub('.*[\\/]([^.]+)[.].*','\\1',Files))
    Files     = Files[order(pre.order)]
  }
  
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
  colnames(dMat)=c('s1','s2','s3','s4','s5','s6','omega','tau')
  dMat = dMat %>% round(3)
  dMat = apply(dMat, 2, function(x){sprintf("%.3f",x)})
  dMat
}


############################################
inF    = "Results/truePar_100.txt"

inDA   = "Results/est100"
pat1   = '5e'
pat2   = '6e'
pat3   = '7e'

inDB   = "Results/nmkb100.1"
inDC   = "Results/est90"
pat0   = ""

#true par
trueP = read.table(inF, header=T, sep="")
sst   = trueP[,5:10]/trueP[,11] 
tMat  = as.matrix(cbind(sst,trueP[,12:11]))
tMat  = tMat %>% round(3)
tMat  = apply(tMat, 2, function(x){sprintf("%.3f",x)})

#est par
dMat1 = datgen(inDA,pat1)
dMat2 = datgen(inDA,pat2)
dMat3 = datgen(inDA,pat3)

nmkb.mat = datgen(inDB,pat0)

dMat90   = datgen(inDC,pat0)
Files    = list.files(inDC,full.names=F)    
spnew  = c()
for (i in 1:90) {
  spnew[i] = str_remove(Files[i],".json")
  spnew[i] = str_remove(spnew[i],"_")
}
rownames(dMat90) = spnew



###output files
ouF      = "Results/supp1/true_par_tab.csv"
ouF1     = "Results/supp1/est_5e_tab.csv"
ouF2     = "Results/supp1/est_6e_tab.csv"
ouF3     = "Results/supp1/est_7e_tab.csv"
ouFnmkb  = "Results/supp1/est_nmkb_tab.csv"
ouF90    = "Results/supp1/est_90_tab.csv"


write.csv(tMat,ouF,row.names=F,quote=F)
write.csv(dMat1,ouF1,row.names=F,quote=F)
write.csv(dMat2,ouF2,row.names=F,quote=F)
write.csv(dMat3,ouF3,row.names=F,quote=F)
write.csv(nmkb.mat,ouFnmkb,row.names=F,quote=F)
write.csv(dMat90,ouF90,row.names=T,quote=F)






