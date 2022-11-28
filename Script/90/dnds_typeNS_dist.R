suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

#setwd("~/Dropbox (ASU)/Indel_project")

readF = function(Files){
  pv = c()
  for(i in 1:length(Files)){
    tmp=read.table(Files[i])
    pv[i]=tmp[1,1]
  }
  return(pv)
}



#inD1 = "test_90_species/Results/dnds"
#inD2 = "test_90_species/Results/typeNS"
#pat  = "max"
main = function(inD1,inD2,ouFig,pat){
  File1 = list.files(inD1,full.names=T,pattern=pat)
  File2 = list.files(inD2,full.names=T,pattern=pat)
  
  dnds_90 = readF(File1)
  tyNS_90 = readF(File2)
  
  df      = matrix(c(dnds_90,tyNS_90),90,2)
  
  cat(sprintf("dnds-mean:%.3f typeNS-mean:%.3f \n dnds-var:%.3f typeNS-var:%.3f",
              mean(dnds_90),mean(tyNS_90),var(dnds_90),var(tyNS_90)),"\n")

  # df.avg = data.frame('dnds'=mean(dnds_90), 'typeN/S'=mean(tyNS_90))
  # write.table(df.avg,ouF)
  
  #draw
  pdf(ouFig)    
  colnames(df) = c('dnds','typeN/typeS')
  Pdf = as.data.frame(as.table(df))
  #ggplot(Pdf, aes(x=Var2, y=Freq)) + geom_boxplot() + xlab("") + ylab("Ratio")
  
  Pdf_sorted = Pdf %>% mutate(para=fct_reorder(Var2,-Freq))
  res        = ggplot(Pdf_sorted, aes(x=Var2,y=Freq,color=para)) + geom_boxplot() + coord_flip() + xlab("") + ylab("Ratio") + 
               stat_summary(fun=mean,geom='point',size=3) 
  print(res)
  
  dev.off()
  
}





##################################
args=commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4])