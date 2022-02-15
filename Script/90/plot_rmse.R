library(utils)
library(ggplot2)
library(dplyr)
#setwd("~/Dropbox (ASU)/Indel_project/Script/90")


#inD = "../../test_90_species/Results"
#######################
main = function(inD,ouF){
  
  rmse.file = list.files(inD,pattern="rmse",full.names=T)
  rmse.data = lapply(rmse.file, function(x){read.table(x,sep='\t')})
  
  rmse.m = cbind(rmse.data[[1]],rmse.data[[2]],rmse.data[[3]])
  t.rmse = t(rmse.m)
  rmse   = as.vector(t.rmse)  #y Axis
  
  simulations = as.numeric(rep(seq(1,100),each=3))  #x Axis
  sample.size = rep(c('5e','6e','7e'),times=100)    #group, one shape per group
  data        =  data.frame(simulations, rmse, sample.size)
  
  # stacked area chart
  pdf(ouF)
  g = ggplot(data, aes(x=simulations, y=rmse, fill=sample.size)) + 
    geom_area()
  print(g)
  dev.off()
  
}





###################
args = commandArgs(trailingOnly=T)
main(args[1],args[2])