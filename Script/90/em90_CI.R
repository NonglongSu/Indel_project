#Bootstrapping of 90 ests.

library(jsonlite)
library(stringr)
library(ggplot2)
library(ggrepel)

#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

CI = function(pm){
  avg    = colMeans(pm)
  se     = apply(pm, 2, sd)/sqrt(10)
  t.crit = qt(0.975,df=nrow(pm)-1)
  CI.low  = avg - t.crit*se
  CI.high = avg + t.crit*se
  return(rbind(CI.low,CI.high))
}

Var = function(pm){
  avg    = colMeans(pm)
  se     = apply(pm, 2, sd)/sqrt(10)
  return(se)
}



# inD   = "../../test_90_species/Results/boot90"
# inF   = "../../test_90_species/Results/specList.txt"
# ouD   = "../../test_90_species/Results/CI90"
############################################
main = function(inD, inF, ouD, ouF){
  #fake name
  #Name = "Results/CI90/01_FcaCaf.json"
  #name = str_extract(basename(Name), "[^.]+")

  #read file
  spec = read.table(inF,header=F,sep="\n")
  n    = nrow(spec)
  pm90 = array(0, c(10,8,n))
  var90= matrix(0,90,8)
  for (i in 1:n) {
    bootlst = list.files(inD, pattern=spec[i,], full.names=TRUE)
    pm      = matrix(0,10,8)
    for (j in 1:length(bootlst)) {
      Jboot   = fromJSON(bootlst[j])
      pm[j,]  = c(Jboot$par,Jboot$p.intermed[nrow(Jboot$p.intermed),10])   
    }
    pm90[,,i] = pm
    #Confidence interval
    p.CI = CI(pm)
    #variance
    var90[i,] = Var(pm)
    #write.table(p.CI, paste0(ouD,'/',spec[i,],'.txt'),
    #           row.names=F, col.names=F)
  }
  
  
  #>>
  Cm = colMeans(var90)
  range(Cm)
  Rm = rowMeans(var90)
  Rm[which(Rm == max(Rm))]
  
  #>>draw w vs tau
  wv = c()
  tv = c()
  for (k in 1:n) {
    wv = c(wv,pm90[,7,k])
    tv = c(tv,pm90[,8,k])
  }
  
  #>>
  #plot(wv~tv, xlab='tau', ylab='omega')
  
  #>>
  pdf(ouF)
  
  df     = data.frame(omega=wv,tau=tv, num1=1:900, num2=rep(1:900, each=10, length.out=900))
  df$col = rep(1:6, each=10, length.out=900) 
  g = ggplot(df, aes(x=tau, y=omega, col=as.factor(col))) 
  plot0 = g + geom_point(size=2, alpha=0.5, shape=1, show.legend=F) + 
              geom_text_repel(aes(label=ifelse(num1%%10==0, num2, '')), col="#000000", size=2.5)
    
  print(plot0)
  dev.off()
  
  
}

# j=1
# dr = data.frame()
# while(j<=l){
#   min.max = c(range(wv[j:(j+9)]), range(tv[j:(j+9)]))
#   dr      = rbind(dr,min.max)
#   j=j+10
# }
# colnames(dr)=c('wmin','wmax','tmin','tmax')

############################
args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2], args[3], args[4])





