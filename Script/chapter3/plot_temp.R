suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))

readFile = function(inD, pat){
  Files= list.files(inD,pattern=pat,full.names=T)
  index = as.numeric(str_extract(basename(Files),'[^.]+')) #rectify the order
  Files = Files[order(index)]
  n    = length(Files)
  res  = list(Files,n)
  return(res)
}

readMatrix = function(inD){
  Files1 = readFile(inD,'est')[[1]]
  Files2 = readFile(inD,'sum')[[1]]
  n      = readFile(inD,'est')[[2]]
  #parameter estimate
  dMat1 = matrix(0,n,17)
  colnames(dMat1)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','omega','tau','r0*t','r1*t','r2*t','ext.I','ext.D')
  for (i in 1:n) {
    Jmp      = fromJSON(Files1[i])
    go       = Jmp$gap.openning
    go.comb  = c(go[1]+go[4],go[2]+go[5],go[3]+go[6])
    rt       = -log(1-go.comb)
    ge       = Jmp$gap.extension
    dMat1[i,] = c(Jmp$nuc.freq, Jmp$sigmas, Jmp$omega, Jmp$branch.length, rt, ge)  
  }
  
  #summary statistics
  dMat2 = matrix(0,n,5)
  colnames(dMat2)=c('num.phase0','num.phase1','num.phase2','avg.gap.size.I','avg.gap.size.D')
  for (i in 1:n) {
    Jmp2      = fromJSON(Files2[i])
    numP      = Jmp2$gap.phases
    num.phase = c(numP[1]+numP[4],numP[2]+numP[5],numP[3]+numP[6])
    dMat2[i,] = c(num.phase,Jmp2$avg.gap.size)  
  }
  
  return(list(dMat1, dMat2))
}

setwd("~/Dropbox (ASU)/Indel_project/chapter3")
setwd("~/Dropbox (ASU)/Indel_project/chapter4")

#read true input
inF  = "trueP.100.txt"
tP   = read.table(inF,header=T)
tMat = as.matrix(tP)


#read simulated input
inD1 = "Results/Gse"
inD2 = "Results/Pise"

dmp     = readMatrix(inD1)
dMat    = dmp[[1]]
dMat1.2 = dmp[[2]]
dmp     = readMatrix(inD2)
dMat    = dmp[[1]]
dMat2.2 = dmp[[2]]


#qqplot
##sigmas
#par(mar=c(5.1,4.1,4.1,2.1))
par(mar=rep(4,4))
par(mfrow=c(2,3))
colorblindP = "#CC6677"
for(j in 5:10){
  qqplot(tMat[,j],dMat[,j],xlab=colnames(tMat)[j],ylab=paste0(colnames(tMat)[j],'_est'),main=NULL)
  abline(0,1,col=colorblindP,lwd=2)
}
dev.off()

##omega,tau
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(1,2))
for(j in 11:12){
  qqplot(tMat[,j],dMat[,j],xlab=colnames(tMat)[j],ylab=paste0(colnames(tMat)[j],'_est'))
  abline(0,1,col=colorblindP,lwd=2)
}
dev.off()

##indel rates
par(mfrow=c(1,3))
for(j in 13:15){
  qqplot(tMat[,j],dMat[,j],xlab=colnames(tMat)[j],ylab=paste0(colnames(tMat)[j],'_est'))
  abline(0,1,col=colorblindP,lwd=2)
}
dev.off()

##gap extensions
par(mfrow=c(1,2))
for(j in 16:17){
  qqplot(tMat[,j],dMat[,j],xlab=colnames(tMat)[j],ylab=paste0(colnames(tMat)[j],'_est'))
  abline(0,1,col=colorblindP,lwd=2)
}
dev.off()


##########################################
#Error percentage distribution plot
name = colnames(tMat)
diff = abs(dMat-tMat)/tMat

gseq=13:17
par(mar=rep(4,4))
par(mfrow=c(2,3))

for (j in gseq){
  hist(diff[,j],prob=T,xlab=name[j],ylim=c(0,25),main=NULL)
  lines(density(diff[,j]),lty='dotted',col=colorblindP,lwd=2)
}
dev.off()

summary(diff[,gseq])
var(diff[,gseq])     #read diag

###########################################
#summary stat
##num of gap phases
colorblindP = "#CC6677"
par(mfrow=c(2,3))
for(j in 1:5){
  qqplot(numg1[,j],numg2[,j],xlab=colnames(numg1)[j])
  abline(0,1,col=colorblindP,lwd=2)
}


name = colnames(dMat1.2)
diff = abs(dMat1.2-dMat2.2)/dMat1.2
for (j in 1:5){
  hist(diff[,j],prob=T,xlab=name[j],ylim=c(0,25),main=NULL)
  lines(density(diff[,j]),lty='dotted',col=colorblindP,lwd=2)
}
dev.off()







