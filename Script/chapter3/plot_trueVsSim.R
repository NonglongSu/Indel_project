#Usage: Rscript --vanilla ../Script/chapter3/plot_trueVsSim.R 1 0


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


#without wz
doMat17 = function(Files1,Files2,n){
  dMat1 = matrix(0,n,nvar)
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
  return(list(dMat1,dMat2))
}

#with wz
doMat18 = function(Files1,Files2,n){
  dMat1 = matrix(0,n,nvar)
  colnames(dMat1)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','omega','tau','r0*t','r1*t','r2*t','ext.I','ext.D','omegaZ')
  for (i in 1:n) {
    Jmp       = fromJSON(Files1[i])
    go        = Jmp$gap.openning
    go.comb   = c(go[1]+go[4],go[2]+go[5],go[3]+go[6])
    rt        = -log(1-go.comb)
    ge        = Jmp$gap.extension
    dMat1[i,] = c(Jmp$nuc.freq, Jmp$sigmas, Jmp$omega, Jmp$branch.length, rt, ge, Jmp$omega.z) 
  }
  
  #summary statistics
  dMat2 = matrix(0,n,8)
  colnames(dMat2)=c('num.zn0','num.zn1','num.zn2','num.zs0','num.zs1','num.zs2','avg.gap.size.I','avg.gap.size.D')
  for (i in 1:n) {
    Jmp2      = fromJSON(Files2[i])
    numN      = Jmp2$zn.phases
    num.zn    = c(numN[1]+numN[4],numN[2]+numN[5],numN[3]+numN[6])
    numS      = Jmp2$zs.phases
    num.zs    = c(numS[1]+numS[4],numS[2]+numS[5],numS[3]+numS[6])
    dMat2[i,] = c(num.zn,num.zs,Jmp2$avg.gap.size)  
  }
  return(list(dMat1,dMat2))
}

####################################
readMatrix = function(inD){
  Files1 = readFile(inD,'est')[[1]]
  Files2 = readFile(inD,'sum')[[1]]
  n      = readFile(inD,'est')[[2]]
  #parameter estimate
 
  if(nvar==17){
    res.mat = doMat17(Files1,Files2,n)
  }else{
    res.mat = doMat18(Files1,Files2,n)
  }
  
  return(res.mat)
}


qqplt = function(tMat,dMat,p,q){
  colorblindP = "#CC6677"
  for(j in p:q){
    qqplot(tMat[,j],dMat[,j],xlab=colnames(tMat)[j],ylab=paste0(colnames(tMat)[j],'_est'),main=NULL)
    abline(0,1,col=colorblindP,lwd=2)
  }
}

hiplt = function(diff,name,p,q){
  colorblindP = "#CC6677"
  for (j in p:q){
    hist(diff[,j],xlab=name[j],main=NULL,freq=FALSE,ylim=c(0,40))
    lines(density(diff[,j]),lty='dotted',col=colorblindP,lwd=2)
  }
}


##true est vs est plot
plot1 = function(tMat,dMat,ouD){
  ouF1=paste0(ouD,"/","qqplot.pdf")
  ouF2=paste0(ouD,"/","error_perc.pdf")
  
  
  #qqplot
  pdf(ouF1,onefile=TRUE)
  
  ##sigmas:par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=rep(4,4))
  par(mfrow=c(2,3))
  qqplt(tMat,dMat,5,10)
  
  ##omega,tau
  par(mar=c(5.1,4.1,4.1,2.1))
  par(mfrow=c(1,2))
  qqplt(tMat,dMat,11,12)
  
  ##indel phase rate
  par(mfrow=c(1,3))
  qqplt(tMat,dMat,13,15)
  
  ##gap extensions
  par(mfrow=c(1,2))
  qqplt(tMat,dMat,16,17)
  
  #omegaZ
  if(nvar==18){
    par(mfrow=c(1,1))
    qqplt(tMat,dMat,18,18)
  }
  
  
  dev.off()
  
  #Error percentage distribution plot
  name = colnames(tMat)
  diff = abs(dMat-tMat)/tMat
  
  pdf(ouF2,onefile=TRUE)
  par(mar=rep(4,4))
  #sigmas
  par(mfrow=c(2,3))
  hiplt(diff,name,5,10)
  
  #omega & tau
  par(mfrow=c(1,2))
  hiplt(diff,name,11,12)

  #indel phase rate
  par(mfrow=c(1,3))
  hiplt(diff,name,13,15)
  
  #gap exts
  par(mfrow=c(1,2))
  hiplt(diff,name,16,17)
 
  if(nvar==18){
    par(mfrow=c(1,1))
    hiplt(diff,name,18,18)
  }
  
  dev.off()
  
  # summary(diff[,1:17])
  # var(diff[,1:17])     
}

##true align vs est plot
plot2 = function(dMat1.2,dMat2.2){
  pdf("Figure/phase.sum.pdf",onefile=TRUE)
  #summary stat: num of gap phases, avergage gap sizes
  
  #colorblindP = "#CC6677"
  # numg1=cbind(dMat1.2[,1]+dMat1.2[,4],dMat1.2[,2]+dMat1.2[,5],dMat1.2[,3]+dMat1.2[,6])
  # numg2=cbind(dMat2.2[,1]+dMat2.2[,4],dMat2.2[,2]+dMat2.2[,5],dMat2.2[,3]+dMat2.2[,6])
  # colnames(numg1)=c('num.phase0','num.phase1','num.phase2')
  # 
  # par(mfrow=c(1,3))
  # for(j in 1:3){
  #   qqplot(numg1[,j],numg2[,j],xlab=colnames(numg1)[j])
  #   abline(0,1,col=colorblindP,lwd=2)
  # }
  
  ##qqplot
  par(mfrow=c(1,3))
  qqplt(dMat1.2,dMat2.2,1,3)
  
  par(mfrow=c(1,2))
  qqplt(dMat1.2,dMat2.2,4,5)
  
  
  ##error perc
  name = colnames(dMat1.2)
  diff = abs(dMat1.2-dMat2.2)/dMat1.2
  
  par(mfrow=c(1,3))
  hiplt(diff,name,1,3)
  
  par(mfrow=c(1,2))
  hiplt(diff,name,4,5)
  
  dev.off()
  
  
  
  # for (j in 1:3){
  #   hist(diff[,j],prob=T,xlab=name[j],ylim=c(0,25),main=NULL)
  #   lines(density(diff[,j]),lty='dotted',col=colorblindP,lwd=2)
  # }
  #dev.off()
}

############################################################
#setwd("~/Dropbox (ASU)/Indel_project/chapter3")
#setwd("~/Dropbox (ASU)/Indel_project/chapter4")



#num1=1;num2=0
#num1=0;num2=1
#num1=1;num2=1

main = function(num1,num2){
  n1=as.numeric(num1)
  n2=as.numeric(num2)
  
  #read true input
  inF  = "trueP.100.txt"
  tP   = read.table(inF,header=T)
  tMat = as.matrix(tP)
  nvar <<- ncol(tMat)
  
  #read simulated input
  inD1 = "Results/Gse"
  inD2 = "Results/Pise"
  ouD1 = "Figure/500_sim"
  ouD2 = "Figure/500_im"
  
  
  
  
  
  if(n1==1){
    dmp1    = readMatrix(inD1)
    dMat1.1 = dmp1[[1]]
    dMat1.2 = dmp1[[2]]
    if(n2==0){
      plot1(tMat,dMat1.1,ouD1)
    }else{
      dmp2    = readMatrix(inD2)
      dMat2.1 = dmp2[[1]]
      dMat2.2 = dmp2[[2]]
      plot2(dMat1.2,dMat2.2)
    }
  }else{
    dmp2    = readMatrix(inD2)
    dMat2.1 = dmp2[[1]]
    dMat2.2 = dmp2[[2]]
    plot1(tMat,dMat2.1,ouD2)
  }
 
}



########################################
args = commandArgs(trailingOnly=TRUE)
main(args[1],args[2])


