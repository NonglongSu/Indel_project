#compare the phase proportions vs the indel rate (phase-importance-sampling)

suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggpubr))

#setwd("~/Dropbox (ASU)/Indel_project/chapter3")


# inD1 = "90/Results/PISE"
# inD2 = "../test_90_species/Results/phases"
# ouFig="90/Figure/r_vs_prop.pdf"

main = function(inD1,inD2,ouFig){
  Files1 = list.files(inD1,full.names=T,pattern='est')
  Files2 = list.files(inD2,full.names=T)
  Files2 = Files2[c(3,1,2)]
  n      = length(Files1)
  
  #indel-rate
  r012 = matrix(0,n,3)
  for(i in 1:n){
    j.tab   = fromJSON(Files1[i])
    g       = j.tab$gap.openning
    g012    = c(g[1]+g[4],g[2]+g[5],g[3]+g[6])
    r012[i,]= -log(1-g012)
  }
  norm.r = r012/rowSums(r012)
  
  #gap-proportion
  gap.lst = list()
  for(j in 1:3){
    gap.phase = read.table(Files2[j],header=T)  
    gap.phase = gap.phase[,-1]
    gap.lst[[j]] = gap.phase       
  }
  
  
  # gmat = matrix(0,90,6)
  # for(i in 1:n){
  #   j.tab   = fromJSON(Files1[i])
  #   g       = j.tab$gap.openning
  #   gmat[i,]    = c(g[1],g[2],g[3],g[4],g[5],g[6])
  # }
  # colnames(gmat)=c('i0','i1','i2','d0','d1','d2')
  # write.csv(gmat,"indel_rates.csv",row.names=F,quote=F)
  
  
  #############################################################plot
  par(mfrow=c(2,1))
  par(mar=c(5,3,2,2)+0.1)
  
  rb = rainbow(3)
  
  gap.phase=gap.lst[[1]]
  
  p0 = gap.phase[,1]
  p1 = gap.phase[,2]
  p2 = gap.phase[,3]
  
  r0 = norm.r[,1]
  r1 = norm.r[,2]
  r2 = norm.r[,3]
  
  hist(p0,col=rb[1],xlim=c(0,1),ylim=c(0,20),freq=F,xlab="phase-proportions",main=NULL)
  hist(p1,add=T,col=rb[2],freq=F)
  hist(p2,add=T,col=rb[3],freq=F)
  lines(sort(p0),dnorm(sort(p0),mean(p0),sd(p0)),lwd=2)
  lines(sort(p1),dnorm(sort(p1),mean(p1),sd(p1)),lwd=2)
  lines(sort(p2),dnorm(sort(p2),mean(p2),sd(p2)),lwd=2)
  
  legend("topright",legend=c('phase0','phase1','phase2'),col=rb,fill=rb, cex=0.8)
  
  hist(r0,col=rb[1],xlim=c(0,1),ylim=c(0,20),freq=F,xlab="normlized-indel-rate",main=NULL)
  hist(r1,add=T,col=rb[2],freq=F)
  hist(r2,add=T,col=rb[3],freq=F)
  lines(sort(r0),dnorm(sort(r0),mean(r0),sd(r0)),lwd=2)
  lines(sort(r1),dnorm(sort(r1),mean(r1),sd(r1)),lwd=2)
  lines(sort(r2),dnorm(sort(r2),mean(r2),sd(r2)),lwd=2)
  
  dev.off()
  
  #linear model  [not working]  
  # lm1   = lm(gap.phase[,1]~norm.r[,1])
  # lm2   = lm(gap.phase[,2]~norm.r[,2])
  # lm3   = lm(gap.phase[,3]~norm.r[,3])
  # 
  # yhat1 = predict(lm1,list(norm.r[,1]))
  # yhat2 = predict(lm2,list(norm.r[,2]))
  # yhat3 = predict(lm3,list(norm.r[,3]))
  # 
  # par(mfrow=c(1,3))
  # plot(gap.phase[,1]~norm.r[,1])
  # lines(norm.r[,1],yhat1,col="#661100",lwd=2)
  # plot(gap.phase[,2]~norm.r[,2])
  # lines(norm.r[,2],yhat2,col="#661100",lwd=2)
  # plot(gap.phase[,3]~norm.r[,3])
  # lines(norm.r[,3],yhat3,col="#661100",lwd=2)
  # 
  # dev.off()
  # 
  # sm1 = summary(lm1)
  # sm2 = summary(lm2)
  # sm3 = summary(lm3)
  # 
  # sm1$r.squared       ##two low
  # sm2$r.squared
  # sm3$r.squared
   
  #grouped boxplot
  pdf(ouFig,onefile=T)
  ##put same phase together
  p0 = c(unlist(lapply(gap.lst,'[[',1)),norm.r[,1])
  p1 = c(unlist(lapply(gap.lst,'[[',2)),norm.r[,2])
  p2 = c(unlist(lapply(gap.lst,'[[',3)),norm.r[,3])
  
  Phases     = rep(c("phase0","phase1","phase2"),each=360)
  Methods    = rep(c("mafft+sw","coati-max","coati-sampling","mle.r"),each=90)
  Ratio      = c(p0,p1,p2)
  df         = data.frame(Phases,Methods,Ratio)
  df$Methods = factor(Methods,levels=unique(Methods))
  g          = ggplot(df, aes(x=Phases,y=Ratio,color=Methods)) + geom_boxplot() + xlab("") +
               geom_point(size=1,position=position_jitterdodge(jitter.width=0.1)) + 
               theme(axis.text.x = element_text(colour=c("red","green","blue")))
    
  #print(g)  
  
  #create a table
  #generate a mean/var table
  mean_each = function(x){
    x.chunk = split(x,ceiling(seq_along(x)/90))
    avg     = unlist(lapply(x.chunk, function(x) round(mean(x),4)))
    return(avg)
  }
  
  tb           = rbind(mean_each(p0),mean_each(p1),mean_each(p2))
  colnames(tb) = c("mafft+sw","coati-max","coati-sampling","mle.r")
  rownames(tb) = c("phase0","phase1","phase2")
  cols       = matrix("black", nrow(tb), ncol(tb))
  col.head   = c("red","green","cyan","purple")
  mytab      = tableGrob(tb,theme=ttheme_default(base_size=7,core=list(fg_params = list(col=cols),
                                                                       bg_params = list(col=NA)),
                                                                       rowhead=list(bg_params = list(col=NA)),
                                                                       colhead=list(bg_params = list(col=col.head))))     
  
  gg = g + annotation_custom(mytab, xmin=0,xmax=Inf,ymin=0.7,ymax=0.8)
  print(gg)
  dev.off()
  
  
}
  



##########################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3])
