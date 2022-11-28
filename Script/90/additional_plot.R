suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(plyr))    ##ddply

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#cal. rolling stats(moving average) default:k=201
ma    = function(x) rollmean(x, k=201)

#extract the species name
naming = function(x) sub("\\..*", "", basename(x))

#cal the GC % 
gc_perc = function(x){
  dat    = read_tsv(x,show_col_types=F)
  datsum = unname(colSums(dat))
  gc.perc= (datsum[1]+datsum[2])/(datsum[3]+datsum[4])
  return(gc.perc)
}

#read the omega
omega_read =function(x){
 w.val = as.vector(unlist(read.table(x,header=F)))
 return(w.val)
}

#dnds
dnds_est = function(x){
  Nd = sum(x$Nd)
  Sd = sum(x$Sd)
  N  = sum(x$N)
  S  = sum(x$S)
  
  PN = (Nd/N) 
  PS = (Sd/S) 
  
  #jukes cantor formula (1986)
  dNdS = log(1-4*PN/3)/log(1-4*PS/3)
  return(dNdS)
}

#generate dnds,ZnZs matrix
dnds_gen = function(x){
  dnds = c()
  for(i in 1:length(x)){
    spec.sum      = read_delim(x[i],col_names=T,show_col_types=F)
    dnds[i]   = dnds_est(spec.sum) 
  }
  dnds
}



##################################################################

# inF  = "Results/ZnZs/coef_max.tsv"
# inD  = "Results/ZD_sum"
# pat   = "max"

main = function(inF,ouFig,pat){
  Files = list.files(inD,full.names=T,pattern=pat)
  dat   = read.table(inF,header=T,sep='\t')
  namev = dat$species  
  n     = nrow(dat)
  
 
  #########################check order
  dnds.max = dnds_gen(Files)
  dat1     = cbind(dat[,1],dnds.max,dat[,2:9])
  eudat    = dat1[1:15,]
  order(eudat$dnds.max,decreasing=T)
  order(eudat$neutral,decreasing=T)
  order(eudat$GC,decreasing=T)
  order(eudat$yintercept,decreasing=T)
  order(eudat$init_slope,decreasing=T)
  order(eudat$ymax,decreasing=T)
  order(eudat$corr,decreasing=T)
  order(eudat$rmse,decreasing=T)
  order(eudat$neutra-eudat$ymax,decreasing=T)
  write.csv(eudat,"Results/supp4/coef_max_euk.csv",row.names=T,col.names=T)
  
  ##################################################nls-coeffficent
  ##slope0
  slope0 = dat$init_slope
  
  #par(mfrow=c(2,1))
  #hist(slope0,freq=F,breaks=90, xlab=("initial slope"),main="")
  #lines(density(slope0),col='#661100',lwd=2)
  
  pdf("Figure/ZnZs/init_slope.pdf", onefile=T)
  dat2 = data.frame("species"=1:90,"init_slope"=slope0)
  dat2$species = as.factor(dat2$species)
  gg = ggplot(dat2,aes(x=species,y=init_slope,fill=species)) + geom_bar(stat="identity") + 
       theme(axis.text.x=element_blank(), legend.position="none")
  print(gg)
  dev.off()
  
  ##check neg species
  which(slope0<0)
  
  #norm.slope0 = scale(slope0)     #(coef-mu)/sd
  #qqnorm(norm.slope0,main=NULL,xlim=c(-3,6),ylim=c(-3,6))
  #abline(0,1,col="#CC6677",lwd=2)
  #dev.off()
  shapiro.test(norm.slope0)    #failed 
  print(summary(slope0))
  
  
  ####################yintercept
  y0 = dat$yintercept
  
  par(mfrow=c(1,2))
  hist(y0,freq=F,breaks=40)
  lines(density(y0),col='#661100',lwd=2)
  
  norm.y0 = scale(y0)     #(coef-mu)/sd
  qqnorm(norm.y0,main=NULL,xlim=c(-4,4),ylim=c(-4,4))
  abline(0,1,col="#CC6677",lwd=2)
  dev.off()
  
  shapiro.test(norm.y0)    #failed 
  print(summary(y0))
  
  ###################asymptote
  ymax = dat$ymax
  
  par(mfrow=c(1,2))
  hist(ymax,freq=F,breaks=40)
  lines(density(ymax),col='#661100',lwd=2)
  
  norm.ymax = scale(ymax)     #(coef-mu)/sd
  qqnorm(norm.ymax,main=NULL,xlim=c(-4,4),ylim=c(-4,4))
  abline(0,1,col="#CC6677",lwd=2)
  dev.off()
  
  shapiro.test(norm.ymax)    #failed 
  print(summary(ymax))
  
  ###########################################asymptote vs neutral
  neuZ = dat$neutral
  ymax = dat$ymax
  
  dat3    = data.frame(neuZ,ymax)
  Species = rep(1:90,each=2)
  Stat    = rep(c('neutral_ZnZs','asympote'),90)
  Ratio   = c(t(dat3))
  df      = data.frame(Species,Stat,Ratio)
  
  
  pdf("Figure/ZnZs/neuZ_vs_ymax.pdf", onefile=T)
  #hist
  mu  = ddply(df, "Stat", summarise, grp.mean=mean(Ratio), grp.sd=sd(Ratio))
  gg1 = ggplot(df,aes(x=Ratio,color=Stat,fill=Stat)) + 
        geom_histogram(position="identity",alpha=0.6,bins=90) + geom_density(alpha=0.6) + 
        geom_vline(data=mu, aes(xintercept=grp.mean, color=Stat),linetype="dashed") +
        labs(x=" ", y = "Density") + theme_classic()
  
  #print(gg1)
  
  ##stack 
  gg2 = ggplot(df, aes(x=Species, y=Ratio, fill=Stat)) +  geom_line(aes(color=Stat)) + geom_point(aes(color=Stat)) +
    geom_vline(xintercept=seq(1,90,by=1),color='gray',size=.5,alpha=.5)
  #print(gg2) 
  
  gg3  = ggarrange(gg1, gg2, labels=c("A","B"), ncol=1,nrow=2)
  print(gg3)
  dev.off()
  
  ##check how many neuz < ymax
  which(neuZ<ymax)
   
   
   ########################################neutral vs GC%
   #neutral Zn/(Zn+Zs)
   gc   = dat$GC
   
   hist(neuZ,freq=F,breaks=40)
   hist(gc,freq=F,breaks=40)
   
   #neutral vs GC% across 90 species
   plot(neuZ~gc,ylab='neutral-ZnZs')
   text(neuZ~gc,labels=1:90,col="#009E73",cex=0.5,pos=1)
   
   df = data.frame(neu=neuZ,gc=gc)
   df = df %>% arrange(gc)
   
   glm.fit= glm(neu~poly(gc,2),data=df)
   lines(sort(gc),glm.fit$fitted.values,col='#661100',lwd=2)
   
   ##how much changed
   (max(neuZ)-min(neuZ))/max(neuZ)   #0.167
   (max(gc)-min(gc))/max(gc)         #0.616
   
   ##check perc of euk Id
   specId = order(gc)
   eukId  = specId[which(specId<=15)]
   
   quat   = summary(df$neu)[3]         #median
   topId  = specId[which(df$neu>quat)]
   
   intersect(topId,eukId)              #10/15
   
   #try observed ZnZs [noise]
   obsZ = rep(0,90)
   for(i in 1:N){
     datmp    = read_tsv(Files[i],show_col_types=F)
     cs     = colSums(datmp)
     obsZ[i]= cs[1]/(cs[1]+cs[2])
   }
   plot(obsZ~gc,ylab='observed-ZnZs')
   text(obsZ~gc,labels=1:90,col="#661100",cex=0.5,pos=1)
   
   
   ################################################omega vs seq-len (noise)
   wv   = rep(0,n)
   lenv = rep(0,n)
   for(i in 1:n){
     datmp  = read_tsv(Files[i],show_col_types=F)
     cs     = colSums(datmp)
     Nd = cs[3]; N = cs[4]; Sd = cs[5]; S = cs[6]
     Pn = Nd/N; Ps = Sd/S
     dn = -0.75*log(1-4*Pn/3)
     ds = -0.75*log(1-4*Ps/3)
     wv[i] = dn/ds
     
     lenA = cs[7]; lenB = cs[8]
     lenv[i] = (lenA+lenB)/2
   }
   
   hist(wv,breaks=40,freq=F,xlab='omega',main=NULL)  
   lines(density(wv),col="#661100",lwd=2)
   
   ##seq_length vs omega [noise]
   plot(wv~lenv,xlab='seq.len', ylab='omega') 
   
   #reverse CDF of omega ???
   
   
   #################################################sample-size vs rmse (noise)
   rmse = dat$rmse
   samp = dat$sample_size
   
   plot(rmse~sample_size,dat)
   
   hist(rmse,breaks=40,freq=F)  
   srmse = sort(rmse)
   lines(srmse,dnorm(srmse,mean(rmse),sd(rmse)),col='#661100',lwd=2)
   
   shapiro.test(rmse)  
   
   
   ################################################correlation hist
   corr = dat$corr
   
   hist(corr,breaks=40,freq=F)  
   lines(density(corr),col='#661100',lwd=2)
   
   length(corr[abs(corr)>0.7])/n     #41% high correlation (>0.7)
   length(corr[abs(corr)<0.5])/n     #38% low 
   
   
   ##corr vs rmse (failed)
   plot(corr~rmse)
   
   
   
   
   
}


##############################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4],args[5],args[6])
