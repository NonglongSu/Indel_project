#Remove the low similarity between pariwise homologs. 
#Remove the high indel counts between homologs. 
#Remove the high indel length between homologs. 

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(car))     #qqPlot
suppressPackageStartupMessages(library(moments)) #agostino test
suppressPackageStartupMessages(library(stats))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Raw_data")

#read files
readF = function(Files){
  sim      = c()
  id.count = c()
  id.len   = list()
  seq.len  = c()
  for(i in 1:length(Files)){
    dna   = readBStringSet(Files[i], format="fasta")
    wid   = width(dna)
    if(length(wid)>2){#MSA
      lw  = length(wid)
      dna = dna[(lw-1):lw]
    }
    specs       = str_split(dna,'')
    sim[i]      = sim_score(specs,wid)
    id.len      = indel_count(specs)
    seq.len[i]  = sum(width(dna))-id.len
  }
  
  res = list(sim,id.len,seq.len)
  return(res)
}

#cal. similiarity
sim_score = function(specs,wid){
  identical_pos = length(which(specs[[1]][specs[[1]] == specs[[2]]] != '-'))
  shared_gaps   = length(which(specs[[1]][specs[[1]] == specs[[2]]] == '-'))
  total_gaps    = length(which(specs[[1]] =='-')) + length(which(specs[[2]] =='-'))
  sim.score     = identical_pos/(wid[1]-(total_gaps-shared_gaps))
  return(sim.score)
}

#sum the indel length of each file
indel_count = function(specs){
  g  = IRangesList(lapply(specs, function(x){IRanges(x=='-')}))
  ug = unlist(g)
  if(length(ug)==0){
    return(list(0,0))
  }else{
    len.g = sum(width(ug))
  }
  return(len.g)
}

#filter out via 1.5*IQR rule
iqr = function(x){
  Q1   = quantile(x,0.25)
  Q3   = quantile(x,0.75)
  IQR  = Q3-Q1
  res  = c(Q1-1.5*IQR,Q3+1.5*IQR) 
  return(res)
}


# inD   = "coati_align"
# ouFig = "QC/sim_indel_dist.pdf"
# ouF2  = "QC/high_indel_length.txt"
main = function(inD,ouF1,ouF3,ouFig){
  
  Files       = list.files(inD,full.names=TRUE)
  n           = length(Files)
  data        = readF(Files)
  sim.score   = data[[1]]
  id.count    = data[[2]]
  id.len      = data[[3]]
  seq.len     = data[[4]]
  
  write.table(id.count,"haha.txt",col.names = F,row.names = F)
  write.table(seq.len,"seq_len.txt",col.names = F,row.names = F)
  
  ##sim dist
  sim.logit   = logit(sim.score,percents=F,adjust=0.025)
  sim.avg     = mean(sim.logit)
  sim.sd      = sd(sim.logit)
  sim.zscore  = (sim.logit-sim.avg)/sim.sd
  
  par(mfrow=c(1,2))
  
  ##similarity check
  hist(sim.zscore,prob=T,main=NULL)
  lines(density(sim.zscore),lty="dotted",col="#CC6677",lwd=2)
  qqPlot(sim.zscore,col.lines="#CC6677")
  
  ago.test = agostino.test(sim.zscore) #test normality
  print(ago.test)
  
  ##count check
  id.count.norm = id.count/(seq.len)
  hist(id.count.norm,breaks=40,prob=T,main=NULL)
  lines(density(id.count.norm),lty="dotted",col="#CC6677",lwd=2)
  
  #poisson 
  par(mfrow=c(1,3))
  set.seed(666)
  pseudo.len    = 1e+3
  obs.count     = id.count.norm*pseudo.len
  lambda        = mean(obs.count) 
  exp.po.count  = rpois(n,lambda)
  qqplot(exp.po.count,obs.count,xlim=c(0,30),ylim=c(0,30))
  abline(0,1,col="#CC6677",lwd=2)
  
  #binomial
  set.seed(666)
  exp.bi.count = rbinom(n,seq.len,prob=mean(id.count.norm))
  qqplot(exp.bi.count,id.count,xlim=c(0,150),ylim=c(0,150))
  abline(0,1,col="#CC6677",lwd=2)
  
  #negative binomial
  set.seed(666)
  exp.nb.len = rnbinom(n, size=id.count+1, prob=mean(id.count.norm)) + id.count
  qqplot(exp.nb.len,seq.len,xlim=c(0,80000),ylim=c(0,80000))
  abline(0,1,col="#CC6677",lwd=2)
  
  
  ##length check
  obs.len  = unlist(id.len)/3
  
  hist(obs.len,breaks=40,prob=T,main=NULL)
  lines(density(obs.len),lty="dotted",col="#CC6677",lwd=2)
  
  p1       = 1-1/(sum(obs.len)/sum(id.count))
  exp.len  = rgeom(n,prob=1-p1)
  qqplot(exp.len,obs.len,xlim=c(0,2000),ylim=c(0,2000))
  abline(0,1,col="#CC6677",lwd=2)
  
  
  dev.off()
  
  #record the low-quality files.
  # id.count1 = id.count[which(id.count!=0)]
  # id.len2   = id.len1[which(id.len1!=0)]
  
  iqr.sim  = iqr(sim.zscore)
  iqr.count= iqr(id.count)
  iqr.len  = iqr(id.len1)
  
  bad.id1 = which(sim.zscore<iqr.sim[1])
  bad.id2 = which(id.count>iqr.count[2])
  
  flag    = unlist(lapply(id.len, function(x){any(x>iqr.len[2])}))
  bad.id3 = which(flag)
  
  bad.ids = unique(c(bad.id1,bad.id2,bad.id3))
    
}



#####################################
args=commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4])

