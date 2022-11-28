#Test the similarity zscore between homologs. 
#Test the normalized indel length between homologs. 

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(moments)) #agostino test
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hexbin))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")

#read files
readF = function(Files){
  sim       = c()
  id.len    = c()
  lenA =lenB= c()
  for(i in 1:length(Files)){
    dna   = readBStringSet(Files[i], format="fasta")
    if(length(width(dna))>2){#MSA
      lw  = length(wid)
      dna = dna[(lw-1):lw]
    }
    wid         = width(dna)
    specs       = str_split(dna,'')
    sim[i]      = sim_score(specs,wid)
    id.len[i]   = indel_len(specs)
    #seq.len[i]  = sum(wid)-id.len[i]
    
    lenA[i] = nchar(str_remove_all(dna[[1]],'-'))
    lenB[i] = nchar(str_remove_all(dna[[2]],'-'))
  }
  
  res = list(sim,id.len,lenA,lenB)
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
indel_len = function(specs){
  g  = IRangesList(lapply(specs, function(x){IRanges(x=='-')}))
  ug = unlist(g)
  if(length(ug)==0){
    return(0)
  }else{
    len.g = sum(width(ug))
  }
  return(len.g)
}

#filter out via 1.5*IQR rule
# iqr = function(x){
#   Q1   = quantile(x,0.25)
#   Q3   = quantile(x,0.75)
#   IQR  = Q3-Q1
#   res  = c(Q1-1.5*IQR,Q3+1.5*IQR) 
#   return(res)
# }




# inD   = "Raw_data/coati_align"
# ouF   = "Raw_data/QC/id.hex.txt"
# ouFig = "Raw_data/Figure/id.hex.pdf"
main = function(inD,ouF,ouFig){
  
  Files       = list.files(inD,full.names=TRUE)
  n           = length(Files)
  data        = readF(Files)
  sim.score   = data[[1]]
  id.len      = data[[2]]
  seqA        = data[[3]]
  seqB        = data[[4]]
  
  dat = data.frame("id.len"=id.len,"lenA"=seqA,"lenB"=seqB)
  
  #sim dist
  sim.logit   = logit(sim.score,percents=F,adjust=.025)
  sim.avg     = mean(sim.logit)
  sim.sd      = sd(sim.logit)
  sim.zscore  = (sim.logit-sim.avg)/sim.sd
  
  #test normality
  ago.test = agostino.test(sim.zscore) 
  print(ago.test)
  
  pdf(ouFig,onefile=T)
  #similarity check
  par(mfrow=c(1,2))
  hist(sim.zscore,prob=T,main=NULL)
  lines(density(sim.zscore),lty="dotted",col="#CC6677",lwd=2)
  qqnorm(sim.zscore)
  abline(0,1,col="#CC6677",lwd=2)
  
  
  ##testing code
  #dat <- readr::read_delim("data.txt")
  # dat1 <- dat %>% mutate(totlen = lenA+lenB,
  #                       diff   = abs(lenA-lenB),
  #                       f      = id.len/totlen,
  #                       fadj   = (id.len-diff)/totlen
  # )
  # gg <- ggplot(dat1, aes(x=diff,y=fadj)) + geom_point()
  # gg <- gg + scale_x_log10()
  # gg <- gg + xlab("Length Difference") + ylab("Fraction Gaps")
  # print(gg)
  
  #hex plot
  # dat2 <- dat %>% mutate(totlen = lenA+lenB,
  #                       diff   = abs(lenA-lenB),
  #                       fadj    = (id.len-diff)/totlen,
  #                       g       = diff/totlen
  # )
  # gg <- ggplot(dat2, aes(x=g,y=fadj)) + geom_hex()
  # gg <- gg + scale_x_log10() + scale_y_log10()
  # gg <- gg + xlab("Fraction Length Diff") + ylab("Fraction Gaps")
  # print(gg)
  
  
  dat3 <- dat %>% mutate(totlen = lenA + lenB,
                         diff    = abs(lenA-lenB),
                         f       = (id.len)/totlen,
                         fdiff   = diff/totlen,
                         fadj    = f-fdiff)
  
  
  gg <- ggplot(dat3, aes(x=fdiff,y=fadj)) + geom_hex()
  gg <- gg + scale_x_log10() + scale_y_log10()
  gg <- gg + coord_fixed()
  gg <- gg + xlab(expression(f[diff])) + ylab(expression(f[adj]))
  gg <- gg + geom_hline(yintercept=0.1,color="#4daf4a",lwd=1)
  gg <- gg + geom_vline(xintercept=0.1,color="#4daf4a",lwd=1)
  print(gg)
  
  dev.off()

  
  #output the illegal files. 
  bad.id  = which((dat3$fdiff>=0.1) | (dat3$fadj>=0.1))
  bad.id  = unique(bad.id)
  rm.file = str_remove(basename(Files[bad.id]),'\\.[^\\.]*$')
  write.table(rm.file,ouF,row.names=F,col.names=F,quote=F)
  
}
#####################################
args=commandArgs(trailingOnly=T)
main(args[1],args[2],args[3])





