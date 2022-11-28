#logistic regression between p(Y=hk|x) vs 
#x can be dnds, ZnZs. 
#Data_6 data

suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(nnet))              #multinom
suppressPackageStartupMessages(library(kknn))              #kknn
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(pROC))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")

traceId = function(x){
  dna = readDNAStringSet(x,format="fasta")
  n   = names(dna)
  gsub("\\..*","",n[1])
}

#cal. rolling stats(moving average) default:k=101
ma = function(x) rollmean(x, k=101)


####################################################

# inF1  = "housekeeping_genes_human.csv"           #https://housekeeping.unicamp.br/
# inF2  = "Results/dnds_znzs_sum.tsv"
# inD1  = "Raw_data/cds_seq"
# inD2  = "Data_6/Mafft/mapped_cds"
# ouFig = "Results/hk_logit.pdf"

main = function(inF1,inF2,inD1,inD2,ouFig){
  
  #read human hk trans ID
  dat  = read.table(inF1,header=T,sep=";")
  hkId = dat[,1]               #2833
  
  #backtrack hk geneID
  File1   = list.files(inD1,pattern='.fa',full.names=T)
  file1   = list.files(inD1,pattern='.fa',full.names=F)
  genes   = sapply(file1, function(x) gsub("\\..*","",x))
  genes   = unname(genes)
  
  trId = sapply(File1, function(x){traceId(x)})
  trId = unname(trId)          #13827
  
  comm   = intersect(hkId,trId)
  id     = which(trId %in% comm)
  geneid = genes[id]           #1282
  
  #%usage
  length(geneid)/length(hkId)  #45.3%
  
  #get sw geneID
  files= list.files(inD2,pattern='.fa',full.names=F)
  swId = sapply(files, function(x) gsub("\\..*","",x))
  swId = unname(swId)          #1521
  
  #interset the swID with hkID
  hk_gene = intersect(geneid,swId)
  
  hkid    = which(swId %in% hk_gene)
  ts_gene = swId[-hkid]
  
  #%hk-gene in sw-gene
  length(hkid)/length(swId)  #5.1%
  
  
  ########################################################test logit regression 
  #lift curve
  source('../../machine-learning-666/sources/robfuns.R')
  source('../../machine-learning-666/sources/rob-utility-funs.R')
  
  #manage data
  HK       = rep(0,length(swId))
  HK[hkid] = 1
  
  dat1    = read_tsv(inF2,show_col_types=F)
  dat1    = cbind(HK,dat1)
  dat1$HK = as.factor(dat1$HK)
  dat1 = dat1 %>% mutate(Pn= Nd/N, Ps=Sd/S)
  dat1 = dat1 %>% mutate(dn= -0.75*log(1-4*Pn/3), 
                         ds= -0.75*log(1-4*Ps/3),
                         w = dn/ds,
                         p = Zn/(Zn+Zs))

  dat2    = dat1 %>% filter(is.finite(w)) %>% arrange(w)              #sort dnds
  dat3    = dat1 %>% filter(is.finite(p)) %>% arrange(p)              #sort p
  
  ##train-split test
  set.seed(99)
  ntrain = round(2/3 * nrow(dat1))
  tr     = sample(1:nrow(dat1), ntrain)
  train  = dat1[tr,]
  test   = dat1[-tr,]
  
  
  #moving average
  tab2  = tibble(w  = ma(dat2$w),
                 HK = ma(as.numeric(dat2$HK))-1)
  
  tab3  = tibble(p  = ma(dat3$p),
                 HK = ma(as.numeric(dat3$HK))-1)
                 
  
  pdf(ouFig,onefile=T)
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>dnds
  m1    = glm(HK~w,family="binomial",data=train)
  ypred = predict(m1,test,type="response")
  
  M1    = glm(HK~w,family="binomial",data=dat2)
  yfit  = predict(M1,tab2,type="response")
  
  ##logit curve
  par(mfrow=c(1,2))
  plot(ypred~test$w,xlab="dnds",ylab="P(gene=HK|x)",col='green', main="train-test-split")
  
  plot(tab2$HK~tab2$w, main="moving average")
  points(tab2$w,yfit,col='blue')
  
  ##fit multinom [not bad!]
  mlogf   = multinom(HK~w,data=train)
  logfit  = predict(mlogf,test)
  print(table(logfit,test$HK))
  mul_tab = table(logfit,test$HK)
  mis1    = (sum(mul_tab)-sum(diag(mul_tab)))/sum(mul_tab)    #misclassfication rate: 0.051
  
  ##Lift curve
  ii    = (1:length(ypred))/length(ypred)
  olift = liftf(test$HK,ypred,dopl=F)
  plot(ii,olift,type='l',lwd=2,xlab='% tried',ylab='% of successes',cex.lab=1,ylim=c(0,1),col='red')
  abline(0,1,lty=2)
  text(0.2,1,paste0("misclassification rate:", round(mis1,3)))
  
  ##ROC curve
  rocR = roc(response=test$HK,predictor=ypred)
  auc1 = auc(rocR)             #0.693
  plot(rocR, col="red")
  text(0.8,1,paste0("auc:", round(auc1,3)))
  mtext("Classfication metrics of dnds vs P(HK|x)", side = 3, line=-2, outer = T, col='red')
  
  
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ZnZs
  m2     = glm(HK ~ p, family=binomial(), data=train)
  ypred2 = predict(m2,test,type="response")
  
  M2     = glm(HK~p,family="binomial",data=dat3)
  yfit2  = predict(M2,tab3,type="response")
  
  ##logit curve
  par(mfrow=c(1,2))
  plot(ypred2~test$p,xlab="Zn/(Zn+Zs)",ylab="P(gene=HK|x)",col='green',main="train-test-split")
  
  plot(tab3$HK~tab3$p)
  points(tab3$p,yfit2,col='blue',main="moving average")
  
  
  ##fit multinom [not bad!]
  mlogf2   = multinom(HK~p,data=train)
  logfit2  = predict(mlogf2,test)
  print(table(logfit2,test$HK))
  mul_tab2 = table(logfit2,test$HK)
  mis2     = (sum(mul_tab2)-sum(diag(mul_tab2)))/sum(mul_tab2) #misclassfication rate 0.052
  
  ##lift curve [bad!]
  ii2    = (1:length(ypred2))/length(ypred2)
  olift2 = liftf(test$HK,ypred2,dopl=F)
  plot(ii2,olift2,type='l',lwd=2,xlab='% tried',ylab='% of successes',cex.lab=1,ylim=c(0,1),col='red')
  abline(0,1,lty=2)
  text(0.2,1,paste0("misclassification rate:", round(mis2,3)))
  
  ##ROC curve [bad!]
  rocR2 = roc(response=test$HK,predictor=ypred2)
  auc2  = auc(rocR2)                 #0.569                         
  plot(rocR2, col="red")
  text(0.8,1,paste0("auc:", round(auc2,3)))
  
  mtext("Classfication metrics of Zn/(Zn+Zs) vs P(HK|x)", side = 3, line=-2, outer = T, col='red')
  dev.off()
  
  
}

###########################################
args = commandArgs(trailingOnly=T)
main(args[1], args[2], args[3], args[4], args[5])

