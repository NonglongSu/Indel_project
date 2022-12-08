#Phase-proportion of insertion and deletion process separately

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))

#generate a table of gap stat
phasing = function(file){
  dna = readBStringSet(file, format="fasta")
  l   = length(dna)
  
  dna.1 = str_split(as.character(dna), '')
  g1 = lapply(dna.1, function(x) {IRanges(x %in% c("-", "+"))})
  g2 = lapply(dna.1, function(x) {IRanges(x == "+")})
  g1 = IRangesList(g1)
  g2 = IRangesList(g2)
  
  mr     = g1[(l-1):l]
  wid.mr = unlist(width(mr))
  pos.mr = unlist(start(mr))
  
  #record "+++"
  flag  = rep(0, length(wid.mr))
  flag1 = unlist(start(g2[(l-1):l]))
  flag[which(pos.mr %in% flag1)] = rep(1, length(flag1))
  
  #record file
  fname = rep(basename(file), length(pos.mr))
  
  #record which species
  spec = c(rep('m', length(mr[[1]])), rep('r', length(mr[[2]])))
  
  data.frame("pos"=pos.mr, "wid"=wid.mr, "flag"=flag, "file" = fname, "spec"=spec)
}  

#max parsimont
indel_watch = function(file, psj) {
  dna     =  readBStringSet(file, format="fasta")
  dna.ref = dna[[1]]
  pos = psj$pos
  wid = psj$wid
  tag = psj$spec
  
  pat = substr(dna.ref, pos, pos+wid-1)
  ins = del = 0

 if(tag == 'm'){
   if(grepl('-', pat) == TRUE){
     ins = 1
   }else{
     del = 1
   }
 }else{
   if(grepl('-', pat) == TRUE){
     ins = 1
   }else{
     del = 1
   }
 }
  
 c(ins, del)
}

#record which phases
Record = function(PST, gapName){
  pos.3  = PST[PST$wid == 3, ]$pos
  pos.6  = PST[PST$wid == 6, ]$pos
  pos.9  = PST[PST$wid == 9, ]$pos
  pos.12 = PST[PST$wid == 12,]$pos
  
  pos.lst  = list(pos.3, pos.6, pos.9, pos.12)
  DF.phase = data.frame() 
  for (i in 1:length(pos.lst)) {
    rem      = pos.lst[[i]] %% 3
    phase.0  = length(which(rem == 1))
    phase.1  = length(which(rem == 2)) 
    phase.2  = length(which(rem == 0)) 
    df.phase = data.frame(Phase_0 = phase.0, Phase_1 = phase.1, Phase_2 = phase.2)
    DF.phase = rbind(DF.phase, df.phase)
  }
  DF.phase = as.table(noquote(t(DF.phase)))
  colnames(DF.phase) = gapName
  return(DF.phase)
}  

#generate a table of phase prop of ins/del separately 
work_one_window = function(dir1,dir2){
  Files = list.files(dir2, full.names=F)
  PST.1 = c()
  PST.2 = c()
  for(i in 1:length(Files)){
    file1 = paste0(dir1, "/", Files[i])
    file2 = paste0(dir2, "/", Files[i])
    
    pscore1 = phasing(file1)
    pscore2 = phasing(file2)
    PST.1   = rbind(PST.1, pscore1)
    PST.2   = rbind(PST.2, pscore2)
    print(i)
  }
  
  #filter out the "+++"
  PST.2f = PST.2[PST.2$flag==0, ]
  PST.1f = PST.1[which(PST.2$flag==0), ]  
  
  #catch the seq state of ref to seperate indels
  Ins = c()
  Del = c()
  for (j in 1:nrow(PST.1f)) {
    file3 = paste0(dir2, "/", PST.1f$file[[j]])
    indel = indel_watch(file3, PST.1f[j, ])
    Ins   = c(Ins, indel[1])
    Del   = c(Del, indel[2])
  }
  PST.tag = cbind(PST.2f, Ins, Del)
  
  PST.i   = PST.tag[PST.tag$Ins == 1, ]
  PST.d   = PST.tag[PST.tag$Del == 1, ]
  gap.len = c(wind1,wind2,wind3,wind4)
  Phase.i = Record(PST.i, gap.len)
  Phase.d = Record(PST.d, gap.len)
  cbind(Phase.i,Phase.d)
}

####
#bar plot of ins/del of each window
bar_plot = function(tab,wind){
  par(mfrow = c(1,2), oma=c(2,0,2,0))
  Phase.i = tab[,1:4]
  Phase.d = tab[,5:8]
  safe_colorblind_palette = c("#0072B2","#009E73","#AA4499")
  barplot(Phase.i, xlab="Insertion", ylab="Count", col=safe_colorblind_palette, ylim=c(0,1000)) 
  barplot(Phase.d, xlab="Deletion", ylab="Count", col=safe_colorblind_palette, ylim=c(0,1000)) 
  #mtext("The proportion of phased-indels across length of gaps", outer = TRUE, cex = 1.5)
  mtext(paste0("window size = ", wind), outer=T, cex=1.5, side=1)
  legend(2.5, 800, bg ="transparent", bty="n", legend = rownames(Phase.i), fill=safe_colorblind_palette, cex=0.75)
}

####
chi_test = function(x1,x2,n1,n2) {prop.test(x=c(x1,x2),n=c(n1,n2),correct=F)}

#non-parameteric chi-squared test for proportion
prop_test = function(tab){
  Itab = tab[,1:4]
  Dtab = tab[,5:8]
  n1   = sum(Itab)
  n2   = sum(Dtab)
  
  I.pha = rowSums(Itab)
  D.pha = rowSums(Dtab)
  
  p.stat = matrix(0,3,3)
  for(k in 1:3){
    res  = chi_test(I.pha[k],D.pha[k],n1,n2)
    prop = res$estimate
    pval = res$p.value
    p.stat[,k] = c(prop,pval)
    #print(res)
  }
  p.stat
  
}

############################################################
#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")

#ouFig = "Results/Phase.mafft.ins.del.pdf"
#ouF = "Results/Phase.prop.test.csv"
main = function(ouFig){
  
  da1  = "Data_3/Mafft/updated_cds"
  da2  = "Data_3/Mafft/mapped_cds"
  db1  = "Data_6/Mafft/updated_cds"
  db2  = "Data_6/Mafft/mapped_cds"
  dc1  = "Data_9/Mafft/updated_cds"
  dc2  = "Data_9/Mafft/mapped_cds"
  dd1  = "Data_12/Mafft/updated_cds"
  dd2  = "Data_12/Mafft/mapped_cds"
  
  wind1<<-3;wind2<<-6;wind3<<-9;wind4<<-12
  
  phase.tab1 = work_one_window(da1,da2)
  phase.tab2 = work_one_window(db1,db2)
  phase.tab3 = work_one_window(dc1,dc2)
  phase.tab4 = work_one_window(dd1,dd2)
  
  
  ##PART II
  pdf(ouFig,onefile=T)
  bar_plot(phase.tab1,wind1)
  bar_plot(phase.tab2,wind2)
  bar_plot(phase.tab3,wind3)
  bar_plot(phase.tab4,wind4)
  dev.off()
  
  ##PART III
  test1 = prop_test(phase.tab1)
  test2 = prop_test(phase.tab2)
  test3 = prop_test(phase.tab3)
  test4 = prop_test(phase.tab4)
  
  testAll = cbind(test1,test2,test3,test4)
  test.df = apply(testAll, 2, function(x){sprintf("%.3f",x)})
  rownames(test.df) = c('Insertion','Deletion','P.value')
  colnames(test.df) = rep(c('phase0','phase1','phase2'),4)
  #write.csv(test.df,ouF,row.names=T,quote=F)
}
  

#################################################
args = commandArgs(trailingOnly=TRUE)
main(args[1])
  

