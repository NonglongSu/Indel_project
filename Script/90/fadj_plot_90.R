#Filter alignments via the fadj/fdiff per-gene-pair.
#from table

suppressWarnings(suppressMessages(library(readr)))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(car))        #logit
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hexbin))
suppressPackageStartupMessages(library(mixtools))

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

#plot the index hex plot [adjust the fadj-boundary]
plot.gap.hex = function(dat,y,pen){
  gg <- ggplot(dat, aes(x=fdiff,y=fadj)) + geom_hex()
  gg <- gg + scale_x_log10() + scale_y_log10()
  gg <- gg + coord_fixed()
  gg <- gg + xlab(expression(f[diff])) + ylab(expression(f[adj]))
  gg <- gg + geom_hline(yintercept=pen,color="#4daf4a",lwd=1) #fadj
  gg <- gg + geom_vline(xintercept=0.1,color="#4daf4a",lwd=1)
  gg <- gg+ggtitle(y)
  print(gg)
}



###############################################
# inD    = "Raw_data/JCdis_sum"
# ouF    = "Raw_data/QC2/"
# ouFig  = "Figure/QC2/id.hex.max.pdf"
# pat    = "max"
# fadj   = '0.5'

main = function(inD,ouF,ouFig,pat,fadj){
  
  Files  = list.files(inD,full.names=TRUE)
  N      = length(Files)
  pen    = as.numeric(fadj)
  
  namev = c()
  for(i in 1:N){
    namev[i] = gsub('\\..*','',basename(Files[i]))
  }
  
  #fadj/fdiff hex 
  pdf(ouFig,onefile=T,pointsize=18,paper='a4',width=8, height=15)
  for (i in 1:N){
    dat = read.table(Files[i],header=T)
    dat <- dat %>% mutate(totlen  = lenA + lenB,
                          diff    = abs(lenA-lenB),
                          f       = (gapA_len+gapB_len)/totlen,
                          fdiff   = diff/totlen,
                          fadj    = f-fdiff)
    plot.gap.hex(dat,namev[i],pen)
    
    ##output the bad genes. 
    bad.id  = which((dat$fdiff>=0.1) | (dat$fadj>=pen))
    bad.id  = unique(bad.id)
    rm.file = dat$locus[bad.id]
    write.table(rm.file,paste0(ouF,namev[i],'.',pat,'.txt'),row.names=F,col.names=F,quote=F)
  }
  dev.off()
  
  
}
#####################################
args=commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4],args[5])



