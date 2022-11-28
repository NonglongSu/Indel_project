suppressWarnings(suppressMessages(library(tidyr)))
suppressPackageStartupMessages(library(stringr))


setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

inF = "MLE90.tab.txt"
ouF = "MLE90.tab.csv"

dat = read.table(inF,header=T,sep="\t",quote=NULL)

spnew  = c()
for (i in 1:90) {
  spnew[i] = str_remove(dat$species[i],"_aligned_cds")
  spnew[i] = str_remove(spnew[i],"_")
}
rownames(dat) = spnew


dat1 = dat[,5:14]
write.csv(dat1,ouF,row.names=T,quote=F)
