#ensemble the mle est of chapter 3 as 1-line file

suppressPackageStartupMessages(library(jsonlite))

#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

cal_e = function(e3){
  x = 1/(1-e3)
  1 - 1/(x*3)
}



inD   = "../chapter3/90/Results/PISE"
Files = list.files(inD, full.names=TRUE, pattern='est')
n     = length(Files)

dmat = matrix(0,n,14)
colnames(dmat)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','omega','tau','g','e')

for(i in 1:n){
  Jtmp = fromJSON(Files[i])
  f    = Jtmp$nuc.freq
  s    = Jtmp$sigmas
  w    = Jtmp$omega 
  t    = Jtmp$branch.length
  g    = mean(Jtmp$gap.openning)*2
  e3   = mean(Jtmp$gap.extension)
  e    = cal_e(e3)
  dmat[i,] = c(f, s, w, t, g, e)
}

spec.name = read.table("Species.txt",header=F)[[1]]
ouD = "MLE/"
for(i in 1:n){
  cmd = paste0(c('-p',dmat[i,1:4],'-x',dmat[i,5:10],'-w',dmat[i,11],'-t',dmat[i,12],'-g',dmat[i,13],'-e',dmat[i,14]), sep='')
  write(cmd,paste0(ouD,spec.name[i],'.txt'),sep=' ',ncolumns=20)
}

#also generate a table
df = data.frame(dmat)
df['species'] = spec.name
write.table(df,"MLE90.tab.txt",row.names=F,sep='\t',quote=F)
