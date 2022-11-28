
#setwd("~/Dropbox (ASU)/Indel_project/test_90_species")

inD   = "Results/phase_sample"
ouF   = "Results/phases/phase_coati_sample.txt"

Files = list.files(inD,full.names=T)
n     = length(Files)
spec.name = unlist(lapply(Files,function(x){basename(x)}))

dmat = matrix(0,90,3)
for(i in 1:n){
  tmp     = read.table(Files[i],header=F)
  dmat[i,]= as.numeric(tmp)
}

df = data.frame(dmat)
df = cbind(spec.name,df)
colnames(df) =  c("spec.name","phase-0","phase-1","phase-2")

write.table(df,ouF,quote=F,row.names=F,sep='\t')
#capture.output(print(df, print.gap=3),file=ouF)
