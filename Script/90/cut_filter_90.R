#measure how much loss of data from after concatenation.

library(Biostrings)
library(jsonlite)
library(stringr)
library(ggplot2)
library(ggrepel)
library(R.utils)
#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

countL = function(f){
  if(isDirectory(f)){
    ff = list.files(f,full.names=T)
    l  = 0
    for (j in 1:length(ff)) {
      Seq  = readBStringSet(ff[j])
      l    = l+sum(width(Seq))
      # Seqs = str_split(Seq,'')
      # g    = lapply(Seqs, function(x){IRanges(x=='-')})
      # g    = IRangesList(g)
      # g    = unique(unlist(g))
      # g.rem= width(Seq)[1]-sum(width(g))
      # l    = l+g.rem
    }
  }else{
    Seq = readBStringSet(f)
    l   = sum(width(Seq))
  }
  return(l)
}


# inD1 = "../../test_90_species/Raw_data/concat"
# inD2 = "../../test_90_species/Raw_data/cds"
# inD3 = "../../test_90_species/Results/est90"
#########################
main = function(inD1,inD2,inD3,ouF){
  File1 = list.files(inD1, full.names=T)
  File2 = list.files(inD2, full.names=T)

  seqM  = matrix(0,90,2)
  sname = str_extract(basename(File1), "[^.]+")
  tname = basename(File2)
  for (i in 1:90) {
    testi = which(grepl(sname[i],tname,ignore.case=T))
    if(length(testi)!=0){
      seqM[i,1] = countL(File1[i])
      seqM[i,2] = countL(File2[testi])
      print(seqM[i,])
    }else{
      print(i)
      next
    }
  }

#perctage of loss
fperc = (seqM[,2]-seqM[,1])/seqM[,2]

# a = str_extract(basename(File1), "[^.]+")
# b = str_extract(basename(File1), "[^aligned]+")
# setdiff(a,b)

#>>reading tau
File3 = list.files(inD3, full.names=T)
tv    = c()
for(i in 1:length(File3)){
  Jtmp    = fromJSON(File3[i])
  tv[i]   = Jtmp$p.intermed[nrow(Jtmp$p.intermed),10]
}

#>>plot tau vs fperc
#plot(fperc~tv)

pdf(ouF)
df  = data.frame(tau=tv, perc=fperc)
lof = ggplot(df, aes(x=tau, y=perc)) + geom_point(alpha=0.6, col="#56B4E9") +
  geom_text_repel(label=rownames(df), size=2.5)
print(lof)
dev.off()



}
##########
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4])
