#Update the Indel events with multiple optimal alignment (---/+++)

library(Biostrings)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(seqinr)
library(stringr)
library(stringi)

#m1 = "AAT===AAACAAAGAATGCTTACTGT+++ATAAGGCTTACTGTTCTAGCG++++++ATCACCGCG===TCATGTCTAGTTATGAACGGC"
#m2 = "AAT===AAACAAAGAATGCTTACTGT---ATAAGGCTTACTGTTCTAGCG------ATCACCGCG===TCATGTCTAGTTATGAACGGC"

# inFile2="~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Data/test/tmp/test1.fa"
# inFile1="~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Data/test/tmp1/test1.fa"

main = function(inFile1,inFile2,ouDir){
  
  #ouDir = "../Data/mapped_cds_plus/"
  data1 = readBStringSet(inFile1,format = "fasta")
  data2 = readBStringSet(inFile2,format = "fasta")
  
  #From mapped_cds_Out/
  dna.1 = str_split(as.character(data1), '')
  g = lapply(dna.1, function(x) { IRanges(x == '+')})
  g = IRangesList(g)
  
  m1 = g[[1]]
  r1 = g[[2]]
  
  wid.m = width(m1)
  wid.r = width(r1)
  
  l.m = length(wid.m)
  l.r = length(wid.r)
  
  pos.m = start(m1)
  pos.r = start(r1)
  
  en.m = end(m1)
  en.r = end(r1)
  
  #mapped_cds_In_Out/
  m2 = toString(data2[[1]],width=NULL)
  r2 = toString(data2[[2]],width=NULL)
  name = names(data2)
  
  if( l.m == 0 & l.r == 0 ){
    m2 = m2
    r2 = r2
  }else if(l.m>0 & l.r == 0){
    for(i in 1:l.m){
      subseq(m2,start=pos.m[i], end=en.m[i]) = stri_rand_strings(1,wid.m[i],'[+]')
    }
  }else if(l.m == 0 & l.r>0 ){
    for(i in 1:l.r){
      subseq(r2,start=pos.r[i], end=en.r[i]) = stri_rand_strings(1,wid.r[i],'[+]')
    }
  }else{
    for(i in 1:l.m){
      subseq(m2,start=pos.m[i], end=en.m[i]) = stri_rand_strings(1,wid.m[i],'[+]')
    }
    for(i in 1:l.r){
      subseq(r2,start=pos.r[i], end=en.r[i]) = stri_rand_strings(1,wid.r[i],'[+]')
    }
  }
  
  new_seq = list(m2,r2)
  write.fasta(sequences = new_seq,names = name, nbchar=80,
              open = "w",as.string = TRUE, file.out =paste0(ouDir,basename(inFile2))  )
  
}
#######################################
args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])

