#Map aa back to DNA.
suppressPackageStartupMessages(library(seqinr))

#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

map_cds2aa = function(file1,file2,ouD){

for (k in 1:length(file1)) {
  ouF = paste0(ouD,'/',basename(file1[k]))
  if(file.exists(ouF)){
    next()
  }else{
    data1 = read.fasta(file1[k], seqtype="AA", as.string=FALSE,
                       set.attributes=FALSE)
    data2 = read.fasta(file2[k], seqtype="DNA", as.string=FALSE,
                       set.attributes=FALSE, forceDNAtolower=FALSE)
    for(i in 1:length(data1)){
      if(length(which(data1[[i]]=="-")) != 0){
        pos = which(data1[[i]]=="-")
        for(j in 1: length(pos)){
          data2[[i]] = append(data2[[i]], c("-","-","-"), after = (pos[j]-1)*3)
        }
      }
    }
    name = names(data1)
    write.fasta(sequences=data2, names=name, nbchar=100,
                open="w", as.string=FALSE, file.out=ouF)
    print(k)
  }
}
  
}


####################################################
# dir1  = "../../test_90_species/Raw_data/align_mafft/07_yeast_aligned_cds"
# dir2  = "../../test_90_species/Raw_data/cds/07_yeast_aligned_cds"
# ouD   = "../../test_90_species/Raw_data/mapped_cds/07_yeast_aligned_cds"

main = function(dir1,dir2,ouD){
  file1     = list.files(dir1, full.names=TRUE)
  file2     = list.files(dir2, full.names=TRUE)
  print(basename(dir1))
  map_cds2aa(file1,file2,ouD)
}


#######################################
args = commandArgs(trailingOnly=TRUE)
main(args[1],args[2],args[3])
