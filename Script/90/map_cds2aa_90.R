#Abnormal alignment: align amino acid to DNA.
library(seqinr)

#setwd("~/Dropbox (ASU)/Indel_project/Script/90")

map_cds2aa = function(file1, file2, subD){

for (n in 1:length(file1)) {
  ouFile = paste0(subD, basename(file1[n]))
  if(file.exists(ouFile)){
    next()
  }else{
    print(basename(file1[n]))
    data1 = read.fasta(file1[n], seqtype = "AA", as.string = FALSE,
                       set.attributes = FALSE)
    data2 = read.fasta(file2[n], seqtype = "DNA", as.string = FALSE,
                       set.attributes = FALSE, forceDNAtolower = FALSE)
    for(i in 1:length(data1)){
      if(length(which(data1[[i]] == "-")) != 0){
        pos = which(data1[[i]] == "-")
        for(j in 1: length(pos)){
          data2[[i]] = append(data2[[i]], c("-", "-", "-"), after = (pos[j] - 1) * 3)
        }
      }
    }
    name = names(data1)
    write.fasta(sequences = data2,names = name, nbchar = 80,
                open = "w", as.string = FALSE, file.out = ouFile)
  }
 }
}


# dir1  = "../../test_90_species/Raw_data/align"
# dir2  = "../../test_90_species/Raw_data/cds"
# ouD   = "../../test_90_species/Raw_data/mapped_cds/"

main = function(dir1, dir2, ouD){
  Dirs1     = list.files(dir1, full.names = TRUE)
  Dirs2     = list.files(dir2, full.names = TRUE)
  for (i in 1:length(Dirs1)) {
    subdir  = paste0(ouD, basename(Dirs1[i]), "/")
    dir.create(subdir)
    Files1  = list.files(Dirs1[i], full.names = TRUE)
    Files2  = list.files(Dirs2[i], full.names = TRUE)
    print(basename(Dirs1[i]))
    map_cds2aa(Files1, Files2, subdir)
  }
}


###########################################
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3])
