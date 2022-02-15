#Abnormal alignment: align amino acid to DNA.

library(seqinr)

map_cds2aa = function(file1, file2, file3){

#file1 = file.path("..","Raw_data",file1)
#file2 = file.path("..","Raw_data",file2)

#file1 = file.path("..","Raw_data_2",file1)
#file2 = file.path("..","Raw_data_2",file2)

# file1 = "../test_human_mouse_rat/Raw_data_2/Alignments/prank_out/ENSG00000155657.fa"
# file2 = "../test_human_mouse_rat/Raw_data_2/cds_seq/ENSG00000155657.fa"

data1 = read.fasta(file1, seqtype="AA", as.string = FALSE)
data2 = read.fasta(file2, seqtype="DNA", as.string = FALSE,forceDNAtolower = FALSE)

for(i in 1:length(data1)){
  if(length(which(data1[[i]]=="-")) != 0){
    pos = which(data1[[i]]=="-")
    for(j in 1: length(pos)){
      data2[[i]] = append(data2[[i]], c("-","-","-"), after=(pos[j]-1)*3)
    }
  }
}



name = names(data1)


#file3 = "../Raw_data/mapped_cds_seq/"
write.fasta(sequences = data2,names = name, nbchar=80,
            open = "w", as.string = FALSE, file.out = paste0(file3, basename(file1)))
}

args = commandArgs(trailingOnly = TRUE)
map_cds2aa(args[1], args[2], args[3])
