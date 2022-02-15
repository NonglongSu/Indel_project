library(Biostrings)
library(stringr)
library(seqinr)

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

#inDir = "../Raw_data.2.outgroups/mapped_cds_mafft_Tmp"

main = function(inDir,ouFile){
	File.holder = list.files(inDir,full=TRUE)
	nameList = c()

	for(i in 1:length(File.holder)){
		dna = readDNAStringSet(File.holder[i],format="fasta")
		len = width(dna)
		if(length(unique(len)) != 1){
			nameList = c(nameList,basename(File.holder[i]))
		}else{
			next
		}
	}

	cat(nameList, file=ouFile, sep='\n', append=FALSE)

}

args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2]) 
