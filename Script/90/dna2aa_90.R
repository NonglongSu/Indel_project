#Translate DNA to Amino acid
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(stringr))

#setwd("~/Dropbox (ASU)/Indel_project/Script/90")


# dir  = "../../test_90_species/Raw_data/cds/01_FcaCaf_aligned_cds"
# ouD = "../../test_90_species/Raw_data/aa/01_FcaCaf_aligned_cds"

dna_to_aa = function(f,ouD){
 for(i in f){
   cds = readDNAStringSet(i,format="fasta")
   aa  = Biostrings::translate(cds, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
   
   writeXStringSet(aa, paste0(ouD,'/',basename(i)), append=F,
                   compress=F, compression_level=NA, format="fasta")
 }  
}

main = function(dir, ouD){
  Files = list.files(dir, full.names=TRUE)
  print(basename(dir))
  dna_to_aa(Files,ouD)
}

##########################################
args = commandArgs(trailingOnly=TRUE)
main(args[1], args[2])

