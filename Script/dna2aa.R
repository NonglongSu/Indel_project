#Translate DNA to Amino acid

#setwd("~/Dropbox (ASU)/Indel_project/Script")
library(Biostrings)
library(scales)

# file_name = "../test_human_mouse_rat/Raw_data/cds_seq/ENSG00000250722.fa"
# output    = "../test_human_mouse_rat/Raw_data/aa_seq/"

dna_to_aa = function(file_name, output, name.spec){
 
 cds     = readDNAStringSet(file_name, format = "fasta")
 protein = Biostrings::translate(cds, genetic.code = GENETIC_CODE, if.fuzzy.codon = "solve")

 names(protein) = name.spec
 writeXStringSet(protein, paste0(output, basename(file_name)), append = FALSE,
                 compress = FALSE, compression_level = NA, format = "fasta")
}

args = commandArgs(trailingOnly = TRUE)
dna_to_aa(args[1], args[2], eval(parse(text = args[3])))

