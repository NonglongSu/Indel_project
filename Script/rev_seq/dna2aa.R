#Translate DNA to Amino acid

#setwd("~/Dropbox (ASU)/Indel_project/Script")
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(scales))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Rev_seq")

# file_name = "Raw_data/ENSG00000250722.rcds.fa"
# ouF       = "Raw_data/ENSG00000250722.aa.fa"

dna_to_aa = function(file_name, ouF, name.spec){
 
 cds     = readDNAStringSet(file_name, format = "fasta")
 protein = Biostrings::translate(cds, genetic.code = GENETIC_CODE, if.fuzzy.codon = "solve")

 names(protein) = name.spec
 writeXStringSet(protein, ouF, append=FALSE,
                 compress=FALSE, compression_level=NA, format = "fasta")
}

#########################################
args = commandArgs(trailingOnly=TRUE)
dna_to_aa(args[1], args[2], eval(parse(text=args[3])))

