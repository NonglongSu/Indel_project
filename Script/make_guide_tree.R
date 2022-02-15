#Create a hash table corresponding species name to the common name.
#Make guide tree for prank alignment.

library(ape)     # read.tree()
library(stringr)
library(readr)

# setwd("~/Dropbox (ASU)/Indel_project/Script")

#Using vim to add the species instead of text editor.
# file = "../Input/homo3.1.txt"
# file = "../Input/homo3.2.txt"
# file = "../Input/homo4.1.txt"

#Create a hash table
hash_table = function(){
        spec.map = list("hamster" = "cricetulus_griseus_chok1gs_hdv1",
                        "mouse" = "mus_musculus_reference_cl57bl6_strain",
                        "rat" = "rattus_norvegicus",
                        "human" = "homo_sapiens",
                        "macaque" = "macaca_mulatta"
                        )
        return(spec.map)
}


main = function(file){
        
        name = read_delim(file, delim = "\t", col_names = FALSE)
        pat  = name[, 2][[1]]
        
        #Url changes every year.
        url = "ftp://ftp.ensembl.org/pub/current_compara/species_trees/49_mammals_EPO_default.nh"
        a           = read.tree(url)
        a$tip.label = str_to_lower(a$tip.label) 
        
        # keep it lower case.
        # keep = c("homo_sapiens", "mus_musculus_reference_cl57bl6_strain", "rattus_norvegicus")
        # keep = c("homo_sapiens", "mus_musculus_strain_reference_cl57bl6", "rattus_norvegicus", "cricetulus_griseus")
          keep = unlist(lapply(pat, function(x){grep(x, a$tip.label, value = TRUE)}))
        
          drop   = !(a$tip.label %in% keep)
          b      = drop.tip(a, a$tip.label[drop])
                
        # b$tip.label = c("hamster", "mouse", "rat", "human")  >>>>>>>
        # b$tip.label = c("mouse", "rat", "human")     
          hash.spec   = hash_table()
          comm.name   = names(hash.spec) 
          new.order   = unlist(lapply(b$tip.lab, function(x){which(hash.spec %in% x)}))
          b$tip.label = comm.name[new.order]
                  
        #write out
        cat(write.tree(b))
}

args = commandArgs(trailingOnly = TRUE)
main(args[1])




 
 
 
 
