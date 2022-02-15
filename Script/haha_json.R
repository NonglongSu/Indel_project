library(rjson)


res    = fromJSON(file = "../../Simul/Input_json/ENSG00000003249.json")

score  = res$dist[[1]]$aln[[1]]$score
cigar  = res$dist[[1]]$aln[[1]]$cigar
seq    = res$dist[[1]]$seq

