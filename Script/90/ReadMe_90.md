## 90 species analysis
Dataset request(from zhengting) \rightarrow  
update_gap                      \rightarrow  
sw_gap                          \rightarrow  
phasing indels                  


##1. Quality control of the data
**.**  remove all gaps
**.**  convert '?' to 'N'
**.**  remove genes with "NNN" only
**.**  remove genes that length of not multiple of 3.
**.**  convert all prestop codons to 'X' in translation.
**.**  check if any non-standard codon species
