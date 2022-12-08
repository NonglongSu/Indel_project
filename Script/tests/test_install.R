install.packages("tidyverse")
install.packages("parallel")
install.packages("S4Vectors")
install.packages("seqinr")
install.packages("stringi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")