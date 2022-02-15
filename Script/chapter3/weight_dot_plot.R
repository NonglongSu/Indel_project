library(tidyverse)

#####################

dot_matrix <- function(a, b, n=1) {
  ab  <- str_split(c(a,b),"",simplify=TRUE)
  ab  <- t(ab != "-")
  len <- nrow(ab) - colSums(!ab)
  m   <- matrix(0, nrow = 2*len[1]+1, ncol = 2*len[2]+1)
  ai  <- 1
  bi  <- 1
  for(i in seq_len(nrow(ab))) {
    m[ai,bi] <- n
    if(ab[i,1] == TRUE && ab[i,2] == TRUE) {#match
      m[ai+1,bi+1] <- n
      ai <- ai+2
      bi <- bi+2
    } else if(ab[i,1] == TRUE) {#deletion
      m[ai+1,bi] <- n
      ai <- ai+2
    } else if(ab[i,2] == TRUE) {#insertion
      m[ai,bi+1] <- n
      bi <- bi+2
    }
  }
  m[ai,bi] <- n
  m
}

count_matrix0 = function(f){
  dat = jsonlite::fromJSON(f)
  dat = dat %>% unpack(aln)
  aln = dat %>% count(Seq1,Seq2)
  mat = 0
  for(i in seq_len(nrow(aln))) {
    mat = mat + dot_matrix(aln$Seq1[i], aln$Seq2[i], aln$n[i])
  }
  return(mat)
}

count_matrixW = function(f,wi){
  dat = jsonlite::fromJSON(f)
  dat = dat %>% unpack(aln)
  aln = dat %>% count(Seq1,Seq2)
  mat = 0
  
  us.wi = unique(sort(wi))
  if(length(aln$n)!=length(us.wi)){
    print("Warning: the number of weight is not equal to number of paths!")
  }
  norm.us.wi = us.wi/sum(us.wi)
  for(i in seq_len(nrow(aln))) {
    mat = mat + dot_matrix(aln$Seq1[i], aln$Seq2[i], aln$n[i])*norm.us.wi[i]
  }
  return(mat)
}


setwd("~/Dropbox/Indel_project/chapter3")
####################################################
#read input 
wtable = read.table("weight_mat.txt")
wm     = as.matrix(wtable)

#Final json
Dir   = "Json"
Files = list.files(Dir, full.names=T)
n     = length(Files)

#original json
Dir0  = "Json0"
File0 = list.files(Dir0, full.names=T)


for (j in 1:n) {
  mat0 = count_matrix0(File0[j])
  
 #mat  = count_matrixW(Files[j], wm[j,])
  mat  = count_matrix0(Files[j])
  
  par(mfrow=c(1,2))
  plot(as.raster(1-mat0/100,interpolate=F))
  plot(as.raster(1-mat/100, interpolate=F))
}
 
 



