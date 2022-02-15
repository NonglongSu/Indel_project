# Update all disqualified gaps
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringi))

#setwd("~/Dropbox (ASU)/Indel_project/Script")

#Initial word size (Window) as 33(15+3+15)
#convert ---/---       to ===/===
#convert ---AAA---     to ===AAA===

# Find gap range
find_gap_range = function(X){
  X.1   = str_split(as.character(X), '')
  g     = lapply(X.1, function(x) { IRanges(x == '-')})
  g     = IRangesList(g)
  return(g)
}

# Remove any large gaps
filter_Long = function(x, X){
  if(length(x) == 0){
    return(X)
  }else{
    st  = start(x)
    en  = end(x)
    wid = width(x) 
    for(i in 1:length(x)){
      subseq(X, start = st[i], end = en[i]) = stri_rand_strings(1, wid[i], '[=]')
    }
  }
  return(X)
}

# Cleanup of start/end position of the sequences
start_stop_test = function(X){
  Start = substr(X, 1, Window)
  IsGap.1 = grepl('-', Start)
  
  Stop = substr(X, (nchar(X) - (Window - 1)), nchar(X))
  IsGap.2 = grepl('-', Stop) 
  
  if(IsGap.1 == TRUE){
    Start = gsub('-', '=', Start)
    substr(X, 1, Window) = Start
  }
  
  if(IsGap.2 == TRUE){
    Stop = gsub('-', '=', Stop)
    substr(X, (nchar(X) - (Window - 1)), nchar(X)) = Stop
  }
  return(X)
}

#Adjust close gap-gap inside the species
ad_gap_test = function(seq, x){
  
  if(length(x) == 0){
    return(seq)
  }else{
   len       = length(x)
   pos.start = start(x)
   pos.end   = end(x)
   wid.all   = width(x) 
   # at least 3 extra nts to separate windows.
   filter = Wall + 3
   
   for(j in 1:len){
     pos = pos.start[j]
     wid = wid.all[j]
     en  = pos.end[j]
     
     left  = substr(seq, start = pos-filter, stop = pos - 1)
     right = substr(seq, start = pos+wid, stop = pos + wid + (filter - 1))
     l1    = grep("[-=]", left)
     r1    = grep("[-=]", right)
    
     if((length(l1)!= 0) || (length(r1) != 0)){
       subseq(seq, start = pos, end = en) = stri_rand_strings(1, wid, '[=]')
     }
   }
  }
  return(seq)
}

#Deal with special scenarios (---    || ------  || AAA---AAA
#                             ------ || ---     || ---------  )
overlap_gap_test = function(Z, G){
  
  for (i in 1:length(G)) {
    if(length(G[[i]]) == 0){
      next
    }else{
      # test_sta1
      start.i    = start(G[i])[[1]]
      end.i      = end(G[i])[[1]]
      start.res  = start(G[-i])
      end.res    = end(G[-i])
      
      if(length(unlist(start.res)) == 0){#only one seq has gaps
        next
      }else{
        startMatch = lapply(start.res, function(x){start.i %in% x})
        endMatch   = lapply(end.res, function(x){end.i %in% x})
        
        # test_sta2
        test.box = list()
        for (k in 1:length(start.i)) {
          test.box[[k]] = start.i[k]:end.i[k]
        }
        
        rest.box = list(unlist(start(G[-i])), unlist(end(G[-i])))
        Rest.box = list()
        for (k in 1:length(rest.box[[1]])){
          Rest.box[[k]] = rest.box[[1]][k]:rest.box[[2]][k] 
        }
        
        for (j in 1:length(startMatch[[1]])){ 
          flag1 = unlist(lapply(startMatch, "[[", j))
          flag2 = unlist(lapply(endMatch, "[[", j))
          if(any(flag1 != flag2)){ #---/------
            subseq(Z[i], start = start.i[j], end = end.i[j]) = stri_rand_strings(1, width(G[[i]])[j], '[=]')
          }else{
            if(any(unlist(lapply(Rest.box, function(x){any(test.box[[j]] %in% x) && any(test.box[[j]] != x)})))){ #AAA---AAA/---------
              subseq(Z[i], start = start.i[j], end = end.i[j]) = stri_rand_strings(1, width(G[[i]])[j], '[=]')
            }else{
              next
            }
          }
        }
      }
    }
  }
  return(Z)
}
  
# Deal with gaps based on Phylogeny-Aware Gap Placement
phylo_gap_test = function(Z, G){
  
  # Update the remained free-gaps within human/hamster species 
  G.mr = G[(length(G) - 1):length(G)]
  
  for (i in 1:(length(G) - 2)) {
    if(length(G[[i]]) == 0){
      next
    }else{
      toMatch = lapply(G.mr, function(x){G[[i]] %in% x})   
      Start   = start(G[[i]])
      End     = end(G[[i]])
      Wid     = width(G[[i]])
      for (j in 1:length(G[[i]])) {
        pat = unlist(lapply(toMatch, `[[`, j))
        if(all(pat == FALSE)){
          subseq(Z[[i]], start = Start[j], end = End[j]) = stri_rand_strings(1, Wid[j], '[=]')
        }
      }
    }
  }
  
  # Use focal-species gaps to match the rest group.
  pattern.1 = c(rep(TRUE, length(G) - 2), FALSE)
  pattern.2 = rep(FALSE, length(G) - 1)
  
  for (i in (length(G) -1):length(G)) {
    if(length(G[[i]]) == 0){
      next
    }else{
      toMatch = lapply(G[-i], function(x){G[[i]] %in% x})   
      Start   = start(G[[i]])
      End     = end(G[[i]])
      Wid     = width(G[[i]])
      for (j in 1:length(G[[i]])) {
        pat = unlist(lapply(toMatch, `[[`, j))
        if(all(pat == pattern.1) || all(pat == pattern.2)){
          for (k in 1:length(G)) {
            lf = substr(Z[k], start = Start[j] - Wall - 3, stop = Start[j] - 1)
            rt = substr(Z[k], start = Start[j] + Wid[j], stop = Start[j] + Wid[j] + (Wall + 2))
            l2 = grep("[-=]", lf)
            r2 = grep("[-=]", rt)
            # Check if gaps are independent
            if((length(l2) != 0) || (length(r2) != 0)){
              for (w in 1:length(G)) {                                   # Do not use subseq() = stri_rand_strings() because it will replace everything even there are no gaps.
                test.gap = substr(Z[[w]], Start[j], End[j])
                test.gap = gsub('-', '=', test.gap)
                substr(Z[[w]], Start[j], End[j]) = test.gap
              }
              break
            }
          }
        }else{ # Updating wrong phylogeny tree 
          for (w in 1:length(G)) {
            test.gap = substr(Z[[w]], Start[j], End[j])
            test.gap = gsub('-', '=', test.gap)
            substr(Z[[w]], Start[j], End[j]) = test.gap
          }
        }
      }
    }
  }
  return(Z)
}


# inFile = "../test_human_mouse_rat/Raw_data_2/mapped_cds_mafft/ENSG00000001036.fa"
# ouDir  = "../test_human_mouse_rat/Data_6.2/Mafft/updated_cds/"

# num1 = 6
# num2 = 12
main = function(inFile, ouDir, num1, num2){
  
  dna   = readBStringSet(inFile, format = "fasta")
  name  = names(dna)

  # Set up vars
  num1 = as.numeric(num1)
  num2 = as.numeric(num2)
  Window  <<- num1
  Wall <<- num2
  
  # String mode
  spec = c()
  for (i in 1:length(dna)) {
    spec = c(spec, toString(dna[[i]]))
  }
  
  # Find gap range
  gap  = find_gap_range(dna)
  gaps = c(3, 6, 9, 12)
  
  for (i in 1:length(dna)) {
    gi       = gap[[i]]
    wid.i    = width(gi)
    g.null   = gi[!(wid.i %in% gaps)]   
    g.null.l = length(g.null)
    
    spec[i] = filter_Long(g.null, spec[i])
    spec[i] = start_stop_test(spec[i])
  }
  
  # Keep the gaps updated
  gap.1 = find_gap_range(spec)
  
  # Cleanup of close gaps inside the species
  for (i in 1:length(dna)) {
    spec[i]  = ad_gap_test(spec[i], gap.1[[i]])
  }
  gap.2 = find_gap_range(spec)
  
  # Cleanup of the overlap gaps in different species
  spec.3  = overlap_gap_test(spec, gap.2)
  gap.3   = find_gap_range(spec.3)

  # Cleanup of phylogeny-aware indels
  spec.4 = phylo_gap_test(spec.3, gap.3) 
  spec.4 = as.list(spec.4)
  
  # Output
  write.fasta(sequences = spec.4, names = name, nbchar = 80,
              open = "w", as.string = TRUE , file.out = paste0(ouDir, basename(inFile)) )
}







# Compatible with run_test.R
# if(interactive()){
args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4])
#}

