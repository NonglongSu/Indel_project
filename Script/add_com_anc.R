#Where the outgroups provide no information will be coated as N.
#This is dna-based MP method.

library(Biostrings)
library(stringr)
library(seqinr)

# setwd("~/Dropbox (ASU)/Indel_project/Script")

# Deal with codon without gaps
Find_com_anc_good = function(seq1, seq2, Ant, site, Tag){     
  
  init = Tag                                                #init cannot overpass Tag
  for(w in 1:length(site)){
    most = max(as.data.frame(table(site[[w]]))$Freq)        #Find how many times the most frequent base occurs
    diff = length(unique(site[[w]]))                        #Find how many different bases in this tree-set 
    df   = as.data.frame(table(site[[w]])) 
    var  = as.character(df[[1]])                       
    if(any(site[[w]] %in% '=') || any(site[[w]] %in% '+')){ #A A T  
      seq1 = seq1[-init:(-init - 2)]                        #= = = 
      seq2 = seq2[-init:(-init - 2)]                        #A T G 
      Ant  = Ant [-init:(-init - 2)]                        #= = = 
      Tag  = init  
      break
    }else{
      if(most == 1){                                        #A/T/G/C
         Ant[Tag] = 'N'
      }else if(most == 2 && diff == 2){                     #A/T/A/T , A/T/T/A , A/A/G/G  | randomly select one of two bases.
         Ant[Tag] = site[[w]][1]
      }else{                                                #A/A/G/T , G/T/A/A , A/A/A/G, A/A/A/A
         base     = var[which(df$Freq == max(df$Freq))]  
         Ant[Tag] = base
      }
      Tag         = Tag + 1 
      }
    }  
  res = list(seq1, seq2, Ant, Tag)
  return(res)
}

# Deal with codon containing gaps
Find_com_anc_bad = function(seq1, seq2, Ant, site, Tag, Wind){
  
  init = Tag
  for(w in 1:Wind){
    most = max(as.data.frame(table(site[[w]]))$Freq)
    diff = length(unique(site[[w]]))
    df   = as.data.frame(table(site[[w]]))
    var  = as.character(df[[1]])
    if(most == 2 && diff == 2){      #A/T/A/T , A/T/T/A , A/A/G/G  
      Ant[Tag] = site[[w]][1]
    }else if(most == 1){             #A/G/C/T , A/G/T/-
      Ant[Tag] = 'N'
    }else{                           #A/A/A/- , -/-/A/- , A/T/-/A , A/A/T/- , A/A/A/A 
      base     = var[which(df$Freq == max(df$Freq))]
      Ant[Tag] = base 
    }
    Tag        = Tag + 1
    }
  res = list(seq1, seq2, Ant, Tag)
  return(res)
}
 

# inFile = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_cds/ENSG00000100344.fa"
# dna    = readBStringSet(inFile, format = "fasta")
# inD    = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_cds"
# ouD    = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_anc/"
# H = "AA---A"
# C = "AA---A"
# M = "TT---T"
# R = "ATGATA"

main = function(inD, ouD){
  
  Files = list.files(inD, full.names = TRUE)
  
  for (n in 1:length(Files)) {
    print(basename(Files[n]))
    
    dna     = readBStringSet(Files[n], format = "fasta")
    name    = names(dna)
    len     = min(width(dna))
    seq.set = list()
    for (k in 1:length(dna)) {
      seq.set[[k]] =  str_split(dna[[k]], "")[[1]]
    }
    
    # Grap all legit gaps
    gaps      = lapply(seq.set, function(x){IRanges(x == '-')})
    freeGaps  = unlist(IRangesList(gaps))
    
    # Find commomn ancestor and update the focal sequences
    # Determine the focal species first
    p1  = seq.set[[3]]
    q1  = seq.set[[4]]
    tag = i = 1 
    anc = c()
    
    while(i < len){                               #accmulation unit = 3
      cSite = list()
      for(j in i:(i+2)){                          #Start with the codon set
        Site  = unlist(lapply(seq.set, `[[`, j))
        cSite = c(cSite, list(Site))
      }
      #Create a decision-maker
      Bell = any(lapply(cSite, function(x){any(x  %in% '-')}) == TRUE)
      if(Bell){                                   #codon frame breaks due to the phase-1,2 indels.
        for (k in 1:3) {
          if(any(cSite[[k]] == '-')){             #Determine the correct gap position
            gap.pos = i + k - 1
            gap.wid = width(freeGaps[which(start(freeGaps) == gap.pos)])[1]
            if(gap.pos %% 3 == 1){                #phase-0 (between-codon gaps)
              Window = i:(i + gap.wid - 1)         
            }else{                                #phase-1,2 (inside-codon gaps+edges)     
              Window = i:(i + gap.wid + 2)         
            }
            break
          }
        }
        #Update the window 
        dSite = list()
        for(m in Window){
          Site  = unlist(lapply(seq.set, `[[`, m))
          dSite = c(dSite, list(Site))
        }
        res = Find_com_anc_bad(p1, q1, anc, dSite, tag, length(Window))
        i   = tail(Window, 1) + 1 
      }else{                                      #codon frame is integrate.
        res = Find_com_anc_good(p1, q1, anc, cSite, tag)
        i   = i + 3
      }
      p1  = res[[1]]
      q1  = res[[2]]
      anc = res[[3]]
      tag = res[[4]] 
    }
    
    # check if the length equally 
    p2      = paste(p1,  collapse = "")
    q2      = paste(q1,  collapse = "")
    anc2    = paste(anc, collapse = "")
    seq.new = list(anc2, p2, q2)
    Name    = c("ancestor", name[3:4])
    write.fasta(seq.new, Name, nbchar = 80, open = "w", as.string = TRUE, file.out = paste0(ouD, basename(Files[n])))
  }
  
}







args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])


  
  
 