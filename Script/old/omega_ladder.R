#Cal the dnds value for every sw alignment
#Create a table (geneId::dnds:gap)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(dplyr))

#setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat")
########################################################

#Find 3-mer gaps of -/=/+
find_gap = function(seq, sign){
  g        = lapply(str_split(seq, ""), function(x){IRanges(x == sign)})
  g        = IRangesList(g)
  gaps     = length(which(unlist(width(g)) == 3))
  return(gaps)
}

#Determine number of non-syn/syn sites
cal_sites = function(x, tag){
  unit.x  = substr(x, tag, tag + 2)
  toMatch = c('-', '\\+', '=', 'N')
  s.site  = 0
  
  if(grepl(paste0(toMatch, collapse = "|"), unit.x)){
    return(c(0, 0))
  }else{
    section = codon[[which(sapply(codon, function(X){unit.x %in% X}))]]
    for (k in 1:3) {
      Base     = DNA_BASES[substr(unit.x, k, k) != DNA_BASES]
      unit.sub = unit.x
      for (w in 1:length(Base)) {
        substr(unit.sub, k, k) = Base[w]
        if(unit.sub %in% section){
          s.site = s.site + 1/3
        }
      }
    }
  }
  n.site = 3 - s.site
  res = c(n.site, s.site)
  return(res)
}

#Determine number of non-syn/syn substitutions
cal_subs = function(x, y, tag){
  
  unit.x  = substr(x, tag, tag + 2)
  unit.y  = substr(y, tag, tag + 2)
  sub.s   = 0   
  sub.n   = 0  
  toMatch = c('-', '\\+', '=', 'N')
  
  if(grepl(paste0(toMatch, collapse = "|"), unit.x) || grepl(paste0(toMatch, collapse = "|"), unit.y)){
    return(c(sub.s, sub.n))
  }else{
    section = codon[[which(sapply(codon, function(X){unit.x %in% X}))]]
    mut.num = mapply(function(X, Y) sum(X != Y), strsplit(unit.x, ""), strsplit(unit.y, ""))
    mut.pos = mapply(function(X, Y) which(X != Y), strsplit(unit.x, ""), strsplit(unit.y, ""))
    if(mut.num == 1){#ATG/TTG
      if(unit.y %in% section){
        sub.s = sub.s + 1
      }else{
        sub.n = sub.n + 1
      }
    }else if(mut.num == 2){#ATG/TCG
      unit.list = collect_allSubs(unit.x, unit.y, mut.pos)
      sub.rate  = cal_subsub(unit.list, sub.s, sub.n, section)
      sub.s = sub.rate[1] / 2
      sub.n = sub.rate[2] / 2
    }else if(mut.num == 3){#ATG/CCT
      unit.list = collect_allSubs(unit.x, unit.y, mut.pos)
      sub.rate  = cal_subsub(unit.list, sub.s, sub.n, section)
      sub.s = sub.rate[1] / 6
      sub.n = sub.rate[2] / 6
    }else{#no substitution ATG/ATG
      return(c(sub.s, sub.n))
    }
  }
  return(c(sub.s, sub.n))
}

#Find all possible substitution pathways.        
# unit.x = "TTT"
# unit.y = "TAC"
# unit.y = "GAC"
# mut.pos = 2:3
# mut.pos = 1:3
collect_allSubs = function(unit.x, unit.y, mut.pos){
  unit.set  = c()
  for (i in 1:length(mut.pos)) {
    unit.sub = unit.x
    substr(unit.sub, mut.pos[i], mut.pos[i]) = substr(unit.y, mut.pos[i], mut.pos[i])
    unit.set = c(unit.set, unit.sub)
    if (length(mut.pos) == 3) {
      for (j in mut.pos[mut.pos != mut.pos[i]]) {
        unit.sb               = unit.sub
        substr(unit.sb, j, j) = substr(unit.y, j, j)
        unit.set              = c(unit.set, unit.sb)
      }
    }
  }
  unit.list = list()
  if(length(unit.set) == 2){
    for (k in 1:length(unit.set)) {
      unit.list = c(unit.list, list(c(unit.x, unit.set[k], unit.y)))
    }
  }else{
    for (k in 1:length(unit.set)) {
      if (k %% 3 == 1) {
        for (w in (k + 1):(k + 2)) {
          unit.list = c(unit.list, list(c(unit.x, unit.set[k], unit.set[w], unit.y)))
        }
      }
    }
  }
  return(unit.list)
}

#Sum the number of substitutions of all possible pathways
# section = c("TTT","TTC")
# 12 , 6
cal_subsub = function(unit.list, sub.s, sub.n, sec){
  for (m in 1:length(unit.list)) {
    sec.var = sec
    for (n in 2:(length(unit.list[[m]])) ) {
      if(unit.list[[m]][n] %in% sec.var){
        sub.s = sub.s + 1
      }else{
        sub.n = sub.n + 1
      }
      sec.var = codon[[which(sapply(codon, function(X){unit.list[[m]][n] %in% X}))]]
    }
  }
  return(c(sub.s, sub.n))
}

#Determine the bin size
binSum = function(gVec, wVec, binMax){
  t.limits = 0
  g.limits = 0
  g.bound  = c()
  w.bound  = c()
  for (i in 1:length(gVec)) {
    t.limits = t.limits + gVec[i]
    g.limits = g.limits + 1
    if(t.limits >= binMax){
      g.bound  = c(g.bound, i)
      w.bound  = c(w.bound, wVec[i])
      t.limits = 0
      g.limits = 0
    }
  }
  g.bound = c(g.bound, i)
  w.bound = c(w.bound, wVec[i])
  res     = list(g.bound, w.bound)
  return(res)
}



######################################################################
# inD   = "Data_6.2/Mafft/mixed_cds"
# ouF1  = "Data_6.2/Results/omega_ladder.txt"
# ouF2  = "Data_6.2/Results/omega_bin.txt"
# ouFig = "Data_6.2/Figure/omega_dis.pdf"

main = function(inD, ouF1, ouF2, ouFig){
  
  codon <<- list (c("TTT","TTC"),
                  c("TTA","TTG","CTT","CTC","CTA","CTG"),
                  c("ATT","ATC","ATA"),
                  c("ATG"),
                  c("GTT","GTC","GTA","GTG"),
                  c("TCT","TCC","TCA","TCG","AGT","AGC"),
                  c("CCT","CCC","CCA","CCG"),
                  c("ACT","ACC","ACA","ACG"),
                  c("GCT","GCC","GCA","GCG"),
                  c("TAT","TAC"),
                  c("CAT","CAC"),
                  c("CAA","CAG"),
                  c("AAT","AAC"),
                  c("AAA","AAG"),
                  c("GAT","GAC"),
                  c("GAA","GAG"),
                  c("TGT","TGC"),
                  c("TGG"),
                  c("CGT","CGC","CGA","CGG","AGA","AGG"),
                  c("GGT","GGC","GGA","GGG"),
                  c("TAA","TGA","TAG")
  )
  
  Files = list.files(inD,full.names=T)
  
  DF = c()
  for (i in 1:length(Files)) {
    dna = readBStringSet(Files[i])
    len = width(dna)[1]
  
    M   = toString(dna[1]) 
    R   = toString(dna[2]) 
    
    #Counting gaps
    Gaps  = find_gap(dna, "-") 
    
    #Cal. the dN/dS of each one 
    n.m   = 0 #number of non-synonymous sites (relative to species)
    s.m   = 0 #number of synonymous sites (relative to species)
    n.r   = 0
    s.r   = 0
    
    N.sub = 0 #number of non-synonymous mutations
    S.sub = 0 #number of synonymouse mutations
    
    j = 1       #reset to 1 each run
    while (j < len) {
      Subs    = cal_subs(M, R, j)
      S.sub   = S.sub + Subs[1]
      N.sub   = N.sub + Subs[2]
      
      #Mouse as ref
      nSite_M  = cal_sites(M, j)
      n.m      = n.m + nSite_M[1]
      s.m      = s.m + nSite_M[2]
      
      #Rat as ref
      nSite_R  = cal_sites(R, j) 
      n.r      = n.r + nSite_R[1]
      s.r      = s.r + nSite_R[2]
      
      j = j + 3  #use codon length as increment unit
    }
    
    #Cal. the W (dN/dS) (Nei and Gojobori 1986)
    avg.nSite = (n.m + n.r) / 2
    avg.sSite = (s.m + s.r) / 2
      
    P.n  = N.sub / avg.nSite
    P.s  = S.sub / avg.sSite
    dNdS = P.n / P.s
    
    #Dataframe
    df = data.frame("geneId" = gsub("\\..*", "", basename(Files[i])), 
                    "dNdS"   = dNdS, 
                    "gap_3"  = Gaps)
    DF = rbind(DF, df)
  }
  
  # Rank the dataframe by dnds value
  DF.full  = DF %>% filter(dNdS != "NAN" & dNdS != "Inf")
  DF.rank  = DF.full[order(DF.full$dNdS), ]
  
  # Determine the bin size then plot
  g.Vec    = DF.rank$gap_3 
  w.Vec    = DF.rank$dNdS  
  binMax   = 96
  repeat{#At least 10 groups
    Threshold = sum(g.Vec)/ binMax
    if(Threshold >= 9){
      break
    }else{
      binMax = binMax - 1
    }
  }
  bounD = binSum(g.Vec, w.Vec, binMax)
  
  # Generate a database
  write.table(DF.rank, file = ouFile1,
              sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE )
  write.table(bounD[[1]], file = ouFile2,
              sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #Plot
  pdf(ouFig)
  safe_palate = "#88CCEE"
  breaDex     = bounD[[2]]
  if(all(breaDex != 0)){
    breaDex = c(0, breaDex)
  }
  hist(DF.rank$dNdS, breaks = breaDex, col = safe_palate, 
       xlim = range(0, max(w.Vec) + 0.2), 
       xlab = "omega", main = "The distribution of dNdS ratio across the genome-wide coding regions")
  dev.off()

  
}



args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4])