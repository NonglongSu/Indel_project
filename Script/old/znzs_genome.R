#One real-speciic insert occurs in all locations.

library(Biostrings)
library(stringr)
library(seqinr)

# setwd("~/Dropbox (ASU)/Indel_project/Script")

#Cal. the number of non-syn/syn sites of deletions.
cal_d_sites = function(anc.seq,wid,len){
  D_mer = len-wid+1
  anc.char = str_split(anc.seq,"")[[1]] 
  
  s.d = 0
  n.d = 0
  i = 1
  while (i <= D_mer) {
    if(i %% 3 == 1){#Phase-0
      unit.1 = substr(anc.seq,i,i+2)          #left-amino acid
      unit.2 = substr(anc.seq,i+wid,i+wid+2)  #right-amino acid
      if(grepl('N',unit.1) || grepl('N',unit.2)){
        i = i + 3
      }else{
        s.d = s.d + 1
        i = i + 1
      }
    }else{
      if(i %% 3 == 0){#phase-2
        sub = paste0(c(anc.char[i-2],anc.char[i-1],anc.char[i+wid]),collapse = "")
        sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
        if(unit.1 %in% sec || unit.2 %in% sec){
          s.d = s.d + 1
        }else{
          n.d = n.d + 1  
        }
      }else{#phase-1
        sub = paste0(c(anc.char[i-1],anc.char[i+wid],anc.char[i+wid+1]),collapse = "")
        sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
        if(unit.1 %in% sec || unit.2 %in% sec){
          s.d = s.d  + 1
        }else{
          n.d = n.d + 1
        }
      }
      i = i + 1
    }
  }
  res = c(s.d,n.d)
  return(res)
}

#Cal. the number of non-syn/syn sites of insertions (use real data)
cal_i_sites = function(anc.seq,wid,len,insert){
  I_mer = len-1  #ignore before start and after stop codon
  ins.char = str_split(insert,"")[[1]] 
  
  s.i = 0
  n.i = 0
  i = 1
  while (i <= I_mer) {
    if(i %% 3 == 1){#Phase-1
      unit = substr(anc.seq,i,i+2)          #original codon
      unit.char = str_split(unit,"")[[1]] 
      sec = codon[[which(sapply(codon, function(X){unit %in% X}))]]
      if(grepl('N',unit)){
        i = i + 3
      }else{
        sub.1 = paste0(c(unit.char[1],ins.char[1],ins.char[2]),collapse = "")
        sub.2 = paste0(c(tail(ins.char,1),unit.char[2],unit.char[3]),collapse = "")
        if(sub.1 %in% sec || sub.2 %in% sec){
          s.i = s.i  + 1
        }else{
          n.i = n.i + 1
        }
        i = i + 1
      }
    }else{
      if(i %% 3 == 2){#phase-2
        sub.1 = paste0(c(unit.char[1],unit.char[2],ins.char[1]),collapse = "")
        sub.2 = paste0(c(tail(ins.char,2),unit.char[3]),collapse = "")
        if(sub.1 %in% sec || sub.2 %in% sec){
          s.i = s.i  + 1
        }else{
          n.i = n.i + 1
        }
      }else{#phase-0
        s.i = s.i + 1
      }
      i = i + 1
    }
  }
  res = c(s.i,n.i)
  return(res)
}


#Cal. the number of observed substitutions
#deletion: x--anc,  y--focal
#inserton: x--focal,y--anc
cal_subs = function(x,y,Gj){
  
  pos.sta = start(Gj)
  pos.end = end(Gj)
  
  y.char = str_split(y,"")[[1]] #y is the sequence with gaps. 
  
  n.indel = s.indel = 0
  
  if(pos.sta %% 3 == 0){#phase-2
    sub = paste0(c(y.char[pos.sta-2],y.char[pos.sta-1],y.char[pos.end+1]),collapse = "")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = substr(x,pos.sta-2,pos.sta)
    unit.2 = substr(x,pos.end-1,pos.end+1)
    if(unit.1 %in% sec || unit.2 %in% sec){
      s.indel = 1
    }else{
      n.indel = 1
    }
  }else if(pos.sta %% 3 == 1){#phase-0
    s.indel = 1
  }else{#phase-1
    sub = paste0(c(y.char[pos.sta-1],y.char[pos.end+1],y.char[pos.end+2]),collapse = "")
    sec = codon[[which(sapply(codon, function(X){sub %in% X}))]]
    unit.1 = substr(x,pos.sta-1,pos.sta+1)
    unit.2 = substr(x,pos.end,pos.end+2)
    if(unit.1 %in% sec || unit.2 %in% sec){
      s.indel = 1
    }else{
      n.indel = 1
    }
  }
  
  res = c(s.indel,n.indel)
  return(res)
}



# inDir  = "../test_human_mouse_rat/Data_6.2/Mafft/mapped_anc"
# ouFile = "../test_human_mouse_rat/Data_6.2/Results/zNzS.txt"

main = function(inDir, ouFile){
  
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
  
  Files = list.files(inDir,full.names = TRUE)
 
  s.md    = 0 #number of synonymous deletion sites in mouse.
  n.md    = 0 #number of non-synonymous deletion sites in mouse.
  S.del_m = 0 #number of synonymous deletions in mouse.
  N.del_m = 0 #number of non-synonymous deletions in mouse. 
  
  s.mi = 0 #number of synonymous insertion sites in mouse.
  n.mi = 0 #number of non-synonymous insertion sites in mouse.
  S.ins_m = 0 #number of synonymous insertions in mouse.
  N.ins_m = 0 #number of non-synonymous insertions in mouse.
  
  s.ri = 0 # number of synonymous insertion sites in rat.
  n.ri = 0 # number of non-synonymous insertion sites in rat. 
  S.ins_r = 0 #number of synonymous insertions in rat.
  N.ins_r = 0 #number of non-synonymous insertions in rat. 

  s.rd = 0 #number of synonymous deletion sites in rat.
  n.rd = 0 #number of non-synonymous deletion sites in rat.
  S.del_r = 0 #number of synonymous deletions in rat. 
  N.del_r = 0 #number of non-synonymous deletions in rat. 
  
  
  for (i in 1:length(Files)) {
    dna = readDNAStringSet(Files[i],format = "fasta")
    len = width(dna)[1]
    
    A = toString(dna[[1]])
    M = toString(dna[[2]])
    R = toString(dna[[3]])
    
    #remove all gaps from ancestor sequence
    A1 = str_remove_all(A,'-')
    M1 = str_remove_all(M,'-')
    R1 = str_remove_all(R,'-')
    
    len.A1  = nchar(A1)
    len.M1  = nchar(M1)
    len.R1  = nchar(R1)
    
    dna.str = str_split(as.character(dna),'')
    g = lapply(dna.str, function(x){IRanges(x=='-')})
    g = IRangesList(g)
    
    a.start = start(g[[1]])
    a.end   = end(g[[1]])
    a.width = width(g[[1]])
    
    m.start = start(g[[2]])
    m.end   = end(g[[2]])
    m.width = width(g[[2]])
    
    r.start = start(g[[3]])
    r.end   = end(g[[3]])
    r.width = width(g[[3]])
    #Cal. the zN/zS from mouse 
    if(length(g[[2]])>0){
      toMatch = lapply(g, function(x){g[[2]] %in% x})
      pat.1 = c(TRUE,TRUE,FALSE)
      pat.2 = c(FALSE,TRUE,FALSE)
      
      for (j in 1:length(g[[2]])) {
        pat.obs = unlist(lapply(toMatch,`[[`,j), use.names = FALSE)
        if(all(pat.obs==pat.1)){#insertion in rat, using gaps from mouse instead of ancestor because this is a mouse loop!
          ins   = substr(R1,m.start[j],m.end[j]) #Find the insert  
          nSite = cal_i_sites(A1,m.width[j],len.A1,ins)
          s.ri  = s.ri + nSite[1]
          n.ri  = n.ri + nSite[2]
          
          Ni       = cal_subs(R,A,g[[2]][j])
          S.ins_r  = S.ins_r + Ni[1]
          N.ins_r  = N.ins_r + Ni[2]
        }
        if(all(pat.obs==pat.2)){#deletion in mouse
          nSite = cal_d_sites(A1,m.width[j],len.A1)
          s.md  = s.md + nSite[1]
          n.md  = n.md + nSite[2]
          
          Nd    = cal_subs(A,M,g[[2]][j])
          S.del_m  = S.del_m + Nd[1]
          N.del_m  = N.del_m + Nd[2]
        }
      }
    }
    #Cal. the dN/dS from rat
    if(length(g[[3]])>0){
      toMatch = lapply(g, function(x){g[[3]] %in% x})
      pat.1 = c(TRUE,FALSE,TRUE)
      pat.2 = c(FALSE,FALSE,TRUE)
      
      for (j in 1:length(g[[3]])) {
        pat.obs = unlist(lapply(toMatch,`[[`,j), use.names = FALSE)
        if(all(pat.obs==pat.1)){#insertion in mouse
          ins   = substr(M1,r.start[j],r.end[j]) #Find the insert
          nSite = cal_i_sites(A1,r.width[j],len.A1,ins)
          s.mi  = s.mi + nSite[1]
          n.mi  = n.mi + nSite[2]
          
          Ni    = cal_subs(M,A,g[[3]][j])
          S.ins_m  = S.ins_m + Ni[1]
          N.ins_m  = N.ins_m + Ni[2]
        }
        if(all(pat.obs==pat.2)){#deletion in rat
          nSite = cal_d_sites(A1,r.width[j],len.A1)
          s.rd  = s.rd + nSite[1]
          n.rd  = n.rd + nSite[2]
          
          Nd    = cal_subs(A,R,g[[3]][j])
          S.del_r  = S.del_r + Nd[1]
          N.del_r  = N.del_r + Nd[2]
        }
      }
    }
  }
  
  #Comparison of observed vs expected. (KaKs = n/s vs N/S) 
  KaKs_M_ins_obs = N.ins_m/S.ins_m 
  KaKs_M_ins_exp = n.mi/s.mi
  
  KaKs_M_del_obs = N.del_m/S.del_m 
  KaKs_M_del_exp = n.md/s.md
  
  KaKs_R_ins_obs = N.ins_r/S.ins_r 
  KaKs_R_ins_exp = n.ri/s.ri
  
  KaKs_R_del_obs = N.del_r/S.del_r 
  KaKs_R_del_exp = n.rd/s.rd
  
  
  #Cal. the zNzS without any model. (to be fixed in the future)
  zN_M_ins = N.ins_m/n.mi 
  zS_M_ins = S.ins_m/s.mi
  
  zN_M_del = N.del_m/n.md
  zS_M_del = S.del_m/s.md
  
  zN_R_ins = N.ins_r/n.ri
  zS_R_ins = S.ins_r/s.ri
  
  zN_R_del = N.del_r/n.rd
  zS_R_del = S.del_r/s.rd
  #Jukes Cantor formula (1969)
  zNzS_M_ins = log(1-4*zN_M_ins/3,exp(1)) / log(1-4*zS_M_ins/3,exp(1))
  zNzS_M_del = log(1-4*zN_M_del/3,exp(1)) / log(1-4*zS_M_del/3,exp(1))
  
  zNzS_R_ins = log(1-4*zN_R_ins/3,exp(1)) / log(1-4*zS_R_ins/3,exp(1))
  zNzS_R_del = log(1-4*zN_R_del/3,exp(1)) / log(1-4*zS_R_del/3,exp(1))
  
  #Generate a database
  df = data.frame("Mouse_zNzS_ins"=zNzS_M_ins, "Mouse_zNzS_del"=zNzS_M_del,
                                           "Rat_dNdS_ins"=zNzS_R_ins, "Rat_zNzS_del"=zNzS_R_del)
  write.table(df, file = ouFile,
              sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE )
}

args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2])
  


