
library(Biostrings)
library(stringr)
library(seqinr)
library(plyr)
library(R.utils)

#setwd("~/Dropbox (ASU)/Indel_project/Script")

#only simulate indels of multiple of 3
#insertion starts to the left(include immortal link)
sim_ins = function(s, pos, k){
  ins = sample(DNA_BASES,k,prob=Pi,replace=TRUE)  #throw a k-sided dice for inserted nucs.

  #type-N:0, type-S:1
  po = pos %% 3
  flag = 0

  if(po == 2){#phase1
    ref    = paste0(s[pos-1], s[pos], s[pos+1], collapse="")
    sub1   = paste0(s[pos-1], ins[1], ins[2],   collapse="")
    sub2   = paste0(ins[k],   s[pos], s[pos+1], collapse="")
    sec    = syn[[ref]]
    if(sub1 %in% sec || sub2 %in% sec){
      flag = 1
    }
  }else if(po == 0){#phase2
    ref    = paste0(s[pos-2],  s[pos-1], s[pos], collapse = "")
    sub1   = paste0(s[pos-2],  s[pos-1], ins[1], collapse = "")
    sub2   = paste0(ins[k-1],  ins[k],   s[pos], collapse = "")
    sec    = syn[[ref]]
    if(sub1 %in% sec || sub2 %in% sec){
      flag = 1
    }
  }else{#phase0
      flag = 1
  }

  #selection check
  if( (flag == 1) || (runif(1L) < Wz) ){
    s = insert(s, pos, ins)
  }
  return(s)
}

#deletion starts to the right(pos included)
sim_del = function(s, pos, k){
  po = pos %% 3
  flag = 0
  if(po == 2){#phase 1
    ref    = paste0(s[pos-1],  s[pos+k], s[pos+k+1],collapse="")
    sub1   = paste0(s[pos-1],  s[pos],   s[pos+1],  collapse="")
    sub2   = paste0(s[pos+k-1],s[pos+k], s[pos+k+1],collapse="")
    sec    = syn[[ref]]
    if(sub1 %in% sec || sub2 %in% sec){
      flag = 1
    }
  }else if(po == 0){#phase2
    ref    = paste0(s[pos-2],  s[pos-1],  s[pos+k],collapse="")
    sub1   = paste0(s[pos-2],  s[pos-1],  s[pos],  collapse="")
    sub2   = paste0(s[pos+k-2],s[pos+k-1],s[pos+k],collapse="")
    sec    = syn[[ref]]
    if(sub1 %in% sec || sub2 %in% sec){
      flag = 1
    }
  }else{#phase0
      flag = 1
  }

  #selection check
  if(flag == 1 || (runif(1L) < Wz)){
    s = s[-pos: -(pos-1+k)]
  }
  return(s)
}

#substition starts on the choosing position.
sim_sub = function(s, pos, gtr){
  mu = sample(DNA_BASES, 1L, prob = gtr[which(DNA_BASES == s[pos]),])  #If nuc A, it must change to C,G,T.

  po = pos %% 3
  flag = 0
  if(po == 2){#phase 1
    flag  = 0
  }else if(po == 0){#phase2
    ref = paste0(s[(pos-2):pos],collapse="")
    sub = paste0(s[pos-2],s[pos-1],mu,collapse="")
    sec = syn[[ref]]
    if(sub %in% sec){
      flag = 1
    }
  }else{#phase0
    ref = paste0(s[pos:(pos+2)],collapse="")
    sub = paste0(mu, s[pos+1],s[pos+2],collapse="")
    sec = syn[[ref]]
    if(sub %in% sec){
      flag = 1
    }
  }

  if(flag == 1 || (runif(1L) < W)){
    s[pos] = mu
  }
  return(s)
}



# One-time simulation
SID.onestep = function(r1, r2, gtr, dna0){
  dna = dna0
  tau = 0
  repeat {
    #total rates
    phase   = (seq_along(dna)+2) %% 3 + 1
    phase.r = cbind(r1[phase]/2, r1[phase]/2, r2[dna])  #ins:del:sub

    total_rate = sum(phase.r) #sum whole seq

    #poisson waiting time
    if((total_rate == 0) || (tau >= brlen)){
      break
    }else{
      tau = tau + rexp(1,total_rate)
    }

    #sample pos & mut type
    Po  = sample(length(phase.r), 1L, prob=phase.r)
    Pos = arrayInd(Po, dim(phase.r))

    #draw geometric dist of indel len(1codon=3nucs)
    k   = (rgeom(1, prob=ins_exit)+1)*3

    #update dna
    if(Pos[2] == 1){#insertion
      dna  = sim_ins(dna, Pos[1], k)
    }else if(Pos[2] == 2){#deletion
      dna  = sim_del(dna, Pos[1], k)
    }else{#substitution
      dna = sim_sub(dna, Pos[1], gtr)
    }
  }

  print(tau)
  return(paste0(dna,collapse=""))
}


main = function(l, ouFile){
  ##################################PART I pars setup
  # Global parameters
  W         <<- 0.2                            #dnds
  Wz        <<- 1                              #znzs
  ins_exit  <<- 1-1/6                          #exit prob
  Pi        <<- c(0.2, 0.3, 0.3, 0.2)          #nuc freqs: A, C, G, T
  brlen     <<- 0.5                            #branch length

  # construct codons and its degeneracy
  codons = cbind(rep(DNA_BASES, each = 16),
                 rep(DNA_BASES, times = 4, each = 4),
                 rep(DNA_BASES, 16))
  codonstrs  = apply(codons, 1, stringr::str_c, collapse = "")
  syn        = syncodons(codonstrs)
  names(syn) = toupper(names(syn))
  syn        <<- lapply(syn, toupper)


  # 0.6,0.1,0.3 are based on the phased-indels across human-mouse-rat coding regions.
  # treate insertion and deletion as the same process: gammaID = gammaI + gammaD

  #assume indel parameters -- https://doi.org/10.1093/molbev/msn275 (12-16 indels per 100 subs)
  g = 14/100 *brlen
  #instantaneous indel rate
  gamma = -log(1-g)/brlen

  #indel rates:
  phase = c(0.6,0.15,0.25)   #use phase to simulate the {g0,g1,g2}
  r1    = gamma*phase*1/2

  #subs rates: GTR model.
  Sigma   = 1:6
  r       = matrix(0,4,4)
  r[lower.tri(r)] = Sigma
  r       = r+t(r)
  r       = t(r*Pi)
  diag(r) = -rowSums(r)
  T   = -sum(diag(r)*Pi)
  gtr = r/T
  print(-sum(diag(gtr)*Pi))

  r2        = -diag(gtr)
  names(r2) = DNA_BASES
  diag(gtr) = 0            #avoid negative prob.


  #########################################PART II Run simulation
  set.seed(8088)
  len = as.numeric(l)
  dna0  = sample(DNA_BASES, len, prob = Pi, replace = TRUE)
  dna00 = paste0(dna0,collapse='')


  #run 5 loops
  sim.IDS = list()
  for(i in 1:5){
    sim.IDS[[i]] = SID.onestep(r1, r2, gtr, dna0)
  }

  print(sim.IDS)

  evo.seq = DNAStringSet(c(dna00,sim.IDS[[1]]))
  names(evo.seq) = c('Seq1','Seq2')
  writeXStringSet(evo.seq, ouFile)
}

####################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2])
