#1: output the prob of ziqi's model
#2: output the summuary stats

#gap length>12 are unlikely to preserve the originial aa.


#Generate the phase matrix [zn,zs,long]
#Generate a list of non-gap-site-indicator-removed positions
#Count the number of indels of {3,6,9,12}
#Sum the length of indels
gap_sum = function(A, syn, maxwid){
  l     = length(A)
  seqc  = str_split(A,'')
  g     = IRangesList(lapply(seqc, function(x){IRanges(x=='-')}))
  
  Imat  = list()
  phase = matrix(0,6,2)            
  lphase= matrix(0,3,2)              #long gap phase
  posEdt= list('ins'=c(),'del'=c(),'l.ins'=c(),'l.del'=c())  
  
  
  if(length(unlist(g)) == 0){#no gaps
    return(list(phase,posEdt,lenG,0))
  }else{
    for (k in seq(2)) {
      wid = width(g[[k]])
      
      if(length(wid)==0){next}
      if(any(wid>maxwid)){
        idx = which(wid>maxwid)
        pos = start(g[[k]][idx])
        stp = end(g[[k]][idx])
        
        rem = pos %% 3
        lphase[1,k] = length(which(rem == 1))
        lphase[2,k] = length(which(rem == 2))
        lphase[3,k] = length(which(rem == 0))
        
        posEdt[[k+2]] = stp+1
        
        print("Warning: long gaps(>12) occured!")
      }
      g369 = g[[k]][which(wid<=12)] 
      
      #if(length(g369)==0){next}
      Imat[[k]] = znzs_cure(seqc,g369,k,syn) 
      if(k==1){
        k0 = k
      }else{
        k0 = k+2
      }
      phase[k0:(k0+2),] = Imat[[k]][[1]] 
      posEdt[[k]]       = Imat[[k]][[2]]
    }
  }
  
  lG = c(length(g[[1]]),length(g[[2]]))
  lT = c(sum(width(g[[1]])), sum(width(g[[2]])))/3
  
  phase           = cbind(phase,c(lphase))
  colnames(phase) = c('zn','zs','long')
  rownames(phase) = c('I0','I1','I2','D0','D1','D2')
  
  res  = list(phase,posEdt,lG,lT)
  return(res)
}

#Judge whether the gap is zn/zs. 
znzs_cure = function(seqc, gk, k, syn){
  if(k==1){#ins
    s=seqc[[2]]
  }else{#del
    s=seqc[[1]]
  }
  
  imat = matrix(0,3,2)
  ps1  = start(gk)
  ps2  = end(gk)
  
  rem = ps1 %% 3
  pos = c()      #record the 0,1-edge phase position
  
  #1A-- -AA AAA  2A-- -AA AAA  3AAA AAA AAA  4AAA AA- --A  >phase1
  #1AAA A-- -AA  2AAA AA- --A  3A-- -A- --A  4A-- -AA AAA
  
  #5AA- --A AAA   >phase2
  #5AAA AA- --A  
  
  #6--- NNN NNN  7--- NNN  8NNN --- >phase0 
  #6NNN --- NNN  7NNN NNN  8NNN NNN
  
  for (i in 1:length(rem)) {
    if(rem[i] == 2){#phase 1 
      if(s[ps2[i]+1]=='-'){#1
        j = ps2[i]+1
        repeat{
          j      = j+1
          tmp.1  = s[j]
          if(tmp.1!='-'){
            break
          }
        }
        tmp.0=s[ps1[i]-1]
        tmp.2=s[j+1]
        pos  = c(pos,j)
      }else if(s[ps2[i]+2]=='-'){#2
        tmp.0=s[ps1[i]-1]
        tmp.1=s[ps2[i]+1]
        j = ps2[i]+2
        repeat{
          j     = j+1
          tmp.2 = s[j]
          if(tmp.2!='-'){
            break
          }
        }
        pos  = c(pos,ps2[i]+1,j)
      }else if(s[ps1[i]-1]=='-'){#1
        j = ps1[i]-1
        repeat{
          j      = j-1
          tmp.0  = s[j]
          if(tmp.0!='-'){
            break
          }
        }
        tmp.1 = s[ps2[i]+1]
        tmp.2 = s[ps2[i]+2]
        pos   = c(pos,ps2[i]+1)
      }else{
        tmp.0 = s[ps1[i]-1]
        tmp.1 = s[ps2[i]+1]
        tmp.2 = s[ps2[i]+2]
        pos   = c(pos,ps2[i]+1)
      }
      ref    = paste0(tmp.0,tmp.1,tmp.2,          collapse="")
      sub1   = paste0(tmp.0,s[ps1[i]],s[ps1[i]+1],collapse="")
      sub2   = paste0(s[ps2[i]],tmp.1,tmp.2,      collapse="")
      sec    = syn[[ref]]
      if(sub1 %in% sec || sub2 %in% sec){ 
        imat[2,2]=imat[2,2]+1 
      }else{
        imat[2,1]=imat[2,1]+1 
      }
    }else if(rem[i] == 0){#phase2
      if(s[ps1[i]-2]=='-'){#2
        j = ps1[i]-2
        repeat{
          j      = j-1
          tmp.0  = s[j]
          if(tmp.0!='-'){
            break
          }
        }
        tmp.1 = s[ps1[i]-1]
        tmp.2 = s[ps2[i]+1]
        pos   = c(pos,ps1[i]-1,ps2[i]+1)
      }else if(s[ps2[i]+1]=='-'){#5
        j = ps2[i]+1
        repeat{
          j      = j+1
          tmp.2  = s[j]
          if(tmp.2!='-'){
            break
          }
        }
        tmp.0=s[ps1[i]-2]
        tmp.1=s[ps1[i]-1]
        pos  = c(pos,j)
      }else if(s[ps1[i]-1]=='-'){#5
        j = ps1[i]-1
        repeat{
          j      = j-1
          tmp.1  = s[j]
          if(tmp.1!='-'){
            break
          }
        }
        tmp.0 = s[j-1]
        tmp.1 = s[j]
        tmp.2 = s[ps2[i]+1]
        pos   = c(pos,ps2[i]+1)
      }else{
          tmp.0 = s[ps1[i]-2]
          tmp.1 = s[ps1[i]-1]
          tmp.2 = s[ps2[i]+1]
          pos   = c(pos,ps2[i]+1)
      }
      ref    = paste0(tmp.0,tmp.1,tmp.2,          collapse="")
      sub1   = paste0(tmp.0,tmp.1,s[ps1[i]],      collapse="")
      sub2   = paste0(s[ps2[i]-1],s[ps2[i]],tmp.2,collapse="")
      sec    = syn[[ref]]
      if(sub1 %in% sec || sub2 %in% sec){ 
        imat[3,2]=imat[3,2]+1  
      }else{
        imat[3,1]=imat[3,1]+1  
      }
    }else{#phase0
      if((s[ps2[i]+1]=='-') || (is.na(s[ps2[i]+1]))){#6,#8
        if(!is.na(s[ps2[i]+1])){
          j = ps2[i]+1
          repeat{
            j    = j+1
            tmp  = s[j]
            if(tmp!='-'){
              break
            }
          }
          pos =c(pos,j)
        }
      }else if((s[ps1[i]-1]=='-') || (ps1[i]==1)){#6,#7
        pos =c(pos,ps2[i]+1)
      }else{
        pos =c(pos,ps2[i]+1)            
      }
      imat[1,2]=imat[1,2]+1 
    }
  }
  
  res = list(imat,pos)
  return(res)
}


#Judge whether zn/zs for non-gap state
znzs_treat = function(rs, rs1, L2k, i, syn){
  rem  = i %% 3
  wid  = c(3,6,9,12)
  
  #indicator of each position
  idcat = rep(0,4)
  if(rem==2){#phase1
    for (j in seq(4)) {
      if((i+wid[j])>(L2k+1)){
        next
      }
      sub1=substring(rs,i-1,i+1)
      sub2=substring(rs,i+wid[j]-1,i+wid[j]+1)
      ref =paste0(rs1[i-1],rs1[i+wid[j]],rs1[i+wid[j]+1],collapse='')
      sec =syn[[ref]]
      if(!(sub1 %in% sec) && !(sub2 %in% sec)){
        idcat[j]=1
      }
    }
  }else if(rem==0){#phase2
    for (j in seq(4)) {
      if((i+wid[j])>(L2k+1)){
        next
      }
      sub1=substring(rs,i-2,i)
      sub2=substring(rs,i+wid[j]-2,i+wid[j])
      ref =paste0(rs1[i-2],rs1[i-1],rs1[i+wid[j]],collapse='')
      sec =syn[[ref]]
      if(!(sub1 %in% sec) && !(sub2 %in% sec)){
        idcat[j]=1
      }
    }
  }else{#phase0
    return(idcat)
  }
  return(idcat)
}

#Generate the indicator matrix of non-gap edge of every position for two seqs separately
nogap_sum = function(A, rA, syn){
  L1 = width(A)[1]
  L2 = width(rA)
  
  A1 = str_split(A,'')
  rA1= str_split(rA,'')
  
  idm= array(NA,c(L1,4,2))
  for (k in seq(2)) {
    s1 = A1[[k]]
    rs = toString(rA[[k]])
    rs1=rA1[[k]]
    
    i=j=1
    while (i <= L2[k]-2) {
      if(s1[j]=='-'){
        repeat{
          j=j+1
          if(s1[j]!='-'){
            break
          }
        }
      }
      idm[j,,k] = znzs_treat(rs,rs1,L2[k],i,syn)
      i=i+1
      j=j+1
    }
  }
  
  #check
  if((L1==j+1) && (L2[2]==i+1)){
    return(idm)
  }else{
    print("Warning:The indicator matrix setup is wrong!")
  }
}

#Remove gaps only
rmgap1 = function(A){
  Av   = str_split(A,'')
  g    = IRangesList(lapply(Av, function(x){IRanges(x=='-')}))
  rm.A1 = stri_sub_replace_all(A[[1]],from=sort(start(g[[1]])),to=sort(end(g[[1]])),replacement='')
  rm.A2 = stri_sub_replace_all(A[[2]],from=sort(start(g[[2]])),to=sort(end(g[[2]])),replacement='')
  rm.A  = DNAStringSet(c(rm.A1,rm.A2))
  if(all(width(rm.A)%%3 == 0)){
    return(rm.A)
  }else{
    print("Warning:removed-gap of each sequence is not multiple of three")
  }
}

#Trimming the non-gap edges 
edge_trim = function(M.012, posEdt){
  
  ##put the non-gap-indicator-removed-site as NA
  for (k in seq(2)) {
    M.012[c(posEdt[[k]],posEdt[[k+2]]),,k]=NA
  }
  ##put the long-gap-indicator-site of the other seq as 0
  M.012[posEdt[[3]],,2]=rep(0,4)
  M.012[posEdt[[4]],,1]=rep(0,4)
  
  return(M.012)
}


#Remove all gap-positioned string
rmgap2 = function(A){
  Av   = str_split(A,'')
  g    = unlist(IRangesList(lapply(Av, function(x){IRanges(x=='-')})))
  rm.A = stri_sub_replace_all(A,from=sort(start(g)),to=sort(end(g)),replacement='')
  rm.A = DNAStringSet(rm.A)
  if(all(width(rm.A)%%3 == 0)){
    return(rm.A)
  }else{
    print("Warning:removed-gap sequences are not multiple of three")
  }
}

#Generate obs-codon matrix
countN = function(A, codonstrs){
  seqs = str_split(A,'') 
  
  #codon trans matrix  
  nmat = matrix(0,64,64)
  i=1
  while(i<width(A)[1]) {
    c1 = paste0(seqs[[1]][i:(i+2)], collapse = '')
    c2 = paste0(seqs[[2]][i:(i+2)], collapse = '')
    coor1 = which(codonstrs %in% c1)
    coor2 = which(codonstrs %in% c2)
    nmat[coor1,coor2] = nmat[coor1,coor2] + 1
    i=i+3
  }
  
  return(nmat)
}

#remove effect of scaling factor [ancestor/descendant]
test_pab = function(ab,f0,cf0){
  dnaAB    = DNAStringSet(gsub('-','',ab))
  ab.count = trinucleotideFrequency(dnaAB,step=3)
  sumPab   = sum(ab.count*log(cf0))
  return(sumPab)
}



#>testing using long gap

########################################
ziqi_prob_wz = function(A,g0,e3,p0,f0){
  
  g012= g0[1:6]
  wz  = g0[7]
  ei  = e3[1]
  ed  = e3[2]
  
  co.res = codon_call()
  codons    = co.res[[1]]
  codonstrs = co.res[[2]]
  syn       = co.res[[3]]
  
  #summary the gap indicator matrix
  #setup a maximum gap length of 12 for omegaz-affected indels.
  sumG  = gap_sum(A,syn,12)       
  N.012 = sumG[[1]]
  posEdt= sumG[[2]]
  lG    = sumG[[3]] 
  lenT  = sumG[[4]] 
  
  #average gap length
  scoreE= (lenT[1]-lG[1])*log(ei)+(lG[1])*log(1-ei) + (lenT[2]-lG[2])*log(ed)+(lG[2])*log(1-ed)
  avg.g = lenT/lG
  
  
  #summary the non-gap indicator matrix
  rA     = rmgap1(A)                         
  M.012  = nogap_sum(A, rA, syn) 
  M.trim = edge_trim(M.012, posEdt)          
  
  
  #>>>>>>>>>>>>>>>>Cal. the score
  ##context-based gap len dist
  # if(sum(Lg)==0){#no gap
  #   scoreE = 0
  #   lk     = rep(0,4)
  # }else{
  #   lk     = Lg/sum(Lg)
  #   eId    = lk[which(lk!=0)]
  #   lg     = Lg[which(Lg!=0)]
  #   scoreE = sum(lg*log(eId))
  # }
  #geometric gap length 
  e3.avg = mean(e3)
  lk     = (1-e3.avg)*e3.avg^c(0,1,2,3)
   
  #gap prob
  scoreG = sum(N.012[,1]*log(g012*wz) + N.012[,2]*log(g012) +  N.012[,3]*log(g012)) 
  
  #non-gap prob
  k1 = t(t(M.trim[,,1])*lk)
  k2 = t(t(M.trim[,,2])*lk)
  scoreM = sum(log(1-g012[1:3] + (1-wz)*g012[1:3]*rowSums(k1)), na.rm=TRUE) +
           sum(log(1-g012[4:6] + (1-wz)*g012[4:6]*rowSums(k2)), na.rm=TRUE) 
          
  
  ##substition prob 
  cf0      = sapply(seq(64), function(x){prod(f0[match(cod[x, ], DNA_BASES)])})
  scal.pab = test_pab(A[[1]],f0,cf0) + test_pab(A[[2]],f0,cf0)
  print(-scal.pab)
  
  rA2    = rmgap2(A)
  dat    = countN(rA2,codonstrs)
  scoreP = LL(codons,codonstrs,syn,p0,f0,cf0,dat) 
  
  ##sum LL
  score_ziqi = scoreE + scoreG + scoreM + scoreP - scal.pab
  
  #summary stat
  res = list(score_ziqi,M.trim,N.012,lG,lenT,dat)
  return(res)
}




