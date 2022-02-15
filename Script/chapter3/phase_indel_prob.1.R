#>>>>>>>>>>
#M0->M0   (1-i0)(1-d0)
#M0->D0   (1-i0)*d0
#M0->I0   i0
#M0->End  1-i0       #conditonal on no insertion.
#D0->M0   1-e
#D0->D0   e
#D0->I0   0
#D0->End  1
#I0->M0   (1-e)(1-d0)
#I0->D0   (1-e)*d0
#I0->I0   e
#I0->End  1-e

#M1->M1   (1-i1)(1-d1)
#M1->D1   (1-i1)*d1
#M1->I1   i1
#D1->M1   1-e
#D1->D1   e
#D1->I1   0
#I1->M1   (1-e)(1-d1)
#I1->D1   (1-e)*d1
#I1->I1   e

#M2->M2   (1-i2)(1-d2)
#M2->D2   (1-i2)*d2
#M2->I2   i2
#D2->M2   1-e
#D2->D2   e
#D2->I2   0
#I2->M2   (1-e)(1-d2)
#I2->D2   (1-e)*d2
#I2->I2   e

decompose_A = function(A){
  ab = str_split(A,'',simplify=T)
  ab = t(ab!='-')
  st = rep(0,width(A)[1]+2)
  st[c(1,length(st))] = c(7,10)       #start/end
  eflag = FALSE

  for (i in 1:width(A)[1]) {
    rem   = i%%3
    j     = i+1
    if(ab[i,1]==TRUE && ab[i,2]==TRUE){#match
      eflag = 'M'
      if(rem==1){#phase0
        st[j] = 1
      }else if(rem==2){#phase1
        st[j] = 4
      }else{#phase2
        st[j] = 7
      }
    }else if(ab[i,1]==TRUE){#del
      if(eflag=='D'){
        st[j] = st[j-1]
      }else{
        if(rem==1){
          st[j] = 2
        }else if(rem==2){
          st[j] = 5
        }else{
          st[j] = 8
        }
        eflag = 'D'
      }
    }else if(ab[i,2]==TRUE){#ins  ###I need a flag
      if(eflag=='I'){
        st[j] = st[j-1]
      }else{
        if(rem==1){
          st[j] = 3
        }else if(rem==2){
          st[j] = 6
        }else{
          st[j] = 9
        }
        eflag = 'I'
      }
    }
  }
  #M0 D0 I0 | M1 D1 I1 | M2 D2 I2
  return(st)
 }

p_matrix = function(i0,i1,i2,d0,d1,d2,e){
  M2M0  = (1-i0)*(1-d0)
  M2D0  = (1-i0)*d0
  M2I0  = i0
  D0M0  = 1-e[2]
  D0D0  = e[2]
  D0I0  = 0
  I0M0  = (1-e[1])*(1-d0)
  I0D0  = (1-e[1])*d0
  I0I0  = e[1]
  
  M0M1  = (1-i1)*(1-d1)
  M0D1  = (1-i1)*d1
  M0I1  = i1
  D1M1  = 1-e[2]
  D1D1  = e[2]
  D1I1  = 0
  I1M1  = (1-e[1])*(1-d1)
  I1D1  = (1-e[1])*d1
  I1I1  = e[1]
  
  M1M2  = (1-i2)*(1-d2)
  M1D2  = (1-i2)*d2
  M1I2  = i2
  D2M2  = 1-e[2]
  D2D2  = e[2]
  D2I2  = 0
  I2M2  = (1-e[1])*(1-d2)
  I2D2  = (1-e[1])*d2
  I2I2  = e[1]
  
  M2End = 1-i0
  D0End = 1
  I0End = 1-e[1]
  
  pm = matrix(NA,10,10)
  pm[1,4:6] = c(M0M1,M0D1,M0I1)
  pm[2,1:3] = c(D0M0,D0D0,D0I0)
  pm[3,1:3] = c(I0M0,I0D0,I0I0)
  pm[4,7:9] = c(M1M2,M1D2,M1I2)
  pm[5,4:6] = c(D1M1,D1D1,D1I1)
  pm[6,4:6] = c(I1M1,I1D1,I1I1)
  pm[7,1:3] = c(M2M0,M2D0,M2I0)
  pm[8,7:9] = c(D2M2,D2D2,D2I2)
  pm[9,7:9] = c(I2M2,I2D2,I2I2)
  pm[c(7,2,3),10] = c(M2End,D0End,I0End)
  
  colnames(pm) = c('M0','D0','I0','M1','D1','I1','M2','D2','I2','End')
  rownames(pm) = c('M0','D0','I0','M1','D1','I1','M2','D2','I2','End')
  return(pm)
}

trans_matrix = function(decA){
  tm = matrix(0,10,10)
  colnames(tm) = c('M0','D0','I0','M1','D1','I1','M2','D2','I2','End')
  rownames(tm) = c('M0','D0','I0','M1','D1','I1','M2','D2','I2','End')

  for (i in 1:(length(decA)-1)) {
    tm[decA[i],decA[i+1]] = tm[decA[i],decA[i+1]] + 1
  }
  return(tm)
}

#################################################
i0=g0[1]
i1=g0[2]
i2=g0[3]
d0=g0[4]
d1=g0[5]
d2=g0[6]

#test
#A=DNAStringSet(c("A---AAAAA","AAAA---AA"))

pm    = p_matrix(i0,i1,i2,d0,d1,d2,e3)
pathA = decompose_A(A)

#we need to adjust the gap length by dividing 3.


trm   = trans_matrix(pathA)


sum(log(pm)*trm, na.rm=TRUE)



