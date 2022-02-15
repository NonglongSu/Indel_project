set.seed(8888)

##-Log-likelihood of indels
# LL_min_z = function(gv){
#   g012 = gv[1:6]
#   wz   = gv[7]
#   #gap length LL
#   if(sum(lgw)==0){#no gap
#     scoreE = 0
#     lk     = rep(0,4)
#   }else{
#     lk     = lgw/sum(lgw)
#     lkreal = lk[which(lk!=0)]
#     lgreal = lgw[which(lgw!=0)]
#     #scoreE = sum(lgreal*log(lkreal))
#   }
#   #gap LL
#   scoreG = sum(phasew[,1]*log(g012*wz) + phasew[,2]*log(g012))
#   #non-gap LL
#   k1 = t(t(matchw[,,1])*lk)
#   k2 = t(t(matchw[,,2])*lk)
#   scoreM = sum(log(1-gv[1:3] + (1-wz)*gv[1:3]*rowSums(k1)), na.rm=TRUE) +
#     sum(log(1-gv[4:6] + (1-wz)*gv[4:6]*rowSums(k2)), na.rm=TRUE)
#   
#   -(scoreG + scoreM)
# }

#>weighted non-gap indicator matrix
# for (h in 1:K) {
#   if(h==1){
#     Izw1 = Iz[[h]][,,1]*Wi[h]
#     Izw2 = Iz[[h]][,,2]*Wi[h]
#   }else{
#     Izw1 = rbind(Izw1,Iz[[h]][,,1]*Wi[h])
#     Izw2 = rbind(Izw2,Iz[[h]][,,2]*Wi[h])
#   }
# }
# if(nrow(Izw1)==nrow(Izw2)){
#   Izw = array(c(Izw1,Izw2),c(nrow(Izw1),4,2))
#   matchw = Izw
# }else{
#   print("Warning:weighted non-gap indicator matrix has incorrect row dimensions!")
# }
#matchw<<-matchw
#matchw <<- Iz[[9]]
#gb  = nmkb(fn=LL_min_z, par=g0, lower=0.0, upper=1.0, control=list(tol=1e-4,trace=T)) 


#del-->subs-->ins.
#Judge whether the gap state is zn/zs. 
#record the one-after-gap-right-end pos.
znzs_cure = function(seqc, gz, k, syn){
  if(k==1){#ins 
    s =seqc[[2]]
    s3=seqc[[1]]
  }else{#del
    s =seqc[[1]]
    s3=seqc[[2]]
  }
  
  imat = matrix(0,3,2)
  ps1  = start(gz)
  ps2  = end(gz)
  
  rem = ps1 %% 3
  pos = c()      #record the 0,1-edge phase position
  
  #1A-- -AA AAA  2A-- -AA AAA   4AAA AA- --A  3AAA AAA AAA  >phase1
  #1AAA A-- -AA  2AAA AA- --A   4A-- -AA AAA  3A-- -A- --A
  
  #5AA- --A AAA   >phase2
  #5AAA AA- --A  
  
  #6--- NNN NNN  7--- NNN  8NNN NNN  9NNN --- end  10NNN NNN end >phase0 
  #6NNN --- NNN  7NNN NNN  8--- NNN  9NNN NNN end  10NNN --- end
  
  for (i in 1:length(rem)) {
    if(rem[i] == 2){#phase 1 
      if(s3[ps2[i]+2]=='-'){#3
        tmp.0=s[ps1[i]-1]
        tmp.1=s[ps2[i]+1]
        j = ps2[i]+2
        repeat{
          j     = j+1
          if(s3[j]!='-'){
            break
          }
        }
        tmp.2  = s[j]
        pos    = c(pos,ps2[i]+1)
        ref    = paste0(tmp.0,tmp.1,tmp.2,          collapse="")
        sub1   = paste0(tmp.0,s[ps1[i]],s[ps1[i]+1],collapse="")
        sub2   = paste0(s[ps2[i]],tmp.1,s[ps2[i]+2],collapse="")
        sub3   = paste0(s[j-2],s[j-1],tmp.2        ,collapse="")
        sec    = syn[[ref]]
        if(sub1 %in% sec || sub2 %in% sec || sub3 %in% sec){ 
          imat[2,2]=imat[2,2]+1 
        }else{
          imat[2,1]=imat[2,1]+1 
        }
      }else{
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
        }else if(s[ps2[i]+2]=='-'){#2,4
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
          pos  = c(pos,ps2[i]+1)
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
        }else{#normal
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
      }
    }else if(rem[i] == 0){#phase2
      if(s3[ps1[i]-2]=='-'){#3
        repeat{
          j      = j-1
          if(s3[j]!='-'){
            break
          }
        }
        tmp.0  = s[j]
        tmp.1  = s[ps1[i]-1]
        tmp.2  = s[ps2[i]+1]
        pos    = c(pos,ps2[i]+1)
        ref    = paste0(tmp.0,tmp.1,tmp.2,          collapse="")
        sub1   = paste0(tmp.0,s[j+1],s[j+2],        collapse="")
        sub2   = paste0(s[ps1[i]-2],tmp.1,s[ps1[i]],collapse="")
        sub3   = paste0(s[ps2[i]-1],s[ps2[i]],tmp.2,collapse="")
        sec    = syn[[ref]]
        if(sub1 %in% sec || sub2 %in% sec || sub3 %in% sec){ 
          imat[3,2]=imat[3,2]+1 
        }else{
          imat[3,1]=imat[3,1]+1 
        }
      }else{
        if(s[ps1[i]-2]=='-'){#2,4
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
          pos   = c(pos,ps2[i]+1)
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
      }
    }else{#phase0
      if(is.na(s[ps2[i]+1])){#9,10
        pos = 0
      }else{
        if(s[ps2[i]+1]=='-'){#6 
          j = ps2[i]+1
          repeat{
            j    = j+1
            tmp  = s[j]
            if(tmp!='-'){
              break
            }
          }
          pos =c(pos,j)
        }else{#7,8
          pos =c(pos,ps2[i]+1) 
        }
      }
      imat[1,2]=imat[1,2]+1   
    }
  }
  res = list(imat,pos)
  return(res)
}






