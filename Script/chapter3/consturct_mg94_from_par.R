library(stringr)
library(Matrix)

sub1 = "../Script/sources/"
source(paste0(sub1,"gtr.R"))
source(paste0(sub1,"mg94.R"))
#Input: assume you have one line of 11 parameters
# A C G T s1 s2 s3 s4 s5 s6 w

#Output: a P94 matrix 

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#setup codon table
cod = cbind(rep(DNA_BASES, each = 16),
              rep(DNA_BASES, times = 4, each = 4),
              rep(DNA_BASES, 16))
codonstrs  = apply(cod, 1, stringr::str_c, collapse = "")
syn        = syncodons(codonstrs)
names(syn) = toupper(names(syn))
syn        = lapply(syn, toupper)


#True parameters, unnormalized
Pi    = c(0.308,0.185,0.199,0.308)         
Sigma = c(0.009489730, 0.039164824, 0.004318182, 0.015438693, 0.038734091, 0.008550000)
brlen = 0.0133
omega = 0.2

#Set up GTR matrix
gtr = GTR(Sigma, Pi)
print(-sum(diag(gtr)*Pi))

T=-sum(diag(gtr)*Pi)
gtr=gtr/T

#Set up mg94 matrix and obtain the codon freq. 
mg94 = MG94(gtr, omega, cod, codonstrs, syn)
Pi2  = sapply(seq(64), function(x){prod(Pi[match(cod[x, ], DNA_BASES)])})
print(sum(Pi2))

#set up P matrix
p94 = expm(mg94*brlen)             ##This is what reed did in the end
P94 = p94*Pi2
print(sum(P94))

###GTR matrix func
GTR = function(Si, pai){
  r1 = matrix(0, 4, 4)
  r1[lower.tri(r1)] = Si
  r  = r1 + t(r1)
  r  = t(r*pai)
  diag(r) = -rowSums(r)
  return(r)
}

###MG94 matrix func
MG94 = function(g, oga, cod, codonstrs, syn, num=64){
  R  = matrix(0, num, num)
  for (i in 1:num) {
    for (j in 1:num) {
      if(i == j){
        R[i, j] = 0
      }else if(sum(cod[i, ] != cod[j, ]) > 1){#more than 1 nucleotide changes
        R[i, j] = 0
      }else{
        if(codonstrs[j] %in% syn[[codonstrs[i]]]){#syn_subs
          w = 1
        }else{#nonsyn_subs
          w = oga
        }
        pos = which(cod[i, ] != cod[j, ])
        x   = which(DNA_BASES == cod[i, pos])
        y   = which(DNA_BASES == cod[j, pos])
        R[i, j] = w * g[x, y]
      }
    }
  }
  diag(R) = -rowSums(R)
  
  return(R)
}
