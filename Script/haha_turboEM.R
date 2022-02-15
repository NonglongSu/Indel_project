library(Biostrings)
library(stringr)
library(seqinr)
library(plyr)
library(purrr)
library(R.utils)
library(Matrix)
library(matlib)
library(expm)
library(doParallel)
library(turboEM)

##Our EM method based on doi:10.1006/jmbi.2001.5405 Holmes&Rubin

################################# PART I Function set up
# Set up 4*4 GTR matrix
GTR = function(Si, pai){
  r1 = matrix(0, 4, 4)
  r1[lower.tri(r1)] = Si
  r  = r1 + t(r1)
  r  = t(r*pai)
  diag(r) = -rowSums(r)
  return(r)
}

# Construct a 64*64 MG94 matrix
MG94 = function(g, oga, cd){
  R  = matrix(0, 64, 64)
  for (i in 1:64) {
    for (j in 1:64) {
      if(i == j){
        R[i, j] = 0
      }else if(sum(cd[i, ] != cd[j, ]) > 1){#more than 1 nucleotide changes
        R[i, j] = 0
      }else{
        if(codonstrs[j] %in% syn[[codonstrs[i]]]){#syn_subs
          w = 1
        }else{#nonsyn_subs
          w = oga
        }
        pos = which(cd[i, ] != cd[j, ])
        x   = which(DNA_BASES == cd[i, pos])
        y   = which(DNA_BASES == cd[j, pos])
        R[i, j] = w * g[x, y]
      }
    }
  }
  diag(R) = -rowSums(R)
  
  return(R)
}


#############################################################Part II Model building

# Construct codons and its degeneracy
codons = cbind(rep(DNA_BASES, each = 16),
               rep(DNA_BASES, times = 4, each = 4),
               rep(DNA_BASES, 16))
codonstrs  <<- apply(codons, 1, stringr::str_c, collapse = "")
syn        = syncodons(codonstrs)
names(syn) = toupper(names(syn))
syn        <<- lapply(syn, toupper)


# True parameters, unnormalized
Pi    = c(0.2, 0.3, 0.3, 0.2) 
Sigma = 1:6
Tau   = 0.1
omega = 0.2                    ###>>>

# Set up GTR matrix
gtr = GTR(Sigma, Pi)
print(sum(gtr))

# Set up mg94 matrix and normalize it.
Pi2  = sapply(seq(64), function(x){prod(Pi[match(codons[x, ], DNA_BASES)])})  
mg94 = MG94(gtr, omega, codons)
T    = -sum(diag(mg94)*Pi2)
mg94 = mg94 / T

# Create Symmetric matrix
o   = outer(sqrt(Pi2), 1/sqrt(Pi2))      
s94 = mg94 * o
print(s94[1, ]) 
print(s94[, 1])

# Create Probability matrix
p94 = expm(s94*Tau)* t(o)
P94 = p94* Pi2        
print(rowSums(p94))              
print(sum(P94))

#Normlize sigma as well
Sigma = Sigma/T


#Build omega list for M-step.
##Locating all non-syn locations in 64*64 R matrix.
ome.id = c()
for (i in 1:64) {
  for (j in 1:64) {
    if((i != j) && 
       (sum(codons[i, ] != codons[j, ]) == 1) && 
       (!(codonstrs[j] %in% syn[[codonstrs[i]]])) ){
      ome.id = c(ome.id, i, j)          
    }
  }
}
omega.id <<- split(ome.id, ceiling(seq_along(ome.id) / 2))

#Build sigma list for M-step
##Locating all 6-sigma locations in 64*64 R matrix.
ii = GTR(1:6, rep(1,4))
diag(ii) = 0
I = matrix(0, 64, 64)
for (i in 1:64) {
  for (j in 1:64) {
    if(i == j){
      I[i, j] = 0
    }else if(sum(codons[i, ] != codons[j, ]) > 1){
      I[i, j] = 0
    }else{
      pos = which(codons[i, ] != codons[j, ])
      x   = which(DNA_BASES == codons[i, pos])
      y   = which(DNA_BASES == codons[j, pos])
      I[i, j] = ii[x, y]
    }
  }
}
sigma.id <<- sapply(1:6, function(x){which(I == x)})


##########################################################Part III simulate data and run EM.
# Create sample data
set.seed(8088)
dat  = sample(64*64, 10000, replace = TRUE, prob = P94)
dat  = table(dat)
dat  = as.data.frame(dat)
id1  = as.numeric(as.vector(dat[[1]]))
id2  = as.numeric(as.vector(dat[[2]]))
dat1 = matrix(0, 64, 64)
for (i in 1:length(id1)) {
  dat1[id1[i]] = id2[i]
}

###>>>>>>>use expected freq 
dat1 = 10000 * P94

##initial parameters
f1     = rowSums(dat1) + colSums(dat1)  
f1.id  = sapply(seq(64), function(x){match(codons[x, ], DNA_BASES)})
base.f = c()
for (i in seq(4)) {
  base.id = sapply(seq(64), function(x){length(which(f1.id[, x] == i))})
  base.f  = c(base.f, sum(base.id * f1))
}
f  <<- base.f/sum(base.f)
cf <<- sapply(seq(64), function(x){prod(f[match(codons[x, ], DNA_BASES)])})
print(sum(f)) 
print(sum(cf))



#Construnct EM func.
emEst = function(p, y){
  #input
  s = p[1:6]
  w = p[7]
  dat1 = y
  
  #construct R
  r  = GTR(s, f)
  R  = MG94(r, w, codons)
  #print(sum(R) %>% round(13))
  
  #Simulate branch length
  T1   = -sum(diag(R)*cf)       
  
  ##make r, R symmetric
  fmat = outer(sqrt(cf), 1/sqrt(cf))
  S    = R * fmat
  
  ##calculate eigenvectors and values.
  eig = eigen(S)
  D   = eig$values
  V   = eig$vectors
  
  ##calculate Prob(b|a)
  pab = V %*% diag(exp(D)) %*% t(V)           
  Pab = pab * t(fmat)
  #print(rowSums(Pab))
  
  # PART II
  ##construct the Jkl matrix       
  J = outer(D/T1, D/T1, function(x,y) {
    ifelse(x-y == 0,
           T1*exp(x*T1),
           exp(y*T1)*(expm1((x-y)*T1))/(x-y))
  })
  
  ##calculate the expected values
  #W[a,b,i,i] is the expected time spent state i on a branch going from a -> b
  #U[a,b,i,j] is the expected number of events going from i -> j on a branch going from a->b
  W = array(0, c(64,64,64,64))     
  U = array(0, c(64,64,64,64))
  
  tm = system.time(
    for(a in 1:64) {
      for(b in 1:64) {
        for(i in 1:64) {
          for(j in 1:64) {
            ff = sqrt(cf[i]*cf[b]/cf[a]/cf[j])
            o  = outer(V[a,]*V[i,], V[j,]*V[b,])   
            W[a,b,i,j] = ff * sum(o*J)
          }
        }
        W[a,b,,] = W[a,b,,] / Pab[a,b]
        U[a,b,,] = R * W[a,b,,]
      }
    }
  ) 
  
  ##calculate expected values by summing over observations --a,b is sumable. 
  Wh = array(0, c(64,64))
  Uh = array(0, c(64,64))
  for(i in 1:64) {
    for(j in 1:64) {
      Wh[i,j] = sum(W[,,i,j] * dat1)
      Uh[i,j] = sum(U[,,i,j] * dat1)
    }
  }
  
  ##M-Step maximize sigmas.
  Wh = diag(Wh)
  sigma.Cij = sapply(seq(6), function(x){sum(Uh[sigma.id[,x]])})
  sigma.Wij = c()
  for (k in 1:6) {
    ichunks   = ceiling(sigma.id[, k] / 64)
    sigma.Wij = c(sigma.Wij, sum(Wh[ichunks]* t(R)[sigma.id[, k]])/s[k] )  
  }
  s = sigma.Cij/sigma.Wij
  
  
  ##M-Step maximize omega
  w.Cij = sum(sapply(omega.id, function(x){Uh[x[1], x[2]]}))
  Rii   = c()
  for (i in 1:64) {
    ith.non = omega.id[which(sapply(omega.id, "[[", 1) == i)]
    ith.sum = sum(sapply(ith.non, function(x){R[x[1], x[2]]})) / w    
    Rii     = c(Rii, ith.sum)
  }
  w.Wij = sum(Wh * Rii)
  w     = w.Cij/w.Wij        
  
  pNew = c(s, w, T1)
  p    = pNew
  return(pNew)
}

#Construnct LL func.
LL = function(p, y){
  #input
  s = p[1:6]
  w = p[7]
  dat1 = y
  
  #construct R
  r  = GTR(s, f)
  R  = MG94(r, w, codons)
  
  ##make r, R symmetric
  fmat = outer(sqrt(cf), 1/sqrt(cf))
  S    = R * fmat
  
  ##calculate eigenvectors and values.
  eig = eigen(S)
  D   = eig$values
  V   = eig$vectors
  
  ##calculate Prob(b|a)
  pab = V %*% diag(exp(D)) %*% t(V)           
  Pab = pab * t(fmat)
  
  ##Log likelihood "the smaller, the better"
  ll = sum(log(Pab*cf) * dat1)
  #cat(sprintf("%i: %.6f\n", counter, ll))
  return(-ll)
}


##>>>>>>>>>>>>>>running turboEM
p0 =  runif(8)    #starting value c(s,w,tau)


##decme method
boundary = function(par, dr){
  lower = rep(0, length(par))
  upper = rep(1, length(par))
  low1 = max(pmin((lower-par)/dr, (upper-par)/dr))
  upp1 = min(pmax((lower-par)/dr, (upper-par)/dr))
  return(c(low1, upp1))
}

tm4 = system.time({
  res2 = turboem(par = p0, y = dat1, fixptfn = emEst, objfn = LL,
                 boundary = boundary, method = "decme")
})

print(res2$convergence)
est2 = res2$pars
print(est2[8])             #tau
print(est2[7])             #omega
print(est2[1:6]/est2[8])   #sigmas


























