
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(expm))

##Our EM method based on doi:10.1006/jmbi.2001.5405 Holmes&Rubin
################################# PART I Function set up
#Set up 4*4 GTR matrix
GTR = function(Si, pai){
  r1 = matrix(0, 4, 4)
  r1[lower.tri(r1)] = Si
  r  = r1 + t(r1)
  r  = t(r*pai)
  diag(r) = -rowSums(r)
  return(r)
}

#Construct a 64*64 MG94 matrix
MG94 = function(g, oga, cd, m){
  R  = matrix(0, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
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

#Return LL
LL = function(theta, f, cf, Dat, m){
  r  = GTR(theta[1:6], f)
  R  = MG94(r, theta[7], codons, m)
  ll = sum(log(expm(R)*cf)*Dat)
  return(ll)
}


####################################Part II Model building
dEM = function(Dat,f0,cf0,rms){
  
  #stop codon position
  stp = c(49,51,57)
  
  #Construct codons and its degeneracy
  codons = cbind(rep(DNA_BASES, each = 16),
                 rep(DNA_BASES, times = 4, each = 4),
                 rep(DNA_BASES, 16))
  codons = codons[-stp,]
  codonstrs  <<- apply(codons, 1, stringr::str_c, collapse = "")
  syn        = syncodons(codonstrs)
  names(syn) = toupper(names(syn))
  syn        <<- lapply(syn, toupper)
  
  ###################################Part III simulate data and run EM.
  #define mg94 dim
  m=61
  
  ##Determine init value
  #distance method for guessing simgas and tau
  #create nuc transition matrix [4-fold-degeneracy-codons only]
  fourD_id = which(sapply(syn, function(x){length(x) == 4}))
  f1.id    = sapply(seq(m), function(x){match(codons[x,], DNA_BASES)})
  k=1
  ndat = matrix(0,4,4)
  while (k<length(fourD_id)) {
    i = j = fourD_id[k:(k+3)]
    ndat = ndat+Dat[i,j]
    k=k+4
  }
  fn = rowSums(ndat) + colSums(ndat)        #neutral freq
  fn = fn/sum(fn)
  print(sum(fn))
  
  #symmetric average of dat
  sdat = matrix(0,4,4)
  diag(sdat) = diag(ndat)
  for (i in 1:4) {
    for (j in 1:4) {
      if((i!=j) && (i<j)){
        sdat[i,j] = sdat[j,i] = (ndat[i,j]+ndat[j,i])/2
      }
    }
  }
  sf = colSums(sdat)/sum(sdat)
  print(sum(sf))
  Dhat  = diag(sf)
  phat  = t(sdat/colSums(sdat))
  Athat = logm(phat)
  that  = -sum(diag(Athat*Dhat))
  print(that)      
  
  Ahat = t(Athat)/that
  print(Ahat)
  
  s = Ahat[lower.tri(t(Ahat))]
  w  = 0.20  #smart guess, I guess!   
  
  p0 = c(s,w)
  r  = GTR(s,f)
  R  = MG94(r, w, codons, m)
  print(sum(R) %>% round(13))
  
  

  #Prepare omega list for M-step.
  ##Locating all non-syn locations in 64*64 R matrix.
  omega.id = c()
  for (i in 1:m) {
    for (j in 1:m) {
      if((i != j) && 
         (sum(codons[i, ] != codons[j, ]) == 1) && 
         (!(codonstrs[j] %in% syn[[codonstrs[i]]])) ){
        omega.id = c(omega.id, i, j)          
      }
    }
  }
  omega.id = split(omega.id, ceiling(seq_along(omega.id) / 2))
  
  #Prepare sigma list for M-step
  ##Locating all 6-sigma locations in 64*64 R matrix.
  ii = GTR(1:6, rep(1,4))
  diag(ii) = 0
  I = matrix(0, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
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
  sigma.id = sapply(1:6, function(x){which(I == x)})
  
  
  #>>>EM
  llv = c()
  
  for(counter in 1:100) {

    #Simulate branch length
    T1 = -sum(diag(r)*f)
    print(T1)
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

    ##-Log likelihood
    ll = sum(log(Pab*cf) * Dat)
    llv[counter] = ll
    cat(sprintf("%i: %.6f\n", counter, ll))
    if(ll>=ll.true){
      break
    }


    #PART II
    ##construct the Jkl matrix
    J = outer(D/T1, D/T1, function(x,y) {
      ifelse(x-y == 0,
             T1*exp(x*T1),
             exp(y*T1)*(expm1((x-y)*T1))/(x-y))
    })

    ##calculate the expected values
    # W[a,b,i,i] is the expected time spent state i on a branch going from a -> b
    # U[a,b,i,j] is the expected number of events going from i -> j on a branch going from a->b
    W = array(0, c(m,m,m,m))
    U = array(0, c(m,m,m,m))

    tm = system.time(
      for(a in 1:m) {
        for(b in 1:m) {
          for(i in 1:m) {
            for(j in 1:m) {
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
    Wh = array(0, c(m,m))
    Uh = array(0, c(m,m))
    for(i in 1:m) {
      for(j in 1:m) {
        Wh[i,j] = sum(W[,,i,j] * Dat)
        Uh[i,j] = sum(U[,,i,j] * Dat)
      }
    }


    ##M-Step maximize sigmas.
    Wh = diag(Wh)
    sigma.Cij = sapply(seq(6), function(x){sum(Uh[sigma.id[[x]]])})
    sigma.Wij = c()
    for (k in 1:6) {
      ichunks   = ceiling(sigma.id[[k]]/m)
      sigma.Wij = c(sigma.Wij, sum(Wh[ichunks]* t(R)[sigma.id[[k]]])/s[k])
    }
    s = sigma.Cij/sigma.Wij


    ##M-Step maximize omega
    w.Cij = sum(sapply(omega.id, function(x){Uh[x[1], x[2]]}))
    Rii   = c()
    for (i in 1:m) {
      ith.non = omega.id[which(sapply(omega.id, "[[", 1) == i)]
      ith.sum = sum(sapply(ith.non, function(x){R[x[1], x[2]]})) / w
      Rii     = c(Rii, ith.sum)
    }
    w.Wij = sum(Wh * Rii)
    w     = w.Cij/w.Wij

    ##reconstruct gtr and mg94.
    pNew = c(s,w)
    print(pNew)
    r    = GTR(pNew[1:6], f)
    R    = MG94(r, pNew[7], codons, m)
    
    ##>>>>>
    if(counter>1){#aitken's method 
      theta.c = pNew
      theta.d = abs((theta.a*theta.c-theta.b^2)/(theta.a+theta.c-2*theta.b))
      lld     = LL(theta.d, f, cf, Dat, m)
      if(lld < lla){#overshooting
        theta.d = theta.c
      }
      #root-mean-square
      delta = (theta.d-theta.a)/theta.a
      rms   = sqrt(mean(delta*delta))
      print(rms)
      #if(rms < tol){                         
      # if(ll>=ll.true){
      #   break
      # }
      theta.a = theta.b
      theta.b = theta.c
      lla     = ll
    }else{#catch 1st round
      theta.a = p0
      theta.b = pNew
      lla     = ll
    }
  }
  
  
  ##################Print EM results
  ##>>
  # T1=0.2
  # s = c(1,2,3,4,5,6)
  # llv = rep(10,10)
  
  
  # Estimate T and normalized Sigmas
  print(T1)     #est of Tau.
  print(s/T1)   #est of normalized sigmas.
  print(w)      #est of omega
  
  est = c(s/T1,T1,w)
  res = list(est,llv)
  
 return(res)
}
