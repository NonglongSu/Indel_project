library(expm)
library(stats)
library(Biostrings)

###################function setup
# Set up 4*4 GTR matrix
GTR = function(si, pai){
  r1 = matrix(0, 4, 4)
  r1[lower.tri(r1)] = si
  r  = r1 + t(r1)
  r  = t(r*pai)
  diag(r) = -rowSums(r)
  return(r)
}

# Construct a 64*64 MG94 matrix
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





set.seed(8088)
#######################PART I: parameter space set up
##pi
##2(A+C)=1; highC-G or highA-T case.
nd  = 100
PiA = runif(nd,0.0,0.50)
PiT = PiA
PiC = 0.50 - PiA
PiG = PiC
Pi.all = list(PiA, PiC, PiG, PiT)
Pi.all = sapply(1:nd, function(x){sapply(Pi.all, "[[", x)})
print(dim(Pi.all))
print(colSums(Pi.all)) #1

##omega
##from paper
wv = runif(nd,0.02,0.5) #sim
#wv = runif(nd,0.02,0.3) #real

##if we want a more general range. 
#wv = runif(nd,0.001,1.2)  

##six sigmas
##set up the mean rate first; then b =1/a ; 1-sigma rule
##cv = 1; a = 1/(cv^2)

##distance Tau (total mutations load--only determined by neutral gtr matrix)
##from paper
#Tau = runif(nd,0.05,5.5)  #(put 80% of the simulations in the area of saturated third-positions,seems too unreal!!)
#we setup the a*b=0.4
cv = 1
a  = 1/(cv^2)
b  = 0.4/a
Sigmas = matrix(0,6,100)
tv = c()
for (i in 1:100) {
  pai = Pi.all[,i]
  si  = rgamma(6, shape=a, scale=b)
  gtr = GTR(si, pai)
  Sigmas[,i] = si 
  tv[i]  = -sum(diag(gtr)*pai)
}

mean(tv)
hist(tv, prob = TRUE, xlim = c(0,2))
stv = sort(tv)
lines(stv,dgamma(stv,shape=a,scale=b),col='magenta',lwd=3)
plot(density(tv))


##normalize sigma. 
nSigmas = sapply(1:100, function(x){Sigmas[,x]/tv[x]})
print(nSigmas)

##output the 100 parameters
tPar = matrix(0,100,12)
for (i in 1:nd) {
  tPar[i,1:4]  = Pi.all[,i]
  tPar[i,5:10] = Sigmas[,i]
  tPar[i,11]   = tv[i]
  tPar[i,12]   = wv[i]
}
colnames(tPar)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','tau','omega')
write.table(tPar,"Results/truePar_100.txt",quote=F,sep="\t",
            row.names=F)



#######################PART II: distance method for init value setup. 
##P204~210 (Inferring phylogenies 2004)
#assume true parameters, unnormalized
Sigma = 1:6
Pi    = c(0.2,0.3,0.3,0.2)
Tau   = 0.5

#setup matrix and normalize
rmat = matrix(0,4,4)
rmat[lower.tri(rmat)] = Sigma
rmat = rmat + t(rmat)
rmat = t(rmat * Pi)
diag(rmat) = -rowSums(rmat)
T    = -sum(diag(rmat)*Pi)                   
rmat = rmat/T

#create Probability matrix
pmat = expm(rmat*Tau) 
Pmat = pmat * Pi 
print(sum(Pmat))

#generate simualted data
dat = sample(16, 30000, replace=TRUE, prob=Pmat)   
dat = table(dat)
dat = matrix(dat,4,4)

f = rowSums(dat) + colSums(dat)
f = f/sum(f)
print(sum(f))

#symmetric average of dat
sdat = matrix(0,4,4)
diag(sdat) = diag(dat)
for (i in 1:4) {
  for (j in 1:4) {
    if((i!=j) && (i<j)){
      sdat[i,j] = sdat[j,i] = (dat[i,j]+dat[j,i])/2
      }
  }
}
sf = colSums(sdat)/sum(sdat)
print(sum(sf))
Dhat  = diag(sf)
phat  = t(sdat/colSums(sdat))
Athat = logm(phat)
that  = -sum(diag(Athat*Dhat))
print(that)   #pretty good

Ahat = t(Athat)/that
print(Ahat)
print(rmat)   #compared with rmat: not too bad. 

shat=t(Ahat)/sf
print(shat)
print(t(rmat)/Pi)



###################distance method partII "use 4-dengeneracy sites"
#stop codon position
stp = c(49,51,57)
#define mg94 dimensions
m = 61

#Construct codons and its degeneracy
codons = cbind(rep(DNA_BASES, each = 16),
               rep(DNA_BASES, times = 4, each = 4),
               rep(DNA_BASES, 16))
codons = codons[-stp,]
codonstrs <<- apply(codons, 1, stringr::str_c, collapse = "")            
syn        = syncodons(codonstrs)
names(syn) = toupper(names(syn))
syn        <<- lapply(syn, toupper)

s    = Sigmas[,1]
Pi   = Pi.all[,1]
omega= wv[1]
tau  = tv[1]

gtr  = GTR(s,Pi) 
gtr  = gtr/tau
mg94 = MG94(gtr, omega, codons, m)


#codon freq
Pi2 = sapply(seq(m), function(x){prod(Pi[match(codons[x, ], DNA_BASES)])})  
Pi2 = Pi2/sum(Pi2)
print(sum(Pi2))

#Create Symmetric matrix
o   = outer(sqrt(Pi2), 1/sqrt(Pi2))      
s94 = mg94*o
print(s94[1,])
print(s94[,1])

#Create Probability matrix
p94 = expm(s94*tau)* t(o)
P94 = p94*Pi2        
print(rowSums(p94))              
print(sum(P94))

#Normlize sigma as well
sigma = s/tau

dat  = sample(m*m, 10000, replace=TRUE, prob=P94)     
dat  = table(dat)
dat  = as.data.frame(dat)
id1  = as.numeric(as.vector(dat[[1]]))
id2  = as.numeric(as.vector(dat[[2]]))
dat1 = matrix(0, m, m)
for (i in 1:length(id1)) {
  dat1[id1[i]] = id2[i]
}
print(diag(dat1))


#create nuc transition matrix [4-fold-degeneracy-codons only]
fourD_id = which(sapply(syn, function(x){length(x) == 4}))
f1.id    = sapply(seq(m), function(x){match(codons[x,], DNA_BASES)})
k=1
ndat = matrix(0,4,4)
while (k<length(fourD_id)) {
  i = j = fourD_id[k:(k+3)]
  ndat = ndat+dat1[i,j]
  k=k+4
}

f = rowSums(ndat) + colSums(ndat)
f = f/sum(f)
print(sum(f))

#distance method for guessing simgas and tau.
#symmetric average of 4*4 matrix
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
print(gtr)   #compared with rmat: not too bad. 

shat=t(Ahat)/sf
print(shat)
print(t(gtr)/Pi)







#########################PART III: define a reasonable dna sample size. 
#set up a gradient 1e+5,1e+6,1e+7.
#each did with 100 simulations.
#plot average accuracy of est of each sample size.

#########################PART IV:  Paired t-test of theta: sample size-100. 
theta1 = rnorm(nd,5,1)
theta2 = rnorm(nd,3,1)
t.test(theta1, theta2, paired = TRUE)


###################PART V:why GY94 is not okay. 
#1) GY style existed some scenario that going against mutational bias: CGC->CTC | ATA->AGA.
#2) it won't put the mcmc approaching the stationary distribution.

######################PART VI:
##print the LL after every E-step






#$$$$$$$$$$$$$$$$$$random walk (metropolis hastings algo)
##not too bad, very fast!
nd      = 1000
wd      = rep(1,nd)
sd      = runif(6)
sigma   = .1

r    = GTR(sd[[1]], f)
R    = MG94(r, w, codons)
Pabt = expm(R)
print(rowSums(expm(R)))
LLold=sum(log(Pabt*cf) * dat1)

print(LLold)                ##init LL
print(sum(log(P94)*dat1))   ##true LL.


#update omega
for (t in 2:nd) {
  repeat{
    w = wd[t-1] + sigma*rnorm(1)
    if(w>0){#parameters has to be positive.
      break
    }
  }
  r    = GTR(s, f)
  R    = MG94(r, w, codons)
  Pab  = expm(R)
  LL   = sum(log(Pab*cf) * dat1)
  rat  = exp(LL-LLold)
  
  #print(LL)
  #print(rat)
  
  alpha = min(1, rat)
  if(rbinom(1,1,alpha)){
    wd[t]   = w
    LLold   = LL
  }else{
    wd[t]   = wd[t-1]
  }
  print(wd[t])
}

mean(wd)
plot(wd, type = "b")  ##plot mcmc

#update theta
w  = mean(wd)
sd = matrix(runif(6), 6, nd, byrow=FALSE)
td = rep(1,nd)
r    = GTR(sd[,1], f)
R    = MG94(r, w, codons)
Pabt = expm(R)
print(rowSums(expm(R)))
LLold=sum(log(Pabt*cf) * dat1)
print(LLold)


  for (t in 2:nd) {
    sd[,t] = sd[,t-1]  #copy t-1 to t
    for (i in 1:6) {
      repeat{
        si = sd[i,t-1] + sigma*rnorm(1)
        if((si>0)){
          s      = sd[,t]
          s[i]   = si
          break
        }
      }
      r    = GTR(s, f)
      R    = MG94(r, w, codons)
      Pab  = expm(R)
      LL   = sum(log(Pab*cf) * dat1)
      rat  = exp(LL-LLold)
      print(LL)
      #print(rat)
      
      alpha = min(1, rat)
      if(rbinom(1,1,alpha)){
        sd[i,t] = s[i]
        LLold   = LL
      }
      #print(sd[i,t])
    }
    td[t] = -sum(diag(R)*cf)
  }

that = mean(td)
print(that)
mean(sd)
shat = sapply(1:6, function(x){mean(sd[x,])})
print(shat/that)
print(Sigma)    #true Sigma


##plot mcmc
plot(td, type = "b")
plot(sd[1,], type = "b")


##################################
#theta = {6-sigmas,omega}
#update one fix the rest
#repeat the loop
thetad = matrix(runif(7), 7, nd, byrow=FALSE)
td     = rep(1,nd)
s0     = thetad[1:6,1]
w0     = thetad[7,1]
r      = GTR(s0, f)
R    = MG94(r, w0, codons)
Pabt = expm(R)
print(rowSums(expm(R)))
LLold=sum(log(Pabt*cf) * dat1)

print(LLold)                ##init LL
print(sum(log(P94)*dat1))   ##true LL.


for (t in 2:nd) {
  thetad[,t] = thetad[,t-1]  #copy [t-1]state to t state
  for (i in 1:7) {
    repeat{
      ti = thetad[i,t-1] + sigma*rnorm(1)
      if((ti>0)){#parameters has to be positive
        theta    = thetad[,t]
        theta[i] = ti
        break
      }
    }
    s    = theta[1:6]
    w    = theta[7]
    r    = GTR(s, f)
    R    = MG94(r, w, codons)
    Pab  = expm(R)
    LL   = sum(log(Pab*cf) * dat1)
    rat  = exp(LL-LLold)
    print(LL)
    
    alpha = min(1, rat)
    if(rbinom(1, 1, alpha)){
      thetad[i,t] = theta[i]
      LLold       = LL
    }
    #print(thetad[i,t])
  }
  td[t] = -sum(diag(R)*cf)    #tau for every column.
}

that = mean(td)
print(that)
sd   = thetad[1:6,]
shat = sapply(1:6, function(x){mean(sd[x,])})
print(shat/that)
print(Sigma)    #true Sigma

wd   = thetad[7,]
what = mean(wd)
print(what)
print(omega)    #true omega


##plot mcmc
par(mfrow=c(1,3))
plot(td, type="b")               #marginal draws of tau
plot(sd[1,], type="b")           #marginal draws of sigma1 
plot(wd, type="b")               #marginal draws of omega

acf(td,col="blue")               #acf draws of tau
acf(sd[1,],col="blue")           #acf draws of sigma1 
acf(td,col="blue")               #acf draws of omega


###################Now, let's think about this problem in a different way: 
##let's give a prior 
#assuming theta beta(1,1) which is uniform distribution. 
#since the prior is conjugate, we can easily write down the posterior Beta(alpha+y, beta+n-y)
#I claim y(# of success) is the # of times passing the min(1,rat) in mh. 
yd = rep(0,7) 
for (t in 2:nd) {
  thetad[,t] = thetad[,t-1]  #copy [t-1]state to t state
  for (i in 1:7) {
    repeat{
      ti = thetad[i,t-1] + sigma*rnorm(1)
      if((ti>0)){#parameters has to be positive
        theta    = thetad[,t]
        theta[i] = ti
        break
      }
    }
    s    = theta[1:6]
    w    = theta[7]
    r    = GTR(s, f)
    R    = MG94(r, w, codons)
    Pab  = expm(R)
    LL   = sum(log(Pab*cf) * dat1)
    rat  = exp(LL-LLold)
    print(LL)
    
    alpha = min(1, rat)
    if(rbinom(1, 1, alpha)){
      thetad[i,t] = theta[i]
      LLold       = LL
      yd[i]       = yd[i] + 1  #>>>>>>add this
    }
    #print(thetad[i,t])
  }
  td[t] = -sum(diag(R)*cf)    #tau for every column.
}

ta   = list()
tahat= c()
for (i in 1:7) {
  ta[[i]]  = rbeta(nd, 1+yd[i], 1+nd-yd[i])
  tahat[i] = mean(ta[[i]])
}
print(tahat[1:6]/that)         #bad haha. 
print(Sigma)
print(tahat[7])
print(omega)

hist(ta[[1]], prob=TRUE, xlim = c(0,0.15))
abline(h=1, col="red")
sta1 = sort(ta[[1]])
lines(sta1,dbeta(sta1,1+yd[1],1+nd-yd[1]),col='magenta',lwd=3)

###it tells me the combination of prior is a problem. 


