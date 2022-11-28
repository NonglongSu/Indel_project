m  = matrix(0,2,2)
pi = c(1/3,2/3)
names(pi)=c('M','D')

m[1,1] = 0.5
m[1,2] = 0.5
m[2,1] = 0.25
m[2,2] = 0.75
rownames(m) = c('M','D')
colnames(m) = c('M','D')


y1 = pi %*% m  # ~pi
y1

#pi = pi*m^t
pi %*% m %*% m %*% m

# if pi is stationary, Pt -> pi as t -> inf
m %*% m %*% m %*% m %*% m %*% m

############################################
library(MASS)
#Note that the equation π*P=π implies that the vector π is a left eigenvector of P with eigenvalue equal to 1
r = eigen(m)
rvec = r$vectors
#left eigenvector are the inverse of the right eigenvectors
lvec = ginv(r$vectors)
#The eigenvalues
lam = r$values
#check the spectral decomposition
rvec %*% diag(lam) %*% lvec

#We see the first eigenvalue is 1 and so the first left eigenvector, suitably normalized, should contain the stationary distribution:
pi_eig = lvec[1,]/sum(lvec[1,])
pi_eig
sum(pi_eig)



