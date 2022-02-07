rm(list=ls())
library(dplyr)

## Bridge network - the length of shortest path ##

set.seed(1)
H = function(X) {
  path <- matrix(0, nrow=nrow(X), ncol=5)
  path[,1] <- X[,1]+X[,5]+X[,8]
  path[,2] <- X[,2]+X[,6]+X[,8]
  path[,3] <- X[,2]+X[,7]
  path[,4] <- X[,3]+X[,8]
  path[,5] <- X[,4]
  apply(path, 1, min)
}

N = 10000
u = c(1, 1, 0.5, 0.5, 2, 2, 1.5, 1.5) # X_i ~ Exp(u_i)



## 1. CMC estimator
U = matrix(u, nrow = N, ncol = length(u), byrow = T)
X = -log(runif(N*length(u))) * U
HX = H(X) 
CMC = HX %>% mean
var.CMC = var(HX)
RE.CMC = sd(HX)/CMC/sqrt(N)



## 2. antithetic variable method
U1 = matrix(u, nrow = N/2, ncol = length(u), byrow = T)
U2 = matrix(u, nrow = N/2, ncol = length(u), byrow = T)
u = runif(N/2*length(u))
X1 = -log(u) * U1
X2 = -log(1-u) * U2
# X1 = matrix(rexp(N * X.length / 2), nrow = N/2, ncol=X.length)
# X2 = matrix(rexp(N * X.length / 2), nrow = N/2, ncol=X.length)

HX1 = H(X1); HX2 = H(X2)
AV = c(HX1, HX2) %>% mean

var.AV = var(HX1) + var(HX2) + 2*cov(HX1, HX2)
RE.AV = sqrt(var.AV) / AV /sqrt(N)




## 3. control variable method 
# U = matrix(u, nrow = N, ncol = length(u), byrow = T)
# X = -log(runif(N*length(u))) * U
X
C1 = X[,1]+X[,5]+X[,8]
C2 = X[,2]+X[,6]+X[,8]
C3 = X[,2]+X[,7]
C4 = X[,3]+X[,8]
C5 = X[,4]

C = matrix(c(C1, C2, C3, C4, C5), ncol = 5)
Sigma.C = cov(C)


Sigma.XC = c()
for(i in 1:ncol(C)) Sigma.XC[i] = cov(HX, C[,i])
alpha = solve(Sigma.C) %*% Sigma.XC # minimizes Var(X.alpha)

w = matrix(c(1,0,0,0,1,0,0,1,
             0,1,0,0,0,1,0,1,
             0,1,0,0,0,0,1,0,
             0,0,1,0,0,0,0,1,
             0,0,0,1,0,0,0,0),
           nrow = 5, ncol = 8, byrow=T)
r = matrix(w %*% (1/u), nrow = N, ncol = 5, byrow = T) # the expectation matrix of C

HX.star = HX - t(alpha) %*% t(C-r)

CV = mean(HX.star)
var.CV = var(HX) - t(Sigma.XC) %*% solve(Sigma.C) %*% Sigma.XC
RE.CV = sqrt(var.CV)/sqrt(N)/CV




## 4. CE method 
u = c(1, 1, 0.5, 0.5, 2, 2, 1.5, 1.5)

# estimate CE optimal v
# U = matrix(u, nrow = N, ncol = length(u), byrow = T)
Y = -log(runif(N*length(u))) * U # random number of Exp(u)
HY = H(Y) * Y
v = apply(HY, 2, mean) / mean(H(Y))

# actual estimation

# likelihood
f = function(X, v) {
  V = matrix(v, nrow = N, ncol = length(v), byrow = T)
  pdf = sapply(1:nrow(X), function(i) exp(- X[i,] / V[i,]) / V[i,])
  apply(pdf, 2, prod)
}

V = matrix(v, nrow = N, ncol = length(u), byrow = T)
X = -log(runif(N*length(v))) * V  # random number of Exp(v)
W = f(X,u) / f(X,v)
HW = H(X) * W
CE = mean(HW)
var.CE = var(HW)
RE.CE = sd(HW)/sqrt(N)/CE


## 5. Table
result = data.frame(estimator = c("CMC", "AV", "CV", "CE"),
                    estimate = c(CMC, AV, CV, CE),
                    Var = c(var.CMC, var.AV, var.CV, var.CE),
                    RE = c(RE.CMC, RE.AV, RE.CV, RE.CE))
result[,-1] %>% round(6)





