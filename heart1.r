library(copent) # Copula Entropy
library(energy) # Distance Correlation
library(dHSIC) # Hilbert-Schmidt Independence Criterion
## for additional tests
library(HHG) # Heller-Heller-Gorfine Tests of Independence
library(independence) # Hoeffding's D test or Bergsma-Dassios T* sign covariance
library(Ball) # Ball correlation
library(qad) # Quantification of Asymmetric Dependence
library(BET) # Binary Expansion Testing
library(MixedIndTests) # Cramer-von Mises statistics
library(NNS) # Nonlinear Nonparametric Statistics

scan_heart_data <-function(filename1, nl = 0){
  data1 = scan(filename1, nlines = nl, what = c(as.list(rep(0,75)),list("")))
  l = length(data1[[1]])
  data1m = matrix(unlist(data1), l, 76)
  matrix(as.numeric(data1m[,1:75]), l, 75)
}
#### load heart disease data (899 samples)
h1 = scan_heart_data("http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/cleveland.data", 282*10)
h2 = scan_heart_data("http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/hungarian.data")
h3 = scan_heart_data("http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/switzerland.data")
h4 = scan_heart_data("http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/long-beach-va.data")

heart1 = as.matrix( rbind(h1,h2,h3,h4) )
m = dim(heart1)[1]
n = dim(heart1)[2]

## statistical dependence with attr #58
# Copula Entropy
l = 50
ce58 = rep(0,n)
for (i in 1:n){
  for (j in 1:l){
    data2 = heart1[,c(i,58)]
    data2[,1] = data2[,1] + max(abs(data2[,1])) * 0.000001 * rnorm(m)
    data2[,2] = data2[,2] + max(abs(data2[,2])) * 0.000001 * rnorm(m)
    ce58[i] = ce58[i] + copent(data2)
  }
}
ce58 = ce58 / l
ce58[c(1,2,58)] = min(ce58)

# Other Independence Measures
dcor58 = rep(0,n) # Distance Correlation
dhsic58 = rep(0,n)  # Hilbert-Schmidt Independence Criterion
hhg58 = rep(0,n)  # Heller-Heller-Gorfine Tests
ind58 = rep(0,n)  # Hoeffding's D test or Bergsma-Dassios T* sign covariance
ball58 = rep(0,n) # Ball correlation
qad58 = rep(0,n) # Quantification of Asymmetric Dependence
bet58 = rep(0,n) # Binary Expansion Testing
mixed58 = rep(0,n) # Cramer-von Mises statistics
nns58 = rep(0,n) # Nonlinear Nonparametric Statistics
for (i in 1:n){
  dcor58[i] = dcor(heart1[,i],heart1[,58])
  dhsic58[i] = dhsic(heart1[,i],heart1[,58])$dHSIC
  Dx = as.matrix(dist((heart1[,i]), diag = TRUE, upper = TRUE))
  Dy = as.matrix(dist((heart1[,58]), diag = TRUE, upper = TRUE))
  hhg58[i] = hhg.test(Dx,Dy, nr.perm = 1000)
  ind58[i] = hoeffding.D.test(heart1[,i],heart1[,58])$Dn
  #ind58[i] = hoeffding.refined.test(heart1[,i],heart1[,58])$Rn
  #ind58[i] = tau.star.test(heart1[,i],heart1[,58])$Tn
  ball58[i] = bcor(heart1[,i],heart1[,58])
  qad58[i] = qad(heart1[,i],heart1[,58])$`q(X,Y)`
  bet58[i] = BEAST(heart1[,c(i,58)], 3, index = list(1,2))$BEAST.Statistic
  mixed58[i] = TestIndCopula(heart1[,c(i,58)])$stat
  nns58[i] = NNS.dep(heart1[,i],heart1[,58])
}
dcor58[c(1,2,58)] = 0
dhsic58[c(1,2,58)] = 0
hhg58 = unlist(hhg58)
hhg58[c(1,2,58)] = 0
ind58[c(1,2,58)] = 0
ball58[c(1,2,58)] = 0
qad58[c(1,2,58)] = 0
bet58[c(1,2,58)] = 0
mixed58[c(1,2,58)] = 0
nns58[c(1,2,58)] = 0

#### plot
# ce
x11(width = 10, height = 5)
plot(ce58, xlab = "Variable", ylab = "Copula Entropy", xaxt = 'n')
lines(ce58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16a = rep(ce58[16],75)
lines(th16a, col = "red")
# dcor
x11(width = 10, height = 5)
plot(dcor58, xlab = "Variable", ylab = "dCor", xaxt = 'n')
lines(dcor58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16b = rep(dcor58[16],75)
lines(th16b, col = "red")
# dhsic
x11(width = 10, height = 5)
plot(dhsic58, xlab = "Variable", ylab = "dHSIC", xaxt = 'n')
lines(dhsic58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16c = rep(dhsic58[16],75)
lines(th16c, col = "red")
# hhg
x11(width = 10, height = 5)
plot(hhg58, xlab = "Variable", ylab = "HHG", xaxt = 'n')
lines(hhg58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16d = rep(hhg58[16],75)
lines(th16d, col = "red")
# independence
x11(width = 10, height = 5)
plot(ind58, xlab = "Variable", ylab = "Hoeffding", xaxt = 'n')
lines(ind58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16e = rep(ind58[16],75)
lines(th16e, col = "red")
# ball
x11(width = 10, height = 5)
plot(ball58, xlab = "Variable", ylab = "Ball", xaxt = 'n')
lines(ball58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16f = rep(ball58[16],75)
lines(th16f, col = "red")
# qad
x11(width = 10, height = 5)
plot(qad58, xlab = "Variable", ylab = "qad", xaxt = 'n')
lines(qad58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16g = rep(qad58[16],75)
lines(th16g, col = "red")
# BET
x11(width = 10, height = 5)
plot(bet58, xlab = "Variable", ylab = "BET", xaxt = 'n')
lines(bet58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16h = rep(bet58[16],75)
lines(th16h, col = "red")
# MixedIndCopula
x11(width = 10, height = 5)
mixed58 = unlist(mixed58)
plot(mixed58, xlab = "Variable", ylab = "Mixed", xaxt = 'n')
lines(mixed58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16i = rep(mixed58[16],75)
lines(th16i, col = "red")
# NNS
x11(width = 10, height = 5)
nns58 = unlist(nns58)
plot(nns58, xlab = "Variable", ylab = "Mixed", xaxt = 'n')
lines(nns58)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16j = rep(nns58[16],75)
lines(th16j, col = "red")


