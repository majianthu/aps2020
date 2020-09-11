library(copent) # Copula Entropy
library(energy) # Distance Correlation
library(dHSIC) # Hilbert-Schmidt Independence Criterion

## for additional tests
library(HHG) # Heller-Heller-Gorfine Tests of Independence
library(independence) # Hoeffding's D test or Bergsma-Dassios T* sign covariance

scan_heart_data <-function(filename1, nl = 0){
  data1 = scan(filename1, nlines = nl, what = list(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ""))
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

## stat dependence with attr #58
# ce
l = 50
ce58 = matrix(0,l,n)
for(j in 1:l){
for (i in 1:n){
	data2 = heart1[,c(i,58)]
	data2[,1] = data2[,1] + max(abs(data2[,1])) * 0.0005 * runif(m)
	data2[,2] = data2[,2] + max(abs(data2[,2])) * 0.0005 * runif(m)
	ce58[j,i] = copent(data2)
}
ce58[j,c(1,2,58)] = min(ce58[j,])
}
ce58m = colMeans(ce58)
# dcor
dcor58 = rep(0,n)
for (i in 1:n){
  dcor58[i] = dcor(heart1[,i],heart1[,58])
}
dcor58[c(1,2,58)] = 0
# dhsic
dhsic58 = rep(0,n)
for (i in 1:n){
  dhsic58[i] = dhsic(heart1[,i],heart1[,58])$dHSIC
}
dhsic58[c(1,2,58)] = 0
# hhg
hhg58 = rep(0,n)
for (i in 1:n){
  Dx = as.matrix(dist((heart1[,i]), diag = TRUE, upper = TRUE))
  Dy = as.matrix(dist((heart1[,58]), diag = TRUE, upper = TRUE))
  hhg58[i] = hhg.test(Dx,Dy, nr.perm = 1000)
}
hhg58 = unlist(hhg58)
hhg58[c(1,2,58)] = 0
# independence
ind58 = rep(0,n)
for (i in 1:n){
  #ind58[i] = hoeffding.D.test(heart1[,i],heart1[,58])$Dn
  #ind58[i] = hoeffding.refined.test(heart1[,i],heart1[,58])$Rn
  ind58[i] = tau.star.test(heart1[,i],heart1[,58])$Tn
}
ind58[c(1,2,58)] = 0


#### plot
# ce
x11(width = 10, height = 5)
plot(ce58m, xlab = "Variable", ylab = "Copula Entropy", xaxt = 'n')
lines(ce58m)
axis(side = 1, at = c(seq(1,75, by = 5)), labels = c(seq(1,75, by = 5)))
th16a = rep(ce58m[16],75)
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
