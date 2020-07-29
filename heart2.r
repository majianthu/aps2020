# Experiments on heart disease data
library(e1071)
library(glmnet)
library(parcor)
library(msgps)

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

#### prepare data
#58 (num) (the predicted attribute) 
cls = heart1[,58]
#### 13 attributes recommended by dataset creator
#3 (age) #4 (sex) #9 (cp) #10 (trestbps) #12 (chol) 
#16 (fbs) #19 (restecg) #32 (thalach) #38 (exang) 
#40 (oldpeak) #41 (slope) #44 (ca) #51 (thal) 
data1 = heart1[,c(3,4,9,10,12,16,19,32,38,40,41,44,51)]
#### attributes selected with Copula entropy, distance correlation, dHSIC
data2 = heart1[,c(3,4,6,7,9,12,16,28:32,38,40,41,44,51,59:68)] # copula entropy
data2b = heart1[,c(3,4,6,7,9,12,13,16,28:33,38,40,41,52,59:68)] # dcor
data2c = heart1[,c(3,4,6,7,9,12,13,16,25,28:32,38,40,41,44,59:68)] # dhsic
#### for GLMs
data3 = heart1[,c(3:57,59:68)]

## training models
# SVMs
svm1 = svm(data1,cls, type = "nu-regression", gamma = 2)
svm2 = svm(data2,cls, type = "nu-regression", gamma = 2)
svm2b = svm(data2b,cls, type = "nu-regression", gamma = 2)
svm2c = svm(data2c,cls, type = "nu-regression", gamma = 2)
# GLM
lasso1 = cv.glmnet(data3, cls, family = "multinomial", type.measure = "class", nfolds = 10)
ridge1 = cv.glmnet(data3, cls, family = "multinomial", alpha = 0, type.measure = "class", nfolds = 10)
elnet1 = cv.glmnet(data3, cls, family = "multinomial", alpha = 0.5, type.measure = "class", nfolds = 10)
#### Adaptive Lasso
ada1 = adalasso(data3,cls)
ada2 = msgps(data3,cls,penalty = "alasso")

## prediction
pred1 = predict(svm1,data1)
pred2 = predict(svm2,data2)
pred2b = predict(svm2b,data2b)
pred2c = predict(svm2c,data2c)
pred3 = predict(lasso1, newx = data3, s = "lambda.min", type = "class")
pred4 = predict(ridge1, newx = data3, s = "lambda.min",  type = "class")
pred5 = predict(elnet1, newx = data3, s = "lambda.min",  type = "class")
# adaptive lasso
## ada1 parcor
intada1 = ada1$intercept.adalasso
coefada1 = ada1$coefficients.adalasso
predada1 = data3 %*% coefada1 + intada1
## ada2 msgps
coefada2 = coef.msgps(ada2)
predada2 = predict(ada2, data3)

## correct numbers
pred1num = sum(abs(cls-pred1)<=0.5) # recommandation
pred2num = sum(abs(cls-pred2)<=0.5) # ce
pred2bnum = sum(abs(cls-pred2b)<=0.5) # dcor
pred2cnum = sum(abs(cls-pred2c)<=0.5) # dhsic
pred3num = sum(cls==pred3) # lasso
pred4num = sum(cls==pred4) # ridge regression
pred5num = sum(cls==pred5) # elastic net
prednum1ada = sum(abs(cls-predada1)<0.5) # ada1 parcor
prednum2ada = sum(abs(cls-predada2[,3])<0.5) # ada2 msgps

#### plotting for GLMs
## coefficients of lasso
coef10 = coef(lasso1, s = "lambda.min") 
coef10a = abs(coef10$'0'[2:66])+abs(coef10$'1'[2:66])+abs(coef10$'2'[2:66])+abs(coef10$'3'[2:66])+abs(coef10$'4'[2:66])
x11();
#pdf("~/Rworks/heart/lasso1.pdf")
barplot(abs(coef10a * colMeans(data3)), xlab = "Variable ID", ylab = "Coefficients Value", main = "LASSO")
axis(side = 1, at = c(1:65)*65/54-0.5, labels = c(3:57,59:68))
#dev.off()

## coefficients of ridge regression
coef40 = coef(ridge1, s = "lambda.min")
coef40a = abs(coef40$'0'[2:66])+abs(coef40$'1'[2:66])+abs(coef40$'2'[2:66])+abs(coef40$'3'[2:66])+abs(coef40$'4'[2:66])
x11();
#pdf("~/Rworks/heart/ridge1.pdf")
barplot(abs(coef40a * colMeans(data3)), xlab = "Variable ID", ylab = "Coefficients Value", main = "Ridge Regression")
axis(side = 1, at = c(1:65)*65/54-0.5, labels = c(3:57,59:68))
#dev.off()

## coefficients of elastic net
coef5 = coef(elnet1, s = "lambda.min")
coef5a = abs(coef5$'0'[2:66])+abs(coef5$'1'[2:66])+abs(coef5$'2'[2:66])+abs(coef5$'3'[2:66])+abs(coef5$'4'[2:66])
x11();
#pdf("~/Rworks/heart/elnet1.pdf")
barplot(abs(coef5a * colMeans(data3)), xlab = "Variable ID", ylab = "Coefficients Value", main = "Elastic Net")
axis(side = 1, at = c(1:65)*65/54-0.5, labels = c(3:57,59:68))
#dev.off()

## coefficients of Adaptive Lasso
# ada1 parcor
names(coefada1) = {} #c(3:57,59:68)
x11(); 
#pdf("~/Rworks/heart/adalasso1.pdf")
barplot(abs(coefada1), xlab = "Variable ID", ylab = "Coefficients Value", main = "Adaptive LASSO")
axis(side = 1, at = c(1:65)*65/54-0.5, labels = c(3:57,59:68))
#dev.off()

# ada2 msgps
coefada23 = abs(coefada2[2:66,3])
names(coefada23) ={}
x11(); barplot(coefada23, xlab = "Variable ID", ylab = "Coefficients Value", main = "Adaptive LASSO2")
axis(side = 1, at = c(1:65)*65/54-0.5, labels = c(3:57,59:68))

######## step_AIC(glm)
d31 = data.frame(data3)
names(d31) = c(3:57,59:68)
d3a = cbind(d31,cls)
glm1 = glm(cls~., data = d3a, family = poisson(link = "log"))
glm1s = step(glm1)
pred1glm1 = predict(glm1s, newdata = d3a, type = "response")
prednum1aic = sum(abs(cls-pred1glm1)<0.5)
cat("Accuracy of stepwise GLM (AIC) :", prednum1aic,"\n")
## plot coefficients
idx1 = c(3,4,5,9,12,16,18,26,28,29,30,32,38,40,44,47,50,53,54,60,61,63,64,65,66,67)
coef1step = rep(0,68) 
coef1step[idx1] = abs(coef(glm1s)[2:27])
x11(); 
#pdf("~/Rworks/heart/step1.pdf")
barplot(coef1step[c(3:57,59:68)], xlab = "Variable ID", ylab = "Coefficients Value", main = "Stepwise GLM(AIC)")
axis(side = 1, at = c(1:65)*65/54-0.5, labels = c(3:57,59:68))
#dev.off()

#### step_BIC(glm)
d31 = data.frame(data3)
names(d31) = c(3:57,59:68)
d3a = cbind(d31,cls)
glm2 = glm(cls~., data = d3a, family = poisson(link = "log"))
glm2s = step(glm2, k = log(length(cls)))
pred2glm2 = predict(glm2s, newdata = d3a, type = "response")
prednum1bic = sum(abs(cls-pred2glm2)<0.5)
cat("Accuracy of stepwise GLM (BIC) :", prednum1bic,"\n")
## plot coefficients
idx2 = c(3,4,5,9,16,18,29,30,40,53,63,66,67)
coef2step = rep(0,68) 
coef2step[idx2] = abs(coef(glm2s)[2:14])
x11(); 
#pdf("~/Rworks/heart/step2.pdf")
barplot(coef2step[c(3:57,59:68)], xlab = "Variable ID", ylab = "Coefficients Value", main = "Stepwise GLM(BIC)")
axis(side = 1, at = c(1:65)*65/54-0.5, labels = c(3:57,59:68))
#dev.off()
