rm(list=objects());graphics.off()



load("cookie.app.RData")
load("cookie.val.RData")

xtrain <- as.matrix(cookie.app)
xtest <- as.matrix(cookie.val)

ytrain <- xtrain[,1]
ytest <- xtest[,1]


xtrain <- xtrain[,-1]
xtest <- xtest[,-1]


plot.new()
library(car)
par(mfrow=c(2,2),mar=c(2,2,2,1))
boxplot(xtrain[,1:50])
boxplot(xtrain[,51:100])
boxplot(xtrain[,101:150])
boxplot(xtrain[,151:200])

boxplot(xtrain[,201:250])
boxplot(xtrain[,251:300])
boxplot(xtrain[,301:350])
boxplot(xtrain[,351:400])

boxplot(xtrain[,401:450])
boxplot(xtrain[,451:500])
boxplot(xtrain[,501:550])
boxplot(xtrain[,551:600])

par(mfrow=c(1,2),mar=c(6,1,6,1))
boxplot(xtrain[,601:650])
boxplot(xtrain[,651:700])

par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
matplot(x=1:700,y=t(xtrain),type="l",xlab="variables (fréquences)", ylab="absorption",lty=c(1,1))

C <- cor(xtrain)
C.eig <- eigen(C)
plot(C.eig$values,xlab="valeur propre",ylab="")
points(C.eig$values,col="red")


library(FactoMineR)
res=PCA(xtrain)
nbr.vp <- nrow(res$eig) #nombre de valeurs propres
#graphe des valeurs propres
par(mfrow=c(1,1))
plot(res$eig[,1],xlab="valeur propre",ylab="")
points(res$eig[,1],pch=20,col="red")

barplot(res$eig[,2],main="% inertie",names=paste("Dim",1:nrow(res$eig)))

par(mfrow=c(1,2))
#premier plan principal
plot(res,axes = c(1,2), choix = "var",graph.type = "classic")
plot(res,axes = c(1,2), choix = "ind",graph.type = "classic")
#deuxième plan principal
plot(res,axes = c(2,3), choix = "var",graph.type = "classic")
plot(res,axes = c(2,3), choix = "ind",graph.type = "classic")
#troisième plan principal
plot(res,axes = c(3,4), choix = "var",graph.type = "classic")
plot(res,axes = c(3,4), choix = "ind",graph.type = "classic")
#quatrième plan principal
plot(res,axes = c(4,5), choix = "var",graph.type = "classic")
plot(res,axes = c(4,5), choix = "ind",graph.type = "classic")
#cinquième plan principal
plot(res,axes = c(5,6), choix = "var",graph.type = "classic")
plot(res,axes = c(5,6), choix = "ind",graph.type = "classic")


reconstruct<-function(res,nr,Xm,Xsd)
{
  coord.var = as.matrix(res$var$coord)[,1:nr,drop=F]
  coord.ind = as.matrix(res$ind$coord)[,1:nr,drop=F]
  x_rec = coord.ind%*%t(sweep(coord.var,2,sqrt(res$eig[1:nr,1]),FUN="/"))
  x_rec = sweep(x_rec,2,Xsd,FUN="*")
  x_rec = sweep(x_rec,2,Xm,FUN="+")
  return(x_rec)
}
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)

res <- PCA(xtrain,ncp=39)
library("Metrics")
mat=as.matrix(xtrain)
Xm <- apply(mat,2,mean)
Xsd <- apply(mat,2,sd)

Xrec1=reconstruct(res,1,Xm,Xsd)
Xrec2=reconstruct(res,2,Xm,Xsd)
Xrec3=reconstruct(res,3,Xm,Xsd)
Xrec4=reconstruct(res,4,Xm,Xsd)
Xrec5=reconstruct(res,5,Xm,Xsd)
Xrec39=reconstruct(res,39,Xm,Xsd)

par(mfrow=c(3,2),mar=c(2,4,4,2)+0.1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 1",type="l",main=paste("RMSE=",round(rmse(mat,Xrec1),5),"MAE=",round(mae(mat,Xrec1),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec2),xlab="",ylab="nr = 2",type="l",main=paste("RMSE=",round(rmse(mat,Xrec2),5),"MAE=",round(mae(mat,Xrec2),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 3",type="l",main=paste("RMSE=",round(rmse(mat,Xrec3),5),"MAE=",round(mae(mat,Xrec3),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 4",type="l",main=paste("RMSE=",round(rmse(mat,Xrec4),5),"MAE=",round(mae(mat,Xrec4),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 5",type="l",main=paste("RMSE=",round(rmse(mat,Xrec5),5),"MAE=",round(mae(mat,Xrec5),5),collapse = " "),cex.main=1)
matplot(x=1:700,y=t(Xrec1),xlab="",ylab="nr = 39",type="l",main=paste("RMSE=",round(rmse(mat,Xrec39),5),"MAE=",round(mae(mat,Xrec39),5),collapse = " "),cex.main=1)

par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
matplot(Xrec1[,24],type="l",col="gray0",ylab="X24")
matlines(Xrec2[,24],col="red")
matlines(Xrec3[,24],col="darkblue")
matlines(Xrec4[,24],col="chartreuse")
matlines(Xrec5[,24],col="gray")
matlines(Xrec39[,24],col="darkmagenta")
legend("topleft",c("Xrec1","Xrec2","Xrec3","Xrec4","Xrec5","Xrec39"),
       col=c("gray0","red","darkblue","chartreuse","gray","darkmagenta"),lty=c(1:6),lwd = 0.5 , xpd = T )
#pas de grande variation entre les courbes de X24 aux six niveaux, le premier axe principale représente bien cette variable



library(glmnet)

grid=10^seq(6,-10,length=100) # la grille de lambda
ridge.glm=glmnet(xtrain,ytrain,alpha=0,lambda=grid)

coef(ridge.glm)[1,] #l'intercept
plot(log(grid),coef(ridge.glm)[1,],xlab="log(k)",ylab="intercept")
points(log(grid),coef(ridge.glm)[1,],col="red",pch=20)
abline(h=mean(ytrain),lwd=2,col="blue")
text(-20, mean(ytrain)+2, "mean(y)", col = "blue")

#recalcule de l'intercept
theta0 <- mean(ytrain)*rep(1,length(grid)) - colMeans(xtrain%*%coef(ridge.glm)[-1,])
plot(log(grid),theta0)
points(log(grid),theta0,lwd=2,col="red")

#si on centre ytrain seulement
ridge.glm=glmnet(xtrain,scale(ytrain,scale=FALSE),alpha=0,lambda=grid)
plot(log(grid),coef(ridge.glm)[1,],xlab="log(k)",ylab="intercept")
points(log(grid),coef(ridge.glm)[1,],col="red",pch=20)
abline(h=0,lwd=2,col="blue")
text(-20, 2, "mean(y)", col = "blue")

#si on centre xtrain seulement
ridge.glm=glmnet(scale(xtrain,scale=FALSE),ytrain,alpha=0,lambda=grid)
plot(log(grid),coef(ridge.glm)[1,],xlab="log(k)",ylab="intercept")
points(log(grid),coef(ridge.glm)[1,],col="red",pch=20)
abline(h=mean(ytrain),lwd=2,col="blue")
text(-20, mean(ytrain)+0.02, "mean(y)", col = "blue")

#si on centre tous les deux
ridge.glm=glmnet(scale(xtrain,scale=FALSE),scale(ytrain,scale=FALSE),alpha=0,lambda=grid)
plot(log(grid),coef(ridge.glm)[1,],xlab="log(k)",ylab="intercept")
points(log(grid),coef(ridge.glm)[1,],col="red",pch=20)
abline(h=0,lwd=2,col="blue")
text(-20, 0.02, "mean(y)", col = "blue")

#estimation de theta quand lambda tend vers 0
xscaled <- scale(xtrain)
yscaled <- scale(ytrain)

C <- eigen(t(xscaled)%*%xscaled)
ind <- which(C$values > 0.0001) #length(ind)=rang de xscaled
eig.values <- C$values[ind]
v <- C$vectors[ind,]
u <- eigen(xscaled%*%t(xscaled))$vectors[ind,]

mat <- matrix(0,ncol(xtrain),nrow(xtrain))
for (i in 1:length(ind)){
  m <- matrix(0,ncol(xtrain),nrow(xtrain))
  for (j in 1:nrow(xtrain)){
    m[,j] <- t(v[i,])*u[i,j]
  }
  m <- m*(1/sqrt(C$values[ind[i]]))
  mat <- mat+m
}
theta_l <- mat%*%yscaled
theta_l

library(MASS)
ridge.lm <- lm.ridge(ytrain~xtrain,lambda=grid)
coef(ridge.lm)[1,]
plot(ridge.lm$a0)
all(coef(ridge.lm)==t(coef(ridge.glm))) # FALSE, pas les mêmes coefficients

coefstheta <- matrix(0,701,length(grid))
coefstheta[1,] <- theta0

for (i in 1:length(grid)){
  coefstheta[-1,i] <- solve(t(xtrain)%*%xtrain + grid[i]*diag(1,700,700))%*%(t(xtrain)%*%ytrain-t(xtrain))
}
plot(log(grid),coef(ridge.glm)[2,])
plot(log(grid),thetas[1,])

library(pls)
set.seed(1)
cv.seg <- cvsegments(nrow(xtrain),4,type="random")
cv.errors = vector()
grid=10^seq(5.5,-5.5,length=100)

x.train <- xtrain
y.train <- ytrain

cv.models  = matrix(0,length(grid),2)

for (i in 1:length(grid)){
  k = grid[i]
  errorsk <- rep(0,length(cv.seg))
  for (j in 1:length(cv.seg)){
    xval <- xtrain[unlist(cv.seg[j]),]
    yval <- ytrain[unlist(cv.seg[j])]
    xt <- xtrain[unlist(cv.seg[-j]),]
    yt <- ytrain[unlist(cv.seg[-j])]
    reg <- glmnet(xt,yt,alpha=0,lambda=k)
    ypred <- predict(reg,xval)
    err <- mean((yval-ypred)^2)
    errorsk[j] = err
  }
  cv.models[i,1] = k
  cv.models[i,2] = mean(errorsk)
}

ind <- which(cv.models[,2]==min(cv.models[,2]))
kappa <- grid[ind]
print(kappa)

par(mfrow=c(1,2))
plot(log(cv.models[,1]),cv.models[,2])
points(log(cv.models[,1]),cv.models[,2],col="red",pch=20)


set.seed(1)
cv.out=cv.glmnet(xtrain,ytrain,alpha=0,nfolds = 4,lambda=grid) 

plot(cv.out)
print(cv.out$lambda.min)
print(cv.out)
# différent de kappa que l'on a trouvé

# on va utiliser le paramètre cv.out$labmda.min 
reg.ridge <- glmnet(x.train,y.train,alpha=0,lambda=cv.out$lambda.min)
pred <- predict(reg.ridge,scale(xtest))
err.gener <- mean((scale(ytest,scale=FALSE)-pred)^2)




par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)

z=sapply(ytrain,function(x){if(x<18){x<-0}
  else{x<-1}})

ztest=sapply(ytest,function(x){if(x<18){x<-0}
  else{x<-1}})
sum(z[which(z==1)])
sum(ztest[which(ztest==1)])

logridge<-cv.glmnet(as.matrix(xtrain),z,alpha=0,lambda=grid,family="binomial")
loglasso<-cv.glmnet(as.matrix(xtrain),z,alpha=1,lambda=grid,family="binomial")

plot(logridge)
plot(loglasso)

ridgemodel=glmnet(as.matrix(xtrain),z,alpha=0,family = "binomial",lambda = logridge$lambda.min)
x.test <- model.matrix(sucres ~., cookie.val)[,-1]
library(magrittr)
probabilities <- ridgemodel %>% predict(newx = x.test)
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
observed.classes <- ztest
mean(predicted.classes == observed.classes)
library(pROC)
res.roc <- roc(observed.classes, probabilities)
plot(1-res.roc$specificities,res.roc$sensitivities,type="l")
lines(c(0,1),c(0,1),col=2,lty=2)  
segments(c(0,0),c(0,1),c(0,1),c(1,1),lty=3,lwd=2,col=2) 



lassomodel=glmnet(as.matrix(xtrain),z,alpha=1,family = "binomial",lambda = loglasso$lambda.min)
x.test <- model.matrix(sucres ~., cookie.val)[,-1]
probabilities <- lassomodel %>% predict(newx = x.test)
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
observed.classes <- ztest
mean(predicted.classes == observed.classes)
res.roc <- roc(observed.classes, probabilities)
plot(1-res.roc$specificities,res.roc$sensitivities,type="l")
lines(c(0,1),c(0,1),col=2,lty=2)  
segments(c(0,0),c(0,1),c(0,1),c(1,1),lty=3,lwd=2,col=2)

