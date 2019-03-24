####Generate The Data For both hetergeneous and homogenour 
library(WhatIf)
library(Matching)
##example
sqdata <- data.frame(t = c(1, 1, 1, 1, 0, 0, 0, 0),x = c(0, 0, 1, 1, .5, .5, 1.5, 1.5),y = c(1, 0, 0, 1, .5, 1.5, .5, 1.5))
summary(whatif(~ x + y, data = sqdata[sqdata$t == 1,], cfact = sqdata[sqdata$t == 0,]))
summary(whatif(~ x + y, data = sqdata[sqdata$t == 0,], cfact = sqdata[sqdata$t == 1,]))



#Generate Data
library(mvtnorm)


#Set Seed and N
set.seed(1298)
N = 10000

#Step 1:  Generate Multivariate normal data.  3 covariates.
mydata = rmvnorm(N,sigma = diag(2,nrow = 3))
names(mydata) = c("X1", "X2", "X3")

#Define Propensity score
trueps = (4 - mydata[,1]^2 - mydata[,2]^2)/3.5
trueps = pmax(trueps,0)
trueps = pmin(trueps,1)

trt = sapply(trueps, function(x){rbinom(1,1,x)})
mydata = data.frame(mydata,trt)

#Response Homog
response = with(mydata, 1/(X1^2 + X2^2 + X3^2) + 2*trt)

potouttrt =  with(mydata, 1/(X1^2 + X2^2 + X3^2) + 2)
potoutcon =  with(mydata, 1/(X1^2 + X2^2 + X3^2) )



#Response Heterog
responsehet = with(mydata, 1/(X1^2 + X2^2 + X3^2) + 2*X1*trt + 3*X2*trt + 2*trt)

potouttrthet =  with(mydata, 1/(X1^2 + X2^2 + X3^2) + 2*X1 + 3*X2 + 2)
potoutconhet =  with(mydata, 1/(X1^2 + X2^2 + X3^2) )


fitmodel = glm(trt~(X1 + X2+ X3)^2+ I(X1^2) + I(X2^2) + I(X3^2), data = mydata, family = "binomial")
estpropscore = fitmodel$fitted.values
mydata=data.frame(cbind(mydata),"response"=response,"tp"=trueps,"ep"=estpropscore,"responsehet"=responsehet,"potouttrthet"=potouttrthet,"potoutconhet"=potoutconhet)


#######LCC
##
## True Propensity
#################
t0.tp=Sys.time()
df.tpop=LLCabspropensityscore(dataframe=mydata,trt=4,criticalvalue = 0.005,propensity=trueps)
coms.tp=LLCgetLargest(dflist=df.tpop,trt=4,bytrtsize =TRUE,order=1,bindit = TRUE)
t1.tp=Sys.time()
t.tp=t1.tp-t0.tp
par(mfrow=c(1,2))
plot(mydata$X1,mydata$X2,xlab="X1",ylab="X2",col=c("blue","red")[mydata$trt+1],main="Original Data")
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=1.05,pch=1)
plot(coms.tp$X1,coms.tp$X2,xlab="X1",ylab="X2",col=c("blue","red")[coms.tp$trt+1],main="Common Support 1",xlim=c(-4,4),ylim=c(-4,4))
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=0.8,pch=1)
llc_hete_tp_tatt=mean(coms.tp[coms.tp$trt==1,]$potouttrthet)-mean(coms.tp[coms.tp$trt==1,]$potoutconhet)

#genmatch
genmat.ep = GenMatch(coms.tp$trt,coms.tp[,1:3],max.generations=100000)
matGenWei.ep = Match(coms.tp$responsehet,coms.tp$trt, coms.tp[,1:3], Weight.matrix = genmat.ep$Weight.matrix)
summary(matGenWei.ep)

##1-1 Propensity match
matchATT.ep = Match(coms.tp$responsehet,coms.tp$trt, coms.tp$ep, estimand = "ATT")
summary(matchATT.ep)

## Estimated Propensity
#################
t0.ep=Sys.time()
df.epop=LLCabspropensityscore(dataframe=mydata,trt=4,criticalvalue = 0.005,propensity=estpropscore)
coms.ep=LLCgetLargest(dflist=df.epop,trt=4,bytrtsize =TRUE,order=1,bindit = TRUE)
t1.ep=Sys.time()
t.ep=t1.ep-t0.ep
par(mfrow=c(1,2))
plot(mydata$X1,mydata$X2,xlab="X1",ylab="X2",col=c("blue","red")[mydata$trt+1],main="Original Data")
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=1.05,pch=1)
plot(coms.ep$X1,coms.ep$X2,xlab="X1",ylab="X2",col=c("blue","red")[coms.ep$trt+1],main="Common Support 2",xlim=c(-4,4),ylim=c(-4,4))
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=0.8,pch=1)
llc_hete_ep_tatt=mean(coms.ep[coms.ep$trt==1,]$potouttrthet)-mean(coms.ep[coms.ep$trt==1,]$potoutconhet)

#genmatch
genmat.ep = GenMatch(coms.ep$trt,coms.ep[,1:3],max.generations=100000)
matGenWei.ep = Match(coms.ep$responsehet,coms.ep$trt, coms.ep[,1:3], Weight.matrix = genmat.ep$Weight.matrix)
summary(matGenWei.ep)

##1-1 Propensity match
matchATT.ep = Match(coms.ep$responsehet,coms.ep$trt, coms.ep$ep, estimand = "ATT")
summary(matchATT.ep)


## General Euclidean distance
#################
t0.eu=Sys.time()
df.eu=MultiLLCoperation(dataframe=mydata,start=1,end=3,trt=4,criticalvalue=0.35)
coms.eu=LLCgetLargest(dflist=df.eu,trt=4,bytrtsize = FALSE,order=1,bindit = TRUE)
t1.eu=Sys.time()
t.eu=t1.eu-t0.eu
par(mfrow=c(1,2))
plot(mydata$X1,mydata$X2,xlab="X1",ylab="X2",col=c("blue","red")[mydata$trt+1],main="Original Data")
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=1.05,pch=1)
plot(coms.eu$X1,coms.eu$X2,xlab="X1",ylab="X2",col=c("blue","red")[coms.eu$trt+1],main="Common Support 3",xlim=c(-4,4),ylim=c(-4,4))
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=0.8,pch=1)
llc_hete_eu_tatt=mean(coms.eu[coms.eu$trt==1,]$potouttrthet)-mean(coms.eu[coms.eu$trt==1,]$potoutconhet)

#genmatch
genmat.ep = GenMatch(coms.eu$trt,coms.eu[,1:3],max.generations=100000)
matGenWei.ep = Match(coms.eu$responsehet,coms.eu$trt, coms.eu[,1:3], Weight.matrix = genmat.ep$Weight.matrix)
summary(matGenWei.ep)

##1-1 Propensity match
matchATT.ep = Match(coms.eu$responsehet,coms.eu$trt, coms.eu$ep, estimand = "ATT")
summary(matchATT.ep)

## Mahalanobis distance
#################
t0.maha=Sys.time()
s=matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,byrow=TRUE)
df.maha=MultiLLCMahalanobis(dataframe=mydata,start=1,end=3,trt=4,criticalvalue=0.09,covs=s,inverted=FALSE)
coms.maha=LLCgetLargest(dflist=df.maha,trt=4,bytrtsize = FALSE,order=1,bindit = TRUE)
t1.maha=Sys.time()
t.maha=t1.maha-t0.maha

par(mfrow=c(1,2))
plot(mydata$X1,mydata$X2,xlab="X1",ylab="X2",col=c("blue","red")[mydata$trt+1],main="Original Data")
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=1.05,pch=1)
plot(coms.maha$X1,coms.maha$X2,xlab="X1",ylab="X2",col=c("blue","red")[coms.maha$trt+1],main="Common Support 4",xlim=c(-4,4),ylim=c(-4,4))
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=0.8,pch=1)

llc_hete_maha_tatt=mean(coms.maha[coms.maha$trt==1,]$potouttrthet)-mean(coms.maha[coms.maha$trt==1,]$potoutconhet)

genmat.ep = GenMatch(coms.maha$trt,coms.maha[,1:3],max.generations=100000)
matGenWei.ep = Match(coms.maha$responsehet,coms.maha$trt, coms.maha[,1:3], Weight.matrix = genmat.ep$Weight.matrix)
summary(matGenWei.ep)

##1-1 Propensity match
matchATT.ep = Match(coms.maha$responsehet,coms.maha$trt, coms.maha$ep, estimand = "ATT")
summary(matchATT.ep)

###Largest Caliper
####################
t0.lcp=Sys.time()
df.lcp=LargestCaliper(dataframe=mydata,start=1,end=3,trt=4,criticalvalue =1, standerdize = FALSE,weights=c(0.3,0.3,0.3))
coms.lcp=LLCgetLargest(dflist=df.lcp,trt=4,bytrtsize = FALSE,order=1,bindit = TRUE)
t1.lcp=Sys.time()
t.lcp=t1.lcp-t0.lcp
par(mfrow=c(1,2))
plot(mydata$X1,mydata$X2,xlab="X1",ylab="X2",col=c("blue","red")[mydata$trt+1],main="Original Data")
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=1.05,pch=1)
plot(coms.lcp$X1,coms.lcp$X2,xlab="X1",ylab="X2",col=c("blue","red")[coms.lcp$trt+1],main="Common Support 5",xlim=c(-4,4),ylim=c(-4,4))
legend("topright",legend=c("Control","Trt"),col=c("blue","red"),text.width=0.8,pch=1)
llc_hete_lcp_tatt=mean(coms.lcp[coms.lcp$trt==1,]$potouttrthet)-mean(coms.lcp[coms.lcp$trt==1,]$potoutconhet)

genmat.ep = GenMatch(coms.lcp$trt,coms.lcp[,1:3],max.generations=100000)
matGenWei.ep = Match(coms.lcp$responsehet,coms.lcp$trt, coms.lcp[,1:3], Weight.matrix = genmat.ep$Weight.matrix)
summary(matGenWei.ep)

##1-1 Propensity match
matchATT.ep = Match(coms.lcp$responsehet,coms.lcp$trt, coms.lcp$ep, estimand = "ATT")
summary(matchATT.ep)
