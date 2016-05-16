setwd("E:/Mstat Courses/Statistical Modelling/Final Project")
# we have been provided with a dataset from medical researchers who wanted to investigate
# the number of irregular heart beats during the fixed time period of a certainmedical test. Irregular
# heart beats may be indicative for certain diseases. Several other clinical variables have been measured
# in this study too. There are repeated observations (1, 2 or 3 measurements per patient) for several
# patients in the dataset. The task is to construct a good model to estimate the expected number of
# irregular heart beats during the fixed time period of the medical test using the data.
#subset to work
fulldata.B = read.table("MedicalData2015.txt",header=T)
set.seed(0479671)
rownumbers.B = sample(unique(fulldata.B$Patient),80,replace=F)
mydata.B = fulldata.B[fulldata.B$Patient%in%rownumbers.B,]

#checking the subset
names(mydata.B)
summary(mydata.B)

#checking the subset
names(mydata.B)
summary(mydata.B)
attach(mydata.B)
length(mydata.B[,1])

library(hglm)

options(scipen=100)
options(digits=2)

#histogram for the response
hist(Irregbeat,density=30,breaks=8,freq=FALSE)
lines(density(Irregbeat, adjust=2),lty="dotted")
# Two probable distribution in this case are Gaussian and Poisson distribution. The idea behind checking
# for normal distribution is the fact that there are 203 observations. Hence the suspicion is that
# several discrete values would have converged to the normal case.

#test for normality 
shapiro.test(Irregbeat)
#The data is not normal
# the null hypothesis can be rejected and we
# can state that the Irregular heartbeat is not normally distributed

#test for poisson 
library (vcd)
temp<-goodfit (Irregbeat , type ="poisson",method ="ML")
summary(temp)
# Chi-sq Goodness of Fit test with the following hypothesis
# H0 : The data are consistent with Poisson distribution

# the null hypothesis can be rejected. But Irregular heart
# beat is a count data with non negative values. Hence Poisson would be a better approximation than
# any other distribution.

boxplot(Irregbeat~Patient, main="Boxplot for Patient", xlab="Patient",ylab="Irregbeat")
boxplot(Irregbeat~Day, main="Boxplot for Day", xlab="Day",ylab="Irregbeat")
# By observing the boxplots, it is clear that Pateint has a lot more variation than Day. Hence Patient
# can be chosen as the variable for random effect. Now we fit a Poisson glm with Irregular heart beat as
# the response with Patient as the random effect along with all other covariates excluding Day as fixed
# effects.

# library(gmodels)
# CrossTable(Patient,Enzyme)

attach(mydata.B)
fit1<-hglm ( fixed = Irregbeat ~ Enzyme +Hunthess + Heartperf + Bloodpress +
           Heartfail + Age + Gender + Hypertn +Smoke+Diabetes+Highchol , random = ~1| Patient ,
     data = mydata.B, calc.lik =T, family = poisson ( link =log ))
summary(fit1)
fit1$likelihood$cAIC

combinations = function(n){
  comb = NULL
  for( i in 1:n) comb = rbind(cbind(1,comb),cbind(0,comb))
  return(comb)
}

X<-mydata.B[,c(3:13)]
ncol(X)
combi<-combinations(ncol(X))
cAIC<-c()
xtemp<-c()
fittemp<-c()
#for(i in 1:1){

  temp <- combi[i,]
  rowvar <- which(temp!=0)
  xtemp<-as.matrix((X[,rowvar]))

  fit1<-hglm ( fixed = Irregbeat ~ paste(colnames(mydata.B[,2:13]),collapse="+") , random = ~1| Patient,data = mydata.B, calc.lik =T, family = poisson ( link =log ))
  #cAIC[i]<-fittemp$likelihood$cAIC
#}

# As we observe the figure on the left, we can spot an observation at the right extreme indicating a bad
# leverage point and other outliers having a bigger residual values or good leverage points. When we
# observe the covariates of the bad leverage point, it is found that the count of Irregualar heartbeat is
# 21 which is very high. Hence the coressponding observation is deleted for further model fitting
#remove heartfail #996
plot(fit1$resid,fit1$hv[1:nrow(mydata.B)], main="Residual Plot",xlab="Residuals",ylab="Fitted Values")
which(fit1$resid==max(fit1$resid))
mydata.B<-mydata.B[-which(fit1$resid==max(fit1$resid)),]
a<-paste(paste(colnames(mydata.B[,2:13]),collapse="+"), "+","(1|Patient)" )
fit1<-hglm2( Irregbeat ~ a ,data = mydata.B, calc.lik =T, family = poisson ( link =log ))
summary(fit1)

cAIC<-c()

#remove highcol #991
plot(fit1$resid,fit1$hv[1:nrow(mydata.B)], main="Residual Plot",xlab="Residuals",ylab="Fitted Values")
which(fit1$resid==max(fit1$resid)

fit1<-hglm ( fixed = Irregbeat ~ Enzyme +Hunthess + Heartperf + Bloodpress +
                     Heartfail + Age + Gender + Hypertn +Smoke+Diabetes+Highchol, random = ~1| Patient ,
                   data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[1]<-fit1$likelihood$cAIC

#Highchol removed
fit2<-hglm ( fixed = Irregbeat ~ Enzyme +Hunthess + Heartperf + Bloodpress +
               Heartfail + Age + Gender + Hypertn +Smoke+Diabetes, random = ~1| Patient ,
             data = mydata.B, calc.lik =T, family = poisson ( link =log ))
cAIC[2]<-fit2$likelihood$cAIC


#diabetes removed
fit3<-hglm ( fixed = Irregbeat ~ Enzyme +Hunthess + Heartperf + Bloodpress +
               Heartfail + Age + Hypertn +Smoke+Diabetes+Highchol, random = ~1| Patient ,
             data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[3]<-fit3$likelihood$cAIC


fit4<-hglm ( fixed = Irregbeat ~ Enzyme +Hunthess + Heartperf + Bloodpress +
               Heartfail + Hypertn +Smoke+Diabetes+Highchol, random = ~1| Patient ,
             data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[4]<-fit4$likelihood$cAIC


fit5<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf + Bloodpress +
               Heartfail + Hypertn +Smoke+Diabetes+Highchol, random = ~1| Patient ,
             data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[5]<-fit5$likelihood$cAIC


fit6<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf  +
               Heartfail + Hypertn +Smoke+Diabetes+Highchol, random = ~1| Patient ,
             data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[6]<-fit6$likelihood$cAIC

fit7<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf  +
               Heartfail + Hypertn +Smoke+Diabetes, random = ~1| Patient ,
             data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[7]<-fit7$likelihood$cAIC
summary(fit7)
library(nlme)
intervals(fit7,which="var-cov")


fit9<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf  +
               Heartfail + Hypertn +Smoke+Diabetes+Hunthess:Smoke+Heartperf:Heartfail, random = ~1| Patient ,
             data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[9]<-fit9$likelihood$cAIC

fit10<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf  +
               Heartfail + Hypertn +Smoke+Diabetes+Hunthess:Smoke+Heartperf:Heartfail+Heartperf:Hypertn+Heartperf:Diabetes, random = ~1| Patient ,
             data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[10]<-fit10$likelihood$cAIC

fit11<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf  +
                Heartfail + Hypertn +Smoke+Diabetes+Hunthess:Smoke+Heartperf:Heartfail+Heartperf:Hypertn+Heartperf:Diabetes+Heartfail:Smoke, random = ~1| Patient ,
              data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[11]<-fit11$likelihood$cAIC

fit12<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf  +
                Heartfail + Hypertn +Smoke+Diabetes+Hunthess:Smoke+Heartperf:Heartfail+Heartperf:Hypertn+Heartperf:Diabetes+Heartfail:Smoke+Hypertn:Smoke+Hypertn:Diabetes+Smoke:Diabetes, random = ~1| Patient ,
              data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[12]<-fit12$likelihood$cAIC

fit13<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf  +
                Heartfail + Hypertn +Smoke+Diabetes+Hunthess:Smoke+Heartperf:Hypertn+Heartperf:Diabetes+Heartfail:Smoke+Hypertn:Smoke+Hypertn:Diabetes+Smoke:Diabetes, random = ~1| Patient ,
              data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[13]<-fit13$likelihood$cAIC

fit14<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf  +
                Heartfail + Hypertn +Smoke+Diabetes+Hunthess:Smoke+Heartperf:Hypertn+Heartfail:Smoke+Hypertn:Smoke+Hypertn:Diabetes+Smoke:Diabetes, random = ~1| Patient ,
              data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[14]<-fit14$likelihood$cAIC

fit15<-hglm ( fixed = Irregbeat ~ Hunthess + Heartperf  +
                Heartfail + Hypertn +Smoke+Diabetes+Heartperf:Hypertn+Heartfail:Smoke+Hypertn:Smoke+Hypertn:Diabetes+Smoke:Diabetes, random = ~1| Patient ,
              data = mydata.B, calc.lik =T, family = poisson ( link =log ))      
cAIC[15]<-fit15$likelihood$cAIC

summary(fit15)
#B2
# By penalizing, we create a case where some of the
# parameter estimates may be exactly zero. The larger the penalty applied, more estimates are shrunk
# towards zero. This is convenient whenwewant to performmodel selection as LASSO performs estimation
# and model selection simeltaneously. In this case, we study systolic blood pressure as an outcome.
 library(lmmlasso)
#checking the subset
names(mydata.B)
summary(mydata.B)
attach(mydata.B)

#predictors
x.matrix=as.matrix(cbind(1,mydata.B[,c(2:5,7:14)]))
# head(x.matrix)
# x.matrix
summary(Irregbeat)
z = matrix(rep(1,nrow(mydata.B)),ncol =1)
 colnames(z)="Intercept"
grp=mydata.B$Patient
#View(mydata.B)
#y<-as.matrix(Bloodpress)
s<-seq(0.01,1,0.01)
aicval<-c()
mylambda<-0.1
fit<-c()
for(i in 1:100){
temp= lmmlasso(y=mydata.B$Bloodpress,x=x.matrix,z=z,grp=grp,lambda=s[i],pdMat="pdIdent")
aicval[i]<-temp$aic
}

which(aicval==min(aicval))
plot(s,aicval,main="AIC values", xlab="lambda", ylab="AIC")
# the penalty parameter around 0.35 would be a suitable one
# 
# The variables Age, Gender, Highchol and Diabetes were set as 0 for the fit.
fitlasso<-lmmlasso(y=mydata.B$Bloodpress,x=x.matrix,z=z,grp=grp,lambda=0.35,pdMat="pdIdent")

summary(fitlasso)
