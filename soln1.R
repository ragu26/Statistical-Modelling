setwd("E:/Mstat Courses/Statistical Modelling/Final Project")
# The data to be analyzed is obtained from the European Social Survey conducted
# in Germany. A model is required to be designed with variable stsfeco, indicting
# the satisfaction with the present state of the economy as the response with all
# other variables as indictors.
fulldata.A = read.table("SatisfactionEconomy2015.txt", header = T)
set.seed(479671)
rownumbers = sample(1:1631, size = 800)

mydata.A = fulldata.A[rownumbers, ]

mydata.A$stsfeco <- factor(mydata.A$stsfeco)
mydata.A$pdwrk <- factor(mydata.A$pdwrk)
mydata.A$chldhm <- factor(mydata.A$chldhm)
mydata.A$gincdif <- factor(mydata.A$gincdif)
mydata.A$gndr <- factor(mydata.A$gndr)

attach(mydata.A)
detach(mydata.A)
# Creating a mean centered data set
mydata.mc = mydata.A
mydata.mc[, 2:8] = data.frame(scale(mydata.A[, 2:8], center = T))
attach(mydata.mc)
detach(mydata.mc)

# #Main effects model The response is binary, hence a logistic regression of the
# following form is fitted
fit0 <- glm(stsfeco ~ ., family = binomial(link = "logit"), data = mydata.A)
summary(fit0)

# it will be appropriate to use stepwise algorithm through which covariates can
# be added or removed by means of their p-value. For stepwise algorithm we have
# to define the scope of the model, i.e. lower and upper limits for models which
# can be used for model fitting. stepwise AIC model
library(MASS)
fitAIC = stepAIC(object = fit0, direction = "both", scope = stsfeco ~ . + .^2, )
summary(fitAIC)

# BIC model
fitBIC = stepAIC(object = fit0, direction = "both", scope = stsfeco ~ . + .^2, k = log(800))
summary(fitBIC)

### eof#

mydata.A = fulldata.A[rownumbers, ]
mydata.A <- data.frame(mydata.A)
# look at the seed data
summary(mydata.A)

# initializing
mydata.A$stsfeco <- factor(mydata.A$stsfeco)
mydata.A$pdwrk <- factor(mydata.A$pdwrk)
mydata.A$chldhm <- factor(mydata.A$chldhm)
mydata.A$gndr <- factor(mydata.A$gndr)
mydata.A$gincdif <- factor(mydata.A$gincdif)
summary(mydata.A)

attach(mydata.A)
# reposne binomial
fitnull <- glm(stsfeco ~ NULL, family = binomial)
fit1 <- glm(stsfeco ~ Stsf + Trust + Pol + agea + eduyrs + wkhtot + wkhsch + pdwrk + 
              chldhm + hhmmb + gincdif + gndr, family = binomial)
temp <- glm(stsfeco ~ poly(Stsf, 2) * poly(Trust, 2), family = binomial)
summary(temp)
options(na.action = "na.omit")
fittemp <- c()


fitfull <- glm(stsfeco ~ .^2, family = binomial, data = , mydata.A)
summary(fitfull)
library(sjPlot)
sjp.int(fitfull)

# model selection
library(MASS)


fitfullalt <- glm(stsfeco ~ Stsf + Trust + Pol + agea + eduyrs + wkhtot + wkhsch + 
                    pdwrk + chldhm + hhmmb + gincdif + gndr + wkhsch:hhmmb + Trust:Stsf + Trust:wkhtot + 
                    eduyrs:Stsf, family = binomial)
library(MASS)

# model 1
modelAICback <- stepAIC(fitfullalt, direction = "backward")
summary(modelAICback)
# #stsfeco ~ Stsf + Trust + eduyrs + wkhtot + wkhsch + hhmmb + wkhsch:hhmmb +
# Stsf:Trust + Trust:wkhtot + Stsf:eduyrs #AIC 86k9

# model2
modelAICforward <- stepAIC(fit1, scope = list(upper = fitfullalt, lower = ~1), direction = "forward")
summary(modelAICforward)
## stsfeco ~ Stsf + Trust + Pol + agea + eduyrs + wkhtot + wkhsch + pdwrk + chldhm
## + hhmmb + gincdif + gndr + Stsf:Trust + Stsf:eduyrs + Trust:wkhtot +
## wkhsch:hhmmb AIC 882

# model3
modelAICstepwisefull <- stepAIC(fitfullalt, direction = "both")
summary(modelAICstepwisefull)
# stsfeco ~ Stsf + Trust + eduyrs + wkhtot + wkhsch + hhmmb + wkhsch:hhmmb +
# Stsf:Trust + Trust:wkhtot + Stsf:eduyrs AIC 869

# model 1
modelBICback <- stepAIC(fitfullalt, direction = "backward", k = log(800))
summary(modelBICback)

# model2
modelBICforward <- stepAIC(fitnull, scope = list(upper = fitfullalt, lower = ~1), 
                           direction = "forward", k = log(800))
summary(modelBICforward)
# stsfeco ~ Stsf + Trust + Stsf:Trust AIC 875

# model3
modelBICstepwisefull <- stepAIC(fitfullalt, direction = "both", k = log(800))
summary(modelBICstepwisefull)
# stsfeco ~ Stsf + Trust + Stsf:Trust AIC 875

# model4
modelBICstepwisenull <- stepAIC(fitnull, scope = list(upper = fitfullalt, lower = ~1), 
                                direction = "both")
summary(modelBICstepwisenull)
# stsfeco ~ Trust + wkhtot + eduyrs + Trust:wkhtot AIC 977

# alternative model with gndr
fitAIC <- glm(stsfeco ~ Stsf + Trust + Pol + agea + eduyrs + wkhtot + wkhsch + pdwrk + 
                chldhm + hhmmb + gincdif + gndr + Stsf:Trust + Stsf:eduyrs + Trust:wkhtot + wkhsch:hhmmb, 
              family = binomial)
summary(fitAIC)

# data for prediction

fitconf <- glm(stsfeco ~ Stsf + Trust + Pol + pdwrk + chldhm + gincdif, family = binomial)


# These two models had lower AIC compared to all othermodels viz. former(0.3)
# with AIC of 882 and the latter(0.4) with AIC of 875 with very less to choose
# between them. But the glaring difference between them is the number of
# covariates. Both the models have their advantages with the former having the
# capability to be more precise for estimation and latter satisfying the
# principle of parsimony.

# we have been provided with three covariates, Tust, Stsf and intercation between
# Trust and Stsf. We have to test in a nonparametric way the null hypothesis that
# an appropriate transformation of the expected satisfaction with the present
# state of the economy is a function of only these three covariates.
# non-parametric testing

# finding critical value for 1% significance i.e. Cn
t_os <- function(x) {
  temp = 0
  for (j in 1:100) {
    temp = (temp + (1 - pchisq(x * j, j))/j)
  }
  print(exp(-temp))
}

t_os(6.745)

# with null hypothesis we use the model with just the three covariates and check
# if the fitted model is appropriate by choosing one additive alternate
# model(having more covariates)
options(scipen = 100)
options(digits = 2)
library(lmtest)
fitlist <- c()
mAIC <- c()
fit0 <- glm(stsfeco ~ Stsf + Trust + Stsf:Trust, family = binomial)
mAIC0 <- (-2) * logLik(fit0) + 6.75 * length(fit0$coefficients)
tos <- NULL
mAIC <- NULL
fitlist <- NULL
library(stats)

tos <- c()
fittemp <- NULL
fitlist <- c()
mAIC <- c()
n <- 1
for (i in 2:4) {
  for (j in 2:4) {
    print(i)
    print(j)
    fittemp <- glm(stsfeco ~ Trust + Stsf + Trust:Stsf + poly(Trust, i) * poly(Stsf, 
                                                                               j), family = "binomial")
    fitlist[[n]] <- fittemp
    tos[[n]] <- (lrtest(fittemp, fit0)$Chisq[2]/length(fittemp$coefficients))
    mAIC[[n]] <- (-2) * logLik(fittemp) + 6.75 * length(fittemp$coefficients)
    n <- n + 1
  }
}
par(mfrow = c(1, 1))
which(mAIC == min(mAIC))

max(tos)
fitlist[which((tos[] == max(tos)))]

# Here we will add succesively the polynomial terms(one at a time) for Trust and
# Stsf till the fourth degree and compute Tn,OS . None of the said fitted models
# had Tn,OS greater than the critical value of 6.75. Hence we can conclude that
# the expected satisfaction with the present state of the economy is a function
# of only the covariates given in themodel represented in null hypothesis.

### Prediction for full model
fitconf <- glm(stsfeco ~ Stsf + Trust + Pol + agea + eduyrs + wkhtot + wkhsch + pdwrk + 
                 chldhm + hhmmb + gincdif + gndr, family = binomial)


# prediction male
datamale <- data.frame(stsfeco = NA, Stsf = mean(Stsf), Trust = mean(Trust), Pol = mean(Pol), 
                       agea = mean(agea), eduyrs = mean(eduyrs), wkhtot = mean(wkhtot), wkhsch = mean(wkhsch), 
                       pdwrk = as.factor(0), chldhm = as.factor(2), hhmmb = mean(hhmmb), gincdif = as.factor(2), 
                       gndr = as.factor(1))

fitpredict <- predict(fitconf, newdata = datamale, se.fit = T)
upper = fitpredict$fit + 1.96 * fitpredict$se.fit
lower = fitpredict$fit - 1.96 * fitpredict$se.fit
plogis(upper)
plogis(lower)
# [0.57,0.75]

# prediction female
datafemale <- data.frame(stsfeco = NA, Stsf = mean(Stsf), Trust = mean(Trust), Pol = mean(Pol), 
                         agea = mean(agea), eduyrs = mean(eduyrs), wkhtot = mean(wkhtot), wkhsch = mean(wkhsch), 
                         pdwrk = as.factor(0), chldhm = as.factor(2), hhmmb = mean(hhmmb), gincdif = as.factor(2), 
                         gndr = as.factor(2))

fitpredict <- predict(fitconf, newdata = datafemale, se.fit = T)
upper = fitpredict$fit + 1.96 * fitpredict$se.fit
lower = fitpredict$fit - 1.96 * fitpredict$se.fit
plogis(upper)
plogis(lower)
# [0.59,0.76]


### Prediction for model through AIC

fitAIC <- glm(stsfeco ~ Trust + agea + gndr + eduyrs + wkhtot + wkhsch + pdwrk + 
                hhmmb + wkhsch:hhmmb + agea:pdwrk + Trust:wkhtot, family = binomial)
summary(fitAIC)

# prediction male

fitpredict <- predict(fitAIC, newdata = datamale, se.fit = T)
upper = fitpredict$fit + 1.96 * fitpredict$se.fit
lower = fitpredict$fit - 1.96 * fitpredict$se.fit
plogis(upper)
plogis(lower)
# [0.61,0.76]


# prediction female

fitpredict <- predict(fitAIC, newdata = datafemale, se.fit = T)
upper = fitpredict$fit + 1.96 * fitpredict$se.fit
lower = fitpredict$fit - 1.96 * fitpredict$se.fit
plogis(upper)
plogis(lower)
# [0.64,0.78]

fittry <- glm(stsfeco ~ Stsf + Trust + Trust:Stsf, family = binomial)
summary(fittry)

plogis(32.99) 