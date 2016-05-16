setwd("E:/Mstat Courses/Statistical Modelling/Final Project")
# We have been provided with data taken from 97 men before a radical
# prostatectomy. The variables in the dataset mainly consists of clinical
# information (the exact medical interpretation is not important for this
# project).

fulldata.C = read.table("prostate.txt", header = T)
set.seed(479671)
rownumbers = sample(1:length(fulldata.C[, 1]), size = 70)
mydatatr.C = fulldata.C[rownumbers, ]
mydatats.C = fulldata.C[-rownumbers, ]
mydatatr.C <- data.frame(mydatatr.C)
# look at the seed data
summary(mydatatr.C)
attach(mydatatr.C)
names(mydatatr.C)
names(mydatats.C)
# The resposne is log of prostate antigen. The fact that the response is log
# transformed makes us suspicious that it will be normally distributed. test for
# normality
shapiro.test(lpsa)
# we do not have sufficient proof to reject the null hypothesis. We can say that
# the variable lpsa is normally distributed. Model Selections with AIC, AICc, BIC


library(MuMIn)
options(na.action = "na.fail")
# construct models selected through various information criteria. For this
# purpose we fit the following model
fit1 <- lm(lpsa ~ ., data = mydatatr.C)
plot(fit1)

options(na.action = "na.fail")
library(pls)
library(MuMIn)

# AIC and mspe
mAIC <- dredge(fit1, evaluate = T, rank = "AIC")
mAIC
a <- get.models(mAIC, subset = 1)[[1]]
a
dataAIC <- mydatats.C[, -c(6:8)]
temp <- predict(a, newdata = dataAIC, se.fit = T)
mean((mydatats.C$lpsa - temp$fit)^2)


# AICc and mspe
mAICc <- dredge(fit1, evaluate = T, rank = "AICc")
model <- get.models(mAICc, subset = 1)[[1]]
model
dataAICc <- mydatats.C[, c(1, 2, 5)]
temp <- predict(model, newdata = dataAICc, se.fit = T)

mean((mydatats.C$lpsa - temp$fit)^2)

# BICc and mspe
mBIC <- dredge(fit1, evaluate = T, rank = "BIC")
model <- get.models(mBIC, subset = 1)[[1]]
model
dataBIC <- mydatats.C[, c(1, 2, 5)]
temp <- predict(model, newdata = dataBIC, se.fit = T)
mean((mydatats.C$lpsa - temp$fit)^2)

# have to derive Takeuchi´s information criterion TIC that is applicable to
# normal response.
TICnormal <- function(Y, X) {
  n <- length(Y)
  fit <- lm(Y ~ X - 1)
  summary(fit)
  
  
  sigma <- summary(fit)$sigma
  coef <- as.matrix(fit$coefficients)
  yhat <- X %*% coef
  
  # partials for j and k matrices
  deriv1.beta <- X * as.vector(Y - yhat)/sigma^2  # preceeding negative sign ignore for easy computation in matrices
  deriv1.sigmasq <- ((Y - yhat)^2/(2 * sigma^4)) - 1
  deriv2.beta <- t(X) %*% X/sigma^2
  deriv2.sigmasq <- (1/(2 * sigma^4) + (t(Y - yhat) %*% (Y - yhat)))/sigma^6  # preceeding negative sign ignore for easy computation in matrices
  deriv.mixed <- (-1/sigma * 4) * t(Y - yhat) %*% X
  
  kmatrix1 <- t(deriv1.beta) %*% deriv1.beta/n
  print("kmatrix1")
  print(dim(kmatrix1))
  kmatrix3 <- t(deriv1.sigmasq) %*% deriv1.beta/n
  print("kmatrix3")
  print(dim(kmatrix3))
  kmatrix2 <- t(kmatrix3)
  print("kmatrix2")
  print(dim(kmatrix2))
  
  kmatrix4 <- t(deriv1.sigmasq) %*% deriv1.sigmasq/n
  print("kmatrix4")
  print(dim(kmatrix4))
  
  
  kmatrix <- cbind(rbind(kmatrix1, kmatrix3), rbind(kmatrix2, kmatrix4))
  dim(kmatrix)
  n <- length(X[1, ])
  jmatrix1 <- deriv2.beta/n
  print("jmatrix1")
  print(dim(jmatrix1))
  
  jmatrix4 <- deriv2.sigmasq/n
  print("jmatrix4")
  print(dim(jmatrix4))
  
  jmatrix3 <- matrix(rep(0, n), ncol = n)
  print("jmatrix3")
  print(dim(jmatrix3))
  jmatrix2 <- t(jmatrix3)
  print("jmatrix2")
  print(dim(jmatrix2))
  jmatrix <- cbind(rbind(jmatrix1, jmatrix3), rbind(jmatrix2, jmatrix4))
  
  logLikfit = c(logLik(fit))
  
  penaltyparameter = sum(diag(solve(jmatrix) %*% kmatrix))
  
  TIC = 2 * (-logLikfit + penaltyparameter)
  
  
  return(list(TIC = TIC, penaltyparameter = penaltyparameter))
}

# all possible combinations given all the variables
combinations = function(n) {
  comb = NULL
  for (i in 1:n) comb = rbind(cbind(1, comb), cbind(0, comb))
  return(comb)
}

# initializing variables to calculate TIC
combi <- NULL
Y <- lpsa
Xcombi <- as.matrix(cbind(rep(1, length(Y)), mydatatr.C[1:8]))
combi <- combinations(length(Xcombi[1, ]))

# calculate TIC values for all the combinations listed as per the function seen
# before
TICvalues <- c()
# combi<-cbind(rep(1,length(combi[,1])), combi)
dim(combi)
for (i in 1:length(combi[1, ])) {
  temp <- combi[i, ]
  rowvar <- which(temp != 0)
  X <- Xcombi[, rowvar]
  tempTIC <- TICnormal(Y, X)
  TICvalues[i] <- tempTIC$TIC
  
}
temp <- which(TICvalues == min(TICvalues))
min(TICvalues)

combi
temp
# we will proceed to check the models with Mean Square Prediction Error (MSPE).
# MSPE is mean of sum of squared difference between the observed and predicted
# values. msep with TIC
dataTIC <- mydatats.C[, -c(6:7)]
attach(dataTIC)
model_TIC <- lm(lpsa ~ lcavol + lweight + age + lbph + svi + pgg45)
summary(model_TIC)
temp <- predict(model_TIC, newdata = dataTIC, se.fit = T)
mean((mydatats.C$lpsa - temp$fit)^2)
detach(dataTIC)
attach(mydatatr.C)
# FIC with modelling


# to calculate FIC
FIC.normal = function(Y, X, Z, XZeval = cbind(X, Z), outfile = "FICoutput.txt", allsubsets = T, 
                      nested = F) {
  # Y = response vector X = regression variables that are present in all models,
  # always included, explicitly include an intercept if you want one to be present.
  # Z = regression variables to be used in the variable selection method XZeval =
  # Values of X,Z variables for evaluating the FIC, default is all values of (X,Z)
  # outfile: name of the file where should the output should be written.
  # allsubsets: variable search amongst all possible subsets (default) OR nested:
  # only nested sequences, in the same order as given in the Z matrix (first column
  # comes first, then the second column is added, etc.
  
  qq = ncol(as.matrix(Z))
  pp = ncol(as.matrix(X))
  nn = nrow(as.matrix(X))
  
  XZeval = as.matrix(XZeval)
  X.eval = XZeval[, 1:pp]
  X.eval = matrix(X.eval, ncol = pp)
  Z.eval = XZeval[, -(1:pp)]
  Z.eval = matrix(Z.eval, ncol = qq)
  
  
  ## form submatrices of the big J matrix 1/n sum_i xi * t(xi)
  
  ## estimate sigma^2 in biggest model (for example)
  
  full.fit = lm(Y ~ X + Z - 1)  # exclude intercept since already in the X-part
  sigma.hat.sq = as.vector(t(full.fit$residuals) %*% full.fit$residuals/(nn - pp - 
                                                                           qq))
  # pp + qq already includes the intercept, no need to subtract 1 in the df for the
  # sigma^2.
  
  J = matrix(0, nrow = (pp + qq + 1), ncol = (pp + qq + 1))
  J[1, 1] = 2 * nn
  J[(2:(pp + 1)), (2:(pp + 1))] = t(X) %*% X
  J[(pp + 2):(pp + qq + 1), (2:(pp + 1))] = t(Z) %*% X
  J[(2:(pp + 1)), (pp + 2):(pp + qq + 1)] = t(X) %*% Z
  J[(pp + 2):(pp + qq + 1), (pp + 2):(pp + qq + 1)] = t(Z) %*% Z
  
  J = J/(nn * sigma.hat.sq)
  
  J11 = J[(pp + 2):(pp + qq + 1), (pp + 2):(pp + qq + 1)]
  J00 = J[1:(pp + 1), 1:(pp + 1)]
  J10 = J[(pp + 2):(pp + qq + 1), 1:(pp + 1)]
  
  invJ = solve(J)
  K <- invJ[-(1:(pp + 1)), -(1:(pp + 1))]
  
  
  ## estimate (beta,gamma) in biggest model
  
  betagamma.hat = full.fit$coefficients
  
  deltahat = nn^{
    1/2
  } * betagamma.hat[-(1:pp)]
  
  if (allsubsets) 
    varmatrix = combinations(qq)
  if (nested) {
    id = diag(qq)
    id[lower.tri(id)] = 1
    varmatrix = rbind(rep(0, qq), id)
  }
  
  
  ## Main part of the program!
  
  # This matrix will contain per subject the FIC value for each of the models.
  
  n.eval = nrow(Z.eval)
  
  FIC.subject.model = risk.subject.model = bias.subject.model = bias2.subject.model = var.subject.model = var.naive.subject.model = matrix(NA, 
                                                                                                                                           nrow = n.eval, ncol = nrow(varmatrix))
  
  for (k in 1:n.eval) {
    print(k)
    x.chosen = as.vector(c(0, X.eval[k, ]))
    z.chosen = Z.eval[k, ]
    
    partmutheta = as.vector(x.chosen)
    tau0sq = t(partmutheta) %*% solve(J00) %*% partmutheta
    
    omega = J10 %*% solve(J00) %*% x.chosen - z.chosen
    psi.full = t(omega) %*% deltahat
    
    for (j in 1:nrow(varmatrix)) {
      FIC.out = FIC(varmatrix[j, ], Kn = K, omega = omega, deltahat = deltahat, 
                    psi.full = psi.full)
      FIC.subject.model[k, j] = FIC.out$FIC
      
      # write(FIC.subject.model[k,],file=outfile,ncol=nrow(varmatrix),append=T)
      
      risk.subject.model[k, j] = tau0sq - t(omega) %*% K %*% omega + FIC.subject.model[k, 
                                                                                       j]
      bias.subject.model[k, j] = FIC.out$bias.S/nn
      bias2.subject.model[k, j] = FIC.out$bias2.S/nn
      var.subject.model[k, j] = (tau0sq + FIC.out$var2.S)/nn
      var.naive.subject.model[k, j] = (tau0sq + FIC.out$var2.naive.S)/nn
    }  # end loop over j
    NULL
  }  # end loop over k
  FIC.minima = apply(FIC.subject.model, 1, which.min)
  which.vars.selected = varmatrix[FIC.minima, ]
  return(list(FIC = FIC.subject.model, Bias = bias.subject.model, Bias2 = bias2.subject.model, 
              Risk = risk.subject.model, Var = var.subject.model, Varnaive = var.naive.subject.model, 
              SelectedVars = which.vars.selected))
}

# Set of functions to be used in the calculation of the FIC values.  These are
# just help files and not to be directly called.


# All subsets search, make an indicator matrix containing all possible
# combinations:

combinations = function(n) {
  comb = NULL
  for (i in 1:n) comb = rbind(cbind(1, comb), cbind(0, comb))
  return(comb)
}
print(n)

# Example, models with extra variables 2 and 3 variables = c(0,1,1,0,0)


FIC = function(variables, Kn, omega, deltahat, psi.full) {
  if (sum(variables) > 0) {
    qq = length(variables)
    indic = variables * (1:qq)
    
    Id = diag(rep(1, qq))
    pi.S = matrix(Id[indic, ], ncol = qq)
    
    K.S = solve(pi.S %*% solve(Kn) %*% t(pi.S))
    M.S = t(pi.S) %*% K.S %*% pi.S %*% solve(Kn)
    
    psi.S = t(omega) %*% (t(pi.S) %*% K.S %*% pi.S) %*% solve(Kn) %*% deltahat
    
    omega.S = pi.S %*% omega
    
    FIC.S = (psi.full - psi.S)^2 + 2 * t(omega.S) %*% K.S %*% omega.S
    
    bias.S = t(omega) %*% deltahat - psi.S
    
    bias2.S = sign(bias.S) * sqrt(max(0, t(omega) %*% (Id - M.S) %*% (deltahat %*% 
                                                                        t(deltahat) - Kn) %*% (Id - M.S) %*% omega))
    
    var2.S = t(omega) %*% M.S %*% omega
    
    var2.naive.S = t(omega.S) %*% omega.S
  } else {
    FIC.S = psi.full^2
    bias.S = bias2.S = t(omega) %*% deltahat
    var2.S = 0
    var2.naive.S = 0
  }
  
  outputFIC = list(FIC.S, bias.S, bias2.S, var2.S, var2.naive.S)
  names(outputFIC) = c("FIC.S", "bias.S", "bias2.S", "var2.S", "var2.naive.S")
  return(outputFIC)
  
}




FIC.summary = function(FIC.subject.model, qq, allsubsets = T, nested = F) {
  
  # Equal weight to all covariate combinations for model estimation
  FIC.eq.weighted = apply(FIC.subject.model, 2, mean)
  
  if (allsubsets) 
    varmatrix = combinations(qq)
  if (nested) {
    id = diag(qq)
    id[lower.tri(id)] = 1
    varmatrix = rbind(rep(0, qq), id)
  }
  
  print("Global model selection via FIC")
  print("Focus = estimation of the mean function")
  print("All subjects equal weight")
  print("")
  print("variables in the selected model:")
  print(varmatrix[which.min(FIC.eq.weighted), ])
  
  print("Average FIC values for all models:")
  print(FIC.eq.weighted)
  
  print("model number:")
  print(which.min(FIC.eq.weighted))
  
}



summary.FIC.latex = function(FIC.output, X, Z, XZ.eval, outfile = "outputfic.out", 
                             allsubsets = T, nested = F) {
  attach(FIC.output)
  FIC.minima = apply(FIC.output$FIC, 1, which.min)
  if (allsubsets) {
    qq = log(ncol(as.matrix(FIC.output$FIC)), base = 2)
    varmatrix = combinations(qq)
  }
  if (nested) {
    qq = ncol(as.matrix(FIC.output$FIC)) - 1
    id = diag(qq)
    id[lower.tri(id)] = 1
    varmatrix = rbind(rep(0, qq), id)
  }
  varnamesZ = names(as.data.frame(Z))
  varnamesX = names(as.data.frame(X))  #dimnames(X)[[2]]
  help1 = varmatrix[FIC.minima, ]
  
  library(xtable)
  help2 = help1
  help3 = (1:qq) * help2
  summary.table = as.data.frame(cbind(as.vector(sqrt(Varnaive)), as.vector(sqrt(Var)), 
                                      as.vector(Bias), as.vector(Bias2), as.vector(sqrt(Risk)), as.vector(FIC.output$FIC)))
  names(summary.table) = c("se.S", "se.f", "bias1", "bias2", "risk", "FIC")
  sum.table = summary.table[sort.list(summary.table[, 6]), ][1:min(20, nrow(FIC.output$FIC)), 
                                                             ]
  xtab = xtable(sum.table, align = rep("c", ncol(summary.table) + 1))
  digits(xtab) = c(0, 3, 3, 3, 3, 3, 3)
  write(paste(c(varnamesX, varnamesZ), XZ.eval), file = outfile, append = T)
  write(c("* FIC selected model:", c(varnamesX, varnamesZ[help3])), file = outfile, 
        append = T)
  print.xtable(xtab, file = outfile, append = T, table.placement = "!")
  # write('------------------------------------------------------------------------------------------')
  detach(FIC.output)
}

X <- rep(1, length(lpsa))
Z <- cbind(lcavol, lweight, age, lbph, svi, lcp, as.ordered(gleason), pgg45)

Y <- matrix(lpsa, ncol = 1)


FIC <- FIC.normal(Y, X, Z)
FICvalues <- FIC$FIC[1, ]
min(FICvalues)
model.no <- which(FICvalues == min(FICvalues))
FIC.formodel <- FIC[model.no]
model.FIC <- combi[model.no, ]
model.FIC
detach(mydatatr.C)
attach(mydatats.C)
model <- lm(lpsa ~ gleason)
dataFIC <- data.frame(mydatats.C[, 7])
temp <- predict(model, newdata = dataFIC, se.fit = T)
mean((mydatats.C$lpsa - temp$fit)^2)

# select the model for which the first observation is in focus
selectedmodel <- which(model.FIC == 1)
selectedmodel[1:length(selectedmodel)]
print("Selected Model through FIC")
colnames(Z)[c(selectedmodel[1:length(selectedmodel)])]