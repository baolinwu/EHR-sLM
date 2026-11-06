# sLM
 - Conduct covariate-adjusted linear regression from stratified summary statistics in aggregate electronic health record databases
 - Reference
   - Wu et. al (2026) Exact covariate-adjusted linear regression from stratified summary statistics in aggregate electronic health record databases. tech rep.
 - R scripts for conducting analyses and simulation studies
```R

################################################
#### sLM() and related functions
################################################
source("slm_functions.R")


################################################
#### Analysis of data
################################################
## ConY.txt: Simulated aggregate summary data
##  one row per stratum defined by
##    Covid status (1/0 for Yes/No) and the covariates: Sex (Male/Female), Race (White/Black/Hispanic), Age (18-30, 31-50, 51-64, 65-100)
##    Avg/SD/N: sample mean, standard deviation, and sample size in each stratum

LOS0 = read.table('ConY.txt', header=TRUE)
## start with full factorial of all covariate interactions.
Fint3 = slm( ~ Covid + Sex*Race*Age, data=LOS0, n = "N", mean = "Avg", sd = "SD")
## lack of fit test, assessing whether the constant Covid effect is reasonable
slm_lack_of_fit(Fint3) 
## test 3-way interactions
A = cbind(matrix(0,6,19), diag(6))
slm_ftest(Fint3, t(A))
##
Fint2 = slm( ~ Covid + Sex*Race + Sex*Age + Race*Age, data=LOS0, n = "N", mean = "Avg", sd = "SD")
slm_lack_of_fit(Fint2)
## test all pairwise interactions
A2 = cbind(matrix(0,11,8), diag(11))
slm_ftest(Fint2, t(A2))
## test individual interactions
## Sex:Race
B1 = rbind(matrix(0,8,2), diag(2), matrix(0,9,2))
slm_ftest(Fint2, B1)
## Sex:Age
B2 = rbind(matrix(0,10,3), diag(3), matrix(0,6,3))
slm_ftest(Fint2, B2)
## Race:Age
B3 = rbind(matrix(0,13,6), diag(6))
slm_ftest(Fint2, B3)

## final fitted model
Srg = slm( ~ Covid + Sex*Race + Sex*Age, data=LOS0, n = "N", mean = "Avg", sd = "SD")
slm_lack_of_fit(Srg)
Srg

## Weighted LS, t-test
naive_wls( ~ Covid + Sex*Race + Sex*Age, data=LOS0, n = "N", mean = "Avg", sd = "SD")
crude_ttest(LOS0, group='Covid', n = "N", mean = "Avg", sd = "SD", level = 0.95) 

################################################
#### Simulation study
################################################
### Simulation based on the stratum data
## set.seed(1)
DF0 = LOS0$N-1
sig = sqrt(Srg$sigma2)
eta = drop(Srg$X%*%Srg$coef)
LOSj = LOS0
S2 = rchisq(48,DF0)/DF0*sig^2
Y1 = eta + rnorm(48)*sig/sqrt(DF0)
LOSj$Avg = Y1; LOSj$SD = sqrt(S2)
slm( ~ Covid + Sex*Race + Sex*Age, data=LOSj, n = "N", mean = "Avg", sd = "SD")

################################################
#### Equivalence of lm() and slm()
################################################
### equivalence of sLM vs lm() on individual-level data
Ni = LOS0$N
NN = sum(Ni)
## generate individual level data
sid = rep(1:48, times=Ni)
etai = rep(eta, times=Ni)
Yi = etai + rnorm(NN)*sig
## lm() using individual level data
LOSi = LOS0[sid,]
LOSi$Yi = Yi
summary( lm(Yi ~ Covid + Sex*Race + Sex*Age, data=LOSi) )
## sLM() using summary level data
S2 = as.vector( tapply(Yi, sid, var) )
Y1 = as.vector( tapply(Yi, sid, mean) )
LOSj = LOS0
LOSj$Avg = Y1; LOSj$SD = sqrt(S2)
slm( ~ Covid + Sex*Race + Sex*Age, data=LOSj, n = "N", mean = "Avg", sd = "SD")
```
