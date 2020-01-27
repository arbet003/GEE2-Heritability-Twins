##########################################################################=
##########################################################################
#Example1: 
#
#-Estimating heritability for a heavy-tailed t-distributed trait
#-No mean-level covariate effects
#-No variance level covariate effects (i.e. we are not letting herititability differ by any covariates)
##########################################################################
##########################################################################


load("example1data.RData") #load dataset with heavy tailed t-distributed trait
#load("final_lgpsim1.RData")
#colnames(simdata)=c("zygosity","famid","twin","outcome")
source('twinACEmodels.R') #load functions to fit GEE2-Falconer, GEE2-NACE, and traditional Falconer method

################
# GEE2-Falconer
################

gee2f=gee2falconer(data=simdata,zygo="zygosity",twin_number = "twin",famID="famid",Y="outcome",cilevel=0.95)
gee2f$heritability


################
#Traditional Falconer's method (no standard errors or confidence intervals are reported)
################
falc=falconer(data=simdata,zygo="zygosity",twin_number = "twin",famID="famid",Y="outcome")
falc

###############
#GEE2-NACE
###############
gee2n=gee2nace(data=simdata,zygo="zygosity",twin_number = "twin",famID="famid",Y="outcome",tol=0.001)
gee2n$heritability

###############
#NACE
#
#The mets R package twinlm() function can be used to fit the NACE model.
###############
library("mets")
NACE=twinlm(outcome~1,data=simdata,id="famid",zyg="zygosity",DZ=0)

summary(NACE)$acde  #h2,c2,e2

summary(NACE)$estimate # mean level covariates effects, and standard deviations of ACE variance parameters (note gee2nace() reports the variances "s2a, s2c, s2e" and not standard deviations)
