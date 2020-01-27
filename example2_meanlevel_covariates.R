##########################################################################=
##########################################################################
#Example2: 
#
#-Estimating heritability for a heavy-tailed t-distributed trait where 
#-the outcome mean differs by zygosity (MZ/DZ status).  The true population mean is 5 for MZ and 0 for DZ
#-No variance level covariate effects (i.e. we are not letting herititability differ by any covariates)
##########################################################################
##########################################################################


load("example2data.RData") #load dataset with heavy tailed t-distributed trait. Population mean is 5 for MZ and 0 for DZ twins.
source('twinACEmodels.R') #load functions to fit GEE2-Falconer, GEE2-NACE, and traditional Falconer method

################
# GEE2-Falconer
################

gee2f=gee2falconer(data=simdata,zygo="zygosity",twin_number = "twin",famID="famid",Y="outcome",meanformula="outcome~zygosity",cilevel=0.95)
gee2f$heritability
gee2f$meaneffects #DZ mean = Int; MZ mean = + zygosity


################
#Traditional Falconer's method (no covariate effects, no standard errors or confidence intervals are reported)
################
falc=falconer(data=simdata,zygo="zygosity",twin_number = "twin",famID="famid",Y="outcome")
falc

###############
#GEE2-NACE
###############
Xmat=data.frame(Int=1,zygosity=simdata$zygosity)
gee2n=gee2nace(data=simdata,zygo="zygosity",twin_number = "twin"
               ,famID="famid",Y="outcome",Xmat=Xmat)
gee2n$heritability
gee2n$meaneffects #DZ mean = Int; MZ mean = + zygosity

###############
#NACE
#
#The mets R package twinlm() function can be used to fit the NACE model.
###############

#twinlm() throws an error if you include the same variable 
#used in the 'zyg' argument also in the mean-level effects formula statement.
#Therefore, to allow the mean to differ by zygosity, you can simply create a new zygosity 
#variable 'zygosity2' to use in mean-level effects formula statement.
library("mets")
simdata$zygosity2=simdata$zygosity 
NACE=twinlm(outcome~1+zygosity2,data=simdata,id="famid",zyg="zygosity",DZ=0)

summary(NACE)$acde  #h2,c2,e2

summary(NACE)$estimate # mean level covariates effects, and standard deviations of ACE variance parameters (note gee2nace() reports the variances "s2a, s2c, s2e" and not standard deviations)
#DZ mean = "outcome" (this is the intercept)
#MZ mean = "outcome" + "outcome~zygosity2" (i.e. intercept + coefficient for zygosity)