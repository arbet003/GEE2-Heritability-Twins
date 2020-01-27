##########################################################################=
##########################################################################
#Example3: allowing heritability to vary as a function of covariates
#
#-In this example, we allow the the ACE variance parameters (and thus the proportions h2,c2,e2)
# to differ by a binary variable, sex (males, females)
#
# We will walk you through a step by step example of how to estimate h2,c2,e2 as a 
# function of covariates using the geese() function from the geepack R package.  
# Currently the gee2falconer() function does not support automating this process for 
# an arbitrary set of covariates, but we are working on developing this.
##########################################################################
##########################################################################


load("example3data.RData") 
library("geepack")
library("msm")

################################################
# GEE2-Falconer
################################################

#GEE2-Falconer can allow the correlation between a pair of twins to differ as a function 
#of covariates (see Arbet et al 2020 Section 2.3.2).  Recall that Falconer's estimators are 
#functions of the MZ correlation and DZ correlation.  Thus allowing the twin correlation to 
#differ by covariates also allows heritability h2, and environmental effects c2, e2 to 
#differ by those covariates.
#
# To do this, you need to create a matrix "cormat" of the variables you want to allow the 
# twin correlations to differ by.  There are a few simple rules for creating this matrix:
# 1) The number of rows in cormat should equal the total number of twin pairs in the dataset.
# 2) The 1st column should be a vector of ones (i.e. an intercept)
# 3) The 2nd column should be the zygosity of the twin pair (1 for MZ, 0 for DZ)
# 4) Additional columns represent other variables you want the twin correlation to differ by.
     # Thus far we have only considered variables that have the same value for 
     # both twins within a given pair (e.g. sex, age, age^2, family-level variables, etc.).  
# 5) all variables must be coded numerically (e.g. 0/1 dummy variables for categorical variables).
# 6) For each variable you add, e.g. sex, you also need to include an interaction with zygosity, e.g. sex*zygosity (as explained in Arbet et al 2020 Section 2.3.2)
# 7) p-value thresholds (e.g. alpha=0.05, or multiple testing adjusted p-values using the 
#    p.adjust R function) can be used to determine which variables to include.  
#    But note, for any variable you include "x", you need to also include the 
#    x*zygosity interaction, even if the later has a large p-value.
    

###############Step 1: create 'cormat' matrix to allow heritability to differ by sex

#make sure dataset is ordered by famID (numeric ID identifying each twin pair)
simdata=simdata[order(simdata$famid),]


#the number of rows in cormat must equal the total number of twin pairs (i.e. nrow(dataset)/2)
#each twin pair is represented by a single row.  the columns give the covariate values for that twin pair.
#In this example, we are allowing the correlation to differ by zygosity, sex, and sex*zygosity.
cormat=vector()
for(i in 1:(nrow(simdata)/2)){
  cormat=rbind(cormat,c(int=1,zyg=simdata$zygosity[2*i],sex=simdata$sex[2*i],sex_zygo=simdata$sex[2*i]*simdata$zygosity[2*i]))
}

###############Step 2: fit GEE2-Falconer model using geese() R function and the user-created cormat
#Note that whatever variables you include in cormat also need to be included in 
# "sformula," which allows the twin variances to differ by these covariates.
#
#for example, here sformula needs to be: sformula~1+zygosity+sex+sex:zygosity
gee2f=geese(outcome~1+sex,family="gaussian",sformula~1+zygosity+sex+sex:zygosity,cor.link = "fisherz",corstr="userdefined",id=famid,zcor=cormat,data=simdata)

#estimates of correlation parameters (we will use these to estimate h2,c2,e2 for males and females)
summary(gee2f)$corr

#sandwich covariance of correlation estimates:
gee2f$valpha

###############Step 3: estimate MZ and DZ twin correlations for males and females
#Basic idea: for a given set of covariate values "X", estimate rmz_x and rdz_x, 
#then plut into Falconer's equations to get h2x, c2x, e2x
#
#####if using cor.link="fisherz", the inverse function for fisherz is: (exp(x)-1)/(exp(x)+1).
#Note the inverse of fisher's z-transformation is usually written as (exp(2*x)-1)/(exp(2*x)+1) but this version is incorrect for geepack. You should not include the 2.
rmz_m= (exp(sum(summary(gee2f)$corr[1:4,1]))-1)/(exp(sum(summary(gee2f)$corr[1:4,1]))+1)
rmz_f=(exp(sum(summary(gee2f)$corr[1:2,1]))-1)/(exp(sum(summary(gee2f)$corr[1:2,1]))+1)
rdz_m= (exp(sum(summary(gee2f)$corr[c(1,3),1]))-1)/(exp(sum(summary(gee2f)$corr[c(1,3),1]))+1)
rdz_f= (exp(sum(summary(gee2f)$corr[c(1),1]))-1)/(exp(sum(summary(gee2f)$corr[c(1),1]))+1)

####if instead using cor.link="identity":
#gee2f_identitylink=geese(outcome~1+sex,family="gaussian",sformula~1+zygosity+sex+sex:zygosity,cor.link = "identity",corstr="userdefined",id=famid,zcor=cormat,data=simdata)
#rmz_m= sum(summary(gee2f_identitylink)$corr[1:4,1])
#rmz_f= sum(summary(gee2f_identitylink)$corr[1:2,1])
#rdz_m= sum(summary(gee2f_identitylink)$corr[c(1,3),1])
#rdz_f= sum(summary(gee2f_identitylink)$corr[c(1),1])

###############Step 4: estimate h2,c2,e2 for males and females
#recall Falconer's formulas: h2=2*(rmz-rdz),  c2=2*rdz-rmz,  e2= 1 - h2 - c2 = 1 - rmz
h2m= 2*(rmz_m-rdz_m)
c2m= 2*rdz_m-rmz_m
e2m= 1-rmz_m
h2f= 2*(rmz_f-rdz_f)
c2f= 2*rdz_f-rmz_f
e2f= 1-rmz_f

###############Step 5: standard errors and CIs for h2m,c2m,e2m, h2f, c2f, e2f using delta-method 
####if using cor.link="fisherz":
h2mse=deltamethod(~ 2*((exp(x1+x2+x3+x4)-1)/(exp(x1+x2+x3+x4)+1)-(exp(x1+x3)-1)/(exp(x1+x3)+1)), mean= summary(gee2f)$cor[,1], gee2f$valpha)
c2mse=deltamethod(~ 2*(exp(x1+x3)-1)/(exp(x1+x3)+1)-(exp(x1+x2+x3+x4)-1)/(exp(x1+x2+x3+x4)+1), mean= summary(gee2f)$cor[,1], gee2f$valpha)
e2mse=deltamethod(~ 1-(exp(x1+x2+x3+x4)-1)/(exp(x1+x2+x3+x4)+1), mean= summary(gee2f)$cor[,1], gee2f$valpha)

h2fse=deltamethod(~ 2*((exp(x1+x2)-1)/(exp(x1+x2)+1)-(exp(x1)-1)/(exp(x1)+1)), mean= summary(gee2f)$cor[,1], gee2f$valpha)
c2fse=deltamethod(~ 2*(exp(x1)-1)/(exp(x1)+1)-(exp(x1+x2)-1)/(exp(x1+x2)+1), mean= summary(gee2f)$cor[,1], gee2f$valpha)
e2fse=deltamethod(~ 1-(exp(x1+x2)-1)/(exp(x1+x2)+1), mean= summary(gee2f)$cor[,1], gee2f$valpha)

####if using cor.link="identity"
#h2mse=deltamethod(~ 2*(x1+x2+x3+x4-(x1+x3)), mean= summary(gee2f_identitylink)$cor[,1], gee2f_identitylink$valpha)
#c2mse=deltamethod(~ 2*(x1+x3)-(x1+x2+x3+x4), mean= summary(gee2f_identitylink)$cor[,1], gee2f_identitylink$valpha)
#e2mse=deltamethod(~ 1-(x1+x2+x3+x4), mean= summary(gee2f_identitylink)$cor[,1], gee2f_identitylink$valpha)
#h2fse=deltamethod(~ 2*(x1+x2-x1), mean= summary(gee2f_identitylink)$cor[,1], gee2f_identitylink$valpha)
#c2fse=deltamethod(~ 2*(x1)-(x1+x2), mean= summary(gee2f_identitylink)$cor[,1], gee2f_identitylink$valpha)
#e2fse=deltamethod(~ 1-(x1+x2), mean= summary(gee2f_identitylink)$cor[,1], gee2f_identitylink$valpha)


###############Step 6: return final estimates, standard errors, confidence intervals 
meaneffects=summary(gee2f)$mean
meaneffects

ci.level=0.95 #confidence interval level
alpha=1-ci.level
est=c(h2m=h2m,c2m=c2m,e2m=e2m,h2f=h2f,c2f=c2f,e2f=e2f)
SE=c(h2mse,c2mse,e2mse,h2fse,c2fse,e2fse)
low=est-abs(qnorm(1-alpha/2))*SE; low=ifelse(low<0,0,low)
up=est+abs(qnorm(1-alpha/2))*SE; up=ifelse(up>1,1,up)

tab=data.frame(est,SE,low,up)
tab
