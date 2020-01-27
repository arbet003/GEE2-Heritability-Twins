library("geepack") #geese() function used for GEE2-Falconer
library("msm") #delta method used to construct confidence intervals 


######################################################################
######################################################################
#GEE2-Falconer
######################################################################
######################################################################


gee2falconer=function(data,zygo,twin_number,famID,Y,meanformula=as.formula(paste0(Y,"~1")),cilevel=0.95,cormat=NULL){
  #gee2falconer() is the proposed GEE2-Falconer model from Arbet et al 2020, 'A Robust and Unified Framework for Estimating Heritability in Twin Studies using Generalized Estimating Equations'
  #
  #data: the dataset with column names used for arguments zygo, twin_number, famID, and Y
  #zygo: column name of variable that contains the zygosity for each twin pair (1=MZ, 0=DZ)
  #twin_number: column name of variable identifying twin 1 and twin 2 of each twin pair.  Should only take on values 1 and 2.
  #famID: column name of numeric variable giving the family ID for each twin pair.  The dataset will be sorted by famID to ensure both twins in a pair have adjacent rows in the dataset.
  #Y: column name of outcome/trait
  #meanformula: a formula for mean level covariate effects, should be written in quotes: e.g. "Y ~ 1 + zygo + sex + age".... etc
  #   by default, meanformula is "Y ~ 1", i.e. no covariate effects, assumes the same population mean for MZ and DZ twins.  "Y~1+zygo" would allow the mean to differ by zygosity, but you could include other covariat effects like sex, age, etc.
  #ci.level: numeric variable between 0 and 1 used to construct cilevel*100% confidence intervals for h2,c2,e2
  
  
  ####################error checks
  
  #check to make sure zygo, twin_number, famID, Y are column names in the dataset
  if(!all(c(zygo,twin_number,famID,Y)%in%colnames(data))){
    stop("arguments zygo, twin_number, famID, and Y should refer to column names in 'data'")
  }
  
  #check famID is numeric
  if(is.numeric(data[,famID])==F){
    stop("famID should be the column name of a numeric variable identifying the family ID for each twin pair.  The dataset will be sorted by famID to ensure both twins in a pair have adjacent rows in the dataset.")
  }
  
  #check twin_number is numeric and only contains values 1 and 2
  if(!is.numeric(data[,twin_number])| !1%in%data[,twin_number] | !2%in%data[,twin_number] | length(table(data[,twin_number]))>2){
    stop("twin_number should be the column name of a numeric variable identifying twin 1 and twin 2 in each twin pair.  Should only take on values 1 and 2.")
  }
  
  #check zygo is numeric and only contains 0 and 1 (1=MZ, 0=DZ)
  if(!is.numeric(data[,zygo])| !1%in%data[,zygo] | !0%in%data[,zygo] | length(table(data[,zygo]))>2){
    stop("zygo should be the column name of a numeric variable identifying the zygosity of each twin pair: 1=MZ, 0=DZ.  Should only contain values 1 and 0.")
  }
  
  #check that each twin pair has 2 twins
  if(all(table(data[,famID])!=2)){
    stop("each twin pair should have 2 twins")
  }
  
  #check missing data
  if(sum(is.na(data[,c(zygo,twin_number,famID,Y)]))>0){
    stop("Currently this function does to support missing values in variables zygo, twin_number, famID, or Y")
  }
  
  #check ci.level is numeric and between 0 and 1
  if(!is.numeric(cilevel)){
    stop("cilevel should be a numeric variable between 0 and 1 used to calcualte cilevel*100 % confidence intervals for h2,c2,e2")
  }
  if(cilevel<0 | cilevel >1){
    stop("cilevel should be a numeric variable between 0 and 1 used to calcualte cilevel*100 % confidence intervals for h2,c2,e2")
  }
  
  #if user supplies cormat to allow twin correlation to differ by covariates, 
  #run checks for the rules specified in example3_variancelevel_covariates.R tutorial:
  if(!is.null(cormat)){
    if(nrow(cormat)!=length(unique(data[,famID]))){
      stop("The number of rows in cormat must equal the total number of twin pairs.")
    }
    if(all(is.numeric(cormat))==F){
      stop("cormat must be a numeric design matrix, use 0/1 dummy variables for categorical variables.")
    }
  }
  
  #finished error checks
  #############################
  
  ########sort dataset by famID and twin_number
  data=data[order(data[,famID],data[,twin_number]),]
  
  
  ###########cormat: let correlation differ as function of zygosity
  #cormat is used in the geepack R package geese() function 'zcor' argument for userdefined correlation structure
  #cormat represents the design matrix for equation "h(rho_z)" from Arbet et al 2020 Section 2.3.2
  #should have 1 row per twin pair
  #column 1 is the "intercept"
  #column 2 is the zygosity of the twin pair, allows the correlation within a twin pair to differ by MZ/DZ status
  #
  #The current version of gee2falconer() does not support h2,c2,e2 covariate effects.
  #See example3_variancelevel_covariates.R for a tutorial on how to allow h2,c2,e2 covariate effects.
  cormat=vector()
  for(i in 1:(nrow(data)/2)){
    cormat=rbind(cormat,c(1,data[,zygo][(2*i)]))
  }
  
  famid_ind=which(colnames(data)==famID)
  fit=geese(formula=as.formula(meanformula),family="gaussian",sformula=as.formula(paste("~1+",zygo)),cor.link = "fisherz",corstr="userdefined",id=data[,famid_ind],zcor=cormat,data=data,sca.link="log")
  
  ######estimated twin correlations
  rmzhat= (exp(sum(summary(fit)$corr[,1]))-1)/(exp(sum(summary(fit)$corr[,1]))+1)
  rdzhat= (exp(summary(fit)$corr[1,1])-1)/(exp(summary(fit)$corr[1,1])+1)
 
  #estimates of h2 and c2
  h2= 2*(rmzhat-rdzhat)
  c2= 2*rdzhat-rmzhat
  e2=1-rmzhat
  h2se=deltamethod(~ 2*((exp(x1+x2)-1)/(exp(x1+x2)+1)-(exp(x1)-1)/(exp(x1)+1)), mean= summary(fit)$cor[,1],cov=fit$valpha)
  c2se=deltamethod(~ 2*((exp(x1)-1)/(exp(x1)+1))-(exp(x1+x2)-1)/(exp(x1+x2)+1), mean= summary(fit)$cor[,1],cov=fit$valpha)
  e2se=deltamethod(~ 1-(exp(x1+x2)-1)/(exp(x1+x2)+1), mean= summary(fit)$cor[,1],cov=fit$valpha)
  
  est=c(h2=h2,c2=c2,e2=e2)
  se=c(h2se=h2se,c2se=c2se,e2se=e2se)
  low=vector()
  up=vector()
  alpha=1-cilevel
  low=est-abs(qnorm(1-alpha/2))*se
  up=est+abs(qnorm(1-alpha/2))*se
  
  low=ifelse(low<0,0,low)
  up=ifelse(up>1,1,up)
 
  heritability=data.frame(est=est,se=se,lowCI=low,upCI=up)
  meaneffects=summary(fit)$mean
  
  
  ####################stats to return
  #heritability: dataframe giving results for h2, c2, e2, i.e. proportion of outcome variance explained by additive genetic effects, shared family environmental effects, and unique individual env. effects)
  #              est: parameter estimates
  #              se: robust sandwhich standard error estimate
  #              lowCI: lower confidence interval limit
  #              upCI: upper confidence interval limit
  #meaneffects: mean-level covariate effects
  return(list(heritability=heritability,meaneffects=meaneffects))
}

######################################################################
######################################################################
# Traditional Falconer method
######################################################################
######################################################################

falconer=function(data,zygo,twin_number,famID,Y){
  #falconer() is the traditional Falconer's method.  No standard errors or confidence intervals are reported.
  #
  #data: the dataset with column names used for arguments zygo, twin_number, famID, and Y
  #zygo: column name of variable that contains the zygosity for each twin pair (1=MZ, 0=DZ)
  #twin_number: column name of variable identifying twin 1 and twin 2 of each twin pair.  Should only take on values 1 and 2.
  #famID: column name of numeric variable giving the family ID for each twin pair.  The dataset will be sorted by famID to ensure both twins in a pair have adjacent rows in the dataset.
  #Y: column name of outcome/trait
  
  ####################error checks
  
  #check to make sure zygo, twin_number, famID, Y are column names in the dataset
  if(!all(c(zygo,twin_number,famID,Y)%in%colnames(data))){
    stop("arguments zygo, twin_number, famID, and Y should refer to column names in 'data'")
  }
  
  #check famID is numeric
  if(is.numeric(data[,famID])==F){
    stop("famID should be the column name of a numeric variable identifying the family ID for each twin pair.  The dataset will be sorted by famID to ensure both twins in a pair have adjacent rows in the dataset.")
  }
  
  #check twin_number is numeric and only contains values 1 and 2
  if(!is.numeric(data[,twin_number])| !1%in%data[,twin_number] | !2%in%data[,twin_number] | length(table(data[,twin_number]))>2){
    stop("twin_number should be the column name of a numeric variable identifying twin 1 and twin 2 in each twin pair.  Should only take on values 1 and 2.")
  }
  
  #check zygo is numeric and only contains 0 and 1 (1=MZ, 0=DZ)
  if(!is.numeric(data[,zygo])| !1%in%data[,zygo] | !0%in%data[,zygo] | length(table(data[,zygo]))>2){
    stop("zygo should be the column name of a numeric variable identifying the zygosity of each twin pair: 1=MZ, 0=DZ.  Should only contain values 1 and 0.")
  }
  
  #check that each twin pair has 2 twins
  if(all(table(data[,famID])!=2)){
    stop("each twin pair should have 2 twins")
  }
  
  #check missing data
  if(sum(is.na(data[,c(zygo,twin_number,famID,Y)]))>0){
    stop("Currently this function does to support missing values in variables zygo, twin_number, famID, or Y")
  }
  
  #finished error checks
  #############################
  
  ########sort dataset by famID and twin_number
  data=data[order(data[,famID],data[,twin_number]),]
  
  
  ###calculate MZ correlation
  MZ=data[data[,zygo]==1,]
  MZ1=MZ[MZ[,twin_number]==1,]
  MZ2=MZ[MZ[,twin_number]==2,]
  rMZ = cor(MZ1[,Y],MZ2[,Y]) #MZ correlation
  
  ###calculate DZ correlation
  DZ=data[data[,zygo]==0,]
  DZ1=DZ[DZ[,twin_number]==1,]
  DZ2=DZ[DZ[,twin_number]==2,]
  rDZ = cor(DZ1[,Y],DZ2[,Y]) # DZ correlation
  
  ###calculate and report h2, c2, e2, i.e. the proportion of outcome variance explained by additive genetic effects (h2), shared family environmental effects (c2), and unique individual environmental effects (e2)
  h2= 2*(rMZ-rDZ)
  c2= 2*rDZ-rMZ
  e2=1-rMZ
  stats=c(h2=h2,c2=c2,e2=e2)
  return(stats)
}

######################################################################
######################################################################
# GEE2-NACE
#
# -Currently does not allow ACE variance parameters to differ by covariates
######################################################################
######################################################################
#theta_update(): used to get point estimates
#theta_covar(): sandwich covariance estimator to get standard errors
#gee2nace(): wrapper function to fit GEE2-NACE model

####Newton-Raphson type algorithm used to obtain point estimates for GEE2-NACE
theta_update=function(famID,Y,zygo,Xmat,beta,A,C,E){
  #working = "ind" or "norm"
  Xmat=as.matrix(Xmat)
  beta=as.vector(beta)
  fam=table(famID)
  Dk12=matrix(0,2,3)
  Dk21=matrix(0,3,ncol(Xmat))
  
  one=0
  two=0
  for(i in 1:length(fam)){
    ind=which(famID==names(fam)[i])
    temp.y=Y[ind]
    temp.xmat=as.matrix(Xmat[ind,])
    Dk11=temp.xmat
    yhat= temp.xmat%*%beta
    fk1= temp.y-yhat
    temp.z=zygo[ind]
    w=ifelse(temp.z[1]==1,1,0.5)
    Dk22= rbind(c(1,1,1),c(1,1,1),c(w,1,0))
    Dk=rbind(cbind(Dk11,Dk12),cbind(Dk21,Dk22))
    Vk11 = rbind(c(A+C+E,w*A+C),c(w*A+C,A+C+E))
    Vk12 = matrix(0,2,3)
    Vk22= cbind(c(2*(A+C+E)^2,2*(w*A+C)^2,2*(A+C+E)*(w*A+C)),c(2*(w*A+C)^2,2*(A+C+E)^2,2*(A+C+E)*(w*A+C)),
                  c(2*(A+C+E)*(w*A+C),2*(A+C+E)*(w*A+C),(A+C+E)^2+(w*A+C)^2))
    Vk= rbind(cbind(Vk11,Vk12),cbind(t(Vk12),Vk22))
    
    sigmak= c(A+C+E,A+C+E,w*A+C)
    sk=c((fk1[1])^2,fk1[2]^2, fk1[1]*fk1[2])
    fk2= sk -sigmak
    fk=t(t(c(fk1,fk2)))
    one= one + t(Dk)%*%solve(Vk)%*%(Dk)
    two= two + t(Dk)%*%solve(Vk)%*%(fk)
  }
  theta_new= c(beta,A,C,E)+solve(one)%*%two
  rownames(theta_new)=c(paste("beta",0:(ncol(Xmat)-1),sep=""),"A","C","E")
  return(t(theta_new))
}

####sandwich covariance estimator to get standard errors for GEE2-NACE
theta_covar=function(famID,Y,zygo,Xmat,beta,A,C,E){
  Xmat=as.matrix(Xmat)
  fam=table(famID)
  K=length(names(fam))
  Dk12=matrix(0,2,3)
  Dk21=matrix(0,3,ncol(Xmat))
  
  Sigma=0
  inside=0
  for(i in 1:length(fam)){
    ind=which(famID==names(fam)[i])
    temp.y=Y[ind]
    temp.xmat=as.matrix(Xmat[ind,])
    Dk11=temp.xmat
    yhat= temp.xmat%*%beta
    fk1= temp.y-yhat
    temp.z=zygo[ind]
    w=ifelse(temp.z[1]==1,1,0.5)
    Dk22= rbind(c(1,1,1),c(1,1,1),c(w,1,0))
    Dk=rbind(cbind(Dk11,Dk12),cbind(Dk21,Dk22))
    Vk11 = rbind(c(A+C+E,w*A+C),c(w*A+C,A+C+E))
    Vk12 = matrix(0,2,3)
    Vk22= cbind(c(2*(A+C+E)^2,2*(w*A+C)^2,2*(A+C+E)*(w*A+C)),c(2*(w*A+C)^2,2*(A+C+E)^2,2*(A+C+E)*(w*A+C)),
                  c(2*(A+C+E)*(w*A+C),2*(A+C+E)*(w*A+C),(A+C+E)^2+(w*A+C)^2))
    Vk= rbind(cbind(Vk11,Vk12),cbind(t(Vk12),Vk22))
    
    sigmak= c(A+C+E,A+C+E,w*A+C)
    sk=c((fk1[1])^2,fk1[2]^2, fk1[1]*fk1[2])
    fk2= sk -sigmak
    fk=t(t(c(fk1,fk2)))
    Sigma= Sigma+(1/K)*t(Dk)%*%solve(Vk)%*%Dk
    inside= inside + t(Dk)%*%solve(Vk)%*%fk%*%t(fk)%*%solve(Vk)%*%Dk
  }
  covar= (1/K)*solve(Sigma)%*%inside%*%solve(Sigma)
  return(covar/K)
}

gee2nace=function(data,zygo,twin_number,famID,Y,Xmat=data.frame(int=rep(1,nrow(data))),cilevel=0.95,tol=0.001){
  #gee2falconer() is the proposed GEE2-NACE model from Arbet et al 2020, 'A Robust and Unified Framework for Estimating Heritability in Twin Studies using Generalized Estimating Equations'
  #
  #data: the dataset with column names used for arguments zygo, twin_number, famID, and Y
  #zygo: column name of variable that contains the zygosity for each twin pair (1=MZ, 0=DZ)
  #twin_number: column name of variable identifying twin 1 and twin 2 of each twin pair.  Should only take on values 1 and 2.
  #famID: column name of numeric variable giving the family ID for each twin pair.  The dataset will be sorted by famID to ensure both twins in a pair have adjacent rows in the dataset.
  #Y: column name of outcome/trait
  #Xmat: design matrix for mean-level covariates.  Make sure the first column is a vector of 1's for the intercept.
      # I recommend adding column names to help make the final output easier to interpret.
  #ci.level: numeric variable between 0 and 1 used to construct cilevel*100% confidence intervals for h2,c2,e2
  #tol: tolerance for newton-raphson iterative algorithm to get parameter estimates.
  
  ####################error checks
  
  #check to make sure zygo, twin_number, famID, Y are column names in the dataset
  if(!all(c(zygo,twin_number,famID,Y)%in%colnames(data))){
    stop("arguments zygo, twin_number, famID, and Y should refer to column names in 'data'")
  }
  
  #check famID is numeric
  if(is.numeric(data[,famID])==F){
    stop("famID should be the column name of a numeric variable identifying the family ID for each twin pair.  The dataset will be sorted by famID to ensure both twins in a pair have adjacent rows in the dataset.")
  }
  
  #check twin_number is numeric and only contains values 1 and 2
  if(!is.numeric(data[,twin_number])| !1%in%data[,twin_number] | !2%in%data[,twin_number] | length(table(data[,twin_number]))>2){
    stop("twin_number should be the column name of a numeric variable identifying twin 1 and twin 2 in each twin pair.  Should only take on values 1 and 2.")
  }
  
  #check zygo is numeric and only contains 0 and 1 (1=MZ, 0=DZ)
  if(!is.numeric(data[,zygo])| !1%in%data[,zygo] | !0%in%data[,zygo] | length(table(data[,zygo]))>2){
    stop("zygo should be the column name of a numeric variable identifying the zygosity of each twin pair: 1=MZ, 0=DZ.  Should only contain values 1 and 0.")
  }
  
  #check that each twin pair has 2 twins
  if(all(table(data[,famID])!=2)){
    stop("each twin pair should have 2 twins")
  }
  
  #check missing data
  if(sum(is.na(data[,c(zygo,twin_number,famID,Y)]))>0){
    stop("Currently this function does to support missing values in variables zygo, twin_number, famID, or Y")
  }
  
  #check ci.level is numeric and between 0 and 1
  if(!is.numeric(cilevel)){
    stop("cilevel should be a numeric variable between 0 and 1 used to calcualte cilevel*100 % confidence intervals for h2,c2,e2")
  }
  if(cilevel<0 | cilevel >1){
    stop("cilevel should be a numeric variable between 0 and 1 used to calcualte cilevel*100 % confidence intervals for h2,c2,e2")
  }
  
  #finished error checks
  #############################
  
  ########sort dataset by famID and twin_number
  data=data[order(data[,famID],data[,twin_number]),]
  
  #starting values: betas for mean-level covariate effects, s2a, s2c, s2e
  startvals=c(rnorm(ncol(Xmat),sd=sd(data[,Y])/4),rep(var(data[,Y])/3,3))
  est=vector()
  est=rbind(est,startvals)
  z=1
  delta=tol*2
  iter=1
  while(delta>=tol){
    est=rbind(est,theta_update(famID=data[,famID], Y=data[,Y], zygo=data[,zygo],Xmat=Xmat,beta=est[nrow(est),1:ncol(Xmat)],A=est[nrow(est),ncol(Xmat)+1],C=est[nrow(est),ncol(Xmat)+2],E=est[nrow(est),ncol(Xmat)+3]))
    
    delta = sum(abs(est[nrow(est),]-est[nrow(est)-1,]))
    if(delta<tol){
      print(paste0("estimation converged (delta<tol), ", "delta=",round(delta,nchar(strsplit(as.character(tol), "\\.")[[1]][2])+1)))
    } else {
      print(paste0("Finished iteration ",iter, ", delta=",round(delta,nchar(strsplit(as.character(tol), "\\.")[[1]][2])+1)))
      iter=iter+1 
    }
  }
  
  #converged estimated
  estc=est[nrow(est),]
  
  #sandwich covariance for parameter estimates
  estc_var=theta_covar(famID=data[,famID], Y=data[,Y], zygo=data[,zygo],Xmat=Xmat,beta=estc[1:ncol(Xmat)],A=estc[ncol(Xmat)+1],C=estc[ncol(Xmat)+2],E=estc[ncol(Xmat)+3])
  
  ACE=estc[(ncol(Xmat)+1):length(estc)]
  h2= ACE[1]/sum(ACE)
  c2= ACE[2]/sum(ACE)
  e2= ACE[3]/sum(ACE)
  
  Aind=ncol(Xmat)+1
  Cind=ncol(Xmat)+2
  Eind=ncol(Xmat)+3
  h2se=deltamethod(as.formula(paste0("~x",Aind,"/(","x",Aind,"+","x",Cind,"+x",Eind,")")),estc,estc_var)
  c2se=deltamethod(as.formula(paste0("~x",Cind,"/(","x",Aind,"+","x",Cind,"+x",Eind,")")),estc,estc_var)
  e2se=deltamethod(as.formula(paste0("~x",Eind,"/(","x",Aind,"+","x",Cind,"+x",Eind,")")),estc,estc_var)
  
  finalest=c(estc,h2,c2,e2)
  names(finalest)[(length(finalest)-2):length(finalest)]=c("h2","c2","e2")
  se=c(sqrt(diag(estc_var)),h2se,c2se,e2se)


  low=vector()
  up=vector()
  alpha=1-cilevel
  low=finalest-abs(qnorm(1-alpha/2))*se
  up=finalest+abs(qnorm(1-alpha/2))*se
  
  low[(length(low)-2):length(low)]=ifelse(low[(length(low)-2):length(low)]<0,0,low[(length(low)-2):length(low)])
  up[(length(up)-2):length(up)]=ifelse(up[(length(up)-2):length(up)]>1,1,up[(length(up)-2):length(up)])
   
temp=data.frame(est=finalest,se=se,lowCI=low,upCI=up)
  meaneffects=temp[1:ncol(Xmat),]
  if(is.null(colnames(Xmat))){
    colnames(Xmat)=paste0("beta",0:(ncol(Xmat)-1))
  } 
  rownames(meaneffects)=colnames(Xmat)
  heritability=temp[(ncol(Xmat)+1):nrow(temp),]
  rownames(heritability)=c("s2a","s2c","s2e","h2","c2","e2")
  
  ####################stats to return
  #heritability: dataframe giving results for s2a, s2c, s2e (i.e. the ACE variance parameters), and h2, c2, e2, i.e. proportion of outcome variance explained by additive genetic effects, shared family environmental effects, and unique individual env. effects)
  #              est: parameter estimates
  #              se: robust sandwhich standard error estimate
  #              lowCI: lower confidence interval limit
  #              upCI: upper confidence interval limit
  #meaneffects: mean-level covariate effects
  return(list(heritability=heritability,meaneffects=meaneffects))
}