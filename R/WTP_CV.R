

WTP_CV <- function(object, 
                   disp.pref =  c("mean", "median"), # identical for linear functional form
                   disp.subj = c("individual", "mean"),
                   newdata = NULL,
                   CI      = "KrinskyRobb", # "none" if don't want CIs
                   probs   = c(0.05, 0.95), 
                   reps    = 10000
){
  
  if((disp.pref=="mean"|disp.pref=="median")==FALSE){
    stop("'disp.pref' must either be 'mean' or 'median'.")
  }
  if((disp.subj=="individual"|disp.subj=="mean")==FALSE){
    stop("'disp.subj' must either be 'individual' or 'mean'.")
  }
  
  ##########
  
  ####  
  if(disp.subj == c("individual")){
    if(is.null(newdata)){    
      z  <- model.matrix(object$formula, data = object$model, rhs = 2) 
    }
    if(!is.null(newdata)){
      #          if(isTRUE(all.equal(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),colnames(newdata) ))==FALSE){
      #            stop(paste("The data.frame in the argument 'newdata' must have the following named columns:", 
      #                       paste(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),sep="", collapse=", "),
      #                       ". Alternatively 'newdata' can be 'NULL', then the data from the model estimation are used." ,  sep="")
      #            )}
      z  <- model.matrix(object$formula, data = newdata, rhs = 2) 
    }      
  }
  ####  
  if(disp.subj == c("mean")){
    if(is.null(newdata)){    
      z  <- colMeans(model.matrix(object$formula, data = object$model, rhs = 2) )
    }
    if(!is.null(newdata)){
      #          if(isTRUE(all.equal(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),colnames(newdata) ))==FALSE){
      #            stop(paste("The data.frame in the argument 'newdata' must have the following named columns:", 
      #                       paste(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),sep="", collapse=", "),
      #                       ". Alternatively 'newdata' can be 'NULL', then the data from the model estimation are used." ,  sep="")
      #            )}
      z  <-  colMeans(model.matrix(object$formula, data = newdata, rhs = 2) )
    }
    z<-matrix(z,nrow=1)
  }
  
  
  
  if(CI=="KrinskyRobb"){
    b<-object$coefficients
    k<-length(b)
    r<-ncol(z)
    if(disp.subj == c("mean")){
      WTPsim <- matrix(NA, nrow = 1      , ncol = reps)}
    if(disp.subj == c("individual")){
      WTPsim <- matrix(NA, nrow =nrow(z) , ncol = reps)}
    C<-t(chol(object$vcov))
  }
  
  
  ########################################################### 
  # linear
  
  if(object$functionalForm=="linear"){
    #if((disp.pref=="mean"|disp.pref=="median")==TRUE){ # for linear mean and median are identical
    result <-  (z%*%object$coefficients[1:(ncol(z))])/object$coefficients[ncol(z)+1]   
    #     result <-  (z%*%object$coefficients[1:(length(z))])/object$coefficients[length(z)+1]
    
    
    if(CI=="KrinskyRobb"){    
      for(i in 1:reps){
        b2<- b + C%*%rnorm(k)
        WTPsim[,i] <-  (z%*%b2[1:r])/b2[r+1]
      }
      quants <- t( apply(X=WTPsim, MARGIN=1 , FUN=quantile, probs=probs ))
    }
  }     

  #
  ########################################################### 
  # loglinearWTP
  
  if(object$functionalForm=="loglinearWTP"){

    # loglinearWTP mean
    if((disp.pref=="mean")==TRUE){ # for loglinearWTP mean and median are different
      result <-  exp((z%*%object$coefficients[1:ncol(z)])/object$coefficients[ncol(z)+1]+0.5*(-1/object$coefficients[ncol(z)+1])^2)
      if(CI=="KrinskyRobb"){    
        for(i in 1:reps){
          b2<- b + C%*%rnorm(k)
          WTPsim[,i] <-  exp((z%*%b2[1:r])/b2[r+1]+0.5*(-1/b2[r+1])^2)
        }
        quants <- t( apply(X=WTPsim, MARGIN=1 , FUN=quantile, probs=probs ))
      }
      
    }
    # loglinearWTP median
    if((disp.pref=="median")==TRUE){ # for loglinearWTP mean and median are different
      result <-  exp((z%*%object$coefficients[1:ncol(z)])/object$coefficients[ncol(z)+1])
      if(CI=="KrinskyRobb"){    
        for(i in 1:reps){
          b2<- b + C%*%rnorm(k)
          WTPsim[,i] <- exp((z%*%b2[1:r])/b2[r+1]) 
        }
        quants <- t( apply(X=WTPsim, MARGIN=1 , FUN=quantile, probs=probs ))
      }
    }
    
  }     
  
  #
  ########################################################### 
  # loglinearRUM
  
  
  if(object$functionalForm=="loglinearRUM"){
    # generate income and meanIncome variable
    if(is.null(newdata)){
      income <- object$model[,ncol(object$model)]
      if(disp.subj=="mean"){
        income <- mean(income)}
    }
    if(!is.null(newdata)){
      income <-newdata[,ncol(newdata)]
      if(disp.subj=="mean"){
        income <- mean(income)}
    }
    
    # loglinearRUM mean
    if((disp.pref=="mean")==TRUE){ # for loglinearRUM mean and median are different
      result <-  income-income*exp(-(z%*%object$coefficients[1:ncol(z)])  /object$coefficients[ncol(z)+1  ]+
                                     0.5*(1/((object$coefficients[ncol(z)+1  ])^2)))   
      if(CI=="KrinskyRobb"){    
        for(i in 1:reps){
          b2<- b + C%*%rnorm(k)
          WTPsim[,i] <-  income-income*exp(-(z%*%b2[1:r])/b2[r+1] + 0.5*(1/((b2[r+1])^2)))   
        }
        quants <- t( apply(X=WTPsim, MARGIN=1 , FUN=quantile, probs=probs ))
      }
      
    }
    # loglinearRUM median
    if((disp.pref=="median")==TRUE){ # for loglinearRUM mean and median are different
      result <- income- income*exp(-(z%*%object$coefficients[1:(ncol(z))])/object$coefficients[ncol(z)+1])   
      if(CI=="KrinskyRobb"){    
        for(i in 1:reps){
          b2<- b + C%*%rnorm(k)
          WTPsim[,i] <-  income-income*exp(-(z%*%b2[1:r])/b2[r+1] + 0.5*(1/((b2[r+1])^2)))   
        }
        quants <- t( apply(X=WTPsim, MARGIN=1 , FUN=quantile, probs=probs ))
      }
    }
    
  }
  
  
  
  
  
  ###############################################   
  
  if(CI=="KrinskyRobb"){result<- cbind(result,quants)}
  return(result)  
}     


