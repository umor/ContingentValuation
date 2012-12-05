


WTP_CV <- function(object, 
                   disp.pref =  c("mean", "median"), # identical for linear functional form
                   disp.subj = c("individual", "mean"),
                   newdata = NULL
                   ){

  if((disp.pref=="mean"|disp.pref=="median")==FALSE){
    stop("'disp.pref' must either be 'mean' or 'median'.")
  }
  if((disp.subj=="individual"|disp.subj=="mean")==FALSE){
    stop("'disp.subj' must either be 'individual' or 'mean'.")
  }
  
  #
  ########################################################### 
  # linear
  
  if(object$functionalForm=="linear"){
     if((disp.pref=="mean"|disp.pref=="median")==TRUE){ # for linear mean and median are identical
         if(disp.subj=="individual"){
            if(is.null(newdata)){
            z  <- model.matrix(object$formula, data = object$model, rhs = 2) 
           }
            if(!is.null(newdata)){
#              if(isTRUE(all.equal(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),colnames(newdata) ))==FALSE){
#                stop(paste("The data.frame in the argument 'newdata' must have the following named columns:", 
#                         paste(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),sep="", collapse=", "),
#                           ". Alternatively 'newdata' can be 'NULL', then the data from the model estimation are used." ,  sep="")
#                )}
             z  <- newdata 
             z  <- as.matrix(z)  
              }
           result <-  -(z%*%object$coefficients[1:(ncol(z))])/object$coefficients[ncol(z)+1]   
         }

         
         if(disp.subj=="mean"){
           if(is.null(newdata)){
            z  <- colMeans(model.matrix(object$formula, data = object$model, rhs = 2))
           }
           if(!is.null(newdata)){
#             if(isTRUE(all.equal(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),colnames(newdata) ))==FALSE){
#               stop(paste("The data.frame in the argument 'newdata' must have the following named columns:", 
#                          paste(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),sep="", collapse=", "),
#                          ". Alternatively 'newdata' can be 'NULL', then the data from the model estimation are used." ,  sep="")
#               )} 
             z  <- colMeans(newdata) }
             z  <- matrix(z, nrow=1)
          result <-  -(z%*%object$coefficients[1:(length(z))])/object$coefficients[length(z)+1]
           
         }
        } 
    }     

  #
  ########################################################### 
  # loglinear
  
  
  if(object$functionalForm=="loglinear"){
    # loglinear mean
    if((disp.pref=="mean")==TRUE){ # for loglinear mean and median are different
      if(disp.subj=="individual"){
        if(is.null(newdata)){
          z  <- model.matrix(object$formula, data = object$model, rhs = 2)
          income <- object$model[,ncol(object$model)]
        }
        if(!is.null(newdata)){
#          if(isTRUE(all.equal(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),colnames(newdata) ))==FALSE){
#            stop(paste("The data.frame in the argument 'newdata' must have the following named columns:", 
#                       paste(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),sep="", collapse=", "),
#                       ". Alternatively 'newdata' can be 'NULL', then the data from the model estimation are used." ,  sep="")
#            )}
          z  <- newdata[,-ncol(newdata)] 
          z  <- as.matrix(z)
          income <-newdata[,ncol(newdata)]
          }
         
        result <-  income-income*exp(-(z%*%object$coefficients[1:(ncol(z))])  /object$coefficients[ncol(z)+1  ]+0.5*(1/((object$coefficients[ncol(z)+1  ])^2)))   
      }
     
      if(disp.subj=="mean"){
        
        if(is.null(newdata)){
        z  <- colMeans(model.matrix(object$formula, data = object$model, rhs = 2))
        incomeMean <- mean(object$model[, ncol(object$model)])
        }       
        if(!is.null(newdata)){
#          if(isTRUE(all.equal(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),colnames(newdata) ))==FALSE){
#            stop(paste("The data.frame in the argument 'newdata' must have the following named columns:", 
#                       paste(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),sep="", collapse=", "),
#                       ". Alternatively 'newdata' can be 'NULL', then the data from the model estimation are used." ,  sep="")
#            )}
          z  <- colMeans(newdata[,-ncol(newdata)])
          z  <- matrix(z, nrow=1)
          incomeMean <- mean(newdata[, ncol(newdata)])
          }
     result <-  incomeMean-incomeMean*exp(-(z%*%object$coefficients[1:(length(z))])/object$coefficients[length(z)+1]+
                                       0.5*(1/((object$coefficients[length(z)+1])^2)))
        
      }
    }
   # loglinear median
    if((disp.pref=="median")==TRUE){ # for loglinear mean and median are different
      if(disp.subj=="individual"){
        if(is.null(newdata)){
          z  <- model.matrix(object$formula, data = object$model, rhs = 2) 
          income <- object$model[,ncol(object$model)]
          }
        if(!is.null(newdata)){
#          if(isTRUE(all.equal(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),colnames(newdata) ))==FALSE){
#            stop(paste("The data.frame in the argument 'newdata' must have the following named columns:", 
#                       paste(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),sep="", collapse=", "),
#                       ". Alternatively 'newdata' can be 'NULL', then the data from the model estimation are used." ,  sep="")
#            )}
          z  <- newdata[,-ncol(newdata)] 
          z  <- as.matrix(z)
          income <-newdata[,ncol(newdata)]
        }

        if(length(income)!=nrow(z)){
          stop(paste("The length of '",income, "' must be as long as the number of rows of  the data.", sep=""))
        }
        result <- income- income*exp(-(z%*%object$coefficients[1:(ncol(z))])/object$coefficients[ncol(z)+1])   
      }
      if(disp.subj=="mean"){
        
        if(is.null(newdata)){
        z  <- colMeans(model.matrix(object$formula, data = object$model, rhs = 2))
        incomeMean <- mean(object$model[, ncol(object$model)])
        }
        if(!is.null(newdata)){
#          if(isTRUE(all.equal(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),colnames(newdata) ))==FALSE){
#            stop(paste("The data.frame in the argument 'newdata' must have the following named columns:", 
#                       paste(colnames(model.matrix(object$formula, data = object$model, rhs = 2)),sep="", collapse=", "),
#                       ". Alternatively 'newdata' can be 'NULL', then the data from the model estimation are used." ,  sep="")
#            )}
          z  <- newdata[,-ncol(newdata)] 
          z  <- as.matrix(z)
          incomeMean <- newdata[,ncol(newdata)]
              
        }        
        
        result <- incomeMean- incomeMean*exp(-(z%*%object$coefficients[1:(length(z))])/object$coefficients[length(z)+1])
        
      }
    }
  
  }
    
 return(result)  
}     
     
     


