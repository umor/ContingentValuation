



SBCV <- function(x, ...) UseMethod("SBCV")

SBCV.default <- function(x,y,z, data, initpar, method, functionalForm, ...)
{
  # some checks:
  
  # check if bid is numeric
  if(!all(is.numeric(x))){
    stop(paste("Error: Check your data. At least one element of '", colnames(x)[1], "' is not numeric.", sep=""))}
  # check if answer is numeric
  if(!all(is.numeric(y))){
    stop(paste("Error: Check your data. At least one element of'", colnames(y)[1], "' is not numeric.", sep=""))}
  # check if answer is only 0 and 1
  if(!all((y==1|y==0))){
    stop(paste("Error: Check your data. '", colnames(y)[1], "' can only have the values 0 or 1.", sep=""))}
  
  
  est<-SBCVest(x,y,z, data, initpar, method, functionalForm)
  #est$fitted.values<-
  #est$residuals <- 
  est$call<-match.call
  class(est)<-"SBCV"
  est
}

print.SBCV<-function(x,...)
{ 
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}


summary.SBCV <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object)/se
  TAB <- cbind (Estimate = coef(object),
                 StdErr   = se,
                t.value  = tval,
                p.value  = 2*pt(-abs(tval), df=object$df)
                )
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.SBCV"
  res
  
}  

print.summary.SBCV <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
} 

#formula<-yuYes2~yuBid2|age+gender+experimenter
#data<-maizeWTP

SBCV.formula <- function(formula, data=list(), initpar=NULL, method=NULL, 
                         functionalForm="linear",  ...)
{ 
  if(is.Formula(formula)==FALSE)
  {formula<-Formula(formula)}
  
### some checks  
  # check if first right hand side of
  # the formula consist of two elements for loglinearRUM and of one for linear
  if(functionalForm=="loglinearRUM"){
    if(length(strsplit(as.character(formula(formula, rhs=1, lhs=0))[2], 
                       "+", fixed=TRUE)[[1]])!=2){
      stop(cat("If functionalForm is 'loglinearRUM' the first part of the right hand 
        side of the equation must have two elements"))
    }}

  if(functionalForm=="linear"){
    if(length(strsplit(as.character(formula(formula, rhs=1, lhs=0))[2], 
                       "+", fixed=TRUE)[[1]])!=1){
      stop(cat("If functionalForm is 'linear' the first part of the right hand 
        side of the equation must have one element"))
    }}

  ### data preparation if linear   
  
  if(functionalForm=="linear"){
  
  mf <- model.frame(formula, data = data, 
                    drop.unused.levels = TRUE, na.action="na.omit")
  #attr(mf, "na.action")
  
  y  <- model.part(formula, data = mf, lhs = 1, drop=TRUE)
  x  <- model.matrix(formula, data = mf, rhs = 1)[,2]
  z  <- model.matrix(formula, data = mf, rhs = 2)
  
  if(!is.null(attr(mf, "na.action"))){
  data <- data[-attr(mf, "na.action"),]}
  
  attr(x,"varName") <-  names(data)[which(names(data)==
                                            as.character(formula(formula, lhs=0 , 
                                                                 rhs=1))[2])] 
  }
  
 #####  data preparation if loglinearWTP
  if(functionalForm=="loglinearWTP"){

    temp<-model.part(formula,lhs=0,rhs=1, data=data)
    if(!all(temp>0, na.rm=TRUE)){
      stop(cat("At least one bid is smaller or equal zero. This is a problem
               if functionalForm is loglinearWTP."))
    }
      
     mf <- model.frame(formula, data = data, 
                      drop.unused.levels = TRUE, na.action="na.omit")
    #attr(mf, "na.action")
    
    y  <- model.part(formula, data = mf, lhs = 1, drop=TRUE)
    x  <- model.matrix(formula, data = mf, rhs = 1)[,2]
    x  <- log(x)
    z  <- model.matrix(formula, data = mf, rhs = 2)
    
    if(!is.null(attr(mf, "na.action"))){
    data <- data[-attr(mf, "na.action"),]}
    # add log to data
    bidVarPos<-which(names(data)==as.character(formula(formula, lhs=0 , rhs=1))[2])
    data[, ncol(data)+1]<-log(data[, bidVarPos])
    names(data)[ncol(data)] <- paste("log(",names(data)[bidVarPos],")" , sep="")
    attr(x,"varName") <-  names(data)[ncol(data)] 
  }
  
  
  #### data preparation if loglinearRUM
  
  if(functionalForm=="loglinearRUM"){
    
    temp<-model.part(formula,lhs=0,rhs=1, data=data)
    data$llinInc <- log((temp[,2]-temp[,1])/temp[,2])
    
    
    formula  <- paste(
      as.character(formula(formula,lhs=1, rhs=0))[2], " ~ llinInc | ",
      as.character(formula(formula,lhs=0, rhs=2))[2], sep="")
    formula<-formula(formula)
    
    if(is.Formula(formula)==FALSE)
    {formula<-Formula(formula)}
    
    mf <-  model.frame(formula, data = data, 
                       drop.unused.levels = TRUE, na.action="na.omit")
    #attr(mf, "na.action")
    
    y  <-   model.part(formula, data = mf, lhs = 1, drop=TRUE)
    x  <- model.matrix(formula, data = mf, rhs = 1)[,2]
    z  <- model.matrix(formula, data = mf, rhs = 2)    
    if(!is.null(attr(mf, "na.action"))){
      temp <- temp[-attr(mf, "na.action"),]
      data <- data[-attr(mf, "na.action"),]
  attr(x,"varName") <- paste("log((", names(temp)[2],
               "-",names(temp)[1],")/", names(temp)[2], ")" , sep="")
    }
    
    data<-data.frame(data, temp)
    rm(temp)  
  }
  
  
  
  est <- SBCV.default(x,y,z, data=data, initpar, method, functionalForm, ...)
  est$call <-match.call()
  est$formula <- formula
  est
}


SBCVest<-function(x,y,z, data, initpar, method, functionalForm)  # y= (yes) x=(bid), z=covariates
{
  
  
  # prepare dummies
  yes1<-y
  #colnames(yes1)<-"yes1"
  
  bid1 <-x
  #z<-t(z)
  
  
  
  #covars <-names(z)
  
  ll<-function(coeff){
    b    <- coeff[-c(ncol(z)+1)]
    #rho  <- coeff[ncol(z)+1]
    a  <- coeff[ncol(z)+1]
    
    #calculate the proabability of each decision sequence, dependent on 
    # the parameters we want to estimate
    
    # variance is set to 1 in the probit model (see. eg. Cameron Trived p. 476)   
    py<- pn <- rep(NA, times=length(yes1))
    
    
    if(functionalForm=="loglinearRUM"){
      py[yes1==1]<-     pnorm(crossprod(b,t(z[yes1==1,])) + (a*bid1[yes1==1]))    
      pn[yes1==0]<- 1 - pnorm(crossprod(b,t(z[yes1==0,])) + (a*bid1[yes1==0]))
    }
    if(functionalForm=="linear" | functionalForm=="loglinearWTP"){
      py[yes1==1]<-     pnorm(crossprod(b,t(z[yes1==1,])) - (a*bid1[yes1==1]))    
      pn[yes1==0]<- 1 - pnorm(crossprod(b,t(z[yes1==0,])) - (a*bid1[yes1==0]))
    }
    
    
    #now tell the function what it should produce as output.
    # we want the probability of an event, given a choice of parameters 
    # (alpha, roh)  as output.
    return(sum(
      log(py[yes1==1]),
      log(pn[yes1==0])
    ))   
  }
  
  if(is.null(initpar)){
    
    equInitPar<-"yes1~1+bid1"
    for(lauf in 2:length(colnames(z))){
      equInitPar <- paste(equInitPar,colnames(z)[lauf], sep="+")  
    }
    
    #startMod<-glm( equInitPar ,  family = binomial(link = "probit"),data=data.frame(z, bid1, yes1))
    startMod<-lm( equInitPar, data=data.frame(z, bid1, yes1) )
    
    #  initpar<-c(startMod$coef[-2]/startMod$coef[2],-1/startMod$coef[2])
    initpar<-c(startMod$coef[-2], - startMod$coef[2])
  }
  
  if(is.null(method)){  method<-"Newton-Raphson"}
  #result<-  optim(initpar, fn=ll, method="CG", control=list(fnscale=-1, trace=0), hessian=TRUE)
  result<-maxLik(ll, start=initpar, method=method)
  
 #  browser()
  
  coef    <- result$estimate
  names(coef)<-c(colnames(z), attr(x, "varName") )
  df      <- length(y)-length(coef)
  LogLik <- result$maximum
  #hessian <- result$hessian
  hessian <- numDeriv::hessian(func=ll,x=result$estimate)
  #vcov    <-solve(-result$hessian)
  vcov    <-solve(-hessian)
  colnames(vcov)<-rownames(vcov)<-names(coef)
  
  list(coefficients = coef, 
       vcov = vcov,
       LogLik = LogLik,
       df= df,
       model = data,
       functionalForm = functionalForm)
  
}

