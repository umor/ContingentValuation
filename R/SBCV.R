


#SBCV<- function(formula, data, ... ){
#
# if(is.Formula(formula)==FALSE){  
# formula  <- Formula(formula)
# }
# 
# formulaglm  <-  paste(
# as.character(formula(formula, lhs=1))[2],  
# paste(as.character(formula(formula, lhs = 0, rhs = 2)), collapse=""), " + ",
# as.character(formula(formula, lhs = 0, rhs = 1))[2], sep="")
# 
# 
# result <-  glm(formulaglm,  family = binomial(link = "probit"),data=data, ...)
# 
# result$formula<-formula
# 
# return(result)
#}


SBCV <- function(x, ...) UseMethod("SBCV")

SBCV.default <- function(x,y,z, data, initpar, method, ...)
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
  
  
  est<-SBCVest(x,y,z, data, initpar, method)
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
  # mf <- model.frame(formula=formula, data=data, drop.unused.levels=TRUE)
  if(is.Formula(formula)==FALSE)
  {formula<-Formula(formula)}
  
  if(functionalForm=="loglinear"){
    
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
    bidinc <- temp[-attr(mf, "na.action"),]}else{
      bidinc <- temp  
    } 
    rm(temp)
    data<-data.frame(data[-attr(mf, "na.action"),],bidinc)
  }else{
  
  mf <- model.frame(formula, data = data, 
                    drop.unused.levels = TRUE, na.action="na.omit")
  #attr(mf, "na.action")
  
  y  <- model.part(formula, data = mf, lhs = 1, drop=TRUE)
  x  <- model.matrix(formula, data = mf, rhs = 1)[,2]
  z  <- model.matrix(formula, data = mf, rhs = 2)
  
  }
    
  est <- SBCV.default(x,y,z, data=mf, initpar, method,  ...)
  est$call <-match.call()
  est$formula <- formula
  est
}


SBCVest<-function(x,y,z, data, initpar, method)  # y= (yes) x=(bid), z=covariates
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
    
    
  #  if(a>0){
  #    py[yes1==1]<-     pnorm(crossprod(b,t(z[yes1==1,])) + (a*bid1[yes1==1]))    
  #    pn[yes1==0]<- 1 - pnorm(crossprod(b,t(z[yes1==0,])) + (a*bid1[yes1==0]))
  #  }else{
      py[yes1==1]<-     pnorm(crossprod(b,t(z[yes1==1,])) - (a*bid1[yes1==1]))    
      pn[yes1==0]<- 1 - pnorm(crossprod(b,t(z[yes1==0,])) - (a*bid1[yes1==0]))
  #  }
    
    
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
    
    startMod<-glm( equInitPar ,  family = binomial(link = "probit"),data=data.frame(z, bid1, yes1))
    #startMod<-lm( equInitPar, data=data.frame(z, bid2, yes2) )
    
    #  initpar<-c(startMod$coef[-2]/startMod$coef[2],-1/startMod$coef[2])
    initpar<-c(startMod$coef[-2], - startMod$coef[2])
  }
  
  if(is.null(method)){  method<-"Newton-Raphson"}
  #result<-  optim(initpar, fn=ll, method="CG", control=list(fnscale=-1, trace=0), hessian=TRUE)
  result<-maxLik(ll, start=initpar, method=method)
  
 #  browser()
  
  coef    <- result$estimate
  names(coef)<-c(colnames(z), colnames(data)[2])
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
       model = data)
  
}

