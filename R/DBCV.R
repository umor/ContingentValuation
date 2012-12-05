

DBCV <- function(x, ...) UseMethod("DBCV")


DBCV.default <- function(x,y,z, data, initpar, method, functionalForm, ...)
{
  # some checks:
  
  # check if bid1 and bid2 are numeric
  if(!all(is.numeric(x[,1]))){
    stop(paste("Error: Check your data. At least one element of '", colnames(x)[1], "' is not numeric.", sep=""))}
  if(!all(is.numeric(x[,2]))){
    stop(paste("Error: Check your data. At least one element of '", colnames(x)[2], "' is not numeric.", sep=""))}
  # check if answer1 and answer2 are numeric
  if(!all(is.numeric(y[,1]))){
    stop(paste("Error: Check your data. At least one element of'", colnames(y)[1], "' is not numeric.", sep=""))}
  if(!all(is.numeric(y[,2]))){
    stop(paste("Error: Check your data. At least one element of '", colnames(y)[2], "' is not numeric.", sep=""))}
  # check if answer1 and answer2 are only 0 and 1
  if(!all((y[,1]==1|y[,1]==0))){
    stop(paste("Error: Check your data. '", colnames(y)[1], "' can only have the values 0 or 1.", sep=""))}
  if(!all((y[,2]==1|y[,2]==0))){
    stop(paste("Error: Check your data. '", colnames(y)[2], "' can only have the values 0 or 1.", sep=""))}
  # check if bid1<bid2 for yn&ny
  if(!all(abs(x[(y[,1]==1),1])<abs(x[(y[,1]==1),2]))){
    stop(paste("Error: Check your data. At least once '", colnames(x)[1],"' > '", colnames(x)[2], "' for '",
               colnames(y)[1],"' = 1.", sep=""))  }
  # and bid1>bid2 for ny&nn
  if(!all(abs(x[(y[,1]==0),1])>abs(x[(y[,1]==0),2]))){
    stop(paste("Error: Check your data. At least onece '", colnames(x)[1],"' < '", colnames(x)[2], "' for '",
               colnames(y)[1],"' = 0.", sep=""))}
  
  
  est<-DBCVest(x,y,z, data, initpar, method, functionalForm)
  #est$fitted.values<-
  #est$residuals <- 
  est$call<-match.call
  class(est)<-"DBCV"
  est
}

print.DBCV<-function(x,...)
{ 
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}


summary.DBCV <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object)/se
  TAB <- cbind (Estimate = coef(object),
                StdErr   = se,
                t.value  = tval,
                p.value  = 2*pt(-abs(tval), df=object$df))
  res <- list(call=object$call, 
              coefficients=TAB)
  class(res) <- "summary.DBCV"
  res
  
}  

print.summary.DBCV <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
} 



DBCV.formula <- function(formula, data=list(), initpar=NULL, method=NULL,
                         functionalForm="linear", ...)
{
  # mf <- model.frame(formula=formula, data=data, drop.unused.levels=TRUE)
  if(is.Formula(formula)==FALSE)
  {formula<-Formula(formula)}

  if(functionalForm=="loglinear"){
    
    temp<-model.part(formula,lhs=0,rhs=1, data=data)
    data$llinInc1 <- log((temp[,3]-temp[,1])/temp[,3])
    data$llinInc2 <- log((temp[,3]-temp[,2])/temp[,3])
    
    
    formula  <- paste(
      as.character(formula(formula,lhs=1, rhs=0))[2], " ~ llinInc1 +  llinInc2| ",
      as.character(formula(formula,lhs=0, rhs=2))[2], sep="")
    formula<-formula(formula)
    
    if(is.Formula(formula)==FALSE)
    {formula<-Formula(formula)}
    
    mf <-  model.frame(formula, data = data, 
                       drop.unused.levels = TRUE, na.action="na.omit")
    #attr(mf, "na.action")
    
    y  <-   model.part(formula, data = mf, lhs = 1, drop=TRUE)
    x  <- model.matrix(formula, data = mf, rhs = 1)[,2:3]
    z  <- model.matrix(formula, data = mf, rhs = 2)    
    if(!is.null(attr(mf, "na.action"))){
      temp <- temp[-attr(mf, "na.action"),]
      data <- data[-attr(mf, "na.action"),]
      }
    
    data<-data.frame(data, temp)
    rm(temp)  
    }else{  
  
  mf <- model.frame(formula, data = data, 
                    drop.unused.levels = TRUE, na.action="na.omit")
  y  <- model.part(formula, data = mf, lhs = 1)
  x  <- model.matrix(formula, data = mf, rhs = 1)[,-1]
  z  <- model.matrix(formula, data = mf, rhs = 2)
    }
  
  est <- DBCV.default(x,y,z, data=data, initpar, method, functionalForm, ...)
  est$call <-match.call()
  est$formula <- formula
  est
}


DBCVest<-function(x,y,z, data, initpar, method, functionalForm)  # y= (yes1, yes2) x=(bid1,bid2), z=covariates
{
  
  
  # prepare dummies
  yes1<-y[,1]
  yes2<-y[,2]
  
  bid1 <-x[,1]
  bid2 <-x[,2]
  #z<-t(z)
  
  
  
  #covars <-names(z)
  
  ll<-function(coeff){
    b    <- coeff[-c(ncol(z)+1)]
    #rho  <- coeff[ncol(z)+1]
    a  <- coeff[ncol(z)+1]
    
    #calculate the proabability of each decision sequence, dependent on 
    # the parameters we want to estimate
    
    # variance is set to 1 in the probit model (see. eg. Cameron Trived p. 476)   
    pyy<- pyn <- pny <- pnn <-   rep(NA, times=length(yes1))
    
    #      pyy[yes1==1&yes2==1]<-     pnorm(crossprod(b,t(z[yes1==1&yes2==1,]))/rho - (1/rho*bid2[yes1==1&yes2==1]))    
    #      pyn[yes1==1&yes2==0]<-     pnorm(crossprod(b,t(z[yes1==1&yes2==0,]))/rho - (1/rho*bid1[yes1==1&yes2==0])) - pnorm(crossprod(b,t(z[yes1==1&yes2==0,]))/rho - (1/rho*bid2[yes1==1&yes2==0]))
    #      pny[yes1==0&yes2==1]<-     pnorm(crossprod(b,t(z[yes1==0&yes2==1,]))/rho - (1/rho*bid2[yes1==0&yes2==1])) - pnorm(crossprod(b,t(z[yes1==0&yes2==1,]))/rho - (1/rho*bid1[yes1==0&yes2==1]))
    #      pnn[yes1==0&yes2==0]<- 1 - pnorm(crossprod(b,t(z[yes1==0&yes2==0,]))/rho - (1/rho*bid2[yes1==0&yes2==0]))
    
     if(functionalForm=="linear"){
      pyy[yes1==1&yes2==1]<-     pnorm(crossprod(b,t(z[yes1==1&yes2==1,])) + (a*bid2[yes1==1&yes2==1]))    
      pyn[yes1==1&yes2==0]<-     pnorm(crossprod(b,t(z[yes1==1&yes2==0,])) + (a*bid1[yes1==1&yes2==0])) - pnorm(crossprod(b,t(z[yes1==1&yes2==0,])) + (a*bid2[yes1==1&yes2==0]))
      pny[yes1==0&yes2==1]<-     pnorm(crossprod(b,t(z[yes1==0&yes2==1,])) + (a*bid2[yes1==0&yes2==1])) - pnorm(crossprod(b,t(z[yes1==0&yes2==1,])) + (a*bid1[yes1==0&yes2==1]))
      pnn[yes1==0&yes2==0]<- 1 - pnorm(crossprod(b,t(z[yes1==0&yes2==0,])) + (a*bid2[yes1==0&yes2==0]))
    }
    
    if(functionalForm=="loglinear"){
      pyy[yes1==1&yes2==1]<-     pnorm(crossprod(b,t(z[yes1==1&yes2==1,])) - (a*bid2[yes1==1&yes2==1]))    
      pyn[yes1==1&yes2==0]<-     pnorm(crossprod(b,t(z[yes1==1&yes2==0,])) - (a*bid1[yes1==1&yes2==0])) - pnorm(crossprod(b,t(z[yes1==1&yes2==0,])) - (a*bid2[yes1==1&yes2==0]))
      pny[yes1==0&yes2==1]<-     pnorm(crossprod(b,t(z[yes1==0&yes2==1,])) - (a*bid2[yes1==0&yes2==1])) - pnorm(crossprod(b,t(z[yes1==0&yes2==1,])) - (a*bid1[yes1==0&yes2==1]))
      pnn[yes1==0&yes2==0]<- 1 - pnorm(crossprod(b,t(z[yes1==0&yes2==0,])) - (a*bid2[yes1==0&yes2==0]))
    }
    
    #now tell the function what it should produce as output.
    # we want the probability of an event, given a choice of parameters 
    # (alpha, roh)  as output.
    return(sum(
      log(pyy[yes1==1&yes2==1]),
      log(pyn[yes1==1&yes2==0]),
      log(pny[yes1==0&yes2==1]),
      log(pnn[yes1==0&yes2==0])
    ))   
  }
  
  if(is.null(initpar)){
    
  equInitPar<-"yes2~1+bid2"
  for(lauf in 2:length(colnames(z))){
    equInitPar <- paste(equInitPar,colnames(z)[lauf], sep="+")  
  }
  
  startMod<-glm( equInitPar ,  family = binomial(link = "probit"),data=data.frame(z, bid2, yes2))
  #startMod<-lm( equInitPar, data=data.frame(z, bid2, yes2) )
  
  #  initpar<-c(startMod$coef[-2]/startMod$coef[2],-1/startMod$coef[2])
  initpar<-c(startMod$coef[-2],startMod$coef[2])
  }
  
  if(is.null(method)){  method<-"Newton-Raphson"}
  #result<-  optim(initpar, fn=ll, method="CG", control=list(fnscale=-1, trace=0), hessian=TRUE)
  result<-maxLik(ll, start=initpar, method=method)
  
  
  coef    <- result$estimate
  names(coef)<-c(colnames(z), paste(colnames(x)[1],"=",colnames(x)[2],sep=""))
  df      <- nrow(y)-length(coef)
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


