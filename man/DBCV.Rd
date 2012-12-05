\name{DBCV}
\alias{DBCV}
\alias{DBCV.default}
\alias{DBCV.formula}
\alias{summary.DBCV}
\alias{print.DBCV}
\alias{print.summary.DBCV}

\title{
Double Bounded Contingent Valuation 
}

\description{
The function DBCV estimates coefficients from a double-bounded contingent valuation survey. It basically is a bivariat probit regression with the restriction that all the coefficients of both questions are identical.
}

\usage{
DBCV(x, ...)
\method{DBCV}{default}(x,y,z, data, initpar, method, functionalForm, ...)
\method{DBCV}{formula}(formula, data=list(), initpar=NULL, method=NULL, functionalForm="linear", ...)
\method{summary}{DBCV}(object, ...)
\method{print}{DBCV}(x, ...)
\method{print}{summary.DBCV}(x, ...)
}
\arguments{
  
  \item{y}{a matrix with two columns. The first column is 1 of the first question was answered with 'yes' and and 0 if it was answered with 'no'. The second column is 1 of the second question was answered with 'yes' and and 0 if it was answered with 'no'.}
  \item{x}{a matrix with two columns. The first column is the offerd bid in the first question, the second column is the offered bid in the second question.}
  \item{z}{additional, subject specific, explaining variables}
  \item{formula}{a symbolic description of the model to estimat. See details below.}
  \item{data}{an otional data.frame containing the variables in the model.}
  \item{initpar}{a vector of inital parameters for logLik() function from the logLik-package. The vector consist of an initial parameter for the bid, the intercept and each of the explaining variables, in this order. If NULL, the initial values are estimated with the glm() function from the stats-package with the second answer as dependent variable and the second bid offered as bid variable.}
  \item{method}{a character determining the optimization method in the function logLik(). If NULL, "Newton-Raphson" method is used. For alternatives and details see the documentation of logLik()}
  \item{functionalForm}{Either "linear" or "loglinear".}
  \item{object}{an object of class \code{"DBCV"}, i.e., a fitted model}
  \item{\dots}{
   currently not used 
  }
}
\details{
The formula must be given in the followin form: For functionalForm = "linear"
\deqn{yes1 + yes2 ~ bid1 + bid2 | z1 + z2 + ... + zn}{yes1 + yes2 ~ bid1 + bid2 | z1 + z2 + ... + zn}
where \emph{yes1} and \emph{yes2} are the reply to the first and the second question. These numeric variables are 1 if the subject answered 'yes' and 0 if the subject answered 'no'. \emph{bid1} and \emph{bid2} are numeric variables which provide the first and the second bid offered to the subject. For subjects which answered 'yes' to the first bid, \eqn{bid1<bid2}{bid1<bid2}. for subjects which answered 'no' to the first bid, \eqn{bid1>bid2}{bid1>bid2}. \emph{z1} to \emph{zn} are addtional explaining, subject specific variables. 
 For functionalForm = "loglinear" For loglinear models 
\deqn{yes1 + yes2 ~ bid1 + bid2 + income | z1 + z2 + ... + zn}{yes ~ bid + income | z1 + z2 + ... + zn}
where bid variables of the form 
\deqn{log(\frac{income-bid}{income})}{log((income-bid)/income)}
are generated.
  
For identification the variance is fixed to one.

}
\value{An object of class \code{"DBCV"}, basically a list containing
  \item{coefficients}{a named vector of coefficients}
  \item{vcov}{covariance matrix of coefficients}
  \item{LogLik}{log-likelihood value}
  \item{df}{Degrees of freedom}
  \item{model}{data used to estimate the function}
  \item{functionalForm}{The functional form used for the estimation}
}
\references{
Haab, T.C. and McConnell, K.E. (2003), \emph{Valuing Environmental and Natural Resources. The Econometrics of non-market Valuation.} Cheltenham, UK: Edward Elgar 

Morawetz, U.B., DeGroote, H. and Kimenju, S.C. (2011) Improving the Use of Experimental Auctions in Africa: Theory and Evidence, \emph{Journal of Agriculturl and Resource Economics}, 36(2):263-279.
}
\author{
Ulrich B. Morawetz
}
%\note{
%%  ~~further notes~~
%}


%\seealso{ 
%  \code{\link{SemiParBIVProbit}}, \code{\link{glm}}. 
%}

\examples{
 data("maizeWTP")
 require(Formula)

 maizeDBCV<- 
  DBCV(yuYes1+yuYes2~yuBid1+yuBid2|age+gender+experimenter,
  data=maizeWTP)
 maizeDBCV
 summary(maizeDBCV)

}

\keyword{ double bounded contingent valuation}
\keyword{ bivariate probit model}
