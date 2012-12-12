

\name{SBCV}
\alias{SBCV}
\alias{SBCV.default}
\alias{SBCV.formula}
\alias{summary.SBCV}
\alias{print.SBCV}
\alias{print.summary.SBCV}

\title{
Single Bounded Contingent Valuation 
}

\description{
The function SBCV estimates coefficients from a single-bounded contingent valuation survey. It basically is a probit regression.
}

\usage{
SBCV(x, ...)
\method{SBCV}{default}(x,y,z, data, initpar, method, functionalForm, ...)
\method{SBCV}{formula}(formula, data=list(), initpar=NULL, method=NULL, functionalForm="linear", ...)
\method{summary}{SBCV}(object, ...)
\method{print}{SBCV}(x, ...)
\method{print}{summary.SBCV}(x, ...)
}
\arguments{
  
  \item{y}{a vector with the answer to the 'yes/no' question. It must be 1 if the answer was 'yes' and and 0 if the answer was 'no'.}
  \item{x}{a vector the the bids.}
  \item{z}{additional, subject specific, explaining variables}
  \item{formula}{a symbolic description of the model to estimat. See details below.}
  \item{data}{an otional data.frame containing the variables in the model.}
  \item{initpar}{a vector of inital parameters for logLik() function from the logLik-package. The vector consist of an initial parameter for the bid, the intercept and each of the explaining variables, in this order. If NULL, the initial values are estimated with the glm() function from the stats-package.}
  \item{method}{a character determining the optimization method in the function logLik(). If NULL, "Newton-Raphson" method is used. For alternatives and details see the documentation of logLik()}
  \item{functionalForm}{Either "linear" or "loglinearRUM".}
  \item{object}{an object of class \code{"SBCV"}, i.e., a fitted model}
  \item{\dots}{
   currently not used 
  }
}
\details{
The formula must be given in the following forms: For linear models
\deqn{yes ~ bid | z1 + z2 + ... + zn}{yes ~ bid | z1 + z2 + ... + zn}
where \emph{yes} is the reply to the survey question. This numeric variable are 1 if the subject answered 'yes' and 0 if the subject answered 'no'.
For loglinearRUM models 
\deqn{yes ~ bid + income | z1 + z2 + ... + zn}{yes ~ bid + income | z1 + z2 + ... + zn}
where a bid variable of the form 
\deqn{log(\frac{income-bid}{income})}{log((income-bid)/income)}
is generated.
\emph{bid} is a numeric variables which is the bid offered to the subject.  \emph{z1} to\emph{zn} are addtional explainig, subject specific, variables. 
 
For identification the variance is fixed to one.

}
\value{An object of class \code{"SBCV"}, basically a list containing
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

data("maizeWTP")

# example
maizeSBCV<- SBCV(yuYes2~yuBid2|age+gender+experimenter , data=maizeWTP)
maizeSBCV


}

\keyword{ single bounded contingent valuation}
