\name{WTP_CV}
\alias{WTP_CV}

\title{
Willingness to Pay for Contingent Valuation 
}

\description{
The function WTP_CV() calculates the coefficients from SBCV or DBCV models.
}

\usage{
WTP_CV (object, 
        disp.pref =  c("mean", "median"), 
        disp.subj = c("individual", "mean"),
        newdata = NULL,
        CI      = "KrinskyRobb", 
         probs   = c(0.05, 0.95), 
        reps    = 10000
        )
}
\arguments{
  
  \item{object}{An object created by the estimation of SBCV or DBCV}
  \item{disp.pref}{Either \code{mean} or \code{median}.}   
  \item{disp.subj}{Either \code{individual} or \code{mean}.} 
  \item{newdata}{Either \code{NULL} or a data.frame}
  \item{CI}{Either "KrinskyRobb" or "none".}
  \item{probs}{A numeric vector of of probabilites with values in [0,1]}
  \item{reps}{Nuber of draws for KrinskyRobb CI}
  \item{\dots}{ currently not used  }
}
\details{
The calculation of the willigness to pay is for probit \code{SBCV} and \code{DBCV}. 


Dispersion of the preferences, \code{disp.pref}, is either \code{mean} or \code{median} depending on which measure of central tendency of the preverences should be used. In case of  \code{functionalForm}  is \code{linear} the resulting willingness to pay is the same for \code{mean} and \code{median}.

\code{disp.subj} is either \code{individual} or \code{mean}. If \code{individual} and \code{newdata} is  \code{NULL}, then the willigness to pay is calculated for each person in the data used to estimate the object. If \code{individual} and in \code{newdata} has a \code{data.frame}, then the willigness to pay is calculated for these data. If \code{mean} the mean of the data provided in  \code{newdata} (or, if \code{NULL}, all the the data used  to estimate the object) is calculated and the willigness to pay is calculated for the mean subject. 

If \code{newdata} is \code{NULL}, all subjects used for the estimation of the object are used, if a \code{data.frame} is given, these data are used. The \code{data.frame} must consist of all the variables used to estimate the object. They must be of the same type (factor, numeric, ...) and have the same names as the data used to estimate the object. Only the bid variable(s) and the intercept are not needed. If \code{functionalForm} is \code{loglinear}, then the last variable of data.frame must be the income.

}
\value{A matrix with calculated willingness to pay and the confidence interval.
}
\references{
Haab, T.C. and McConnell, K.E. (2003), \emph{Valuing Environmental and Natural Resources. The Econometrics of non-market Valuation.} Cheltenham, UK: Edward Elgar 

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
# example linear
data("maizeWTP")

maizeDBCV<- 
  DBCV(yuYes1 + yuYes2~yuBid1 + yuBid2|age+gender+experimenter, 
       data=maizeWTP, functionalForm="linear")
maizeDBCV


## willigness to pay for all individually
WTP_maizeDBCV.all.individual  <- WTP_CV(object=maizeDBCV, 
                                        disp.pref   = "mean",
                                        disp.subj = "individual")

# willigness to pay for mean individual (this means that also 
# the mean of dummy variables is used
# which might not always make sense)
WTP_maizeDBCV.all.mean  <- WTP_CV(object=maizeDBCV, 
                                  disp.pref   = "mean",
                                  disp.subj = "mean")                   

## willigness to pay for three individuals

mySubjects<-data.frame("age"=c(47,21,29),
                       "gender"=relevel(factor(c("male","female","female")),
                                        ref=levels(maizeWTP$gender)[1]), 
                       "experimenter"=factor(c(
                        "experimenter3", "experimenter1", "experimenter6"))
                         )


# indivdual willigness to pay 
WTP_maizeDBCV.mySubjects.individual    <- WTP_CV(object=maizeDBCV, 
                                                 disp.pref   = "mean",
                                                 disp.subj = "individual",  
                                                 newdata=mySubjects
                                            )
# willingeness to pay of the mean of the three mySubjects
WTP_maizeDBCV.mySubjects.mean    <- WTP_CV(object=maizeDBCV, 
                                           disp.pref   = "mean",
                                           disp.subj = "mean",  
                                           newdata=mySubjects
                                    )

# plot something similar to a demand curve with confidence interval

plot(sort(WTP_maizeDBCV.all.individual[,1], decreasing=TRUE), type="l", 
     lwd=2, ylab="willingness to pay", xlab="Consumers", ylim=c(0,140))
lines(sort(WTP_maizeDBCV.all.individual[,2], decreasing=TRUE), type="l", 
      lwd=2, col="grey")
lines(sort(WTP_maizeDBCV.all.individual[,3], decreasing=TRUE), type="l", 
      lwd=2, col="grey")

# though, I am not sure you can use the Krinsky-Robb procdure to calculate
# individual confidence intervals. Most likely the confidence interval
# is bigger.
}

\keyword{willingness to pay}
\keyword{contingent valuation method}
