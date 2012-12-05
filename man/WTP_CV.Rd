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
        disp.pref   =  c("mean", "median"),  
        disp.subj = c("individual", "mean"),
        newdata = NULL
                   )
}
\arguments{
  
  \item{object}{An object created by the estimation of SBCV or DBCV}
  \item{disp.pref}{Either \code{mean} or \code{median}.}   
  \item{disp.subj}{Either \code{individual} or \code{mean}.} 
  \item{newdata}{Either \code{NULL} or a data.frame}
  
  
  \item{\dots}{ currently not used  }
}
\details{
The calculation of the willigness to pay is for probit \code{SBCV} and \code{DBCV}. 


Dispersion of the preferences, \code{disp.pref}, is either \code{mean} or \code{median} depending on where the variance of the preverences is to be evaluated. In case of  \code{functionalForm}  is \code{linear} the resulting willingness to pay is the same for \code{mean} and \code{median}.

\code{disp.subj} is either \code{individual} or \code{mean}. If \code{individual} and \code{newdata} is  \code{NULL}, then the willigness to pay is calculated for each person in the data used to estimate the object. If \code{individual} and in \code{newdata} has a \code{data.frame}, then the willigness to pay is calculated for these data. If \code{mean} the mean of the data provided in  \code{newdata} (or, if \code{NULL}, all the the data used  to estimate the object) is calculated and the willigness to pay is calculated for the mean subject. 

If \code{newdata} is \code{NULL}, all subjects are used, if a \code{data.frame} is given, these data are used. The \code{data.frame} must have \code{(length(object$coefficients)-1)}  columns. The dropped column is the bid variable.

}
\value{A vector with calculated willingness to pay.
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
# example linear

data("maizeWTP")
require(Formula)

maizeDBCV<- 
  DBCV(yuYes1 + yuYes2~yuBid1 + yuBid2|age+gender+experimenter, 
       data=maizeWTP, functionalForm="linear")
maizeDBCV


## willigness to pay for all individually
WTP_maizeDBCV.all.individual  <- WTP_CV(object=maizeDBCV, 
                                        disp.pref   = "mean",
                                        disp.subj = "individual"  
)
# willigness to pay for mean individual (this means that also 
# the mean of dummy variables is used
# which might not always make sense)
WTP_maizeDBCV.all.mean  <- WTP_CV(object=maizeDBCV, 
                                  disp.pref   = "mean",
                                  disp.subj = "mean"  
)                   

## willigness to pay for three individuals

mySubjects<-data.frame("(Intercept)"=c(1,1,1), 
                       "age"=c(47,21,29),
                       "genderfemale"=c(0,1,0), 
                       "experimenterexperimenter2"=c(1,0,0), 
                       "experimenterexperimenter3"=c(0,0,0), 
                       "experimenterexperimenter4"=c(0,0,0), 
                       "experimenterexperimenter5"=c(0,0,0),
                       "experimenterexperimenter6"=c(0,0,0))

# colnames are not necessary


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

# plot something similar to a demand curve

plot(sort(WTP_maizeDBCV.all.individual, decreasing=TRUE), type="l", 
     lwd=2, ylab="willingness to pay", xlab="Consumers")

}

\keyword{willingness to pay}
\keyword{contingent valuation method}
