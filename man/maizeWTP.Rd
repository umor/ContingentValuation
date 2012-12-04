\name{maizeWTP}
\alias{maizeWTP}
\docType{data}
\title{
Dataset 'maizeWTP'
}
\description{
151 Observations from double bounded contingent survey about WTP for 2kg packages of maize meal in 2009 in Machakos, Kenya.
}
\usage{data(maizeWTP)}

\format{
  This data frame contains the following columns:
   \describe{
      \item{yuBid1}{A numeric variable giving the first bid for a 2kg package of yellow unfortified maize offered to the respondent (in Kenyan Shillings).}
      \item{yuBid2}{A numeric variable giving the second bid for a 2kg package of yellow unfortified maize offered to the respondent (in Kenyan Shillings).}
      \item{yuYes1}{A numeric variable which is 1 of the respondent said 'yes' to the first bid for a 2k package of yellow unfortified maize and 0 if the reply was 'no'.}
      \item{yuYes2}{A numeric variable which is 1 of the respondent said 'yes' to the second bid for a 2k package of yellow unfortified maize and 0 if the reply was 'no'.}
      \item{gender}{A factor with the levels 'femal' and 'male' giving the gender of the respondent.}
      \item{age}{A numeric variable giving the age of the respondent}
      \item{experimenter}{A factor with six levels depending which experimenter conducted the survey}
     
}
}



\source{
   Previously unblished data which have been collected along with data used in Morawetz, DeGroote and Kimenju (2011).
}

\references{
Morawetz, U.B., DeGroote, H. and Kimenju, S.C. (2011) Improving the Use of Experimental Auctions in Africa: Theory and Evidence, \emph{Journal of Agriculturl and Resource Economics}, 36(2):263-279.

Kimenju, S.C, Morawetz, U.B. and DeGroote, H.(2005)
Comparing Contingent Valuation Method, Choice Experiments and
Experimental Auctions in soliciting Consumer preference for maize
in Western Kenya: Preliminary results.  \emph{African Econometric Society: 10th Annual Conference on Econometric Modeling in Africa, 6-8 July 2005, Nairobi, Kenya}

}

\examples{
   data(maizeWTP)
   summary(maizeWTP)
}

\keyword{datasets}
