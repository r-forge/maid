\name{ICAR60}
\alias{ICAR60}
\docType{data}
\title{ICAR 60 Data}
\description{Aggregated data from some subsamples of the ICAR 60 intelligence test.}
\usage{data("ICAR60")}
\format{
  A data frame with 7 observations on the following 7 variables.
  \describe{
    \item{\code{subsample}}{a factor with levels \code{CollegeDegree}, \code{CurrentInUniv}, \code{GradOrProDegree}, \code{HSgrad}, \code{InGradOrProSchool}, \code{less12yrs} and \code{SomeCollege}.}
    \item{\code{N}}{a numeric vector. The size of the subsamples.}
    \item{\code{test_length}}{a numeric vector. The number of dichotomous items answered by each testee (16 in all subsamples).}
    \item{\code{obsmean}}{a numeric vector. Observed mean.}
    \item{\code{obssd}}{a numeric vector. Observed standard deviation.}
    \item{\code{alpha}}{a numeric vector. The values of Cronbach's alpha calculated from the complete data before aggregation.}
    \item{\code{gradSample}}{a factor with levels \code{yes} and \code{no}. A \code{yes} indicates that the subsample consists of graduates.}
  }
}
\details{This dataset has been compiled from raw test data by exluding persons who answered less than 16 items and then calculating the means, standard deviations and Cronbach alpha values for each subset.}
\source{The original dataset has been made available by David M. Condon and can be obtained from Harvard Dataverse (\doi{10.7910/DVN/AD9RVY}). Also see \url{https://icar-project.com/}.
}
\references{Condon, D., & Revelle, W. (2014). \dQuote{The international cognitive ability resource: development and initial validation of a public-domain measure}. \emph{Intelligence}, \bold{43}, 52--64.
}
\keyword{datasets}
