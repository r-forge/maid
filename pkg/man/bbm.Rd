\name{bbm}
\alias{bbm}
\title{Moment estimates of Beta-binomial parameters}
\description{Moment estimates are computed from raw data or from observed mean, standard deviation and sample size. The asymptotic covariance matrix of the estimates is computed following Moore (1986).}
\usage{
bbm(obsmean = NULL, obssd = NULL, 
    obsskewness = NULL, obskurtosis = NULL, 
    size, x = NULL, N = NULL, 
    calcmoments = TRUE, vcov = TRUE,
    logitrho = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsmean}{numeric. Observed mean. Can be omitted if \code{x} is not \code{NULL}.}
  \item{obssd}{numeric. Observed standard deviation. Can be omitted if \code{x} is not \code{NULL}.}
  \item{obsskewness}{numeric. Observed skewness. Can be omitted if \code{x} is not \code{NULL}. Only used for the estimation of the asymptotic covariance matrix if \code{calcmoments = FALSE}.}
  \item{obskurtosis}{numeric. Observed kurtosis. Can be omitted if \code{x} is not \code{NULL}. Only used for the estimation of the asymptotic covariance matrix if \code{calcmoments = FALSE}.}
  \item{size}{integer. Number of trials.}
  \item{x}{numeric vector of raw data with at least four observations. Can be omitted if \code{obsmean}, \code{obssd} and \code{N} are supplied.}
  \item{N}{Number of observations.}
  \item{calcmoments}{logical. Should the third and fourth moments be calculated from the moment estimates? Defaults to \code{TRUE}.}
  \item{vcov}{logical. Should the asymptotic covariance matrix be computed? Defaults to \code{TRUE}.}
  \item{logitrho}{logical. Should the moment estimates be computed for logit(rho)? Defaults to \code{FALSE}, that is, the moment estimates are computed for rho without any transformation.}
}
\details{The Beta-Binomial distribution is parametrized using the parameters \eqn{\pi}{pi} and \eqn{\rho}{rho} and assumes \eqn{n}{n} (= \code{size}) trials. The mean is \eqn{n\pi}{n*pi} and the variance \eqn{n\pi(1-\pi) + n(n-1)\rho \pi(1-\pi)}{n*pi*(1-pi) + n*(n-1)*rho*pi*(1-pi)} in this parametrization. Following Moore (1986), \eqn{logit(\pi)}{logit(pi)}, \eqn{\rho}{rho} (or \eqn{logit(\rho)}{logit(rho)} if \code{logitrho} is \code{TRUE}) and their asymptotic covariance matrix are estimated.

If the user specifies observed skewness and kurtosis or supplies raw data (and sets \code{calcoments = FALSE} in either case), observed skewness and kurtotis are used in the estimation of the asymptotic covariance matrix. Otherwise theoretical skewness and kurtosis are calculated from the moment estimates using standard formulae for the skewness and kurtosis of the Beta-binomial parameters. Calculating moments is recommended and hence the default.}
\value{Returns a list with components
\item{coef}{Estimates of parameters \eqn{logit(\pi)}{logit(pi)} and \eqn{\rho}{rho} (or \eqn{logit(\rho)}{logit(rho)} if \code{logitrho} is \code{TRUE}).}
\item{vcov}{Asymptotic covariance matrix of the estimates.}
}
\references{Moore, D. F. (1986). \dQuote{Asymptotic properties of moment estimators for overdispersed counts and proportions}. \emph{Biometrika}, \bold{73} (3), 583--588.}
\author{Philipp Doebler <philipp.doebler@googlemail.com>}

\seealso{\code{\link{betabinmoments}}}
\examples{
if (require("VGAM")) {
    # simulate 100 data with 12 trials each, mean 6 = 12*0.5
    # and some overdispersion
    x <- VGAM::rbetabinom(100, 12, 0.5, 0.1)
    # recover parameters with moment estimation
    pars <- bbm(x = x, size = 12)
    pars$coef[1] # should be close to 0 = logit(0.5)            
    1/(1+exp(-pars$coef[1])) # inverse logit
    pars$coef[2] # estimate of rho 
    
    # use observed moments for calculation
    # of asymptotic covariance matrix now:
    pars2 <- bbm(x=x, size = 12, calcmoments = FALSE)
    
    # compare calculated and observed moments
    round(pars$vcov, 5)
    round(pars2$vcov, 5)
    } else {
    # assume 100 data with 12 trials, mean 6 and observed standard deviation of 2.5
    bbm(obsmean = 6, obssd = 2.5, size = 12, N = 100)
  }

}
\concept{Beta-binomial distribution}
\concept{moment estimation}
\keyword{distribution}