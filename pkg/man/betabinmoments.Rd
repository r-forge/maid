\name{betabinmoments}
\alias{betabinmoments}
\title{The first four central moments of the Beta-binomial distribution}
\description{Calculates the first four central moments as well as skewness and kurtosis from two different parametrizations}
\usage{
betabinmoments(alpha = NULL, beta = NULL, 
               pi = NULL, rho = NULL, 
               size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{numeric > 0. Parameter (see Details).}
  \item{beta}{numeric > 0. Parameter (see Details).}
  \item{pi}{numeric. Parameter between 0 and 1 (see Details and \code{\link{bbm}}).}
  \item{rho}{numeric. Parameter between 0 and 1 (see Details and \code{\link{bbm}}).}
  \item{size}{integer. Parameter (number of trials).}
}
\details{The user needs to supply either \code{alpha} and \code{beta} or \code{pi} and \code{rho}. In either case \eqn{n}{n} (= \code{size}) is the number of trials.  In the parametrization using \code{alpha} and \code{beta}, the mean and variance are given by \eqn{n\frac{\alpha}{\alpha + \beta}}{n*alpha/(alpha+beta)} and \eqn{n\alpha\beta(\alpha+\beta+n)/((\alpha+\beta)^2(\alpha+\beta+1))}{n*alpha*beta*(alpha+beta+n)/((alpha+beta)^2*(alpha+beta+1))} respectively. If \code{rho} and \code{pi} are supplied, the mean and variance of the Beta-binomial distributions are given by the expressions stated in the documentation of \code{\link{bbm}}.}
\value{A list with components
\item{mean}{mean}
\item{variance}{variance}
\item{skewness}{skewness}
\item{kurtosis}{kurtosis}
\item{mu3}{third central moment}
\item{mu4}{fourth central moment}
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}
\note{Both parametrizations are related by the following relations:
\deqn{\pi = \frac{\alpha}{\alpha + \beta}}{pi=alpha/(alpha+beta)}
and
\deqn{\rho = \frac{1}{\alpha+\beta+1}.}{rho=1/(alpha+beta+1).}
}
\seealso{\code{\link{bbm}}}
\examples{
# Example 1: alpha and beta as parameters
alpha <- 4.5
beta <- 4.5
betabinmoments(alpha, beta, size = 12)

# Example 2: calculte pi and rho:
pi <- alpha/(alpha + beta) # 0.5
rho <- 1/(alpha + beta + 1) # 0.1
betabinmoments(pi = pi, rho = rho, size = 12)
# the same moments result
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\concept{Beta-binomial distribution}
\concept{central moments}% __ONLY ONE__ keyword per line
\concept{skewness}
\concept{kurtosis}
