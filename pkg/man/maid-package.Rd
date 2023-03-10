\name{maid-package}
\alias{maid-package}
\alias{maid}
\docType{package}
\title{
\packageTitle{maid}
}
\description{
\packageDescription{maid}
}
\details{The reliability of a test or questionnaire is an important measure of its psychometric quality. However, reliability is sample specific, and thus care needs to be taken when extrapolating from reliability estimates to a new application. Researchers needing to pick an instrument can meta-analyze existing reports of reliability. This technique is also known as reliability generalization (Vacha-Haase, 1998). 

However, reliability is also underreported and hence reliability generalization is not always possible. For the case of dichotomous items a meta-analysis of individual differences (MAID) is an alternative to reliability generalization that does not assume that reliabilities are reported. Assuming that test scores follow a Beta-binomial distribution, only the observed means and standard deviations suffice for MAID. 

In a first step the parameters of the Beta-binomial distribution are estimated from the observed means and standard deviations using moment estimation. The Beta-binomial distribution can be parametrized in a fashion, so that one parameter describes item difficulty (\eqn{\pi}{pi}) and the other overdispersion (\eqn{\rho}{rho}). The parameter \eqn{\rho}{rho} can also be interpreted as a measure of the amount of individual differences the test can detect, hence the name meta-analysis of individual differences. Another interpretation is that \eqn{\rho}{rho} is the correlation of two answers from the same person. Since \eqn{\rho}{rho} is also closely linked to KR21 reliability a meta-analysis of \eqn{\rho}{rho} produces estimates of reliability, too. Moment estimation in the current package is performed by the function \code{\link{bbm}}.

Here a bivariate approach that integrates \eqn{\pi}{pi} and \eqn{\rho}{rho} simultaneously is implemented. The central function is \code{\link{mabb}} which calculates the moment estimates and performs bivariate meta-analysis. The KR21 reliability estimates can then be extracted from the output of \code{\link{mabb}} with the \code{\link{alphaKR21}} method. Many further methods are available, including plots. More details are provided in the documentation of the individual functions and in Doebler & Doebler (submitted).}
\author{
\packageAuthor{maid}

Maintainer: \packageMaintainer{maid}
}
\references{
Doebler, P. and Doebler, A. (submitted).  Meta-analysis of individual differences based on sample means and standard deviations: An approach using the parameters of the Beta-binomial distribution.

Gasparrini A., Armstrong, B., Kenward M. G. (2012). Multivariate  meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29). 3821--3839.

Moore, D. F. (1986). Asymptotic properties of moment estimators for overdispersed counts and proportions. \emph{Biometrika}, \bold{73} (3), 583--588.

Vacha-Haase, T. (1998). Reliability generalization: exploring variance in
measurement error affecting score reliability across studies. \emph{Educational
and Psychological Measurement}, \bold{58} (1), 6--20.}
\keyword{ package }
\seealso{\code{\link{mabb}}}
\examples{
#################
# Basic example #
#################
data(ICAR60)
with(ICAR60, {
     mabb(N = N, size = test_length, 
          obsmean = obsmean, obssd = obssd)
})
}
