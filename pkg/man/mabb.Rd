\name{mabb}
\alias{mabb}
\alias{mabb2}
\title{Fit bivariate meta-analysis model to Beta-binomial parameters}
\description{The Beta-binomial parameters and their asymptotic covariances are calculated from the observed means and standard deviations. Then a bivariate meta-analysis is perfomed. Several methods are available. See \code{\link{mvmeta.fit}}. For convenience, the \code{mabb2} function is supplied which simplifies calling \code{mabb} if the columns of the supplied \code{data.frame} follow naming conventions.}
\usage{
mabb(N, size, obsmean, obssd, 
     X = model.matrix(~1, data.frame(N)), 
     method = "reml", bscov = "unstr",
     control = list(), logitrho = FALSE)
mabb2(data, N = "N", size = "size", 
      obsmean = "obsmean", obssd = "obssd", ...)     
}
\arguments{
  \item{N}{numeric. Vector of sample sizes. For \code{mabb2}: character. Name of column in \code{data} with sample sizes.}
  \item{size}{numeric. Vector of number of trials. For \code{mabb2}: character. Name of column in \code{data} with number of trials.}
  \item{obsmean}{numeric. Vector of observed means .For \code{mabb2}: character. Name of column in \code{data} with observed means.}
  \item{obssd}{numeric. Vector of observed standard deviations. For \code{mabb2}: character. Name of column in \code{data} with observed standard deviations.}
  \item{X}{model matrix. Assumes that an intercept is included in the model. See examples below for a convenient way to specify X for moderator analysis or meta-regression}
  \item{method}{Fitting method. Defaults to "reml". See \code{\link{mvmeta}} for other options.}
  \item{bscov}{a string defining the between-study (co)variance structure in likelihood based models. See \code{\link{mvmeta.fit}}.}
  \item{control}{list of parameters for controlling the fitting process. These are passed to \code{\link{mvmeta.fit}} to replace otherwise selected default values.}
  \item{logitrho}{logical. Should the parameter rho be logit-transformed? Defaults to \code{FALSE}, that is, the parameter rho is not transformed.}
  \item{data}{A \code{data.frame} with appropriate column names.}
  \item{...}{Arguments passed on from \code{mabb2} to \code{mabb}.}
}
\details{The \code{mabb} function calls \code{\link{bbm}} to calculate moment estimates from observed means and standard deviations. The workhorse function \code{\link{mvmeta.fit}} is then used to perform bivariate meta-analysis.

Meta-regression is supported via user supplied design matrices.}
\value{An object of class \code{c("maid")} for which additional methods are available; see \code{\link{maid-Methods}}.}
\references{Doebler, P. and Doebler, A. (submitted).  Meta-analysis of individual differences based on sample means and standard deviations: An approach using the parameters of the Beta-binomial distribution.

Gasparrini A., Armstrong, B., Kenward M. G. (2012). Multivariate  meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29). 3821--3839.
}
\author{Philipp Doebler <philipp.doebler@googlemail.com>}
\seealso{\code{\link{maid-Methods}},\code{\link{mvmeta}}, \code{\link{mvmeta.fit}}}
\examples{
# load example data set 
data(ICAR60)

#######################################################
# First example: intercept only model (no moderators) #
#######################################################
with(ICAR60, {
     fit1 <- mabb(N = N, size = test_length, 
                  obsmean = obsmean, obssd = obssd)
     summary(fit1)
}) # with() is used for convenience here (attach() is
   # an alternative)
   
# The same result can also be produced with
fit1 <- mabb(N = ICAR60$N, size = ICAR60$test_length, 
            obsmean = ICAR60$obsmean, obssd = ICAR60$obssd)
summary(fit1)

# the mabb2 function is a bit more convenient if
# the columns of the data frame match the
# argument names at least partially, and is yet
# another way to produce the above result:
fit1_2 <- mabb2(ICAR60, size = "test_length")
summary(fit1_2)

# It is also possible to use a fixed effects model (not 
# recommended) or moment estimation by using the method 
# argument (consult mvmeta and mvmeta.fit for more 
# information and additional estimators).

# Individual moment estimates and alpha KR21 estimates
# are extracted with the alphaKR21 method:
alphaKR21(fit1)
   
###########################################
# Second example: gradSample as moderator #
###########################################

# First step: create model matrix
X <- model.matrix( ~ gradSample, ICAR60)
  # note that the first argument of model.matrix is
  # a formula with no left hand side

# Second step: fit model
fit2 <- with(ICAR60, {
     mabb(N = N, size = test_length, 
                  obsmean = obsmean, obssd = obssd,
                  X = X)
})
summary(fit2)

# with mabb2, this can be achieved like this:
fit2_2 <- mabb2(ICAR60, size = "test_length", X = X)

# many standard methods are available, e.g.
AIC(fit2) # Akaike information criterion
coef(fit2) # fixed effect coefficients
# See maid-Methods for more information on methods.

###########################################
# Third example: sample size as moderator #
###########################################

# This is a toy example for a (more or less)
# continuous moderator.

# First step: create model matrix
X <- model.matrix( ~ N, ICAR60)
  # note that the first argument of model.matrix is
  # a formula with no left hand side

# Second step: fit model
with(ICAR60, {
     fit3 <- mabb(N = N, size = test_length, 
                  obsmean = obsmean, obssd = obssd,
                  X = X)
     summary(fit3)
     # no influence of sample size
})

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{meta-analysis}
\concept{Beta-binomial distribution}
\keyword{models}% __ONLY ONE__ keyword per line
