\name{maid-Methods}
\alias{maid-Methods}
\alias{maid-Object}
\alias{print.maid}
\alias{summary.maid}
\alias{print.summary.maid}
\alias{qtest.maid}
\alias{logLik.maid}
\alias{plot.maid}
\alias{alphaKR21.maid}
\alias{alphaKR21}
\alias{print.alphaKR21.maid}

\title{Methods for \code{maid} objects}
\description{Standard and plot methods are documented here for \code{maid} objects that are the output of the \code{mabb} function}
\usage{
\method{print}{maid}(x, digits = 4, ...)
\method{summary}{maid}(object, ci.level = 0.95, ...)
\method{print}{summary.maid}(x, digits = 4, ...)
\method{logLik}{maid}(object, ...)
\method{qtest}{maid}(object, ...)
\method{plot}{maid}(x, pooled = TRUE, 
    main = "Meta-analysis of Beta-binomial parameters", 
    ellipsecol = "grey50", pool_pch = 19, pool_cex = 1, 
    predellipse_lty = 2, predellipse_col = "black", 
    backtransform = TRUE, xlim = "auto", ylim = "auto", ...)
\method{alphaKR21}{maid}(x, ci.level = 0.95, ...)
\method{print}{alphaKR21.maid}(x, digits = 4, ...)
}

\arguments{
  \item{x}{a \code{maid} object (these are the output of the \code{\link{mabb}} function).}
  \item{object}{a \code{maid} object (these are the output of the \code{\link{mabb}} function).}  
  \item{digits}{an integer specifying the number of digits to which printed results must be rounded.}
  \item{ci.level}{a numerical value between 0 and 1, specifying the confidence level for the computation of confidence intervals.}
  \item{pooled}{logical. Should the pooled value and its prediction ellipse be plotted?}
  \item{main}{caption of plot.}
  \item{ellipsecol}{color of confidence ellipses for study level estimates.}
  \item{pool_pch}{symbol for pooled value.}
  \item{pool_cex}{size of pooled value.}
  \item{predellipse_lty}{style of line of prediction ellipse.}
  \item{predellipse_col}{color of prediction ellipse.}
  \item{backtransform}{logical. Should logit(pi) and, if applicable, logit(rho) values be backtransformed prior to plotting? Defaults to TRUE}
  \item{xlim}{Either \code{"auto"} or a numeric vector of length 2. If the former, the function determines the plotting range for the x-axis.}             
  \item{ylim}{Either \code{"auto"} or a numeric vector of length 2. If the former, the function determines the plotting range for the y-axis.}             
  \item{...}{further arguments passed to or from other methods.}
}
\details{All methods are identical to or derived from their counterparts in the \code{\link{mvmeta-package}}.}
\value{The summary method function for \code{maid} objects produces a list of class \code{summary.maid}. The plot method function produces a plot and returns an invisible \code{NULL}. See the examples for plotting examples. 

The \code{alphaKR21} method for \code{maid} objects extracts the moment estimates of \eqn{logit(\pi)}{logit(pi)} and \eqn{\rho}{rho} (or \eqn{logit(\rho)}{logit(rho)} if \code{logitrho} is \code{TRUE}) on the study level and computes Wald confidence intervals. In addition, the function calculates the Kuder & Richardson (1937) formula 21 reliability estimates. The confidence interval for these estimates is based on the standard errors of the rho estimates. The \code{alphaKR21} produces an object of class \code{c("alphaKR21.maid", "data.frame")} which is essentially a \code{data.frame} with columns \code{logitpi}, \code{se.logitpi}, \code{ci.logitpi.lower}, \code{ci.logitpi.upper}, \code{rho}, \code{se.rho}, \code{ci.rho.lower}, \code{ci.rho.upper}, \code{alphaKR21}, \code{ci.alphaKR21.lower} and \code{ci.alphaKR21.upper}. If \code{logitrho} is \code{TRUE}, the estimates are given for \eqn{logit(\rho)}{logit(rho)} and the respective columns are labeled \code{logitrho}, \code{se.logitrho}, \code{ci.logitrho.lower}, \code{ci.logitrho.upper}.
}

\author{Philipp Doebler <philipp.doebler@googlemail.com>}
\seealso{\code{\link{mabb}}, \code{\link{summary.mvmeta}}}
\examples{
## Plotting a meta-analysis of individual differences

## Fit model
data(ICAR60) # load example data file
fit <- mabb(N = ICAR60$N, size = ICAR60$test_length, 
            obsmean = ICAR60$obsmean, obssd = ICAR60$obssd)
            
## print moment estimates and alpha KR21 estimates
alphaKR21(fit)

## change the default ci.level to .8 and print more digits
print(alphaKR21(fit, ci.level = .8), digits = 5)

## Produce default plot
plot(fit)
## Add a legend
legend("bottomright", 
       legend = c("study level estimate", "95\% confidence ellipse", 
       "pooled value", "95\% prediction ellipse"),
       lty = c(NA, 1, NA, 2), col = c("black", "grey50", "black"),  pch = c(1,NA, 19, NA))

## Variants of this plot (using all arguments and some ... methods)
plot(fit,
     pooled = TRUE, # TRUE is the default value
     main = "My meta-analysis",
     ellipsecol = "red", # color of the small ellipses for each study,
     pool_pch = 13, # symbol for pooled value
     pool_cex = 2, # size of symbol for pooled value
     predellipse_lty = 3, # style of line
     predellipse_col = "blue", # default is "black"
     backtransform = FALSE, # logit(pi) on x-axis
     xlim = "auto", # function computes xlim (default)
     ylim = c(.11, .15), # only plot values between .11 and .15 (+eps)
     # now some further argument passed on to plot.default() :
     pch = 2,  # symbol for study estimates
     cex = 0.5, # size of symbols for study estimates
     col = "green") # color of symbols for study estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{methods}
\keyword{htest}% __ONLY ONE__ keyword per line
