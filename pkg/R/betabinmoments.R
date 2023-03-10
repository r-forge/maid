betabinmoments <-
function(alpha=NULL, beta=NULL, pi = NULL, rho = NULL, size){
  n <- size
  stopifnot(size > 0, size == round(size))
  if(!((!is.null(alpha)&!is.null(beta))|(!is.null(pi)&!is.null(rho)))){
    stop("Either alpha and beta or pi and rho must be specified!")
  }
  if(!is.null(alpha)){
    stopifnot(alpha > 0, beta > 0)
  }
  if(is.null(alpha) & is.null(beta)){
    stopifnot(0 < pi, pi < 1)
    if(rho < 0){
      warning("rho is < 0! Check your data. Setting rho to 0.0001 for now.")
      rho <- 0.0001
    }
    if(rho > 1){
      warning("rho is > 1! Check your data. Setting rho to 0.9999 for now.")
      rho <- 0.9999
    }
    alpha_plus_beta <- 1/rho - 1
    alpha <- pi*alpha_plus_beta
    beta <- alpha_plus_beta - alpha
  }
  mean = n*alpha/(alpha + beta)
  variance <- n*alpha*beta*(alpha+beta+n)/{(alpha+beta)^2*(alpha+beta+1)}
  skewness <- sqrt((1+alpha+beta)/(n*alpha*beta*(n+alpha+beta))) * 
    (alpha + beta + 2*n)*(beta-alpha)/(alpha+beta+2)
  kurtosis <- {(alpha+beta)*(alpha + beta - 1 + 6*n) + 3*alpha*beta*(n-2) + 6*n^2 - 
                (3*alpha*beta*n*(6-n))/(alpha + beta) - (18*alpha*beta*n^2)/((alpha+beta)^2)} * 
{(alpha+beta)^2*(1+alpha+beta)}/{n*alpha*beta*(alpha+beta+2)*(alpha+beta+3)*(alpha+beta+n)}
  mu3 <- skewness*variance^(3/2)
  mu4 <- kurtosis*variance^2
  return(list(mean = mean, variance = variance, skewness = skewness,
              kurtosis = kurtosis, mu3 = mu3, mu4 = mu4))
}
