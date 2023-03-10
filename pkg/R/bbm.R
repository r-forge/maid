bbm <-
function(obsmean = NULL, obssd = NULL, obsskewness = NULL, obskurtosis = NULL,
                size, x = NULL, N = NULL, calcmoments = TRUE,
                vcov = TRUE, logitrho = FALSE){
  stopifnot(!is.null(x) | (!is.null(obsmean) & !is.null(obssd) & !is.null(N)))
  stopifnot(all.equal(round(size), size))
  if(!is.null(N))
    stopifnot(all.equal(round(N), N))
  if(!is.null(x)){
    stopifnot(length(x) > 3, all.equal(round(x), x))
    obssd <- sd(x)
    obsmean <- mean(x)
    N <- length(x)
  }
  pi <- obsmean/size
  
  ## some values needed for later calculations
  s2 <- obssd^2
  Delta <- pi*(1-pi)
  nu1 <- Delta/size
  nu2 <- (size-1)*nu1
  s2_z <- s2/size^2
  
  # the following estimate of rho is identical to the classic moment estimator
  # formula follows from Moore (1986), see latex'ed notes of Philipp Doebler
  rho <- (s2_z*(N-1)/N - nu1)/nu2  
  
  if(!vcov)
    return(list(coef = c(log(pi/(1-pi)), rho)))
  
  if(is.null(x) | calcmoments){
    if(calcmoments | (is.null(obsskewness) | is.null(obskurtosis)))
      bbmoments <- betabinmoments(pi = pi, rho = rho, size = size)
    
    ## use information from observed skewness to calculate third moment if available
    if(calcmoments | is.null(obsskewness)){
      mu3 <- bbmoments$mu3/size^3
    }else{
      mu3 <- (obsskewness*obssd^3)/size^3
    }
    ## use information from observed kurtosis to calculate third moment if available
    if(calcmoments | is.null(obskurtosis)){
      mu4 <- bbmoments$mu4/size^4
    }else{
      mu4 <- (obskurtosis*obssd^4)/size^4
    }
  }else{
    mu3 <- moments::moment(x, 3, central = TRUE)/size^3
    mu4 <- moments::moment(x, 4, central = TRUE)/size^4
  }
  
  ## asymptotic variance covariance matrix
  nu <- nu1 + rho*nu2
  Gamma <- matrix(c(-Delta^2/nu, -(Delta*(1+(size-1)*rho)*(1-2*pi))/(size*nu),
                    0, -nu2/nu), ncol = 2, nrow = 2)
  Lambda <- matrix(c(Delta^2/nu, mu3*Delta/nu^2, 
                     mu3*Delta/(nu^2), mu4/nu^2 - 1),
                   ncol = 2, nrow = 2)
  invGamma <- solve(Gamma)
  Psi <- invGamma %*% Lambda %*% t(invGamma)
  if(logitrho) {
    nablag <- diag(c(1, 1/(rho*(1-rho))))
    Psi <- t(nablag) %*% Psi %*% nablag
    coefs <-  c(log(pi/(1-pi)), log(rho/(1-rho)))
    colnames(Psi) <- rownames(Psi) <- names(coefs) <- c("logit(pi)", "logit(rho)")
  } else {
    coefs <-  c(log(pi/(1-pi)), rho)
    colnames(Psi) <- rownames(Psi) <- names(coefs) <- c("logit(pi)", "rho")
  }
  return(list(coef = coefs,
         vcov = Psi/N))
}
