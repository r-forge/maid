mabb <-
function(N, size, obsmean, obssd, 
         X = model.matrix(~1, data.frame(N)),
         method = "reml", bscov = "unstr",
         control = list(), logitrho = FALSE){
  call <- match.call()  
  k <- length(N)
  # check input
  stopifnot(k > 2,
            length(size) == k,
            length(obsmean) == k,
            length(obssd) == k,
            method %in% c("fixed", "ml", "reml", "mm", "vc"),
            all(size > 1),
            all(N > 1),
            all(obsmean > 0),
            all(obssd > 0),
            all(abs(size - round(size)) < 0.00001),
            all(abs(N - round(N)) < 0.00001),
            is.logical(logitrho))
  S <- array(NA, dim = c(2,2,k))
  y <- matrix(NA, ncol = 2, nrow = k)
  if(logitrho) {
    colnames(y) <- c("logitpi", "logitrho")
  } else {
    colnames(y) <- c("logitpi", "rho")
  }
  for(i in 1:k){
    bbmests <- bbm(obsmean = obsmean[i], obssd[i],
                   size = size[i], N = N[i], logitrho=logitrho)
    S[,,i] <- bbmests$vcov
    y[i,] <- bbmests$coef
  }
  S <- t(apply(S, 3, vechMat))
  fit <- mvmeta.fit(X = X, y = y, S = S,
                    method = method, bscov = bscov,
                    control = control)
  fit$N <- N
  fit$size <- size
  fit$X <- X
  fit$S <- S
  fit$call <- call
  fit$y <- y
  fit$formula <- NA
  fit$terms <- NA
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- NA
  fit$logitrho <- logitrho
  class(fit) <- c("maid")
  fit
}

mabb2 <- function(data, N = "N", size = "size", 
                  obsmean = "obsmean", obssd = "obssd",
                  ...){
  return(mabb(N = data[, N], size = data[, size], 
              obsmean = data[,obsmean], obssd = data[, obssd],
              ...))
}
