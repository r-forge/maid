## methods for maid objects

## can directly reuse some mvmeta methods
print.maid <- print.mvmeta
logLik.maid <- logLik.mvmeta
vcov.maid <- vcov.mvmeta


## new plot method
plot.maid <- function(x, pooled = TRUE,
                      main = "Meta-analysis of Beta-binomial parameters",
                      ellipsecol = "grey50",
                      pool_pch = 19, pool_cex = 1,
                      predellipse_lty = 2, predellipse_col = "black",
                      backtransform = TRUE,
                      xlim = "auto", ylim = "auto",
                      ...){
  if(nrow(coef(x)) > 1L){
    stop("Plotting only available when no covariates are used.")
  }
  
  stopifnot(xlim == "auto" | (is.numeric(xlim) & length(xlim) == 2))
  stopifnot(ylim == "auto" | (is.numeric(ylim) & length(ylim) == 2))
  stopifnot(is.logical(backtransform))
  
  xx <- x$y[,1]
  yy <- x$y[,2]
  N <- length(xx)
  S <- x$S
  logitrho <- x$logitrho
  
  if(pooled){
    coefs <- coef(x)
    pool_ellipse <- ellipse(x$Psi + x$vcov, centre = coefs)
  }
  ellipses <- array(NA, dim = c(length(xx), 100, 2))
  for(i in 1:N){
    ellipses[i,,] <- ellipse(matrix(c(S[i,1], S[i,2], S[i,2], S[i,3]), ncol = 2),
                          centre = c(xx[i], yy[i]))
  }
  maxs <- apply(ellipses, 3, max)
  mins <- apply(ellipses, 3, min)  
  if(is.character(xlim)){
    if(xlim == "auto"){
    xlim <- c(mins[1], maxs[1])
    if(pooled){
      xlim <- c(min(pool_ellipse[,1], xlim[1]), max(pool_ellipse[,1], xlim[2]))
    }
    xlim <- if(backtransform) plogis(xlim) else xlim 
   }
  }
  if(is.character(ylim)){
    if(ylim == "auto"){
    ylim <- c(mins[2], maxs[2])
    if(pooled){
      ylim <- c(min(pool_ellipse[,2], ylim[1]), max(pool_ellipse[,2], ylim[2]))
    } 
    ylim <- if(backtransform & logitrho) plogis(ylim) else ylim
   }
  }
  
  plot(if(backtransform) plogis(xx) else xx, 
       if(backtransform & logitrho) plogis(yy) else yy,
       ylim = ylim, xlim = xlim,
       main = main,
       xlab = if(backtransform) expression(pi) else expression(logit(pi)), 
       ylab = if(!backtransform & logitrho) expression(logit(rho)) else expression(rho), ...)

  for(i in 1:N){
    lines(cbind(if(backtransform) plogis(ellipses[i,,1]) else ellipses[i,,1], 
                if(backtransform & logitrho) plogis(ellipses[i,,2]) else ellipses[i,,2]), 
                col = ellipsecol)
  }
  
  if(pooled){
    points(if(backtransform) plogis(coefs[1,1]) else coefs[1,1], 
           if(backtransform & logitrho) plogis(coefs[1,2]) else coefs[1,2], 
           pch = pool_pch, cex = pool_cex)
    lines(cbind(if(backtransform) plogis(pool_ellipse[,1]) else pool_ellipse[,1], 
                if(backtransform & logitrho) plogis(pool_ellipse[,2]) else pool_ellipse[,2]), 
          lty = predellipse_lty, col = predellipse_col)
  }
  return(invisible(NULL))
}


## minor modifications necessary compared to qtest.mvmeta
qtest.maid <- 
function (object, ...) 
{
  int <- TRUE
  y <- object$y
  X <- object$X
  S <- object$S
  nay <- is.na(y)
  dim <- object$dim
  Xlist <- lapply(seq(dim$m), function(i) diag(1, dim$k)[!nay[i, ], , drop = FALSE] %x% X[i, , drop = FALSE])
  ylist <- lapply(seq(dim$m), function(i) y[i, ][!nay[i, ]])
  if (dim(S)[2] == ncol(y)) 
    S <- inputcov(sqrt(S), object$control$Scor)
  Slist <- lapply(seq(dim$m), function(i) xpndMat(S[i, ])[!nay[i, 
                                                               ], !nay[i, ], drop = FALSE])
  nalist <- lapply(seq(dim$m), function(i) nay[i, ])
  Psi <- diag(0, dim$k)
  gls <- glsfit(Xlist, ylist, Slist, nalist, Psi, onlycoef = FALSE)
  Q <- drop(crossprod(gls$invtUy - gls$invtUX %*% gls$coef))
  df <- with(object$df, nall - fixed)
  if (dim$k > 1L) {
    Q <- c(Q, colSums(do.call("rbind", mapply(function(y, 
                                                       S, X, na) {
      comp <- rep(0, dim$k)
      comp[!na] <- as.vector((y - X %*% gls$coef)^2/diag(S))
      return(comp)
    }, ylist, Slist, Xlist, nalist, SIMPLIFY = FALSE))))
    df <- c(df, colSums(!nay, na.rm = TRUE) - dim$p)
  }
  pvalue <- sapply(seq(length(Q)), function(i) 1 - pchisq(Q[i], 
                                                          df[i]))
  names(Q) <- names(df) <- names(pvalue) <- if (dim$k > 1L) 
    c(".all", object$lab$k)
  else object$lab$k
  qstat <- list(Q = Q, df = df, pvalue = pvalue, residual = object$dim$p - 
                  int > 0L, k = dim$k)
  class(qstat) <- "qtest.mvmeta"
  qstat
}


summary.maid <- function(object, ci.level = 0.95, ...){
  s <- summary.mvmeta(object, ci.level = ci.level, ...)
  class(s) <- "summary.maid"
  s
}

print.summary.maid <- 
function (x, digits = 4, ...) 
{
  methodname <- c("reml", "ml", "fixed", "mm", "vc")
  methodlabel <- c("REML", "ML", "Fixed", "Method of moments", 
                   "Variance components")
  bscovname <- c("unstr", "diag", "id", "cs", "hcs", "ar1", 
                 "prop", "cor", "fixed")
  bscovlabel <- c("General positive-definite", "Diagonal", 
                  "Multiple of identity", "Compound symmetry", "Heterogeneous compound symmetry", 
                  "Autoregressive of first order", "Proportional to fixed matrix", 
                  "Fixed correlation", "Fixed")
  cat("Call:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat(if (x$dim$k == 1L) 
    "Uni"
    else "Multi", "variate ", ifelse(x$method == "fixed", "fixed", 
                                     "random"), "-effects meta-", ifelse(x$dim$p - 1 > 0, 
                                                                         "regression", "analysis"), "\n", sep = "")
  cat("Dimension: ", x$dim$k, "\n", sep = "")
  if (x$method != "fixed") {
    cat("Estimation method: ", methodlabel[which(x$method == 
                                                   methodname)], "\n", sep = "")
  }
  cat("\n")
  cat("Fixed-effects coefficients", "\n", sep = "")
  signif <- symnum(x$coefficients[, "Pr(>|z|)"], corr = FALSE, 
                   na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 
                                             1), symbols = c("***", "**", "*", ".", " "))
  tabletot <- formatC(x$coefficients, digits = digits, format = "f")
  tabletot <- cbind(tabletot, signif)
  colnames(tabletot)[7] <- ""
  if (x$dim$p == 1L) {
    print(tabletot, quote = FALSE, right = TRUE, print.gap = 2)
  }
  if (x$dim$p > 1L) {
    p <- x$dim$p
    for (i in seq(x$dim$k)) {
      ind <- seq((i - 1) * p + 1, (i - 1) * p + p)
      table <- tabletot[ind, , drop = FALSE]
      rownames(table) <- x$lab$p
      if (x$dim$k > 1) 
        cat(x$lab$k[i], ":", "\n")
      print(table, quote = FALSE, right = TRUE, print.gap = 2)
    }
  }
  cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
  if (!x$method == "fixed") {
    cat("Between-study random-effects (co)variance components", 
        "\n", sep = "")
    if (x$dim$k > 1L) 
      cat("\t", "Structure: ", bscovlabel[which(x$bscov == 
                                                  bscovname)], "\n", sep = "")
    sd <- formatC(cbind(`Std. Dev` = sqrt(diag(x$Psi))), 
                  digits = digits, format = "f")
    if (x$dim$k == 1L) 
      rownames(sd) <- ""
    if (x$dim$k > 1L) {
      corRan <- x$corRandom
      corRan[upper.tri(x$corRan, diag = TRUE)] <- NA
      dimnames(corRan) <- NULL
      corRan <- format(corRan[-1, -ncol(corRan), drop = FALSE], 
                       digits = digits, format = "f")
      corRan <- rbind(x$lab$k[-x$dim$k], corRan)
      colnames(corRan) <- c("Corr", rep("", x$dim$k - 2))
      corRan[grep("NA", corRan)] <- ""
    }
    else corRan <- NULL
    print(cbind(sd, corRan), quote = FALSE, right = TRUE, 
          na.print = "", print.gap = 2)
    if (!is.null(x$negeigen) && x$negeigen > 0) {
      cat("(Note: Truncated estimate - ", x$negeigen, " negative eigenvalues set to 0)", 
          "\n", sep = "")
    }
    cat("\n")
  }
  Q <- formatC(x$qstat$Q, digits = digits, format = "f")
  pvalue <- formatC(x$qstat$pvalue, digits = digits, format = "f")
  i2 <- formatC(pmax((x$qstat$Q - x$qstat$df)/x$qstat$Q * 100, 
                     1), digits = 1, format = "f")
  cat(if (x$qstat$k == 1) 
    "Uni"
    else "Multi", "variate ", "Cochran Q-test for ", if (x$qstat$residual) 
      "residual ", "heterogeneity:", "\n", sep = "")
  cat("Q = ", Q[1], " (df = ", x$qstat$df[1], "), p-value = ", 
      pvalue[1], "\n", sep = "")
  cat("I-square statistic = ", i2[1], "%", "\n\n", sep = "")
  cat(x$dim$m, " studies, ", x$df$nall, " observations, ", 
      x$df$fixed, " fixed and ", x$df$random, " random-effects parameters", 
      "\n", sep = "")
  if (na <- length(x$na.action)) 
    cat("(", na, " stud", ifelse(na > 1L, "ies", "y"), " removed due to missingness", 
        ")\n", sep = "")
  if (!x$method %in% c("mm", "vc")) {
    table <- c(x$logLik, x$AIC, x$BIC)
    names(table) <- c("logLik", "AIC", "BIC")
    table <- formatC(table, digits = digits, format = "f")
    print(table, quote = FALSE, right = TRUE, print.gap = 2)
  }
  cat("\n")
}

# identical to mvmeta:::glsfit (but since not exported from
# mvmeta, it is included here)
glsfit <- 
function (Xlist, ylist, Slist, nalist, Psi, onlycoef = TRUE) 
{
  Sigmalist <- mapply(function(S, na) S + Psi[!na, !na, drop = FALSE], 
                      Slist, nalist, SIMPLIFY = FALSE)
  Ulist <- lapply(Sigmalist, chol)
  invUlist <- lapply(Ulist, function(U) backsolve(U, diag(ncol(U))))
  invtUXlist <- mapply(function(invU, X) crossprod(invU, X), 
                       invUlist, Xlist, SIMPLIFY = FALSE)
  invtUylist <- mapply(function(invU, y) crossprod(invU, y), 
                       invUlist, ylist, SIMPLIFY = FALSE)
  invtUX <- do.call("rbind", invtUXlist)
  invtUy <- do.call("rbind", invtUylist)
  coef <- as.numeric(qr.solve(invtUX, invtUy))
  if (onlycoef) 
    return(coef)
  list(coef = coef, Sigmalist = Sigmalist, Ulist = Ulist, invUlist = invUlist, 
       invtUXlist = invtUXlist, invtUX = invtUX, invtUy = invtUy)
}


confint.maid <- function (object, parm, level = 0.95, ...) 
{
  cf <- as.numeric(coef(object))
  if(missing(parm)){parm <- colnames(vcov(object))}
  
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3),
               "%")
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                             pct))
  ses <- sqrt(diag(vcov(object)))
  ci[] <- cf + ses %o% fac
  ci
}


alphaKR21 <- function(x, ...) UseMethod("alphaKR21")
alphaKR21.default <- function(x, ...) x

alphaKR21.maid <- function(x, ci.level = .95, ...){
  q <- qnorm(1-(1-ci.level)/2)
  out <- data.frame(logitpi = x$y[,1])
  out$se.logitpi <- sqrt(x$S[,1])
  out$ci.logitpi.lower <- x$y[,1] - q*out$se.logitpi
  out$ci.logitpi.upper <- x$y[,1] + q*out$se.logitpi

  if(x$logitrho) {
    out$logitrho <- x$y[,2]
    out$se.logitrho <- sqrt(x$S[,3])
    out$ci.logitrho.lower <- x$y[,2] - q*out$se.logitrho
    out$ci.logitrho.upper <- x$y[,2] + q*out$se.logitrho
    
    out$alphaKR21 <- x$size*plogis(x$y[,2])/(1+(x$size-1)*plogis(x$y[,2]))
    out$ci.alphaKR21.lower <- x$size*plogis(out$ci.logitrho.lower)/(1+(x$size-1)*plogis(out$ci.logitrho.lower))
    out$ci.alphaKR21.upper <- x$size*plogis(out$ci.logitrho.upper)/(1+(x$size-1)*plogis(out$ci.logitrho.upper))
  } else {
    out$rho <- x$y[,2]
    out$se.rho <- sqrt(x$S[,3])
    out$ci.rho.lower <- x$y[,2] - q*out$se.rho
    out$ci.rho.upper <- x$y[,2] + q*out$se.rho
    
    out$alphaKR21 <- x$size*x$y[,2]/(1+(x$size-1)*x$y[,2])
    out$ci.alphaKR21.lower <- x$size*out$ci.rho.lower/(1+(x$size-1)*out$ci.rho.lower)
    out$ci.alphaKR21.upper <- x$size*out$ci.rho.upper/(1+(x$size-1)*out$ci.rho.upper)
  }
  
  class(out) <- c("alphaKR21.maid", "data.frame")
  
  attributes(out)$ci.level <- ci.level
  
  return(out)
}


print.alphaKR21.maid <- function(x, digits = 4, ...){
  level_by_2 <- (1-attributes(x)$ci.level)/2
  
  xx <- round(x, digits)
  
  cat("Moment estimates of logit(pi) with standard error and \n ", 
      attributes(x)$ci.level*100, "% confidence intervals \n\n")
  xx1 <- xx[,1:4]
  colnames(xx1) <- c("logit(pi)", "se", paste(100*level_by_2, "%", sep = ""),
                     paste(100*(1-level_by_2), "%", sep = ""))
  print.data.frame(xx1)

  rho.name <- ifelse("logitrho" %in% colnames(x), "logitrho", "rho")
  cat("\n Moment estimates of", rho.name, "with standard error and \n ", 
      attributes(x)$ci.level*100, "% confidence intervals \n\n")
  xx2 <- xx[,5:8]
  colnames(xx2) <- c(rho.name, "se", paste(100*level_by_2, "%", sep = ""),
                     paste(100*(1-level_by_2), "%", sep = ""))
  print.data.frame(xx2)
  

  cat("\n Estimate of alpha (based on Kuder-Richardson formula 21 and \n  Beta-binomial assumption) with ", 
      attributes(x)$ci.level*100, "% confidence interval \n\n")
  xx3 <- xx[,9:11]
  colnames(xx3) <- c("alphaKR21", paste(100*level_by_2, "%", sep = ""),
                     paste(100*(1-level_by_2), "%", sep = ""))
  print.data.frame(xx3)

  return(invisible(NULL))
}