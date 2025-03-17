summary.BCSreg <- function(object, type = "quantile", ...) {
  ## residuals
  type <- match.arg(type, c("quantile"))
  object$residuals <- residuals(object, type = type)
  object$residuals.type <- type

  ## extend coefficient table
  p <- length(object$coefficients$median)
  q <- length(object$coefficients$dispersion)
  cf <- as.vector(do.call("c", object$coefficients))
  cf <- if (is.null(object$lambda)) cf else cf[-(p + q + 1)]
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf / se, 2 * pnorm(-abs(cf / se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  cf <- if (is.null(object$lambda)) {
    list(
      median = cf[seq.int(length.out = p), , drop = FALSE],
      dispersion = cf[seq.int(length.out = q) + p, , drop = FALSE],
      shape = cf[seq.int(length.out = 1) + p + q, 1:2, drop = FALSE]
    )
  } else {
    list(
      median = cf[seq.int(length.out = p), , drop = FALSE],
      dispersion = cf[seq.int(length.out = q) + p, , drop = FALSE]
    )
  }

  rownames(cf$median) <- names(object$coefficients$median)
  rownames(cf$dispersion) <- names(object$coefficients$dispersion)
  if (is.null(object$lambda)) rownames(cf$shape) <- names(object$coefficients$shape)
  object$coefficients <- cf

  ## number of iterations
  mytail <- function(x) x[length(x)]
  object$iterations <- c("optim" = as.vector(mytail(na.omit(object$optim$count))))

  ## AIC
  object$AIC <- -2 * object$loglik + 2 * (p + q + ifelse(is.null(object$lambda), 1, 0))

  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.BCSreg"
  object
}

print.BCSreg <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  digits <- 4
  if (!x$converged) {
    cat("Model did not converge\n")
  } else {
    if (length(x$coefficients$median)) {
      cat(paste("Coefficients (median model with ", x$link$median$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$median, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in median model)\n\n")
    }
    if (length(x$coefficients$dispersion)) {
      cat(paste("Sigma coefficients (relative dispersion model with ", x$link$dispersion$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$dispersion, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in dispersion model)\n\n")
    }
    if (is.null(x$lambda)) {
      cat(paste("Lambda coefficients:\n", sep = ""))
      print.default(format(x$coefficients$shape, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    }
  }
  invisible(x)
}

print.summary.BCSreg <- function(x, ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  digits <- 4
  if (!x$converged) {
    cat("Model did not converge\n")
  } else {
    types <- c("quantile")
    Types <- c("Quantile residuals")
    cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")
    ))

    if (NROW(x$coefficients$median)) {
      cat(paste("\nCoefficients (median model with ", x$link$median$name, " link):\n", sep = ""))
      printCoefmat(x$coefficients$median, digits = digits, signif.legend = FALSE)
    } else {
      cat("\nNo coefficients (in median model)\n")
    }

    if (NROW(x$coefficients$dispersion)) {
      cat(paste("\nSigma coefficients (relative dispersion model with ", x$link$dispersion$name, " link):\n", sep = ""))
      printCoefmat(x$coefficients$dispersion, digits = digits, signif.legend = FALSE)
    } else {
      cat("\nNo coefficients (in relative dispersion model)\n")
    }

    if (is.null(x$lambda)) {
      cat(paste("\nLambda coefficient:\n", sep = ""))
      printCoefmat(x$coefficients$shape, digits = digits, signif.legend = FALSE)
    } else {
      cat(paste("\nFixed shape parameter (lambda = ", x$lambda, ").\n", sep = ""))
    }

    aux <- x$coefficients[c("median", "dispersion")]

    if (getOption("show.signif.stars") & any(do.call("rbind", aux)[, 4L] < 0.1)) {
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    }

    if (!is.null(x$lambda) && x$lambda == 0) {
      if (x$family == "NO" | x$family == "LO") {
        cat("\nFamily: log -", x$family, switch(x$family,
          "NO" = "(log normal)",
          "LO" = "(log type II logistic)"
        ))
      } else {
        cat("\nFamily: log -", x$family, "(", x$zeta, ")", switch(x$family,
          "TF"    = "(log Student-t)",
          "PE"    = "(log power exponential)",
          "Hyp"   = "(log hyperbolic)",
          "SLASH" = "(log slash)",
          "SN"    = "(log sinh-normal)"
        ))
      }
    } else {
      if (x$family == "NO" | x$family == "LO") {
        cat("\nFamily: BCS -", x$family, switch(x$family,
          "NO" = "(Box-Cox normal)",
          "LO" = "(Box-Cox type II logistic)"
        ))
      } else {
        cat("\nFamily: BCS -", x$family, "(", x$zeta, ")", switch(x$family,
          "TF"    = "(Box-Cox Student-t)",
          "PE"    = "(Box-Cox power exponential)",
          "Hyp"   = "(Box-Cox hyperbolic)",
          "SLASH" = "(Box-Cox slash)",
          "SN"    = "(Box-Cox sinh-normal)"
        ))
      }
    }


    cat(
      "\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df"
    )
    if (!is.na(x$pseudo.r.squared)) cat("\nPseudo R-squared:", formatC(x$pseudo.r.squared, digits = digits))
    if (!is.na(x$Upsilon.zeta)) cat("\nUpsilon statistic:", formatC(x$Upsilon.zeta, digits = digits))
    if (!is.na(x$AIC)) cat("\nAIC:", formatC(x$AIC, digits = digits))
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations[1L], "\n"))
  }

  invisible(x)
}


coef.BCSreg <- function(object, ...) {
  cf <- object$coefficients
  name1 <- names(cf$median)
  name2 <- names(cf$dispersion)
  name3 <- names(cf$shape)
  cf <- c(cf$median, cf$dispersion, cf$shape)
  names(cf) <- c(name1, paste("(sigma)", name2, sep = "_"), name3)
  if (is.null(object$lambda)) cf else cf[-length(cf)]
}


vcov.BCSreg <- function(object, ...) {
  object$vcov
}

logLik.BCSreg <- function(object, ...) {
  structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")
}


model.matrix.BCSreg <- function(object, model = c("median", "dispersion"), ...) {
  model <- match.arg(model)
  val <- if (!is.null(object$x[[model]])) {
    object$x[[model]]
  } else {
    model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  }
  return(val)
}


residuals.BCSreg <- function(object,
                             type = c("quantile"), ...) {
  type <- match.arg(type)
  y <- if (is.null(object$y)) model.response(model.frame(object)) else object$y
  X <- if (is.null(object$x$median)) model.matrix(object, model = "median") else object$x$median
  S <- if (is.null(object$x$dispersion)) model.matrix(object, model = "dispersion") else object$x$dispersion
  mu <- object$fitted.values
  lambda <- object$coefficients$shape
  sigma <- object$link$dispersion$linkinv(as.vector(S %*% object$coefficients$dispersion))
  family <- object$family
  zeta <- object$zeta

  res <- switch(type,
    "quantile" = {
      qnorm(pbcs(y, mu, sigma, lambda, zeta = zeta, family = family))
    }
  )

  return(res)
}
