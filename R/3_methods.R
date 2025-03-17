#' @name BCSreg-methods
#' @title Methods for "\code{BCSreg}" Objects
#'
#' @param x,object an object of class "\code{BCSreg}".
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param formula a model formula or terms object or an R object.
#' @param model a character indicating which regression structure should be used.
#'     It can be \code{"mu"} for the scale regression structure, \code{"sigma"} for
#'     the relative dispersion regression structure, or \code{"full"} (when applicable)
#'     for both regression structures.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
NULL

## Model frame
#' @export
#' @rdname BCSreg-methods
model.frame.BCSreg <- function(formula, ...) {
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

## Model matrix
#' @export
#' @rdname BCSreg-methods
model.matrix.BCSreg <- function(object, model = c("mu", "sigma"), ...) {
  model <- match.arg(model, c("mu", "sigma"))
  val <- if (!is.null(object$x[[model]])) {
    object$x[[model]]
  } else {
    stats::model.matrix(object$terms[[model]], stats::model.frame(object),
                        contrasts = object$contrasts[[model]])
  }
  val
}

# Regression coefficients
#' @rdname BCSreg-methods
#' @export
coef.BCSreg <- function(object, model = c("mu", "sigma", "full"), ...) {
  model <- match.arg(model, c("mu", "sigma", "full"))
  cf <- object$coefficients
  switch(model,
         "full"  = list(mu = cf$mu, sigma = cf$sigma),
         "mu"    = cf$mu,
         "sigma" = cf$sigma)
}

#  Variance-covariance matrix
#' @rdname BCSreg-methods
#' @export
vcov.BCSreg <- function(object, model = c("mu", "sigma", "full"), ...) {
  model <- match.arg(model, c("mu", "sigma", "full"))
  covm <- object$vcov
  p <- length(object$coeff$mu)
  q <- length(object$coeff$sigma)
  switch(model,
         "mu" = {
           covm[seq.int(length.out = p), seq.int(length.out = p), drop = FALSE]
         },
         "sigma" = {
           covm[seq.int(length.out = q) + p, seq.int(length.out = q) + p, drop = FALSE]
         },
         "full" = {
           covm
         },
         )

}

# Log-likelihood
#' @rdname BCSreg-methods
#' @export
logLik.BCSreg <- function(object, ...) {
  p <- object$nobs - object$df.residual + as.numeric(!is.null(object$zeta))
  structure(object$loglik, df = p, class = "logLik")
}

# AIC
#' @export
#' @rdname BCSreg-methods
AIC.ugrpl <- function(object, ..., k = 2) {
  p <- object$nobs - object$df.residual + as.numeric(!is.null(object$zeta))
  AIC <- - 2 * object$logLik + k * p
  class(AIC) <- "AIC"
  return(AIC)
}

# Residuals
#' @export
#' @rdname BCSreg-methods
residuals.BCSreg <- function(object, ...) {
#residuals.BCSreg <- function(object, type = c("quantile"), ...) {
  # type <- match.arg(type)
  # y <- if (is.null(object$y)) stats::model.response(model.frame(object)) else object$y
  # X <- if (is.null(object$x$median)) stats::model.matrix(object, model = "mu") else object$x$mu
  # S <- if (is.null(object$x$dispersion)) stats::model.matrix(object, model = "sigma") else object$x$sigma
  # mu <- object$mu
  # lambda <- object$lambda
  # sigma <- object$sigma
  # family <- object$family
  # zeta <- object$zeta
  #
  # res <- switch(type,
  #               "quantile" = {
  #                 qnorm(pbcs(y, mu, sigma, lambda, zeta = zeta, family = family))
  #               }
  # )
  #
  # return(res)
  object$residuals
}

# Print
#' @rdname BCSreg-methods
#' @param digits a non-null value for digits specifies the minimum number of significant digits to
#'     be printed in values. The default, \code{getOption("digits")}.
#' @export
print.BCSreg <- function(x, digits = getOption("digits"), ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if (!x$converged) {
    cat("Model did not converge\n")
  } else {
    cat(paste("Scale submodel with ", x$link$mu, " link:\n", sep = ""))
    print.default(format(x$coefficients$mu, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")

    cat(paste("Relative dispersion submodel with ", x$link$sigma, " link:\n", sep = ""))
    print.default(format(x$coefficients$sigma, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")

    lambda_id <- is.null(x$control$lambda)
    zeta_id <- !is.null(x$zeta)

    if (lambda_id) {
      cat("Skewness parameter:\n", sep = "")
      print.default(format(x$lambda, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    }


    cat("---\nGenerating family:", x$family, if (zeta_id) paste0("(zeta = ", x$zeta, ")"),
        "\nLog-lik value:", x$loglik,
        "\nAIC:", stats::AIC(x),
        "and BIC:", stats::AIC(x, k = log(x$nobs)), "\n")
  }

  invisible(x)
}

# Summary
#' @rdname BCSreg-methods
#' @export
summary.BCSreg <- function(object, ...) {

  ## residuals
  residuals <- stats::residuals(object)

  ## extend coefficient tables
  est.beta <- stats::coef(object, "mu")
  se.beta <- sqrt(diag(stats::vcov(object, "mu")))
  est.tau <- stats::coef(object, "sigma")
  se.tau <- sqrt(diag(stats::vcov(object, "sigma")))

  mu <- cbind(Estimate = est.beta,
              `Std. error` = se.beta,
              `z value` = est.beta / se.beta,
              `Pr(>|z|)` = 2 * pnorm(-abs(est.beta / se.beta)))

  sigma <- cbind(Estimate = est.tau,
                 `Std. error` = se.tau,
                 `z value` = est.tau / se.tau,
                 `Pr(>|z|)` = 2 * pnorm(-abs(est.tau / se.tau)))


  if (is.null(object$control$lambda)) {
    est.lambda <- object$lambda
    se.lambda <- sqrt(utils::tail(diag(object$vcov), 1L))
    lambda <- cbind(Estimate = est.lambda,
                   `Std. error` = se.lambda,
                   `z value` = est.lambda / se.lambda,
                   `Pr(>|z|)` = 2 * pnorm(-abs(est.lambda / se.lambda)))
  } else {
    lambda <- object$lambda
  }

  ## number of iterations
  mytail <- function(x) x[length(x)]
  iterations <- c("optim" = as.vector(mytail(stats::na.omit(object$optim$count))))

  out <- list(call = object$call,
              mu = mu,
              sigma = sigma,
              lambda = lambda,
              zeta = object$zeta,
              link = object$link,
              family = object$family,
              converged = object$converged,
              iterations = iterations,
              loglik = object$loglik,
              df = object$nobs - object$df.residual + as.numeric(!is.null(object$zeta)),
              residuals = residuals,
              Upsilon.zeta = object$Upsilon.zeta,
              pseudo.r.squared = object$pseudo.r.squared,
              AIC = stats::AIC(object),
              BIC = stats::AIC(object, k = object$nobs))

  ## return
  class(out) <- "summary.BCSreg"
  out
}

# Print summary
#' @rdname BCSreg-methods
#' @export
print.summary.BCSreg <- function(x, digits = getOption("digits"), ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if (!x$converged) {
    cat("Model did not converge\n")
  } else {
    cat("Quantile residuals:\n")
    print(structure(round(as.vector(stats::quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")
    ))
    cat("\n")

    cat(paste("Scale submodel with ", x$link$mu, " link:\n", sep = ""))
    stats::printCoefmat(x$mu, digits = digits, signif.legend = FALSE)
    cat("\n")

    cat(paste("Relative dispersion submodel with ", x$link$sigma, " link:\n", sep = ""))
    stats::printCoefmat(x$sigma, digits = digits, signif.legend = FALSE)
    cat("\n")

    if (is.matrix(x$lambda)) {
      cat("Skewness parameter:\n", sep = "")
      stats::printCoefmat(x$lambda, digits = digits, signif.legend = FALSE)

      aux <- rbind(x$mu, x$sigma, x$lambda)[, 4]
      if (getOption("show.signif.stars") & any(aux < 0.1)) {
        cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
      }

      if (x$family == "NO" | x$family == "LOI" | x$family == "LOII") {
        cat("\nFamily: ", x$family, switch(x$family,
                                          "NO"   = " (Box-Cox normal)",
                                          "LOI"  = " (Box-Cox type I logistic)",
                                          "LOII" = " (Box-Cox type II logistic)"
        ), sep = "")
      } else {
        cat("\nFamily: ", x$family, "(", x$zeta, ")", switch(x$family,
                                           "ST" = paste0(" (Box-Cox t with zeta = ", x$zeta, ")"),
                                           "PE" = paste0(" (Box-Cox power exponential with zeta = ", x$zeta, ")"),
                                           "HP" = paste0(" (Box-Cox hyperbolic with zeta = ", x$zeta, ")"),
                                           "SL" = paste0(" (Box-Cox slash with zeta = ", x$zeta, ")"),
                                           "SN" = paste0(" (Box-Cox sinh-normal with zeta = ", x$zeta, ")")
        ), sep = "")
      }

    } else {

      aux <- rbind(x$mu, x$sigma)[, 4]
      if (getOption("show.signif.stars") & any(aux < 0.1)) {
        cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
      }

      if (x$family == "NO" | x$family == "LOI" | x$family == "LOII") {
        cat("\nFamily: ", x$family, switch(x$family,
                                           "NO"   = " (log-normal)",
                                           "LOI"  = " (log-type I logistic)",
                                           "LOII" = " (log-type II logistic)"
        ), sep = "")
      } else {
        cat("\nFamily: ", x$family, "(", x$zeta, ")", switch(x$family,
                                                             "ST" = paste0(" (log-t with zeta = ", x$zeta, ")"),
                                                             "PE" = paste0(" (log-power exponential with zeta = ", x$zeta, ")"),
                                                             "HP" = paste0(" (log-hyperbolic with zeta = ", x$zeta, ")"),
                                                             "SL" = paste0(" (log-slash with zeta = ", x$zeta, ")"),
                                                             "SN" = paste0(" (log-sinh-normal with zeta = ", x$zeta, ")")
        ), sep = "")
      }

    }

    cat(
      "\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", x$df, "Df"
    )
    if (!is.na(x$pseudo.r.squared)) cat("\nPseudo R-squared:", formatC(x$pseudo.r.squared, digits = digits))
    if (!is.na(x$Upsilon.zeta)) cat("\nUpsilon statistic:", formatC(x$Upsilon.zeta, digits = digits))
    cat("\nAIC:", formatC(x$AIC, digits = digits), "and BIC:", formatC(x$BIC, digits = digits))
    cat(paste0("\nNumber of iterations in ", x$method, "optimization: ", x$iterations[1L], "\n"))
  }

  invisible(x)
}







# summary.BCSreg <- function(object, type = "quantile", ...) {
#   ## residuals
#   type <- match.arg(type, c("quantile"))
#   object$residuals <- residuals(object, type = type)
#   object$residuals.type <- type
#
#   ## extend coefficient table
#   p <- length(object$coefficients$median)
#   q <- length(object$coefficients$dispersion)
#   cf <- as.vector(do.call("c", object$coefficients))
#   cf <- if (is.null(object$lambda)) cf else cf[-(p + q + 1)]
#   se <- sqrt(diag(object$vcov))
#   cf <- cbind(cf, se, cf / se, 2 * pnorm(-abs(cf / se)))
#   colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
#
#   cf <- if (is.null(object$lambda)) {
#     list(
#       median = cf[seq.int(length.out = p), , drop = FALSE],
#       dispersion = cf[seq.int(length.out = q) + p, , drop = FALSE],
#       shape = cf[seq.int(length.out = 1) + p + q, 1:2, drop = FALSE]
#     )
#   } else {
#     list(
#       median = cf[seq.int(length.out = p), , drop = FALSE],
#       dispersion = cf[seq.int(length.out = q) + p, , drop = FALSE]
#     )
#   }
#
#   rownames(cf$median) <- names(object$coefficients$median)
#   rownames(cf$dispersion) <- names(object$coefficients$dispersion)
#   if (is.null(object$lambda)) rownames(cf$shape) <- names(object$coefficients$shape)
#   object$coefficients <- cf
#
#   ## number of iterations
#   mytail <- function(x) x[length(x)]
#   object$iterations <- c("optim" = as.vector(mytail(na.omit(object$optim$count))))
#
#   ## AIC
#   object$AIC <- -2 * object$loglik + 2 * (p + q + ifelse(is.null(object$lambda), 1, 0))
#
#   ## delete some slots
#   object$fitted.values <- object$terms <- object$model <- object$y <-
#     object$x <- object$levels <- object$contrasts <- object$start <- NULL
#
#   ## return
#   class(object) <- "summary.BCSreg"
#   object
# }



# print.summary.BCSreg <- function(x, ...) {
#   cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
#   digits <- 4
#   if (!x$converged) {
#     cat("Model did not converge\n")
#   } else {
#     types <- c("quantile")
#     Types <- c("Quantile residuals")
#     cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
#     print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
#                     .Names = c("Min", "1Q", "Median", "3Q", "Max")
#     ))
#
#     if (NROW(x$coefficients$median)) {
#       cat(paste("\nCoefficients (median model with ", x$link$median$name, " link):\n", sep = ""))
#       printCoefmat(x$coefficients$median, digits = digits, signif.legend = FALSE)
#     } else {
#       cat("\nNo coefficients (in median model)\n")
#     }
#
#     if (NROW(x$coefficients$dispersion)) {
#       cat(paste("\nSigma coefficients (relative dispersion model with ", x$link$dispersion$name, " link):\n", sep = ""))
#       printCoefmat(x$coefficients$dispersion, digits = digits, signif.legend = FALSE)
#     } else {
#       cat("\nNo coefficients (in relative dispersion model)\n")
#     }
#
#     if (is.null(x$lambda)) {
#       cat(paste("\nLambda coefficient:\n", sep = ""))
#       printCoefmat(x$coefficients$shape, digits = digits, signif.legend = FALSE)
#     } else {
#       cat(paste("\nFixed shape parameter (lambda = ", x$lambda, ").\n", sep = ""))
#     }
#
#     aux <- x$coefficients[c("median", "dispersion")]
#
#     if (getOption("show.signif.stars") & any(do.call("rbind", aux)[, 4L] < 0.1)) {
#       cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
#     }
#
#     if (!is.null(x$lambda) && x$lambda == 0) {
#       if (x$family == "NO" | x$family == "LO") {
#         cat("\nFamily: log -", x$family, switch(x$family,
#                                                 "NO" = "(log normal)",
#                                                 "LO" = "(log type II logistic)"
#         ))
#       } else {
#         cat("\nFamily: log -", x$family, "(", x$zeta, ")", switch(x$family,
#                                                                   "TF"    = "(log Student-t)",
#                                                                   "PE"    = "(log power exponential)",
#                                                                   "Hyp"   = "(log hyperbolic)",
#                                                                   "SLASH" = "(log slash)",
#                                                                   "SN"    = "(log sinh-normal)"
#         ))
#       }
#     } else {
#       if (x$family == "NO" | x$family == "LO") {
#         cat("\nFamily: BCS -", x$family, switch(x$family,
#                                                 "NO" = "(Box-Cox normal)",
#                                                 "LO" = "(Box-Cox type II logistic)"
#         ))
#       } else {
#         cat("\nFamily: BCS -", x$family, "(", x$zeta, ")", switch(x$family,
#                                                                   "TF"    = "(Box-Cox Student-t)",
#                                                                   "PE"    = "(Box-Cox power exponential)",
#                                                                   "Hyp"   = "(Box-Cox hyperbolic)",
#                                                                   "SLASH" = "(Box-Cox slash)",
#                                                                   "SN"    = "(Box-Cox sinh-normal)"
#         ))
#       }
#     }
#
#
#     cat(
#       "\nLog-likelihood:", formatC(x$loglik, digits = digits),
#       "on", sum(sapply(x$coefficients, NROW)), "Df"
#     )
#     if (!is.na(x$pseudo.r.squared)) cat("\nPseudo R-squared:", formatC(x$pseudo.r.squared, digits = digits))
#     if (!is.na(x$Upsilon.zeta)) cat("\nUpsilon statistic:", formatC(x$Upsilon.zeta, digits = digits))
#     if (!is.na(x$AIC)) cat("\nAIC:", formatC(x$AIC, digits = digits))
#     cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations[1L], "\n"))
#   }
#
#   invisible(x)
# }

