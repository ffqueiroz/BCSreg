#' @name BCSreg_methods
#'
#' @title Extract Information From a Box-Cox Symmetric Regression Fit
#'
#' @description Methods for \code{"BCSreg"} objects.
#'
#' @param object an object of class \code{"BCSreg"}, a result of a call to \link{BCSreg}.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical Akaike information criteria (AIC).
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

# Model frame
#' @export
#' @rdname BCSreg_methods
model.frame.BCSreg <- function(formula, ...) {
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

## Model matrix
#' @export
#' @rdname BCSreg_methods
model.matrix.BCSreg <- function(object, model = c("mu", "sigma", "alpha"), ...) {
  model <- match.arg(model, c("mu", "sigma", "alpha"))
  if (model == "alpha" & is.null(object$alpha)) {
    stop("The model is not zero-adjusted")
  }
  val <- if (!is.null(object$x[[model]])) {
    object$x[[model]]
  } else {
    stats::model.matrix(object$terms[[model]], stats::model.frame(object),
                        contrasts = object$contrasts[[model]])
  }
  val
}

# Regression coefficients
#' @rdname BCSreg_methods
#' @export
coef.BCSreg <- function(object, model = c("mu", "sigma", "alpha", "full"), ...) {
  model <- match.arg(model, c("mu", "sigma", "alpha", "full"))
  cf <- object$coefficients
  if (model == "alpha" & is.null(object$alpha)) {
    stop("The model is not zero-adjusted")
  }
  switch(model,
         "full"  = list(mu = cf$mu, sigma = cf$sigma),
         "mu"    = cf$mu,
         "sigma" = cf$sigma,
         "alpha" = cf$alpha)
}

#  Variance-covariance matrix
#' @rdname BCSreg_methods
#' @export
vcov.BCSreg <- function(object, model = c("mu", "sigma", "alpha", "full"), ...) {
  model <- match.arg(model, c("mu", "sigma", "alpha", "full"))
  if (model == "alpha" & is.null(object$alpha)) {
    stop("The model is not zero-adjusted")
  }
  covm <- object$vcov
  p <- length(object$coeff$mu)
  q <- length(object$coeff$sigma)
  m <- length(object$coefficients$alpha)
  lambda_id <- is.null(object$control$lambda)
  switch(model,
         "mu" = {
           covm[seq.int(length.out = p), seq.int(length.out = p), drop = FALSE]
         },
         "sigma" = {
           covm[seq.int(length.out = q) + p, seq.int(length.out = q) + p, drop = FALSE]
         },
         "alpha" = {
           covm[seq.int(length.out = m) + p + q + as.numeric(lambda_id),
                seq.int(length.out = m) + p + q + as.numeric(lambda_id), drop = FALSE]
         },
         "full" = {
           covm
         },
  )

}

# Log-likelihood
#' @rdname BCSreg_methods
#' @export
logLik.BCSreg <- function(object, ...) {
  p <- object$nobs - object$df.residual
  structure(object$loglik, df = p, class = "logLik")
}

# AIC
#' @export
#' @rdname BCSreg_methods
AIC.ugrpl <- function(object, ..., k = 2) {
  p <- object$nobs - object$df.residual
  AIC <- - 2 * object$logLik + k * p
  class(AIC) <- "AIC"
  return(AIC)
}


# Residuals -----------------------------------------------------------------------------------
#' @name residuals
#'
#' @title Extract Residuals for a Box-Cox Symmetric Regression Fit
#'
#' @param object an object of class \code{"BCSreg"}, a result of a call to \link{BCSreg}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Function \code{residuals} returns a vector with the quantile residuals
#'     resulting from a Box-Cox symmetric regression fit.
#'
#' @export
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' ## Examples
residuals.BCSreg <- function(object, ...) {

  y <- as.numeric(stats::model.response(stats::model.frame(object)))
  ind <- ifelse(y == 0, 1, 0)
  n <- object$nobs

  mu <- object$mu
  sigma <- object$sigma
  lambda <- object$lambda
  zeta <- object$zeta
  alpha <- if (is.null(object$alpha)) rep(0L, n) else object$alpha
  family <- object$family

  residuals <- cdf <- rep(NA, n)
  cdf <- alpha
  cdf[ind == 0] <- alpha[ind == 0] + (1 - alpha[ind == 0]) *
    pbcs(y[ind == 0], mu = mu[ind == 0], sigma = sigma[ind == 0],
         lambda = lambda, zeta = zeta, family = family)

  u <- rep(NA, sum(ind))
  j <- 1L
  for (i in which(ind == 1)){
    u[j] <- stats::runif(1, 0, alpha[i])
    j <- j + 1
  }

  residuals[ind == 1] <- stats::qnorm(u)
  residuals[ind == 0] <- qnorm(cdf[ind == 0])
  residuals

}


# Summary method ------------------------------------------------------------------------------
#'
#' @name summary.BCSreg
#'
#' @title Summarizing a Box-Cox Symmetric Regression Fit
#'
#' @description \code{summary} method for class \code{"BCSreg"}.
#'
#' @param object an object of class \code{"BCSreg"}, a result of a call to \link{BCSreg}.
#' @param x an object of class \code{"summary.BCSreg"}, a result of a call to \link{summary.BCSreg}.
#' @param digits a non-null value for digits specifies the minimum number of significant digits to
#'     be printed in values.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The function \code{summary.BCSreg} returns an object of class \code{"summary.BCSreg"},
#'     which consists of a list with the following components:
#'  \describe{
#'     \item{call}{the original function call, given in \code{object}.}
#'     \item{mu}{summary statistics for the scale submodel.}
#'     \item{sigma}{summary statistics for the relative dispersion submodel.}
#'     \item{lambda}{summary statistics for the skewness parameter. If this parameter is not
#'         statistically different from zero, the fitted Box-Cox symmetric (BCS) distribution
#'         can be reduced to a more parsimonious log-symmetric distribution.}
#'     \item{alpha}{summary statistics for the zero-adjustment submodel when a zero-adjusted
#'         model is considered; and \code{NULL}, otherwise.}
#'     \item{zeta}{the specified value for the extra parameter of the fitted
#'         BCS distribution, if applicable.}
#'     \item{family}{the generating family of the fitted BCS distribution.}
#'     \item{link}{a list with elements \code{"mu"} and \code{"sigma"} with the
#'         specified link functions for the \code{mu} and \code{sigma} regression
#'         structures, respectively. If the model is zero-adjusted, the element
#'         \code{"alpha"} will also be returned with the link function for
#'         the regression structure of the zero-adjustment parameter.}
#'     \item{converged}{logical indicating successful convergence of the iterative
#'         process.}
#'     \item{iterations}{number of iterations reached until the optimization algorithm converges.}
#'     \item{logLik}{log-likelihood of the fitted model.}
#'     \item{df}{number of model parameters}
#'     \item{residuals}{a vector of quantile residuals.}
#'     \item{pseudo.r.squared}{pseudo R-squared value.}
#'     \item{Upsilon.zeta}{an overall goodness-of-fit measure.}
#'     \item{v}{a vector with the \eqn{v(z)} values for all the observations.
#'         The \eqn{v} function is a weighting function that is involved in the
#'         parameter estimation process. It depends on the generating density function
#'         of the corresponding BCS distribution.}
#'     \item{AIC, BIC}{Akaike and Bayesian information criteria.}
#'  }
#'
#' @export
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' ## Examples
summary.BCSreg <- function(object, ...) {

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

  if (!is.null(object$alpha)) {
    est.alpha <- stats::coef(object, "alpha")
    se.alpha <- sqrt(diag(stats::vcov(object, "alpha")))

    alpha <- cbind(Estimate = est.alpha,
                `Std. error` = se.alpha,
                `z value` = est.alpha / se.alpha,
                `Pr(>|z|)` = 2 * pnorm(-abs(est.alpha / se.alpha)))
  }

  ## Goodness-of-fit quantities
  y <- as.numeric(stats::model.response(stats::model.frame(object)))
  ind <- ifelse(y == 0, 1, 0)
  n <- object$nobs

  ## Upsilon statistic
  Upsilon <- function(zeta) {
    cdf <- sort(pbcs(y[ind == 0], mu = object$mu, sigma = object$sigma, lambda = object$lambda,
                     zeta = object$zeta, family = object$family))
    Upsilon_zeta <- mean(abs(qnorm(cdf) - EnvStats::evNormOrdStats(n = length(y[ind == 0]))))
    Upsilon_zeta
  }
  Upsilon.zeta <-  Upsilon(object$zeta)

  ## Pseudo-R2
  eta.1 <- stats::make.link(object$link$mu)$linkfun(object$mu)
  y.1 <- stats::make.link(object$link$mu)$linkfun(y[ind == 0])
  pseudo.r.squared <- ifelse(stats::var(eta.1) * stats::var(y.1) <= 0, NA, stats::cor(eta.1, y.1)^2)

  ## v-function
  v <- v.function(y[ind == 0], mu = object$mu, sigma = object$sigma, lambda = object$lambda,
                  zeta = object$zeta, family = object$family)$v

  ## Quantile residuals
  residuals <- stats::residuals(object)

  ## number of iterations
  mytail <- function(x) x[length(x)]
  iterations <- c("optim" = as.vector(mytail(stats::na.omit(object$optim$count))))

  out <- list(call = object$call,
              mu = mu,
              sigma = sigma,
              lambda = lambda,
              alpha = if (!is.null(object$alpha)) alpha else NULL,
              zeta = object$zeta,
              family = object$family,
              link = object$link,
              converged = object$converged,
              iterations = iterations,
              loglik = object$loglik,
              df = object$nobs - object$df.residual,
              residuals = residuals,
              pseudo.r.squared = pseudo.r.squared,
              Upsilon.zeta = Upsilon.zeta,
              v = v,
              AIC = stats::AIC(object),
              BIC = stats::AIC(object, k = object$nobs))

  ## return
  class(out) <- "summary.BCSreg"
  out
}

# Print summary
#' @rdname summary.BCSreg
#' @export
print.summary.BCSreg <- function(x, digits = getOption("digits"), ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),  sep = "\n")
  if (!x$converged) {
    cat("Model did not converge\n")
  } else {

    cat("\nQuantile residuals:\n")
    print(structure(round(as.vector(stats::quantile(x$residuals)), digits = digits),
                    .Names = c("Min", "1Q", "Median", "3Q", "Max")
    ))

    # Discrete component (if applicable)
    if (!is.null(x$alpha)) {
      cat("\n--- Fit for the discrete component ---\n")
      cat(paste("\nZero-adjustment submodel with ", x$link$alpha, " link:\n", sep = ""))
      stats::printCoefmat(x$alpha, digits = digits, signif.legend = FALSE)
      cat("\n--- Fit for the continuous component ---\n")
    }

    # Scale submodel
    cat(paste("\nScale submodel with ", x$link$mu, " link:\n", sep = ""))
    stats::printCoefmat(x$mu, digits = digits, signif.legend = FALSE)

    # Relative dispersion submodel
    cat(paste("\nRelative dispersion submodel with ", x$link$sigma, " link:\n", sep = ""))
    stats::printCoefmat(x$sigma, digits = digits, signif.legend = FALSE)

    # Skewness parameter
    if (is.matrix(x$lambda)) {
      cat("\nSkewness parameter:\n", sep = "")
      stats::printCoefmat(x$lambda, digits = digits, signif.legend = FALSE)
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

      if (x$family == "NO" | x$family == "LOI" | x$family == "LOII") {
        cat("\nGenerating family: ", x$family, switch(x$family,
                                                      "NO"   = " (Box-Cox normal)",
                                                      "LOI"  = " (Box-Cox type I logistic)",
                                                      "LOII" = " (Box-Cox type II logistic)"
        ), sep = "")
      } else {
        cat("\nGenerating family: ", x$family, "(", x$zeta, ")", switch(x$family,
                                                                        "ST" = paste0(" (Box-Cox t with zeta = ", x$zeta, ")"),
                                                                        "PE" = paste0(" (Box-Cox power exponential with zeta = ", x$zeta, ")"),
                                                                        "HP" = paste0(" (Box-Cox hyperbolic with zeta = ", x$zeta, ")"),
                                                                        "SL" = paste0(" (Box-Cox slash with zeta = ", x$zeta, ")"),
                                                                        "SN" = paste0(" (Box-Cox sinh-normal with zeta = ", x$zeta, ")")
        ), sep = "")
      }

    } else {

      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

      if (x$family == "NO" | x$family == "LOI" | x$family == "LOII") {
        cat("\nGenerating family: ", x$family, switch(x$family,
                                                      "NO"   = " (log-normal)",
                                                      "LOI"  = " (log-type I logistic)",
                                                      "LOII" = " (log-type II logistic)"
        ), sep = "")
      } else {
        cat("\nGenerating family: ", x$family, "(", x$zeta, ")", switch(x$family,
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



# Plot ----------------------------------------------------------------------------------------
#' Diagnostic Plots for a Box-Cox Symmetric Regression Fit
#'
#' This function provides plots for diagnostic analysis of a Box-Cox symmetric
#'     regression fit.
#'
#' @param x an object of class \code{"sdlrm"}.
#' @param which numeric; if a subset of the plots is required, specify a subset
#'     of the numbers \code{1:6}.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param pch,las,cex,lwd,... graphical parameters (see \code{\link[graphics]{par}})
#'
#' @details The \code{plot} method for \code{\link{BCSreg}} objects provides six types
#'     of diagnostic plots in the following order:
#'     \describe{
#'         \item{Residuals vs fitted values}{a plot of the residuals
#'             versus the fitted medians.}
#'         \item{Residuals vs observation indices.}{an index plot of the residuals
#'             versus the observation indices.}
#'         \item{Density plot}{a graph that compares the empirical density of the residuals
#'             with the density of the standard normal distribution.}
#'         \item{Normal probability plot}{a normal probability plot of the residuals with a
#'             confidence region constructed according to Fox (2016) using the
#'             \code{\link[car]{qqPlot}} function.}
#'         \item{Fitted vs observed values}{a dispersion diagram of the fitted values
#'             versus the observed values.}
#'         \item{Residuals vs v(z) function}{a dispersion diagram of the \eqn{v(z)} function
#'             versus the residuals. For someBCS models, the \eqn{v(z)} function
#'             may be interpreted as weights in the estimation process. If \code{family = "NO"},
#'             the \eqn{v(z)} function is constant.}
#'      }
#'
#'      The \code{which} argument can be used to select a subset of the implemented plots.
#'      Default is \code{which = 1:4}.
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @return \code{plot} method for \code{"\link{BCSreg}"} objects returns six types
#'     of diagnostic plots.
#'
#' @export
#' @importFrom graphics abline identify mtext par points text title curve rug legend
#'
#' @examples
#' ## Examples
plot.BCSreg <- function(x, which = 1:4,
                        ask = prod(graphics::par("mfcol")) < length(which) &&
                        grDevices::dev.interactive(),
                        pch = "+", las = 1, cex = 0.8, lwd = 2, ...)
{

  if(!is.numeric(which) || any(which < 1) || any(which > 6))
    stop("`which' must be in 1:6")

  ## Reading
  res <- stats::residuals(x)
  y <- stats::model.response(stats::model.frame(x))

  ## Graphical parameters setting
  if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op))
  }

  ## Plots to shown
  show <- rep(FALSE, 6)
  show[which] <- TRUE

  ## Residuals versus Fitted values
  if (show[1]){
    plot(stats::fitted(x), res, xlab = "Fitted values", ylab = "Quantile residuals",
         pch = pch, las = las, cex = cex, ...)
    abline(h = c(stats::qnorm(0.025), 0, stats::qnorm(0.975)),
           lty = c(2, 1, 2), lwd = lwd, col = "dodgerblue")
  }

  ## Residuals versus observation indices
  if (show[2]){
    n <- x$nobs
    plot(1:n, res, xlab = "Observation indices", ylab = "Quantile residuals",
         pch = pch, las = las, cex = cex, ...)
    abline(h = c(stats::qnorm(0.025), 0, stats::qnorm(0.975)),
           lty = c(2, 1, 2), lwd = lwd, col = "dodgerblue")
  }

  ## Density
  if(show[3]) {
    plot(stats::density(res), main = "", lwd = 2, ...)
    curve(stats::dnorm(x), col = "dodgerblue", add = TRUE, lwd = 2)
    rug(res)
  }

  ## Normal probability plot
  if(show[4]) {
    car::qqPlot(res, col.lines = "dodgerblue", grid = FALSE,
                pch = pch, las = las, cex = cex, ...,
                xlab = "Theoretical quantiles", ylab = "Sample quantiles", ...)
  }

  ## Fitted vs observed values
  if(show[5]) {
    plot(y, stats::fitted(x), xlab = "Observed values", ylab = "Fitted values",
         pch = pch, las = las, cex = cex, ...)
    abline(0, 1, lty = lwd, lwd = lwd, col = "dodgerblue")
  }

  ## ACF of residuals
  if(show[6]) {
    v <- v.function(y[y > 0], x$mu[y > 0], x$sigma[y > 0], x$lambda, x$zeta, x$family)$v
    if(x$family == "NO")
      warning("The v(z) function is constant for this family.")
    plot(res[y > 0], v, xlab = "Quantile residuals", ylab = expression(v(z)),
           pch = pch, las = las, cex = cex, ...)
  }

  invisible(x)

}


