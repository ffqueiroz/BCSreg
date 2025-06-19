#' @name BCSreg-methods
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
#' @returns
#' \itemize{
#' \item \code{model.frame} returns a \code{data.frame} containing the variables required
#'     by \code{formula} and any additional arguments provided via \code{...}.
#' \item \code{model.matrix} returns the design matrix used in the regression structure,
#'     as specified by the \code{model} argument.
#' \item \code{coef} returns a numeric vector of estimated regression coefficients, based
#'     on the \code{model} argument. If \code{parm = "full"}, it returns a list with the
#'     components \code{"mu"} and \code{"sigma"}, each containing the corresponding
#'     coefficient estimates. If the model is zero-adjusted, it will also have a
#'     \code{"alpha"} component with the estimates of the associated regression coefficients.
#' \item \code{vcov} returns the asymptotic covariance matrix of the regression coefficients,
#'     based on the \code{model} argument.
#' \item \code{logLik} returns the log-likelihood value of the fitted model.
#' \item \code{AIC} returns a numeric value representing the Akaike Information Criterion
#'     (AIC), Bayesian Information Criterion, or another criterion, depending on \code{k}.
#' }
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' ## Data set: fishery (for description, run ?fishery)
#' hist(fishery$cpue, xlab = "Catch per unit effort")
#' plot(cpue ~ tide_phase, fishery, pch = 16,
#'    xlab = "Tide phase", ylab = "Catch per unit effort")
#' plot(cpue ~ location, fishery, pch = 16,
#'    xlab = "Location", ylab = "Catch per unit effort")
#' plot(cpue ~ max_temp, fishery, pch = 16,
#'    xlab = "Maximum temperature", ylab = "Catch per unit effort")
#'
#' ## Fit the Box-Cox normal regression as a reference model
#' fit <- BCSreg(cpue ~ location + tide_phase |
#'                 location + tide_phase + max_temp, fishery)
#'
#' ## coef
#' coef(fit)
#' coef(fit, model = "sigma")
#' coef(fit, model = "full")
#'
#' ## vcov
#' vcov(fit)
#' vcov(fit, model = "sigma")
#' vcov(fit, model = "full")
#'
#' ## Log-likelihood value
#' logLik(fit)
#'
#' ## AIC and BIC
#' AIC(fit)
#' AIC(fit, k = log(fit$nobs))
#'
#' ## Model matrices
#' model.matrix(fit)
#' model.matrix(fit, model = "sigma")
NULL

# Model frame
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
#' @rdname BCSreg-methods
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
#' @rdname BCSreg-methods
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
#' @rdname BCSreg-methods
#' @export
logLik.BCSreg <- function(object, ...) {
  p <- object$nobs - object$df.residual
  structure(object$loglik, df = p, class = "logLik")
}

# AIC
#' @export
#' @rdname BCSreg-methods
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
#' @description Residuals resulting from fitting a Box-Cox symmetric or a zero-adjusted
#'     Box-Cox symmetric regression.
#'
#' @param object an object of class \code{"BCSreg"}, a result of a call to \link{BCSreg}.
#' @param approach a character string indicating the approach for calculating residuals
#'     when a zero-adjusted regression is fitted. Should be either \code{"combined"} (default)
#'     for combined residuals or \code{"separated"} for separate residuals. Ignored if
#'     the model is not zero-adjusted.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' If a Box-Cox symmetric regression is fitted to the data, it returns a numeric vector
#' containing the quantile residuals (Dunn and Smyth, 1996).
#'
#' If the model is a zero-adjusted Box-Cox symmetric regression:
#' \itemize{
#'   \item{For \code{approach = "combined"}, it returns a numeric vector with "combined" quantile residuals. See details}
#'   \item{For \code{approach = "separated"}, it returns a list with two components:
#'   \code{continuous} (quantile residuals for strictly positive responses) and
#'   \code{discrete} (standardized Pearson residuals for the discrete component).}
#' }
#'
#' @details
#'
#' For a Box-Cox symmetric regression fit, the residuals are the quantile residuals
#' (Dunn and Smyth, 1996), defined by \eqn{r_i^q = \Phi^{-1}(\widehat{F}(y_i))},
#' where \eqn{\widehat{F}(\cdot)} is the fitted cumulative distribution function and
#' \eqn{\Phi(\cdot)} is cumulative distribution function of the standard normal distribution.
#'
#' For zero-adjusted Box-Cox symmetric regressions, two approaches are available:
#' \itemize{
#'   \item{\strong{Combined approach}: Returns a single vector of residuals defined as
#'
#'   \eqn{
#'   r_i^q =
#'   \begin{cases}
#'     \Phi^{-1}(u_i), & y_i = 0, \\
#'     \Phi^{-1}\left[\widehat{F}^{(0)}(y_i)\right], & y_i > 0,
#'   \end{cases}
#'   }
#'
#'   where \eqn{u_i} is a random variable uniformly distributed in \eqn{(0, \widehat{\alpha}_i]}
#'   and \eqn{F^{(0)}} is the fitted cumulative distribution function of the mixed response.}
#'
#'   \item{\strong{Separated approach}: Returns a list containing:
#'   \itemize{
#'     \item{Quantile residuals for the positive (continuous) component.}
#'     \item{Standardized Pearson residuals for the discrete component, defined by
#'
#'     \eqn{
#'     r_i^p = \frac{\mathbb{I}(y_i = 0) - \widehat{\alpha}_i}
#'     {\sqrt{\widehat{\alpha}_i(1-\widehat{\alpha}_i)(1-\widehat{h}_{ii})}},
#'     }
#'
#'     where \eqn{\widehat{h}_{ii}} is the \eqn{i}th diagonal element of the
#'     "hat matrix" resulting from a fit of a generalized linear model with a
#'     binary response given by \eqn{\mathbb{I}(y_i = 0)}, being \eqn{\mathbb{I}} the indicator
#'     function.}
#'   }
#'   }
#' }
#'
#' See more details in Medeiros and Queiroz (2025).
#'
#' @export
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @references
#'     Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#'     \emph{Journal of Computational and Graphical Statistics}, \bold{5}, 236---244.
#'
#'     Medeiros, R. M. R., and Queiroz, F. F. (2025). Modeling positive continuous data:
#'     Box-Cox symmetric regression models and their extensions
#'
#' @examples
#' ## Data set: fishery (for description, run ?fishery)
#' hist(fishery$cpue, xlab = "Catch per unit effort")
#' plot(cpue ~ tide_phase, fishery, pch = 16,
#'     xlab = "Tide phase", ylab = "Catch per unit effort")
#' plot(cpue ~ location, fishery, pch = 16,
#'     xlab = "Location", ylab = "Catch per unit effort")
#' plot(cpue ~ max_temp, fishery, pch = 16,
#'     xlab = "Maximum temperature", ylab = "Catch per unit effort")
#'
#' ## BCS fit
#' fit <- BCSreg(cpue ~ location + tide_phase + max_temp |
#'                 location + tide_phase + max_temp, fishery)
#'
#' ## Quantile residuals
#' rq <- residuals(fit)
#' rq
#'
#' ## Normal probability plot
#' qqnorm(rq, pch = "+", cex = 0.8)
#' qqline(rq, col = "dodgerblue", lwd = 2)
residuals.BCSreg <- function(object, approach = c("combined", "separated"), ...) {

  approach <- match.arg(approach, c("combined", "separated"))

  y <- as.numeric(stats::model.response(stats::model.frame(object)))
  n <- object$nobs
  mu <- object$mu
  sigma <- object$sigma
  lambda <- object$lambda
  zeta <- object$zeta
  family <- object$family
  alpha <- if (is.null(object$alpha)) rep(0L, n) else object$alpha
  is_zero_adjusted <- any(alpha > 0)

  ind <- as.numeric(y == 0)

  if (!is_zero_adjusted) {
    # Quantile residuals
    cdf <- pBCS(y, mu = mu, sigma = sigma, lambda = lambda,
                zeta = zeta, family = family)
    residuals <- stats::qnorm(cdf)
    return(residuals)
  }

  # Zero-adjusted case
  if (approach == "combined") {

    cdf <- pZABCS(y, alpha = alpha, mu = mu, sigma = sigma,
                  lambda = lambda, zeta = zeta, family = family)

    residuals <- numeric(n)
    u <- numeric(sum(ind))
    j <- 1L
    for (i in which(ind == 1)){
      u[j] <- stats::runif(1, 0, alpha[i])
      j <- j + 1
    }

    residuals[ind == 1] <- stats::qnorm(u)
    residuals[ind == 0] <- stats::qnorm(cdf[ind == 0])
    return(residuals)

  } else if (approach == "separated") {

    # Residuals for the continuous component (positive responses)
    cdf_continuous <- pBCS(y[ind == 0], mu = mu[ind == 0], sigma = sigma[ind == 0],
                           lambda = lambda, zeta = zeta, family = family)
    residuals_continuous <- stats::qnorm(cdf_continuous)

    # Residuals for the discrete part
    Z <- stats::model.matrix(object, model = "alpha")
    dalpha.deta <- c(stats::make.link(object$link$alpha)$mu.eta(Z %*% object$coefficients$alpha))
    Q <- diag(dalpha.deta^2 / (alpha * (1 - alpha)))
    h_diag <- diag(sqrt(Q)%*%Z%*%solve(t(Z)%*%Q%*%Z)%*%t(Z)%*%sqrt(Q))

    residuals_discrete <- (ind - alpha) /
      sqrt(alpha * (1 - alpha) * (1 - h_diag))

    return(list(
      continuous = residuals_continuous,
      discrete = residuals_discrete
    ))
  }
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
#'     \item{v}{a vector with the \eqn{v(z)} values for all the observations.}
#'     \item{AIC, BIC}{Akaike and Bayesian information criteria.}
#'  }
#'
#' @details
#' An object of class \code{"summary.BCSreg"} provides additional information from a Box-Cox
#'     symmetric or zero-adjusted Box-Cox symmetric regression fit. In addition to summary
#'     tables with estimates, standard errors, and individual significance tests of the
#'     regression coefficients, it also provides the quantile residuals (Dunn and Smyth, 1996)
#'     and goodness-of-fit measures.
#'
#' The goodness-of-fit measures include:
#'
#'   \itemize{
#'     \item{A pseudo-\eqn{R_p^2} defined based on Ferrari and Cribari-Neto (2004) as the
#'         square of the sample correlation coefficient between \eqn{d_1(y)} and
#'         \eqn{d_1(\widehat{y})}, where \eqn{y} is the response variable,
#'         \eqn{\widehat{y}} is the adjusted median, and \eqn{d_1()} is the link
#'         function (\code{link}) used in the regression structure for the scale parameter.
#'         By definition, \eqn{0 \leq R_p^2 \leq 1}, and perfect agreement between the fitted
#'         values and the response yields \eqn{R_p^2 = 1}.}
#'
#'     \item{A general goodness-of-fit measure based on Vanegas and Paula (2016), defined as
#'
#'       \deqn{
#'         \Upsilon_\zeta = \frac{1}{n} \sum_{i=1}^n \left| \Phi^{-1}[\widehat{F}(y_{(i)})] - \nu_{(i)} \right|,
#'       }
#'
#'       where \eqn{y_{(i)}} is the \eqn{i}th order statistic of the response, \eqn{\nu_{(i)}}
#'       is the expected value of the \eqn{i}th order statistic from a standard normal sample
#'       of size \eqn{n}, \eqn{\Phi(\cdot)} denotes the cumulative distribution function
#'       of the standard normal distribution, and \eqn{\widehat{F}(\cdot)} denotes the fitted cumulative
#'       distribution function of the assumed distribution for the response.}
#'
#'     \item{A weighting function \eqn{v(\cdot)} that plays a role in the parameter estimation p
#'         rocess, providing information about the individual contribution of each
#'         observation. It depends on the generating density function and is constant for
#'         distributions generated by the normal family (\code{family = "NO"}).}
#'   }
#' @export
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @references
#'     Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#'     \emph{Journal of Computational and Graphical Statistics}, \bold{5}, 236---244.
#'
#'     Ferrari, S., and Cribari-Neto, F. (2004). Beta regression for modelling
#'     rates and proportions. \emph{Journal of Applied Statistics}, \bold{31}, 799---815.
#'
#'     Medeiros, R. M. R., and Queiroz, F. F. (2025). Modeling positive continuous data:
#'     Box-Cox symmetric regression models and their extensions
#'
#'     Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions:
#'     statistical properties and parameter estimation. \emph{Brazilian Journal of
#'     Probability and Statistics}, \bold{30},196---220
#'
#' @examples
#' ## Data set: fishery (for description, run ?fishery)
#' hist(fishery$cpue, xlab = "Catch per unit effort")
#' plot(cpue ~ tide_phase, fishery, pch = 16,
#'     xlab = "Tide phase", ylab = "Catch per unit effort")
#' plot(cpue ~ location, fishery, pch = 16,
#'     xlab = "Location", ylab = "Catch per unit effort")
#' plot(cpue ~ max_temp, fishery, pch = 16,
#'     xlab = "Maximum temperature", ylab = "Catch per unit effort")
#'
#' ## Fit examples
#'
#' ### Fit a single Box-Cox normal regression model:
#' fit_bcno1 <- BCSreg(cpue ~ location + tide_phase + max_temp, fishery)
#' summary(fit_bcno1)
#' # Other quantities can be obtained from a summary.BCSreg object
#' aux <- summary(fit_bcno1)
#' class(aux)
#' str(aux)
#'
#' ### Fit a double Box-Cox normal regression model:
#' fit_bcno2 <- BCSreg(cpue ~ location + tide_phase |
#'                      location + tide_phase + max_temp, fishery)
#' summary(fit_bcno2)
#'
#'
#' ### Fit a double Box-Cox power exponential regression model:
#' fit_bcpe <- BCSreg(cpue ~ location + tide_phase + max_temp |
#'                     location + tide_phase + max_temp, fishery, family = "PE", zeta = 4)
#' summary(fit_bcpe)
#'
#' ### Fit a double log-power exponential regression model:
#' fit_lpe <- BCSreg(cpue ~ location + tide_phase + max_temp |
#'                    location + tide_phase + max_temp, fishery, family = "PE",
#'                  zeta = 4, lambda = 0)
#' summary(fit_lpe)
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
    cdf <- sort(pBCS(y[ind == 0], mu = object$mu[ind == 0],
                     sigma = object$sigma[ind == 0], lambda = object$lambda,
                     zeta = object$zeta, family = object$family))
    Upsilon_zeta <- mean(abs(qnorm(cdf) - EnvStats::evNormOrdStats(n = length(y[ind == 0]))))
    Upsilon_zeta
  }
  Upsilon.zeta <-  Upsilon(object$zeta)

  ## Pseudo-R2
  y.1 <- stats::make.link(object$link$mu)$linkfun(y[ind == 0])
  if (is.null(object$alpha)) {
    yhat.1 <- stats::make.link(object$link$mu)$linkfun(qBCS(0.5, object$mu, object$sigma,
                                                            object$lambda, object$zeta, object$family))[ind == 0]
  } else {
    yhat.1 <- stats::make.link(object$link$mu)$linkfun(qZABCS(0.5, object$alpha, object$mu, object$sigma,
                                                            object$lambda, object$zeta, object$family))[ind == 0]
  }

  pseudo.r.squared <- ifelse(stats::var(yhat.1) * stats::var(y.1) <= 0, NA, stats::cor(yhat.1, y.1)^2)

  ## v-function
  v <- v.function(y[ind == 0], mu = object$mu[ind == 0], sigma = object$sigma[ind == 0], lambda = object$lambda,
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
#'     or a zero-adjusted regression fit.
#'
#' @param x an object of class \code{"BCSreg"}.
#' @param which numeric; if a subset of the plots is required, specify a subset
#'     of the numbers \code{1:7}.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param pch,las,cex,lwd,... graphical parameters (see \code{\link[graphics]{par}})
#'
#' @details The \code{plot} method for \code{\link{BCSreg}} objects provides seven types
#'     of diagnostic plots in the following order:
#'     \describe{
#'         \item{Residuals vs fitted values}{a plot of the residuals
#'             versus the fitted medians.}
#'         \item{Residuals vs observation indices.}{an index plot of the residuals
#'             versus the observation indices.}
#'         \item{Normal probability plot}{a normal probability plot of the residuals with a
#'             confidence region constructed according to Fox (2016) using the
#'             \code{\link[car]{qqPlot}} function.}
#'         \item{Case-weight perturbation}{An index plot of local influence based on the
#'             case-weight perturbation scheme.}
#'         \item{Density plot}{a graph that compares the empirical density of the residuals
#'             with the density of the standard normal distribution.}
#'         \item{Fitted vs observed values}{a dispersion diagram of the fitted values
#'             versus the observed values.}
#'         \item{Residuals vs v(z) function}{a dispersion diagram of the \eqn{v(z)} function
#'             versus the residuals. For some BCS models, the \eqn{v(z)} function
#'             may be interpreted as weights in the estimation process. If \code{family = "NO"},
#'             the \eqn{v(z)} function is constant.}
#'      }
#'
#'      The \code{which} argument can be used to select a subset of the implemented plots.
#'      Default is \code{which = 1:4}. See \code{\link{residuals.BCSreg}} for details on
#'      the residuals.
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @return \code{plot} method for \code{"\link{BCSreg}"} objects returns seven types
#'     of diagnostic plots.
#'
#' @export
#' @importFrom graphics abline identify mtext par points text title curve rug legend
#'
#' @examples
#' ## Data set: fishery (for description, run ?fishery)
#' hist(fishery$cpue, xlab = "Catch per unit effort")
#' plot(cpue ~ tide_phase, fishery, pch = 16,
#'     xlab = "Tide phase", ylab = "Catch per unit effort")
#' plot(cpue ~ location, fishery, pch = 16,
#'     xlab = "Location", ylab = "Catch per unit effort")
#' plot(cpue ~ max_temp, fishery, pch = 16,
#'     xlab = "Maximum temperature", ylab = "Catch per unit effort")
#'
#' ## Fit a double Box-Cox normal regression model:
#' fit <- BCSreg(cpue ~ location + tide_phase |
#'                 location + tide_phase + max_temp, fishery)
#'
#' ## Available plots:
#'
#' ### Residuals vs fitted values (fitted medians)
#' plot(fit, which = 1)
#'
#' ### Residuals vs observation indices
#' plot(fit, which = 2)
#'
#' ### Normal probability plot
#' plot(fit, which = 3)
#'
#' ### Local influence
#' plot(fit, which = 4)
#'
#' ### Density plot
#' plot(fit, which = 5)
#'
#' ### Fitted medians vs response
#' plot(fit, which = 6)
#'
#' ### v(z) function
#' plot(fit, which = 7)
plot.BCSreg <- function(x, which = 1:4,
                        ask = prod(graphics::par("mfcol")) < length(which) &&
                        grDevices::dev.interactive(),
                        pch = "+", las = 1, cex = 0.8, lwd = 2, ...)
{

  dots <- list(...)

  if(!is.numeric(which) || any(which < 1) || any(which > 7))
    stop("`which' must be in 1:7")

  ## Reading
  res <- stats::residuals(x)
  y <- stats::model.response(stats::model.frame(x))
  n <- x$nobs

  ## Graphical parameters setting
  if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op))
  }

  ## Plots to shown
  show <- rep(FALSE, 7)
  show[which] <- TRUE

  ## Residuals versus Fitted values
  if (show[1]){

    xlab <- if (is.null(dots$xlab)) "Fitted values" else dots$xlab
    ylab <- if (is.null(dots$ylab)) "Quantile residuals" else dots$ylab

    plot(stats::fitted(x), res, xlab = xlab, ylab = ylab,
         pch = pch, las = las, cex = cex, ...)
    abline(h = c(-2.5, 0, 2.5), lty = c(2, 1, 2), lwd = lwd, col = "dodgerblue")
  }

  ## Residuals versus observation indices
  if (show[2]){

    xlab <- if (is.null(dots$xlab)) "Observation indices" else dots$xlab
    ylab <- if (is.null(dots$ylab)) "Quantile residuals" else dots$ylab

    plot(1:n, res, xlab = xlab, ylab = ylab,
         pch = pch, las = las, cex = cex, ...)
    abline(h = c(-2.5, 0, 2.5), lty = c(2, 1, 2), lwd = lwd, col = "dodgerblue")
  }

  ## Normal probability plot
  if(show[3]) {

    ## Normal probability plot
    Pi <- (1:n - 0.5) / n
    zi <- stats::qnorm(Pi)
    s.h <- stats::IQR(res) / 1.349
    m.h <- stats::median(res)
    xi <- m.h + s.h * zi
    se <- s.h * sqrt(Pi * (1 - Pi) / n) / stats::dnorm(zi)

    ## qqline
    qqx <- stats::qnorm(c(0.25, 0.75))
    qqy <- stats::quantile(res, probs = c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
    slope <- diff(qqy) / diff(qqx)
    int <- qqy[[1L]] - slope * qqx[[1L]]

    ## Plot
    xlab <- if (is.null(dots$xlab)) "Theoretical quantiles" else dots$xlab
    ylab <- if (is.null(dots$ylab)) "Sample quantiles" else dots$ylab
    ylim <- if (is.null(dots$xlim)) c(min(xi - stats::qnorm(0.975) * se),
                                      max(xi + stats::qnorm(0.975) * se)) else dots$ylim

    plot(zi, sort(res), xlab = xlab, ylab = ylab, type = "n", ylim = ylim, ...)
    graphics::polygon(x = c(zi, rev(zi)),
            y = c(xi - stats::qnorm(0.975) * se,
                  rev(xi + stats::qnorm(0.975) * se)),
            border = NA,
            col = grDevices::rgb(30, 144, 255, alpha = 0.2 * 255, maxColorValue = 255))
    graphics::points(zi, sort(res), pch = pch, las = las, cex = cex, ...)
    graphics::segments(min(zi), int + slope * min(zi),
                       max(zi), int + slope * max(zi),
                       col = "dodgerblue", lwd = 2)
  }

  ## Local influence
  cw <- influence(x, plot = FALSE)$case.weights
  if(show[4]) {

    xlab <- if (is.null(dots$xlab)) "Observation indices" else dots$xlab
    ylab <- if (is.null(dots$ylab)) "Local influence" else dots$ylab

    plot(1:n, cw, xlab = xlab, ylab = ylab, type = "h", las = las, ...)
  }

  ## Density
  if(show[5]) {

    main <- if (is.null(dots$main)) " " else dots$main
    xlab <- if (is.null(dots$xlab)) "Quantile residuals" else dots$xlab

    plot(stats::density(res), main = main, xlab = xlab, lwd = 2, ...)
    curve(stats::dnorm(x), col = "dodgerblue", add = TRUE, lwd = 2)
    legend("topright", c("Sample density", "N(0, 1)"),
           lty = 1, lwd = 2, col = c(1, "dodgerblue"), cex = 0.9, bty = "n")
    rug(res)
  }

  ## Fitted vs observed values
  if(show[6]) {
    plot(y, stats::fitted(x), xlab = "Observed values", ylab = "Fitted values",
         pch = pch, las = las, cex = cex, ...)
    abline(0, 1, lty = lwd, lwd = lwd, col = "dodgerblue")
  }

  ## v(z) function
  if(show[7]) {

    xlab <- if (is.null(dots$xlab)) "Quantile residuals" else dots$xlab
    ylab <- if (is.null(dots$ylab)) expression(v(z)) else dots$ylab

    v <- v.function(y[y > 0], x$mu[y > 0], x$sigma[y > 0], x$lambda, x$zeta, x$family)$v
    if(x$family == "NO")
      warning("The v(z) function is constant for this family.")
    plot(res[y > 0], v, xlab = xlab, ylab = ylab,
           pch = pch, las = las, cex = cex, ...)
  }

  invisible(x)

}


