#' Influence Diagnostics for BCSreg Objects
#'
#' The \code{influence} function provides two influence measures for a Box-Cox symmetric or a
#'     zero-adusted Box-Cox symmetric regression fit.
#'
#' @param object an object of class \code{"BCSreg"}.
#' @param plot logical. If \code{plot = TRUE} (default), the plots are shown.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot, if \code{plot = TRUE}.
#' @param ... currently not used.
#'
#' @return \code{influence} returns a list with two objects:
#'     \item{case.weights}{The values of \eqn{d_{max}} eigenvector based on case
#'     weights perturbation scheme (see Medeiros and Queiroz (2025)).}
#'     \item{totalLI}{The total local influence (see Lesaffre and Verbeke (1998)).}
#'
#' @seealso \code{\link{BCSreg}} for parameter estimation in the class of the Box-Cox
#'     symmetric or zero-adjusted Box-Cox symmetric regression models,
#'     \code{\link{residuals.BCSreg}} for extracting residuals for a fitted model, and
#'     \code{\link{plot.BCSreg}} for diagnostic plots.
#'
#' @references
#'  Lesaffre, E., and Verbeke, G. (1998). Local influence in linear mixed models.
#'  \emph{Biometrics}, 570--582.
#'
#'  Medeiros, R. M. R., and Queiroz, F. F. (2025). Flexible modeling of non-negative continuous
#'     data: Box-Cox symmetric regression and its zero-adjusted extension.
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' ## Data set: raycatch (for description, run ?raycatch)
#' hist(raycatch$cpue, xlab = "Catch per unit effort")
#' plot(cpue ~ tide_phase, raycatch, pch = 16,
#'    xlab = "Tide phase", ylab = "Catch per unit effort")
#' plot(cpue ~ location, raycatch, pch = 16,
#'    xlab = "Location", ylab = "Catch per unit effort")
#' plot(cpue ~ max_temp, raycatch, pch = 16,
#'    xlab = "Maximum temperature", ylab = "Catch per unit effort")
#'
#' ## Fit a double Box-Cox normal regression model:
#' fit <- BCSreg(cpue ~ location + tide_phase |
#'                location + tide_phase + max_temp, raycatch)
#'
#' ## Influence measures under case-weights perturbation scheme:
#' cw <- influence(fit) ## two index plots are shown
#' str(cw)
#' @export
influence <- function(object, plot = TRUE, ask = grDevices::dev.interactive(), ...){

  y <- if(is.null(object$y)) stats::model.response(stats::model.frame(object)) else object$y
  X <- if(is.null(object$x$mu)) stats::model.matrix(object, model = "mu") else object$x$mu
  S <- if(is.null(object$x$sigma)) stats::model.matrix(object, model = "sigma") else object$x$sigma

  alpha <- object$alpha
  mu <- object$mu
  sigma <- object$sigma
  lambda <- object$lambda
  zeta <- object$zeta
  family <- object$family
  mu.link <- object$link$mu
  sigma.link <- object$link$sigma
  J      <- solve(object$vcov)
  beta   <- object$coefficients$mu
  tau    <- object$coefficients$sigma

  p <- ncol(X)
  q <- ncol(S)
  n <- nrow(X)

  T1 <- diag(c(stats::make.link(mu.link)$mu.eta(X%*%beta)))
  T2 <- diag(c(stats::make.link(sigma.link)$mu.eta(S%*%tau)))

  ## Necessary quantities
  vdv <- v.function(y, mu, sigma, lambda, zeta, family)
  xi <- xi.function(mu, sigma, lambda, zeta, family)$xi
  z <- vdv$z
  v <- vdv$v
  dv <- vdv$dv

  mu_star <- c(-lambda / mu + (z * sigma * lambda + 1) * z * v / (mu * sigma))
  if (lambda == 0) {
    dzdlambda <- log(y / mu)^2 / (2 * sigma)
    sigma_star <- c(-1 / sigma + v * z^2 / sigma)
    lambda_star <- as.vector(log(y / mu) - z * v * dzdlambda)
  } else {
    dzdlambda <- (1 / (sigma * lambda^2)) * (((y / mu)^lambda) * (lambda * log(y / mu) - 1) + 1)
    sigma_star <- c((-1 / sigma) + (v * z^2 / sigma) + (xi / (abs(lambda) * sigma^2)))
    lambda_star <- c(log(y / mu) - z * v * dzdlambda + sign(lambda) * xi / (sigma * lambda^2))
  }


  ## Local influence
  if (is.null(alpha)) {
    Delta <- rbind(t(X)%*%T1%*%diag(mu_star), t(S)%*%T2%*%diag(sigma_star), lambda_star)
    case.weights <- abs(eigen(-t(Delta)%*%solve(J)%*%Delta)$vec[,1])
  } else {
    alpha.link <- object$link$alpha
    kappa <- object$coefficients$alpha
    Z <- if(is.null(object$x$alpha)) stats::model.matrix(object, model = "alpha") else object$x$alpha
    ind <- as.numeric(y == 0)
    A <- diag(ind - alpha)
    T0 <- diag(c(stats::make.link(alpha.link)$mu.eta(Z%*%kappa)))
    T1[which(y == 0), ] <- 0
    T2[which(y == 0), ] <- 0
    alpha_dagger <- 1 / (alpha * (1 - alpha))
    mu_dagger <- mu_star
    mu_dagger[which(y == 0)] <- 0L
    sigma_dagger <- sigma_star
    sigma_dagger[which(y == 0)] <- 0L
    lambda_dagger <- lambda_star
    lambda_dagger[which(y == 0)] <- 0L

    Delta <- rbind(t(Z)%*%A%*%T0%*%diag(alpha_dagger),
                   t(X)%*%T1%*%diag(mu_dagger),
                   t(S)%*%T2%*%diag(sigma_dagger),
                   lambda_dagger)
    case.weights <- abs(eigen(-t(Delta)%*%solve(J)%*%Delta)$vec[,1])
  }



  ## Total LI
  totalLI <- rep(NA, n)
  for(i in 1:n){
    totalLI[i] <- 2 * abs(t(Delta[, i])%*%solve(-J)%*%Delta[, i])
  }
  totalLI <- abs((totalLI - mean(totalLI)) / stats::sd(totalLI))

  if(plot == TRUE){

    if (ask) {
      op <- graphics::par(ask = TRUE)
      on.exit(graphics::par(op))
    }

    plot(case.weights, type = "h", main = "Case-weight perturbation", las = 1, xlab = "Index", ylab = "Local influence")
    plot(totalLI, type = "h", main = "Case-weight perturbation", las = 1, xlab = "Index", ylab = "Total local influence")

  }

  list(case.weights = case.weights, totalLI = totalLI)
}


#' Normal Probability Plots with Simulated Envelope for a Box-Cox Symmetric Regression Fit
#'
#' Produce the normal probability plot with simulated envelope of the quantile
#'     residuals obtained from a Box-Cox symmetric regression fit.
#'
#' @param object a fitted model object of class "\code{BCSreg}".
#' @param rep a positive integer representing the number of iterations to calculate
#'     the simulated envelopes. Default is \code{rep = 60}.
#' @param conf a numeric value in the interval (0,1) that represents the confidence
#'     level of the simulated envelope. Default is \code{conf = 0.95}.
#' @param envcol character specifying the color of the envelope.
#' @param ... additional graphical parameters (see par).
#'
#' @details The \code{envelope} uses the idea of Atkinson (1985) to create normal
#'     probability plots with simulated envelope. Under the correct model,
#'     approximately 100 * \code{conf} of the residuals are expected to be inside
#'     the envelope.
#'
#'     Currently, the \code{envelope()} function, when used in a zero-adjusted Box-Cox
#'     symmetric regression fit, returns only the quantile plot for the quantile
#'     residuals obtained under a combined approach (see \link{residuals.BCSreg}).
#'
#' @return \code{envelope} returns normal probability plot with simulated envelopes
#'     for the quantile residuals.
#'
#' @export
#'
#' @references
#'   Atkinson, A. C. (1985). \emph{Plots, Transformations and Regression: An Introduction
#'      to Graphical Methods of Diagnostic Regression Analysis}.
#'      Oxford Science Publications, Oxford.
#'
#'   Medeiros, R. M. R., and Queiroz, F. F. (2025). Flexible modeling of non-negative continuous
#'      data: Box-Cox symmetric regression and its zero-adjusted extension.
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats median model.frame model.response qqnorm quantile residuals as.formula
#' @seealso \code{\link{BCSreg}}, \code{\link{residuals.BCSreg}}
#'
#' @examples
#' ## Data set: raycatch (for description, run ?raycatch)
#' hist(raycatch$cpue, xlab = "Catch per unit effort")
#' plot(cpue ~ tide_phase, raycatch, pch = 16,
#'    xlab = "Tide phase", ylab = "Catch per unit effort")
#' plot(cpue ~ location, raycatch, pch = 16,
#'    xlab = "Location", ylab = "Catch per unit effort")
#' plot(cpue ~ max_temp, raycatch, pch = 16,
#'    xlab = "Maximum temperature", ylab = "Catch per unit effort")
#'
#' ## Fit a double Box-Cox normal regression model:
#' fit <- BCSreg(cpue ~ location + tide_phase |
#'                location + tide_phase + max_temp, raycatch)
#' envelope(fit)
envelope <- function(object, rep = 60, conf = 0.95, envcol, ...){

  dots <- list(...)
  rep <- max(30, floor(rep))

  if(rep != floor(rep) | rep <= 0) stop("The rep argument must be a positive integer.", call. = FALSE)
  if(conf <= 0  | conf >= 1) stop("The conf argument must be within the interval (0, 1).", call. = FALSE)

  family  <- object$family
  n       <- object$nobs
  muhat     <- object$mu
  sigmahat  <- object$sigma
  lambdahat <- object$lambda
  alphahat <- object$alpha
  zetahat <- object$zeta

  alpha_id <- !is.null(alphahat)
  lambda_id <- !is.null(lambdahat)
  zeta_id <- !is.null(zetahat)

  mf <- model.frame(object)
  y <- as.numeric(model.response(mf))
  ind <- as.numeric(y == 0)
  formula <- Formula::as.Formula(object$formula)
  link <- object$link$mu
  sigma.link <- object$link$sigma
  alpha.link <- object$link$alpha

  resRid <- residuals(object)

  resid_env <- matrix(0, n, rep)

  i <- 1
  bar <- txtProgressBar(min = 0, max = rep, initial = 0, width = 50, char = "+", style = 3)
  while(i <= rep){
    tryCatch({

      if (alpha_id) {

        y_env <- rZABCS(n, alphahat, muhat, sigmahat, lambdahat, zetahat, family = family)
        data_env <- cbind(mf, y_env = y_env)
        formula_env <- as.formula(paste("y_env ~",  paste(deparse(formula[[3]]), collapse = "")))

        val <- suppressWarnings(BCSreg(formula_env, data = data_env,
                                       family = family,
                                       zeta = zetahat, link = link,
                                       sigma.link = sigma.link,
                                       alpha.link = alpha.link,
                                       control = object$control))
      } else {

        y_env <- rBCS(n, muhat, sigmahat, lambdahat, zetahat, family = family)
        data_env <- cbind(mf, y_env = y_env)
        formula_env <- as.formula(paste("y_env ~",  paste(deparse(formula[[3]]), collapse = "")))

        val <- suppressWarnings(BCSreg(formula_env, data = data_env,
                                       family = family,
                                       zeta = zetahat, link = link,
                                       sigma.link = sigma.link,
                                       alpha.link = alpha.link,
                                       control = object$control))
      }

      resid_env[,i] <- sort(residuals(val))
      setTxtProgressBar(bar,i)
      i = i + 1
    }, error=function(e){cat("ERROR: ",conditionMessage(e), "\n")})
  }

  if(missing(envcol)) envcol <- "black"
  close(bar)
  cat("\n")

    qqnormInt <- function(y, resid_env, IDENTIFY = TRUE){

      liml <- apply(resid_env, 1, quantile, prob = (1 - conf)/2)
      limu <- apply(resid_env, 1, quantile, prob = (1 - (1 - conf)/2))
      mm   <- apply(resid_env, 1, median)

      if(is.null(dots$xlab)) xlab <- "Theoretical quantiles"
      if(is.null(dots$ylab)) ylab <- "Sample quantiles"
      if(is.null(dots$main)) main <- ""
      if(is.null(dots$ylim)) ylim <- range(resRid , liml , limu)

      qqnorm(y, pch = "+", las = 1, ylim = ylim, xlab = xlab, ylab = ylab, main = main, cex = 0.8,
             lwd = 3, las = 1, ...) -> X
      par(new = TRUE)
      qqnorm(liml, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, lty = 1, main = "",
             col = envcol)
      par(new = TRUE)
      qqnorm(limu, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, lty = 1, main = "",
             col = envcol)
      par(new = TRUE)
      qqnorm(mm, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, lty = 2, main = main,
             col = envcol)

      invisible(X)
      cat(paste("\nClick on points you would like to identify and press Esc."), "\n")
      if(IDENTIFY) return(identify(X, cex = 0.8))

    }

    qqnormInt(resRid, resid_env)

}




