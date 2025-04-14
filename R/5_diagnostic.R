#' Influence Diagnostics for BCSreg Objects
#'
#' The \code{influence} function provides two influence measures for Box-Cox symmetric regression
#'     models.
#'
#' @param object an object of class \code{"BCSreg"}.
#' @param plot logical. If \code{plot = TRUE} (default), the plots are shown.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot, if \code{plot = TRUE}.
#' @param ... currently not used.
#'
#' @return \code{influence} returns a list with three objects:
#'     \item{case.weights}{The values of \eqn{d_{max}} eigenvector based on case
#'     weights perturbation scheme (see Medeiros and Queiroz (2025)).}
#'     \item{totalLI}{The total local influence (see Lesaffre and Verbeke (1998)).}
#'
#' @seealso \code{\link{BCSreg}}, \code{\link{plot.BCSreg}}
#'
#' @references
#'  Lesaffre, E., and Verbeke, G. (1998). Local influence in linear mixed models.
#'  \emph{Biometrics}, 570--582.
#'
#'  Medeiros, R. M. R., and Queiroz, F. F. (2025). Modeling positive continuous data:
#'     Box-Cox symmetric regression models and their extensions
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' ## Examples
#'
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

  ## Generalized leverage
  # L1 <- diag(y^(lambda - 1) * ((2 * z * sigma * lambda + 1) * v + (z * sigma * lambda + 1) * z * dv) / (sigma^2 * mu^(lambda + 1)))
  # L2 <- diag(z * y^(lambda - 1) * (2 * v + z * dv) / (sigma^2 * mu^lambda))
  # L3 <- 1 / y - (y / mu)^lambda * (v + z * dv * dzdlambda + z * v * log(y / mu)) / (sigma * y)
  # Ldot <- cbind(T1%*%X, matrix(0, nrow = n, ncol = q + 1))
  # Lddot <- rbind(t(X)%*%T1%*%L1, t(S)%*%T2%*%L2, L3)
  # GL <- Ldot%*%solve(J)%*%Lddot

  # zp <- rep(NA, n)
  # l0 <- which(rep(lambda, n) <= 0)
  # l1 <- which(rep(lambda, n) > 0)
  # switch(family,
  #        NO = {
  #          zp[l0] <- qnorm(0.5 * pnorm(1 / (sigma[l0] * abs(lambda))))
  #          zp[l1] <- qnorm(1 - 0.5 * pnorm(1 / (sigma[l1] * abs(lambda))))
  #          dzp_dsigma <- 0.5 * dnorm(sqrt(1 / (sigma[l1] * abs(lambda)))) /
  #            (sigma[l1]^2 * lambda * dnorm(sqrt(zp)))
  #          dzp_dlambda <- 0.5 * dnorm(sqrt(1 / (sigma[l1] * abs(lambda)))) /
  #            (sigma[l1] * lambda^2 * dnorm(sqrt(zp)))
  #        },
  #        ST = {
  #          zp[l0] <- qt(0.5 * pt(1 / (sigma[l0] * abs(lambda)), zeta), zeta)
  #          zp[l1] <- qt(1 - 0.5 * pt(1 / (sigma * abs(lambda)), zeta), zeta)
  #          dzp_dsigma <- 0.5 * dt(sqrt(1 / (sigma[l1] * abs(lambda))), zeta) /
  #            (sigma[l1]^2 * lambda * dt(sqrt(zp), zeta))
  #          dzp_dlambda <- 0.5 * dt(sqrt(1 / (sigma[l1] * abs(lambda))), zeta) /
  #            (sigma[l1] * lambda^2 * dt(sqrt(zp), zeta))
  #        },
  #        LOI = {
  #          zp[l0] <- qlogisI(0.5 * plogisI(1 / (sigma[l0] * abs(lambda))))
  #          zp[l1] <- qlogisI(1 - 0.5 * plogisI(1 / (sigma[l1] * abs(lambda))))
  #          dzp_dsigma <- 0.5 * dlogisI(sqrt(1 / (sigma[l1] * abs(lambda)))) /
  #            (sigma[l1]^2 * lambda * dlogisI(sqrt(zp)))
  #          dzp_dlambda <- 0.5 * dlogisI(sqrt(1 / (sigma[l1] * abs(lambda)))) /
  #            (sigma[l1] * lambda^2 * dlogisI(sqrt(zp)))
  #        },
  #        LOII = {
  #          zp[l0] <- qlogis(0.5 * plogis(1 / (sigma[l0] * abs(lambda))))
  #          zp[l1] <- qlogis(1 - 0.5 * plogis(1 / (sigma[l1] * abs(lambda))))
  #          dzp_dsigma <- 0.5 * dlogis(sqrt(1 / (sigma[l1] * abs(lambda)))) /
  #            (sigma[l1]^2 * lambda * dlogis(sqrt(zp)))
  #          dzp_dlambda <- 0.5 * dlogis(sqrt(1 / (sigma[l1] * abs(lambda)))) /
  #            (sigma[l1] * lambda^2 * dlogis(sqrt(zp)))
  #        },
  #        SN = {
  #          zp[l0] <- asinh((zeta / 2) * qnorm(0.5 * pnorm((2 / zeta) * sinh(1 / (sigma[l0] * abs(lambda))))))
  #          zp[l1] <- asinh((zeta / 2) * qnorm(1 - 0.5 * pnorm((2 / zeta) * sinh(1 / (sigma[l1] * abs(lambda))))))
  #        },
  #        PE = {
  #          zp[l0] <- gamlss.dist::qPE(0.5 * gamlss.dist::pPE(1 / (sigma[l0] * abs(lambda)),
  #                                                            mu = 0, sigma = 1, nu = zeta
  #          ), mu = 0, sigma = 1, nu = zeta)
  #          zp[l1] <- gamlss.dist::qPE(1 - 0.5 * gamlss.dist::pPE(1 / (sigma[l1] * abs(lambda)),
  #                                                                mu = 0, sigma = 1, nu = zeta
  #          ), mu = 0, sigma = 1, nu = zeta)
  #          dzp_dsigma <- 0.5 * gamlss.dist::dPE(sqrt(1 / (sigma[l1] * abs(lambda))), nu = zeta) /
  #            (sigma[l1]^2 * lambda * gamlss.dist::dPE(sqrt(zp), nu = zeta))
  #          dzp_dlambda <- 0.5 * gamlss.dist::dPE(sqrt(1 / (sigma[l1] * abs(lambda))), nu = zeta) /
  #            (sigma[l1] * lambda^2 * gamlss.dist::dPE(sqrt(zp), nu = zeta))
  #        },
  #        HP = {
  #          zp[l0] <- GeneralizedHyperbolic::qhyperb(
  #            0.5 * GeneralizedHyperbolic::phyperb(1 / (sigma[l0] * abs(lambda)),
  #                                                 mu = 0, delta = 1, alpha = zeta, beta = 0
  #            ),
  #            mu = 0, delta = 1, alpha = zeta, beta = 0
  #          )
  #          zp[l1] <- GeneralizedHyperbolic::qhyperb(
  #            1 - 0.5 * GeneralizedHyperbolic::phyperb(1 / (sigma[l1] * abs(lambda)),
  #                                                     mu = 0, delta = 1, alpha = zeta, beta = 0
  #            ),
  #            mu = 0, delta = 1, alpha = zeta, beta = 0
  #          )
  #          dzp_dsigma <- 0.5 * GeneralizedHyperbolic::dhyperb(sqrt(1 / (sigma[l1] * abs(lambda))), alpha = zeta) /
  #            (sigma[l1]^2 * lambda * GeneralizedHyperbolic::dhyperb(sqrt(zp), alpha = zeta))
  #          dzp_dlambda <- 0.5 * GeneralizedHyperbolic::dhyperb(sqrt(1 / (sigma[l1] * abs(lambda))), alpha = zeta) /
  #            (sigma[l1] * lambda^2 * GeneralizedHyperbolic::dhyperb(sqrt(zp), alpha = zeta))
  #        },
  #        SL = {
  #          zp[l0] <- qslash(0.5 * pslash(1 / (sigma[l0] * abs(lambda)), zeta = zeta), zeta = zeta)
  #          zp[l1] <- qslash(1 - 0.5 * pslash(1 / (sigma[l1] * abs(lambda)), zeta = zeta), zeta = zeta)
  #          dzp_dsigma <- 0.5 * dslash(sqrt(1 / (sigma[l1] * abs(lambda))), zeta = zeta) /
  #            (sigma[l1]^2 * lambda * dslash(sqrt(zp), zeta = zeta))
  #          dzp_dlambda <- 0.5 * dslash(sqrt(1 / (sigma[l1] * abs(lambda))), zeta = zeta) /
  #            (sigma[l1] * lambda^2 * dslash(sqrt(zp), zeta = zeta))
  #        },
  #        stop(gettextf("%s family not recognised", sQuote(family)), domain = NA)
  # )
  #
  # # Ldot
  # id <- which(rep(lambda, n) != 0)
  #
  # dy_dbeta <- T1%*%X
  # dy_dbeta[id, ] <- T1[id, ]%*%diag((1 + sigma[id] * lambda * zp[id])^(1 / lambda))%*%X[id, ]
  #
  # dy_dtau <- matrix(0, nrow = n, ncol = q)
  # dy_dtau[id, ] <- T2[id, ]%*%diag(mu[id] * (1 + sigma[id] * lambda * zp[id])^(1 / (lambda - 1)) *
  #                                    (zp[id] + sigma[id] * dzp_dsigma))%*%S[id, ]
  #
  #
  # dy_dlambda <- matrix(0, nrow = n, ncol = 1)
  # dy_dlambda[id, ] <- mu[id] * (1 + sigma[id] * lambda * zp[id])^(1 / lambda) *
  #   (- log(1 + sigma[id] * lambda * zp[id]) / lambda^2 + sigma[id] * (zp[id] + lambda * dzp_dlambda) /
  #      (lambda * (1 + sigma[id] * lambda * zp[id])))
  #
  # Ldot <- cbind(dy_dbeta, dy_dtau, dy_dlambda)
  #
  # # Lddot
  # L1 <- diag(y^(lambda - 1) * ((2 * z * sigma * lambda + 1) * v + (z * sigma * lambda + 1) * z * dv) / (sigma^2 * mu^(lambda + 1)))
  # L2 <- diag(z * y^(lambda - 1) * (2 * v + z * dv) / (sigma^2 * mu^lambda))
  # L3 <- 1 / y - (y / mu)^lambda * (v + z * dv * dzdlambda + z * v * log(y / mu)) / (sigma * y)
  # Lddot <- rbind(t(X)%*%T1%*%L1, t(S)%*%T2%*%L2, L3)
  #
  # GL <- Ldot%*%solve(J)%*%Lddot

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
  #totalLI <- abs((totalLI - mean(totalLI)) / sd(totalLI))

  if(plot == TRUE){

    if (ask) {
      op <- graphics::par(ask = TRUE)
      on.exit(graphics::par(op))
    }

    plot(case.weights, type = "h", main = "Case-weight perturbation", las = 1, xlab = "Index", ylab = "Local influence")
    # cat(paste("\nClick on points you would like to identify and press Esc."))
    # identify(seq_len(n), case.weights, cex = 0.8)

    plot(totalLI, type = "h", main = "Case-weight perturbation", las = 1, xlab = "Index", ylab = "Total local influence")
    # identify(seq_len(n), totalLI, cex = 0.8)

    #plot(diag(GL), pch = "+", las = 1, cex = 0.8, main = "Generalized leverage", xlab = "Index", ylab = expression(GL[ii]))
    # identify(seq_len(n), diag(GL), cex = 0.8)
  }

  list(case.weights = case.weights, totalLI = totalLI)#, GL = diag(GL))
}
