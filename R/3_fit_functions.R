# Control the optimization process ------------------------------------------------------------

#' Auxiliary for Controlling BCS Fitting
#'
#' Optimization parameters that control the fitting of Box-Cox symmetric regression
#'     models using the \code{\link{BCSreg}} function.
#'
#' @name BCSregcontrol
#' @param lambda numeric indicating the value of lambda (if \code{NULL}, lambda
#'      will be estimated).
#' @param method character specifying the \code{method} argument passed to \code{\link{optim}}.
#' @param hessian logical. Should the numerical Hessian matrix from the \code{optim} output be
#'      used for estimation of the covariance matrix? By default the analytical solution is employed.
#' @param maxit,trace,... arguments passed to \code{\link{optim}}.
#' @param start an optional vector with starting values for the regression coefficients
#'      associated with the \code{mu} and \code{sigma} submodels (starting value for the
#'      \code{lambda} parameter must not be included).
#'
#' @details The \code{BCSreg.control} controls the fitting process of Box-Cox symmetric models.
#'     Almost all the arguments are passed directly to \code{\link{optim}}, which
#'     is used to estimate the parameters. Starting values for the regression coefficients
#'     associated with the \code{mu} and \code{sigma} submodels can be provided via \code{start}
#'     argument. The starting value for the skewness parameter \code{lambda} is zero; that is,
#'     the fitting process starts from the corresponding log-symmetric regression model. If the
#'     estimation process is to be performed with a fixed \code{lambda}, a value must be specified
#'     for the \code{lambda} argument. A natural value is \code{lambda = 0}, where a log-symmetric
#'     regression model will be estimated.
#'
#' @return A list with components named as the arguments.
#' @seealso \code{\link{BCSreg}}
#' @export
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#' @examples
#' # Set an example
#'
#'
BCSreg.control <- function(lambda = NULL, method = "BFGS", maxit = 2000, hessian = FALSE,
                           trace = FALSE, start = NULL, ...) {
  val <- list(
    lambda = lambda, method = method, maxit = maxit, hessian = hessian,
    trace = trace, start = start
  )

  val <- c(val, list(...))

  if (!is.null(val$fnscale)) {
    warning("fnscale must not be modified")
  }

  val$fnscale <- -1

  if (is.null(val$reltol)) {
    val$reltol <- .Machine$double.eps^(1 / 2)
  }

  val
}

#' @keywords internal
make.dmu.deta <- function(linkstr) {
  switch(linkstr,
         "log" = function(eta) pmax(exp(eta), .Machine$double.eps),
         "sqrt" = function(eta) rep.int(2, length(eta)),
         "1/mu^2" = function(eta) 3 / (4 * eta^2.5),
         "inverse" = function(eta) 2 / (eta^3),
         "identity" = function(eta) rep.int(0, length(eta))
  )
}


# Auxiliary special functions -----------------------------------------------------------------

## v() function
v.function <- function(y, mu, sigma, lambda, zeta, family) {
  n <- length(y)
  if (lambda == 0) {
    z <- log(y / mu) / sigma
  } else {
    z <- ((y / mu)^lambda - 1) / (sigma * lambda)
  }

  if (family == "NO") {
    vz <- rep(1, n)
    dvz <- rep(0, n)
  }
  if (family == "ST") {
    vz <- (zeta + 1) / (zeta + z^2)
    dvz <- -2 * (zeta + 1) * z / ((zeta + z^2)^2)
  }
  if (family == "LOI") {
    vz <- -2 * (exp(-z^2) - 1) / (1 + exp(-z^2))
    dvz <- 8 * z * exp(-z^2) / ((1 + exp(-z^2))^2)
  }
  if (family == "LOII") {
    vz <- (1 - exp(-abs(z))) / (abs(z) * (1 + exp(-abs(z))))
    dvz <- (2 * abs(z) * exp(-abs(z)) + exp(-2 * abs(z)) - 1) / (((z^2)^3 / 2) * ((1 + exp(-abs(z)))^2))
  }
  if (family == "SN") {
    vz <- 4 * sinh(z) * cosh(z) / (zeta^2 * z) - tanh(z) / z
    dvz <- ((cosh(z) * z - sinh(z)) / z^2) * (4 * cosh(z) / zeta^2 - 1 / cosh(z)) +
      (sinh(z) / z) * (4 * sinh(z) / zeta^2 + sinh(z) / (cosh(z)^2))
  }
  if (family == "PE") {
    pzeta <- sqrt(2^(-2 / zeta) * gamma(1 / zeta) * (gamma(3 / zeta)^(-1)))
    vz <- (zeta * (z^2)^(zeta / 2 - 1)) / (2 * pzeta^zeta)
    dvz <- ((zeta^2 / 2 - zeta) * (pzeta^(-zeta)) * (z^2)^(zeta / 2)) / (z^3)
  }
  if (family == "HP") {
    vz <- zeta / sqrt(1 + z^2)
    dvz <- -(z * zeta) / (1 + z^2)^(3 / 2)
  }
  if (family == "SL") {
    s_aux <- z^2 / 2
    beta_aux <- zeta + (1 / 2)
    gama_aux <- ig(beta_aux, s_aux)
    vz <- (2 / z^2) * ig(zeta + (3 / 2), s_aux) / gama_aux
    dvz <- (-2 / z) * vz + (2 / z) * (1 / gama_aux^2) * exp(-s_aux) * (s_aux^beta_aux) * (gama_aux * (1 - beta_aux / s_aux)
                                                                                          + exp(-s_aux) * (s_aux^(beta_aux - 1)))
  }

  list(z = z, v = vz, dv = dvz)
}

## xi function
xi.function <- function(mu, sigma, lambda, zeta, family) {
  n <- max(length(mu), length(sigma))
  if (any(mu <= 0)) {
    stop(paste("mu must be positive ", "\n", ""))
  }
  if (any(sigma <= 0)) {
    stop(paste("sigma must be positive", "\n", ""))
  }
  if (any(zeta <= 0)) {
    stop(paste("zeta must be positive ", "\n", ""))
  }

  if (lambda == 0) {
    xi <- 0
    dxidsigma <- 0
    dxidlambda <- 0
  } else {
    if (family == "NO") {
      vz <- rep(1, n)
      r <- dnorm(1 / ((sigma * lambda)))
      R <- pnorm(1 / (sigma * abs(lambda)))
    }
    if (family == "ST") {
      vz <- (zeta + 1) / (zeta + (1 / ((sigma * lambda)))^2)
      r <- dt(1 / ((sigma * lambda)), zeta)
      R <- pt(1 / (sigma * abs(lambda)), zeta)
    }
    if (family == "LOI") {
      vz <- -2 * (exp(-(1 / ((sigma * lambda)))^2) - 1) / (1 + exp(-(1 / ((sigma * lambda)))^2))
      r <- dlogisI(1 / ((sigma * lambda)))
      R <- plogisI(1 / (sigma * abs(lambda)))
    }
    if (family == "LOII") {
      vz <- (1 - exp(-abs((1 / ((sigma * lambda)))))) / (abs((1 / ((sigma * lambda)))) * (1 + exp(-abs((1 / ((sigma * lambda)))))))
      r <- dlogis(1 / ((sigma * lambda)))
      R <- plogis(1 / (sigma * abs(lambda)))
    }
    if (family == "SN") {
      vz <- 4 * sinh((1 / ((sigma * lambda)))) * cosh((1 / ((sigma * lambda)))) / (zeta^2 * (1 / ((sigma * lambda)))) - tanh((1 / ((sigma * lambda)))) / (1 / ((sigma * lambda)))
      r <- (2 * cosh(sqrt(1 / ((sigma * lambda))^2)) * exp(-2 * sinh(sqrt(1 / ((sigma * lambda))^2)) * sinh(sqrt(1 / ((sigma * lambda))^2)) / zeta^2) / (sqrt(2 * pi) * zeta))
      R <- pnorm((2 / zeta) * sinh(1 / (sigma * abs(lambda))))
    }
    if (family == "PE") {
      pzeta <- sqrt(2^(-2 / zeta) * gamma(1 / zeta) * (gamma(3 / zeta)^(-1)))
      vz <- (zeta * ((1 / ((sigma * lambda)))^2)^(zeta / 2 - 1)) / (2 * pzeta^zeta)
      r <- gamlss.dist::dPE(1 / ((sigma * lambda)), mu = 0, sigma = 1, nu = zeta)
      R <- gamlss.dist::pPE(1 / (sigma * abs(lambda)), mu = 0, sigma = 1, nu = zeta)
    }
    if (family == "HP") {
      vz <- zeta / sqrt(1 + (1 / ((sigma * lambda)))^2)
      r <- GeneralizedHyperbolic::dhyperb(1 / ((sigma * lambda)), mu = 0, delta = 1, alpha = zeta, beta = 0)
      R <- GeneralizedHyperbolic::phyperb(1 / (sigma * abs(lambda)),
                                          mu = 0, delta = 1,
                                          alpha = zeta, beta = 0
      )
    }
    if (family == "SL") {
      s_aux <- (1 / ((sigma * lambda)))^2 / 2
      beta_aux <- zeta + (1 / 2)
      if (any(s_aux == 0)) {
        s_aux[s_aux == 0] <- 0.0001
      }
      gama_aux <- ig(beta_aux, s_aux)
      vz <- (2 / (1 / ((sigma * lambda)))^2) * ig(zeta + (3 / 2), s_aux) / gama_aux
      r <- dslash(1 / ((sigma * lambda)), zeta = zeta)
      R <- pslash(1 / (sigma * abs(lambda)), zeta = zeta)
    }
    xi <- r / R
    dxidsigma <- (xi^2 / (abs(lambda) * sigma^2)) - (2 / ((lambda^2) * (sigma^3))) * (1 / R) * (-r * vz / 2)
    dxidlambda <- (lambda * (xi^2) / (sigma * (abs(lambda)^3))) - (2 / ((lambda^3) * (sigma^2) * R)) * (-r * vz / 2)
  }

  list(xi = xi, dxidsigma = dxidsigma, dxidlambda = dxidlambda)
}


# Work-horse fitting function (not exported) -----------------------------------------------------------------
BCSreg.fit <- function(X, y, S = NULL, family, zeta = NULL, link = "log",
                       sigma.link = "log", control = BCSreg.control()) {
  n <- length(y)
  p <- dim(X)[2]

  ## Model matrices
  # if (is.null(colnames(X))) {
  #   if (p == 1) {
  #     colnames(X) <- "(Intercept)"
  #   } else {
  #     colnames(X) <- c("(Intercept)", paste0("X", 1:(p-1)))
  #   }
  # }

  if (is.null(S)) {
    # q <- 1
    S <- matrix(1, nrow = n)
    #colnames(S) <- "(Intercept)"
    #rownames(S) <- rownames(X)
  } #else {
  #   q <- dim(S)[2]
  #   if (is.null(colnames(S))) {
  #     if (q == 1) {
  #       colnames(S) <- "(Intercept)"
  #     } else {
  #       colnames(S) <- c("(Intercept)", paste0("S", 1:(p-1)))
  #     }
  #   }
  # }
  q <- ncol(S)
  #sigma_const <- (q == 1) && (sigma.link == "identity")

  ## Link functions
  if (is.character(link)) {
    linkstr <- link
    linkobj <- stats::make.link(linkstr)
    linkobj$dmu.deta <- make.dmu.deta(linkstr)
  } else {
    linkobj <- link
    linkstr <- link$name
    if (is.null(linkobj$dmu.deta)) {
      warning("link needs to provide dmu.deta component", call. = FALSE)
    }
  }

  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  dmu.deta <- linkobj$dmu.deta

  if (is.character(sigma.link)) {
    sigma_linkstr <- sigma.link
    sigma_linkobj <- stats::make.link(sigma_linkstr)
    sigma_linkobj$dmu.deta <- make.dmu.deta(sigma_linkstr)
  } else {
    sigma_linkobj <- sigma.link
    sigma_linkstr <- sigma.link$name
    if (is.null(sigma_linkobj$dmu.deta)) {
      warning("sigma.link needs to provide dmu.deta component.", call. = FALSE)
    }
  }

  sigma_linkfun <- sigma_linkobj$linkfun
  sigma_linkinv <- sigma_linkobj$linkinv
  sigma_mu.eta <- sigma_linkobj$mu.eta
  sigma_dmu.deta <- sigma_linkobj$dmu.deta

  ## Optimization control parameters
  ocontrol <- control
  lambda_id <- is.null(control$lambda) # TRUE means that lambda will be estimated
  zeta_id <- !is.null(zeta) # TRUE means that the family has an extra parameter
  lambda_fix <- control$lambda # NULL means that lambda will be estimated
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  control$lambda <- control$method <- control$hessian <- control$start <- NULL

  ## Starting values
  if (is.null(start)) {
    beta <- solve(t(X) %*% X) %*% t(X) %*% log(y)
    tau <- rep.int(0, q)

    # sigma[1L] <- sd(logy)
    CVy <- 0.75 * diff(stats::quantile(y, c(0.25, 0.75))) / stats::median(y)
    tau[1L] <- asinh(CVy / 1.5) / stats::qnorm(0.75)
    if (!isTRUE(sigma_linkinv(tau[1]) > 0)) {
      warning("No valid starting value for dispersion parameter found, using 1 instead", call. = FALSE)
      tau[1L] <- 1
    }

    start <- list(mu = beta, sigma = tau)
  }

  if (is.list(start)) start <- do.call("c", start)
  start <- c(start, if (lambda_id) 0L else NULL)

  ## Log-likelihood function
  logL <- function(theta, lambda_id) {
    beta <- theta[seq.int(length.out = p)]
    tau <- theta[seq.int(length.out = q) + p]

    eta.1 <- c(X %*% beta)
    eta.2 <- c(S %*% tau)

    mu <- linkinv(eta.1)
    sigma <- sigma_linkinv(eta.2)

    if (lambda_id) {
      lambda <- theta[1L + p + q]
    } else {
      lambda <- lambda_fix
    }

    if (lambda == 0) {
      z <- log(y / mu) / sigma
    } else {
      z <- ((y / mu)^lambda - 1) / (sigma * lambda)
    }

    if (any(!is.finite(z))) {
      NaN
    } else {
      ll <- suppressWarnings(dBCS(y,
                                  mu = mu, sigma = sigma, lambda = lambda,
                                  zeta = zeta, family = family, log = TRUE
      ))
      if (any(!is.finite(ll))) {
        NaN
      } else {
        sum(ll)
      }
    }
  }

  # Não precisa definir duas log-verossimilhanças, basta adaptar uma delas
  # para receber o lambda fixado

  # logLfixed <- function(theta, lambda) {
  #   beta <- theta[seq.int(length.out = p)]
  #   tau <- theta[seq.int(length.out = q) + p]
  #
  #   eta.1 <- as.vector(X %*% beta)
  #   eta.2 <- as.vector(S %*% tau)
  #
  #   mu <- linkinv(eta.1)
  #   sigma <- sigma_linkinv(eta.2)
  #
  #   if (lambda == 0) {
  #     z <- (1 / sigma) * (log(y / mu))
  #   } else {
  #     z <- (1 / (sigma * lambda)) * ((y / mu)^lambda - 1)
  #   }
  #
  #   if (any(!is.finite(z))) {
  #     NaN
  #   } else {
  #     ll <- suppressWarnings(dBCS(y,
  #       mu = mu, sigma = sigma, lambda = lambda,
  #       zeta = zeta, family = family, log = TRUE
  #     ))
  #     if (any(!is.finite(ll))) {
  #       NaN
  #     } else {
  #       sum(ll)
  #     }
  #   }
  # }

  ## Score function
  U <- function(theta, lambda_id) {
    beta <- theta[seq.int(length.out = p)]
    tau <- theta[seq.int(length.out = q) + p]

    eta.1 <- as.vector(X %*% beta)
    eta.2 <- as.vector(S %*% tau)

    mu <- linkinv(eta.1)
    sigma <- sigma_linkinv(eta.2)

    if (lambda_id) {
      lambda <- theta[1L + p + q]
    } else {
      lambda <- lambda_fix
    }

    vdv <- v.function(y, mu, sigma, lambda, zeta, family)
    v <- vdv$v
    z <- vdv$z
    xi.f <- xi.function(mu, sigma, lambda, zeta, family)
    xi <- xi.f$xi

    d1dot <- as.vector(1 / mu.eta(eta.1))
    d2dot <- as.vector(1 / sigma_mu.eta(eta.2))

    T1 <- diag(1 / d1dot)
    T2 <- diag(1 / d2dot)

    mu_star <- as.vector(-(lambda / mu) + (1 / (mu * sigma)) * (z * sigma * lambda + 1) * z * v)
    if (lambda == 0) {
      dzdlambda <- (1 / (2 * sigma)) * (log(y / mu)^2)
      sigma_star <- as.vector((-1 / sigma) + (v * z^2 / sigma))
      lambda_star <- as.vector(log(y / mu) - z * v * dzdlambda)
    } else {
      dzdlambda <- (1 / (sigma * lambda^2)) * (((y / mu)^lambda) * (lambda * log(y / mu) - 1) + 1)
      sigma_star <- as.vector((-1 / sigma) + (v * z^2 / sigma) + (xi / (abs(lambda) * sigma^2)))
      lambda_star <- as.vector(log(y / mu) - z * v * dzdlambda + sign(lambda) * xi / (sigma * lambda^2))
    }

    if (lambda_id) {
      U <- c(t(X) %*% T1 %*% mu_star, t(S) %*% T2 %*% sigma_star, t(rep.int(1, n)) %*% lambda_star)
    } else {
      U <- c(t(X) %*% T1 %*% mu_star, t(S) %*% T2 %*% sigma_star)
    }


    U
  }

  # Ufixed <- function(theta, lambda) {
  #   beta <- theta[seq.int(length.out = p)]
  #   tau <- theta[seq.int(length.out = q) + p]
  #
  #   eta.1 <- as.vector(X %*% beta)
  #   eta.2 <- as.vector(S %*% tau)
  #
  #   mu <- linkinv(eta.1)
  #   sigma <- sigma_linkinv(eta.2)
  #
  #   vdv <- v.function(y, mu, sigma, lambda, zeta, family)
  #   v <- vdv$v
  #   z <- vdv$z
  #   xi.f <- xi.function(mu, sigma, lambda, zeta, family)
  #   xi <- xi.f$xi
  #
  #   d1dot <- as.vector(1 / mu.eta(eta.1))
  #   d2dot <- as.vector(1 / sigma_mu.eta(eta.2))
  #
  #   T1 <- diag(1 / d1dot)
  #   T2 <- diag(1 / d2dot)
  #
  #   mu_star <- as.vector(-(lambda / mu) + (1 / (mu * sigma)) * (z * sigma * lambda + 1) * z * v)
  #   if (lambda == 0) {
  #     sigma_star <- as.vector((-1 / sigma) + (v * z^2 / sigma))
  #   } else {
  #     sigma_star <- as.vector((-1 / sigma) + (v * z^2 / sigma) + (xi / (abs(lambda) * sigma^2)))
  #   }
  #
  #   Ufixed <- c(t(X) %*% T1 %*% mu_star, t(S) %*% T2 %*% sigma_star)
  #   Ufixed
  # }

  ## Observed information matrix
  Jfunction <- function(beta, tau, lambda) {
    eta.1 <- as.vector(X %*% beta)
    eta.2 <- as.vector(S %*% tau)

    mu <- linkinv(eta.1)
    sigma <- sigma_linkinv(eta.2)

    d1dot <- as.vector(1 / mu.eta(eta.1))
    dd1dot <- as.vector(-dmu.deta(eta.1) * d1dot^3)

    d2dot <- as.vector(1 / sigma_mu.eta(eta.2))
    dd2dot <- as.vector(-sigma_dmu.deta(eta.2) * d2dot^3)

    vdv <- v.function(y, mu, sigma, lambda, zeta, family)
    v <- vdv$v
    z <- vdv$z
    dv <- vdv$dv
    xi.f <- xi.function(mu, sigma, lambda, zeta, family)
    xi <- xi.f$xi
    dxidsigma <- xi.f$dxidsigma
    dxidlambda <- xi.f$dxidlambda

    dzdmu <- -(1 / (sigma * mu)) * (y / mu)^lambda
    dzdsigma <- -z / sigma

    dz2dmu2 <- ((lambda + 1) / (sigma * mu^2)) * (y / mu)^lambda

    if (lambda == 0) {
      dzdlambda <- (1 / (2 * sigma)) * (log(y / mu)^2)
      dz2dlambda2 <- (1 / (3 * sigma)) * (log(y / mu)^3)
      sigma_star <- as.vector((-1 / sigma) + (v * z^2 / sigma))
      lambda_star <- as.vector(log(y / mu) - z * v * dzdlambda)
      dl2dsigma2 <- (1 / sigma^2) * (-(z^3) * dv - 3 * (z^2) * v + 1)
      dl2dlambda2 <- -(z * dv + v) * (dzdlambda^2) - v * z * dz2dlambda2
    } else {
      dzdlambda <- (1 / (sigma * (lambda^2))) * (((y / mu)^lambda) * (lambda * log(y / mu) - 1) + 1)
      dz2dlambda2 <- (1 / (sigma * (lambda^3))) * (-2 + ((y / mu)^lambda) * (2 - 2 * lambda * log(y / mu) + (lambda^2) * (log(y / mu))^2))
      sigma_star <- as.vector((-1 / sigma) + (v * z^2 / sigma) + (xi / (abs(lambda) * sigma^2)))
      lambda_star <- as.vector(log(y / mu) - z * v * dzdlambda + sign(lambda) * xi / (sigma * lambda^2))
      dl2dsigma2 <- (1 / sigma^2) * (-(z^3) * dv - 3 * (z^2) * v + 1) + (dxidsigma / (abs(lambda) * sigma^2)) - (2 * xi / (abs(lambda) * sigma^3))
      dl2dlambda2 <- -(z * dv + v) * (dzdlambda^2) - v * z * dz2dlambda2 +
        sign(lambda) * ((1 / (sigma * lambda^2)) * dxidlambda - 2 * xi / (sigma * lambda^3))
    }
    dldmu <- -(lambda / mu) - v * z * dzdmu
    dl2dmu2 <- (lambda / mu^2) - (z * dv + v) * dzdmu^2 - v * z * dz2dmu2


    w1 <- dl2dmu2 * (1 / d1dot) - dldmu * dd1dot / d1dot^2
    w2 <- dl2dsigma2 * (1 / d2dot) - sigma_star * (dd2dot / d2dot^2)
    w3 <- dl2dlambda2
    w4 <- (z / sigma) * dzdmu * (z * dv + 2 * v)
    w5 <- (1 / (mu * sigma)) * (-sigma + z * (lambda * sigma * z + 1) * dzdlambda * dv +
                                  v * (sigma * z * (2 * lambda * dzdlambda + z) + dzdlambda))
    w6 <- (1 / ((sigma^2) * abs(lambda)^3)) * (sigma * (abs(lambda)^3) * z * dzdlambda * (z * dv + 2 * v) +
                                                 (lambda^2) * dxidlambda - lambda * xi)

    T1 <- diag(1 / d1dot)
    T2 <- diag(1 / d2dot)

    W1 <- diag(as.vector(w1))
    W2 <- diag(as.vector(w2))
    W3 <- diag(as.vector(w3))
    W4 <- diag(as.vector(w4))
    W5 <- diag(as.vector(w5))
    W6 <- diag(as.vector(w6))

    Jbeta.beta <- t(X) %*% W1 %*% T1 %*% X
    Jtau.tau <- t(S) %*% W2 %*% T2 %*% S
    Jlambda.lambda <- t(rep.int(1, n)) %*% W3 %*% rep.int(1, n)
    Jbeta.tau <- t(X) %*% W4 %*% T1 %*% T2 %*% S
    Jbeta.lambda <- t(X) %*% W5 %*% T1 %*% rep.int(1, n)
    Jtau.lambda <- t(S) %*% W6 %*% T2 %*% rep.int(1, n)

    J <- rbind(
      cbind(Jbeta.beta, Jbeta.tau, Jbeta.lambda),
      cbind(t(Jbeta.tau), Jtau.tau, Jtau.lambda),
      cbind(t(Jbeta.lambda), t(Jtau.lambda), Jlambda.lambda)
    )
    -J # negative
  }

  ## Maximum likelihood estimation
  theta.opt <- stats::optim(
    par = start, fn = logL, gr = U, lambda_id = lambda_id,
    method = method, control = control, hessian = hessian
  )

  beta <- theta.opt$par[seq.int(length.out = p)]
  tau <- theta.opt$par[seq.int(length.out = q) + p]
  if (lambda_id) {
    lambda <- theta.opt$par[p + q + 1L]
  } else {
    lambda <- lambda_fix
  }

  ## Covariance matrix
  if (hessian) {
    vcov <- solve(theta.opt$hessian)
  } else {
    if (lambda_id) {
      vcov <- solve(Jfunction(beta, tau, lambda))
    } else {
      vcov <- solve(Jfunction(beta, tau, lambda)[seq.int(length.out = p + q),
                                                 seq.int(length.out = p + q)])
    }

  }

  # if (lambda_id) {
  #   rownames(vcov) <- colnames(vcov) <- c(colnames(X), if (sigma_const) {
  #     "(sigma)"
  #   } else {
  #     paste("(sigma)",
  #           colnames(S),
  #           sep = "_"
  #     )
  #   }, "(lambda)")
  # } else {
  #   rownames(vcov) <- colnames(vcov) <- c(colnames(X), if (sigma_const) {
  #     "(sigma)"
  #   } else {
  #     paste("(sigma)",
  #           colnames(S),
  #           sep = "_"
  #     )
  #   })
  # }

  theta.opt$vcov <- vcov
  theta.opt$start <- start


  # if (lambda_id) {
  #   theta.opt <- optim(
  #     par = c(start, 0), fn = logL, gr = U, method = method,
  #     control = control, hessian = TRUE
  #   )
  #   beta <- theta.opt$par[seq.int(length.out = p)]
  #   tau <- theta.opt$par[seq.int(length.out = q) + p]
  #   lambda <- theta.opt$par[seq.int(length.out = 1) + p + q]
  #   if (theta.opt$convergence == 0) {
  #     converged <- TRUE
  #   } else {
  #     converged <- FALSE
  #     warning("Optimization failed to converge.", call. = FALSE)
  #   }
  #   vcov <- solve(Jfunction(beta, tau, lambda))
  # } else {
  #   lambda <- lambda_fix
  #   theta.opt <- optim(start[seq.int(length.out = p + q)], logLfixed,
  #     gr = Ufixed,
  #     method = method, lambda = lambda, control = control, hessian = TRUE
  #   )
  #   beta <- theta.opt$par[seq.int(length.out = p)]
  #   tau <- theta.opt$par[seq.int(length.out = q) + p]
  #   if (theta.opt$convergence == 0) {
  #     converged <- TRUE
  #   } else {
  #     converged <- FALSE
  #     warning("optimization failed to converge.", call. = FALSE)
  #   }
  #   vcov <- solve(Jfunction(beta, tau, lambda)[seq.int(length.out = p + q), seq.int(length.out = p + q)])
  # }

  theta.opt

}

# Main fit function ---------------------------------------------------------------------------
#'
#' @name BCSreg
#'
#' @title Box-Cox Symmetric Regression for Positive Data
#'
#' @description Fit the Box-Cox symmetric (BCS) or the zero-adjusted BCS regression models
#'     using maximum likelihood estimation, providing a flexible approach
#'     for modeling positive data.
#'
#' @param formula a symbolic description of the model, allowing the specification of
#' different regression structures for model parameters using the \code{\link[Formula]{Formula}} package.
#' By default, the formula defines the regression structure for the scale parameter.
#' It can include up to three parts, separated by the `|` operator:
#'
#' \itemize{
#'   \item \bold{First part:} specifies the model for the scale parameter.
#'   \item \bold{Second part (optional):} defines a regression structure for the relative dispersion parameter.
#'   \item \bold{Third part (only applicable for zero-inflated positive data):} models the zero-adjustment parameter.
#' }
#'
#' See "Details" for further explanation.
#' @param data,subset,na.action arguments controlling formula processing via
#'     \code{\link[stats]{model.frame}}.
#' @inheritParams BCS
#' @param zeta strictly positive extra parameter. It must be specified with only one value
#'     in cases where the BCS distribution has an extra parameter.
#' @param link,sigma.link character specification of the link functions for
#'     the scale and relative dispersion regression structures, respectively.
#'     Currently, \code{"log"} (default), \code{"sqrt"}, \code{"1/mu^2"},
#'     \code{"inverse"}, and \code{"identity"} are supported.
#' @param alpha.link character specification of the link function for the
#'     zero-adjustment regression structure (only applicable for zero-inflated positive data). Currently, \code{"logit"} (default),
#'     \code{"probit"}, \code{"cloglog"}, \code{"cauchit"}, and \code{"identity"}
#'     are supported.
#' @param control a list of control parameters passed as arguments for
#'     the \code{\link[stats]{optim}} function specified via \code{\link{BCSreg.control}}.
#' @param model,y,x logicals. If \code{TRUE}, the corresponding components of the fit
#'     (the model frame, the response, and the model matrices, respectively) are returned. For
#'     \code{print()}, \code{x} is a fitted model object of class \code{"BCSreg"}.
#' @param ... arguments passed to \code{\link{BCSreg.control}}.
#'
#' @details The \code{BCSreg} function implements maximum likelihood estimation in the
#'     class of the BCS regression models for the analysis of positive
#'     data (Medeiros and Queiroz, 2025). The BCS distributions (Ferrari and Fumes, 2017)
#'     are a broad class of flexible distributions that achieve different levels of
#'     skewness and tail-heaviness. See details in \link{BCS}.
#'
#'
#'     The distributions currently implemented in the \code{BCSreg} package, along
#'     with their abbreviations used in the \code{family} argument, are listed below:
#'     \tabular{llc}{
#'     \bold{Distribution}  \tab \bold{Family abbreviation} \tab \bold{N. of extra parameters}\cr
#'     Box-Cox Hyperbolic  \tab \code{"HP"}      \tab  1  \cr
#'     Box-Cox Type I Logistic  \tab \code{"LOI"}      \tab  0  \cr
#'     Box-Cox Type II Logistic  \tab \code{"LOII"}      \tab  0  \cr
#'     Box-Cox Normal  \tab \code{"NO"}      \tab  0  \cr
#'     Box-Cox Power Exponential  \tab \code{"PE"}      \tab  1  \cr
#'     Box-Cox Sinh-Normal  \tab \code{"SN"}      \tab  1  \cr
#'     Box-Cox Slash  \tab \code{"SL"}      \tab  1  \cr
#'     Box-Cox \emph{t}  \tab \code{"ST"}      \tab  1  \cr
#'     }
#'
#'     The BCS distributions have at least three parameters: scale (\code{mu}),
#'     relative dispersion (\code{sigma}), and skewness (\code{lambda}) parameters.
#'     Some distributions may also depend on an additional parameter (\code{zeta}),
#'     such as the Box-Cox \emph{t} and Box-Cox power exponential distributions. The
#'     BCS distributions reduce to the log-symmetric distributions
#'     (Vanegas and Paula, 2016) when \code{lambda} is fixed at zero. The Log-symmetric
#'     distributions are an important class of probability models for positive data,
#'     which includes well-known distributions such as the log-normal and log-\emph{t}
#'     distributions. The \code{BCSreg} function allows fitting a log-symmetric
#'     regression using \code{lambda = 0} as an argument
#'     (see \code{\link{BCSreg.control}}).
#'
#'     The \code{formula} argument defines the regression structures for different model
#'     parameters using the \code{\link[Formula]{Formula}} package (Zeileis and Croissant, 2010).
#'     It can have up to three parts, separated by the "\code{|}" operator:
#'     \itemize{
#'      \item \bold{First part:} specifies the model for the scale parameter.
#'      \item \bold{Second part (optional):} defines a regression structure for the relative dispersion parameter.
#'      \item \bold{Third part (only applicable for zero-inflated positive data):} models the zero-adjustment parameter.
#'     }
#'
#'     If only the first part is provided, the model applies only to the scale parameter.
#'     When a second part is included, a regression structure is defined for the relative
#'     dispersion parameter. If the data contains zero inflation, a third part can be
#'     specified to model the zero-adjustment parameter.
#'
#'     For instance, consider a dataset where \code{y} is the zero-inflated
#'     dependent variable, and \code{x}, \code{s}, and \code{z} are explanatory
#'     variables associated with the scale, relative dispersion, and zero-adjustment
#'     parameters, respectively. The following formulas illustrate different model structures:
#'     \itemize{
#'      \item \code{y ~ x}  # Scale parameter only
#'      \item \code{y ~ x | s} # Scale and relative dispersion parameters
#'      \item \code{y ~ x | s | z} # Scale, relative dispersion, and zero-adjustment parameters
#'      \item \code{y ~ x | 1 | z} # Scale and zero-adjustment parameters
#'      \item \code{y ~ 1 | s | z} # Relative dispersion and zero-adjustment parameters
#'     }
#'
#' @return The \code{BCSreg} function returns an object of class \code{"BCSreg"},
#'      which consists of a list with the following components:
#'  \describe{
#'     \item{coefficients}{a list containing the elements \code{"mu"} and
#'         \code{"sigma"} that consist of the estimates of the coefficients
#'          associated with the scale and relative dispersion regression structures,
#'          respectively. If the model is zero-adjusted, the element \code{"alpha"}
#'          will also be returned with the estimated coefficients for the regression
#'          structure of the zero-adjustment parameter.}
#'     \item{fitted.values}{a vector with the fitted median responses. Not to be
#'          confused with the fitted values for \code{mu}.}
#'     \item{mu}{a vector with the fitted scale parameters.}
#'     \item{sigma}{a vector with the fitted relative dispersion parameters.}
#'     \item{lambda}{ the maximum likelihood estimate of the skewness parameter (\code{lambda}), or its fixed value
#'         specified in the \code{BCSreg} function.}
#'     \item{zeta}{the specified value for the extra parameter of the corresponding
#'         BCS distribution, if applicable.}
#'     \item{family}{the generating family of the fitted BCS distribution.}
#'     \item{link}{a list with elements \code{"mu"} and \code{"sigma"} with the
#'         specified link functions for the \code{mu} and \code{sigma} regression
#'         structures, respectively. If the model is zero-adjusted, the element
#'         \code{"alpha"} will also be returned with the link function for
#'         the regression structure of the zero-adjustment parameter.}
#'     \item{logLik}{log-likelihood of the fitted model.}
#'     \item{vcov}{asymptotic covariance matrix of the estimators. By default, the asymptotic
#'         covariance matrix is based on a analytical expression of the observed information matrix.
#'         It can be obtained numerically based on the Hessian matrix via \code{\link[stats]{optim}})
#'         if the argument \code{hessian = TRUE} is used in the \code{BCSreg} function.}
#'     \item{nobs}{number of observations.}
#'     \item{df.null}{residual degrees of freedom in the null model (a model without
#'         any regression structure).}
#'     \item{df.residual}{residual degrees of freedom in the fitted model, that is,
#'         the sample size minus the number of model parameters.}
#'     \item{control}{the control arguments passed to the optim call.}
#'     \item{start}{a vector with the starting values used in the iterative process.}
#'     \item{optim}{a list with the output from \code{\link[stats]{optim}}.}
#'     \item{converged}{logical indicating successful convergence of the iterative
#'         process.}
#'     \item{call}{the original function call.}
#'     \item{formula}{the formula used.}
#'     \item{terms}{a list with elements "\code{mu}", "\code{sigma}", and "\code{full}" containing
#'         the term objects for the respective models. If the model is zero-adjusted,
#'         the element \code{"alpha"} will also be returned with the term object for the
#'         zero-adjustment model.}
#'     \item{levels}{a list with elements "\code{mu}", "\code{sigma}", and "\code{full}" containing
#'         the levels of the categorical regressors. If the model is zero-adjusted, the element
#'         \code{"alpha"} will also be returned.}
#'     \item{contrasts}{a list with elements "\code{mu}" and "\code{sigma}"
#'         containing the contrasts corresponding to levels from the respective models. If the model is zero-adjusted, the element
#'         \code{"alpha"} will also be returned. }
#'     \item{model}{the full model frame (if \code{y = TRUE}).}
#'     \item{y}{the response variable (if \code{y = TRUE}).}
#'     \item{x}{a list with elements "\code{mu}" and "\code{sigma}" with the model matrices from
#'         the \code{mu} and \code{sigma} submodels (if \code{x = TRUE}). If the model is zero-adjusted, the element
#'         \code{"alpha"} will also be returned with the model matrix for the
#'         zero-adjustment submodel.}
#'     \item{alpha}{a vector with the fitted zero-adjustment parameters when a zero-adjusted
#'         model is considered; and \code{NULL}, otherwise.}
#'    }
#'
#' @references Cribari-Neto F, Zeileis A (2010). Beta Regression in R. \emph{Journal of Statistical
#'     Software}, \bold{34}, 1---24
#'
#'     Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'     applications to nutritional data. \emph{AStA Advances in Statistical Analysis},
#'     \bold{101}, 321---344.
#'
#'     Medeiros, R. M. R., and Queiroz, F. F. (2025). Modeling positive continuous data:
#'     Box-Cox symmetric regression models and their extensions
#'
#'     Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions:
#'     statistical properties and parameter estimation. \emph{Brazilian Journal of
#'     Probability and Statistics}, \bold{30},196---220
#'
#'     Zeileis A, Croissant Y (2010). Extended Model Formulas in R: Multiple Parts and Multiple
#'     Responses. \emph{Journal of Statistical Software}, \bold{34}, 1---13.
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#' @export
#'
#'
BCSreg <- function(formula, data, subset, na.action,
                    family = "NO", zeta,
                    link = "log", sigma.link = "log", alpha.link,
                    control = BCSreg.control(...), model = FALSE,
                    y = FALSE, x = FALSE, ...) {
  ## Model call
  cl <- match.call()

  ## Model description
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  oformula <- stats::as.formula(formula)
  formula <- Formula::as.Formula(formula)
  formula_len <- length(formula)[2L]
  if (formula_len == 1L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1, ~1)
  } else if (formula_len == 2L) {
    formula <- Formula::as.Formula(formula(formula), ~1)
  } else if (formula_len > 3L) {
    formula <- Formula::Formula(formula(formula, rhs = 1:3))
    warning("formula must not have more than three RHS parts")
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- stats::terms(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1)
  mtS <- stats::delete.response(stats::terms(formula, data = data, rhs = 2))
  mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 3))
  Y <- stats::model.response(mf, "numeric")
  ind <- ifelse(Y == 0, 1, 0)
  p0 <- mean(ind)
  X <- stats::model.matrix(mtX, mf)
  S <- stats::model.matrix(mtS, mf)
  Z <- stats::model.matrix(mtZ, mf)
  if (length(Y) < 1) {
    stop("Empty model", call. = FALSE)
  }
  if (min(Y) < 0) {
    stop("Invalid dependent variable, all observations must be non-negative.", call. = FALSE)
  }

  ## Lenghts
  n <- length(Y)
  p <- ncol(X)
  q <- ncol(S)
  m <- ncol(Z)

  ## Symmetric generating family
  family <- match.arg(family, c("NO", "HP", "LOI", "LOII", "PE", "SN", "SL", "ST"))
  if (family %in% c("HP", "PE", "SN", "SL", "ST")) {
    if (missing(zeta))
      stop(gettextf("%s family has an extra parameter that must be specified via the 'zeta' argument",
                    sQuote(family)), domain = NA)

    if (zeta < 0)
      warning(paste("zeta must be positive", "\n", ""))
  } else {
    zeta <- NULL
  }

  lambda_id <- is.null(control$lambda) # TRUE means that lambda will be estimated
  zeta_id <- !is.null(zeta)            # TRUE means that the family has an extra parameter

  ## Link functions
  if (is.character(link)) {
    link <- match.arg(link, c("log", "sqrt", "inverse", "1/mu^2", "identity"))
  } else {
    stop("The link argument must be a character", call. = FALSE)
  }
  if (is.character(sigma.link)) {
    sigma.link <- match.arg(sigma.link, c("log", "sqrt", "inverse", "1/mu^2", "identity"))
  } else {
    stop("The sigma.link argument must be a character", call. = FALSE)
  }

  ## Model fit
  opt_fit <- BCSreg.fit(
    X = as.matrix(X[ind == 0, ]), y = Y[ind == 0], S = as.matrix(S[ind == 0, ]),
    zeta = zeta, family = family,
    link = link, sigma.link = sigma.link, control = control
  )

  if (opt_fit$convergence == 0) {
    converged <- TRUE
  } else {
    converged <- FALSE
    warning("Optimization failed to converge.", call. = FALSE)
  }

  ## Model parameters
  start <- opt_fit$start
  opt_fit$start <- NULL

  beta <- opt_fit$par[seq.int(length.out = p)]
  names(beta) <- colnames(X)
  if (is.null(names(beta))) names(beta) <- paste0("X", 1:p)

  tau <- opt_fit$par[seq.int(length.out = q) + p]
  names(tau) <- colnames(S)
  if (is.null(names(tau))) names(tau) <- paste0("S", 1:p)

  mu <- stats::make.link(link)$linkinv(as.vector(X %*% beta))
  sigma <- stats::make.link(sigma.link)$linkinv(as.vector(S %*% tau))
  if (lambda_id) {
    lambda <- opt_fit$par[p + q + 1]
    names(lambda) <- "(lambda)"
  } else {
    lambda <- control$lambda
  }

  ## Covariance matrix
  vcov <- opt_fit$vcov
  opt_fit$vcov <- NULL
  rownames(vcov) <- colnames(vcov) <- c(names(beta),
                                        paste0("(sigma)_", names(tau)),
                                        names(lambda))

  ## Fitted values
  fitted.values <- structure(qBCS(0.5, mu, sigma, lambda, zeta, family), .Names = names(y))

  val <- list(
    coefficients = list(mu = beta, sigma = tau),
    fitted.values = fitted.values,
    mu = mu,
    sigma = sigma,
    lambda = lambda,
    zeta = zeta,
    family = family,
    link = list(mu = link, sigma = sigma.link),
    loglik = opt_fit$value,
    vcov = vcov,
    nobs = n,
    df.null = n - 2 - as.numeric(lambda_id) - as.numeric(zeta_id),
    df.residual = n - p - q - as.numeric(lambda_id) - as.numeric(zeta_id),
    control = control,
    start = start,
    optim = opt_fit,
    converged = converged
  )

  val$call <- cl
  val$formula <- oformula
  val$terms <- list(mu = mtX, sigma = mtS, full = mt)
  val$levels <- list(
    mu = stats::.getXlevels(mtX, mf),
    sigma = stats::.getXlevels(mtS, mf),
    full = stats::.getXlevels(mt, mf)
  )
  val$contrasts <- list(mu = attr(X, "contrasts"), sigma = attr(S, "contrasts"))
  if (model) {
    val$model <- mf
  }
  if (y) {
    val$y <- Y
  }
  if (x) {
    val$x <- list(mu = X, sigma = S)
  }

  ## Zero-adjusted model
  if (sum(ind) == 0 & formula_len == 3) {
    warning("The third part of the RHS formula will be ignored: the dependent variable is not zero-inflated")
  }

  if (sum(ind) > 0) {

    if (missing(alpha.link)){
      alpha.link <- "logit"
    }

    glm_fit <- stats::glm.fit(Z, ind, family = stats::binomial(link = alpha.link))

    val$terms <- list(mu = mtX, sigma = mtS, alpha = mtZ, full = mt)
    val$levels <- list(
      mu = stats::.getXlevels(mtX, mf),
      sigma = stats::.getXlevels(mtS, mf),
      alpha = stats::.getXlevels(mtZ, mf),
      full = stats::.getXlevels(mt, mf)
    )
    val$contrasts$alpha <- attr(Z, "contrasts")
    if (x) {
      val$x$alpha <- Z
    }

    kappa <- glm_fit$coefficients
    names(kappa) <- colnames(Z)
    if (is.null(names(kappa))) names(kappa) <- paste0("Z", 1:m)
    val$coefficients$alpha <- kappa
    val$alpha <- as.numeric(glm_fit$fitted.values)
    val$link$alpha <- alpha.link

    ## Fitted values
    fitted.values <- rep(NA, n)
    fitted.values[val$alpha > 0.5] <- 0L
    fitted.values[val$alpha <= 0.5] <- qBCS((0.5 - val$alpha[val$alpha <= 0.5]) /
                                              (1 - val$alpha[val$alpha <= 0.5]),
                                            mu = mu[val$alpha <= 0.5],
                                            sigma = sigma[val$alpha <= 0.5],
                                            lambda = lambda, zeta = zeta,
                                            family = family)
    val$fitted.values <- fitted.values

    ## Covariance matrix
    npar_BCS <- p + q + as.numeric(lambda_id)
    vcov <- matrix(0L, npar_BCS + m, npar_BCS + m)
    vcov[1:npar_BCS, 1:npar_BCS] <- val$vcov
    vcov[1:m + npar_BCS, 1:m + npar_BCS] <- chol2inv(chol(t(Z)%*%diag(glm_fit$weights)%*%Z))
    rownames(vcov) <- colnames(vcov) <- c(colnames(val$vcov), paste0("(alpha)_", names(kappa)))
    val$vcov <- vcov
    val$df.null <- n - 3 - as.numeric(lambda_id) - as.numeric(zeta_id)
    val$df.residual <- n - npar_BCS - as.numeric(zeta_id) - m

    ## Log-likelihood value
    ll <- rep(NA, n)
    ll[ind == 1L] <- log(val$alpha[ind == 1])
    ll[ind == 0L] <- log((1 - val$alpha[ind == 0]) * dBCS(Y[ind == 0],
                                                          mu = mu[ind == 0],
                                                          sigma = sigma[ind == 0],
                                                          lambda = lambda, zeta = zeta,
                                                          family = family))
    val$loglik <- sum(ll)
  }

  class(val) <- "BCSreg"
  val
}

# print() method ------------------------------------------------------------------------------

# Print
#' @rdname BCSreg
#' @param digits a non-null value for digits specifies the minimum number of significant digits to
#'     be printed in values.
#' @export
print.BCSreg <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  lambda_id <- is.null(x$control$lambda)
  zeta_id <- !is.null(x$zeta)
  alpha_id <- !is.null(x$alpha)

  cat("Call:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), sep = "\n")
  if (!x$converged) {
    cat("Model did not converge\n")
  } else {

    # Discrete component (if applicable)
    if (alpha_id) {
      cat("\n--- Fit for the discrete component ---\n")
      cat("\nZero-adjustment submodel with ", x$link$alpha, " link:\n", sep = "")
      print.default(format(x$coefficients$alpha, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n--- Fit for the continuous component ---\n")
    }

    # Scale submodel
    cat("\nScale submodel with ", x$link$mu, " link:\n", sep = "")
    print.default(format(x$coefficients$mu, digits = digits), print.gap = 2, quote = FALSE)

    cat(paste("\nRelative dispersion submodel with ", x$link$sigma, " link:\n", sep = ""))
    print.default(format(x$coefficients$sigma, digits = digits), print.gap = 2, quote = FALSE)
    if (lambda_id) {
      cat("\nSkewness parameter:\n", sep = "")
      print.default(format(x$lambda, digits = digits), print.gap = 2, quote = FALSE)
    }

    cat("\n---\nGenerating family:", x$family, if (zeta_id) paste0("(zeta = ", x$zeta, ")"),
        "\nLog-lik value:", x$loglik, "\n")
  }

  invisible(x)
}

