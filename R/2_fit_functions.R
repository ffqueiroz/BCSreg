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

BCSreg.fit <- function(X, y, S = NULL, family, zeta = zeta, link = "log",
                       sigma.link = "log", control = BCSreg.control()) {
  n <- length(y)
  p <- dim(X)[2]

  ## Model matrices
  if (is.null(colnames(X))) {
    if (p == 1) {
      colnames(X) <- "(Intercept)"
    } else {
      colnames(X)[1] <- "(Intercept)"
      for (i in 2:p) {
        colnames(X)[1] <- paste("V", i, sep = "")
      }
    }
  }

  if (is.null(S)) {
    q <- 1
    S <- matrix(1, ncol = q, nrow = n)
    colnames(S) <- "(Intercept)"
    rownames(S) <- rownames(X)
  } else {
    q <- dim(S)[2]
    if (q < 1L) stop("relative dispersion regression needs to have at least one parameter", call. = FALSE)
  }

  sigma_const <- (q == 1) && (sigma.link == "identity")

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
  lambda_id <- is.null(control$lambda) # TRUE means that lambda is not fixed
  lambda_fix <- control$lambda # NULL means that lambda is not fixed
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

  ## v() function
  v.function <- function(mu, sigma, lambda) {
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
      vz <- -2 * (1 - exp(-z^2)) / (1 + exp(-z^2))
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
      ll <- suppressWarnings(dbcs(y,
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
  #     ll <- suppressWarnings(dbcs(y,
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

  ## xi function
  xi.funcion <- function(mu, sigma, lambda) {
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
        vz <- -2 * (1 - exp(-(1 / ((sigma * lambda)))^2)) / (1 + exp(-(1 / ((sigma * lambda)))^2))
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

    vdv <- v.function(mu, sigma, lambda)
    v <- vdv$v
    z <- vdv$z
    xi.f <- xi.funcion(mu, sigma, lambda)
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
  #   vdv <- v.function(mu, sigma, lambda)
  #   v <- vdv$v
  #   z <- vdv$z
  #   xi.f <- xi.funcion(mu, sigma, lambda)
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

    vdv <- v.function(mu, sigma, lambda)
    v <- vdv$v
    z <- vdv$z
    dv <- vdv$dv
    xi.f <- xi.funcion(mu, sigma, lambda)
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
  lambda <- if (lambda_id) theta.opt$par[1L + p + q] else lambda_fix

  if (theta.opt$convergence == 0) {
    converged <- TRUE
  } else {
    converged <- FALSE
    warning("Optimization failed to converge.", call. = FALSE)
  }

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

  eta.1 <- as.vector(X %*% beta)
  eta.2 <- as.vector(S %*% tau)
  mu <- linkinv(eta.1)
  sigma <- sigma_linkinv(eta.2)

  ## Goodness of fit Upsilon statistic
  Upsilon <- function(zeta) {
    fda <- sort(pbcs(y, mu = mu, sigma = sigma, lambda = lambda, zeta = zeta, family = family))
    Upsilon_zeta <- mean(abs(qnorm(fda) - EnvStats::evNormOrdStats(n = n)))
    Upsilon_zeta
  }

  ## Quantile residuals
  residuals <- stats::qnorm(pbcs(y, mu, sigma, lambda, zeta, family))

  optim.fit <- theta.opt
  ll <- logL(c(beta, tau, lambda), lambda_id = TRUE)
  Ups.zeta <- Upsilon(zeta)

  pseudor2 <- ifelse(stats::var(eta.1) * stats::var(linkfun(y)) <= 0, NA, stats::cor(eta.1, linkfun(y))^2)
  v <- v.function(mu, sigma, lambda)$v

  names(beta) <- colnames(X)
  names(tau) <- if (sigma_const) "(sigma)" else colnames(S)
  names(lambda) <- "(lambda)"

  if (lambda_id) {
    rownames(vcov) <- colnames(vcov) <- c(colnames(X), if (sigma_const) {
      "(sigma)"
    } else {
      paste("(sigma)",
        colnames(S),
        sep = "_"
      )
    }, "(lambda)")
  } else {
    rownames(vcov) <- colnames(vcov) <- c(colnames(X), if (sigma_const) {
      "(sigma)"
    } else {
      paste("(sigma)",
        colnames(S),
        sep = "_"
      )
    })
  }

  ## Out
  val <- list(
    coefficients = list(mu = beta, sigma = tau),
    mu = mu,
    sigma = sigma,
    lambda = lambda,
    zeta = zeta,
    fitted.values = structure(mu, .Names = names(y)),
    family = family,
    link = list(mu = link, sigma = sigma.link),
    loglik = ll,
    vcov = vcov,
    residuals = residuals, ## Quantile residuals
    pseudo.r.squared = pseudor2,
    Upsilon.zeta = Ups.zeta,
    v = v,
    nobs = n,
    df.null = n - ifelse(is.null(lambda_fix), 3, 2),
    df.residual = n - p - q - ifelse(is.null(lambda_fix), 1, 0),
    control = ocontrol,
    start = start,
    optim = optim.fit,
    converged = converged
  )

  val
}



#' Box-Cox Symmetric Regression for Positive Data
#'
#' Maximum likelihood estimation in the class of the Box-Cox symmetric (BCS) regression models
#'     for the analysis of positive data.
#'
#'
#' @param formula a symbolic description of the model, of type \code{y ~ x} to specify
#'     explanatory variables in the scale (\code{mu}) regression structure only or \code{y ~ x | z}
#'     to also specify explanatory variables in the relative dispersion (\code{sigma}) regression
#'     structure. See details below.
#' @param data,subset,na.action arguments controlling formula processing via
#'     \code{\link[stats]{model.frame}}.
#' @inheritParams bcs
#' @inheritParams bcs
#' @param link,sigma.link character specification of the link function for the
#'     scale (\code{mu}) and relative dispersion (\code{sigma}) regression structures, respectively.
#'     Currently, \code{"log"} (default), \code{"sqrt"}, \code{"1/mu^2"},
#'     \code{"inverse"}, and \code{"identity"} are supported.
#' @param control a list of control parameters passed as arguments for
#'     the \code{\link[stats]{optim}} function specified via \code{\link{BCSreg.control}}.
#' @param model,y,x logicals. If \code{TRUE}, the corresponding components of the fit
#'     (model frame, response, and model matrices) are returned.
#' @param ... arguments passed to \code{\link{BCSreg.control}}.
#'
#' @details The \code{BCSreg} function implements maximum likelihood estimation in the
#'     class of the BCS regression models for the analysis of positive
#'     data (Medeiros and Queiroz, 2025). The BCS distributions (Ferrari and Fumes, 2017)
#'     are a broad class of flexible distributions that achieve different levels of
#'     skewness and tail-heaviness. See details in \link{bcs}.
#'
#'     The distributions currently implemented in the \code{BCSreg} package, along
#'     with their abbreviations used in the \code{family} argument, are listed below:
#'     \tabular{llc}{
#'      \bold{Distribution}  \tab \bold{Family abbreviation} \tab \bold{Number of parameters}\cr
#'      Box-Cox Hyperbolic  \tab \code{"HP"}      \tab  4  \cr
#'      Box-Cox Type I Logistic  \tab \code{"LOI"}      \tab  3  \cr
#'      Box-Cox Type II Logistic  \tab \code{"LOII"}      \tab  3  \cr
#'      Box-Cox Normal  \tab \code{"NO"}      \tab  3  \cr
#'      Box-Cox Power Exponential  \tab \code{"PE"}      \tab  4  \cr
#'      Box-Cox Sinh-Normal  \tab \code{"SN"}      \tab  4  \cr
#'      Box-Cox Slash  \tab \code{"SL"}      \tab  4  \cr
#'      Box-Cox \emph{t}  \tab \code{"ST"}      \tab  4  \cr
#'      }
#'
#'     The BCS distributions have at least three parameters: scale (\code{mu}),
#'     relative dispersion (\code{sigma}), and skewness (\code{lambda}) parameters.
#'     Some distributions may also depend on an additional parameter (\code{zeta}),
#'     such as the Box-Cox \emph{t} and Box-Cox power exponential distributions. The
#'     BCS distributions reduce to the log-symmetric distributions when \code{lambda}
#'     is fixed at zero. The Log-symmetric distributions are an important class of
#'     probability models for positive data, which includes well-known distributions such
#'     as the log-normal and log-\emph{t} distributions. The \code{BCSreg} function
#'     allows fitting a log-symmetric regression using \code{lambda = 0} as an
#'     argument (see \code{\link{BCSreg.control}}).
#'
#'
#'     The \code{BCSreg} function allows defining a regression structure in the
#'     \code{mu} and \code{sigma} parameters via the \code{formula} argument. The
#'     basic formula is of type \code{y ~ x1 + x2 + ... + xp} which specifies the
#'     regression structure for \code{mu} only. Following the syntax of the
#'     \code{betareg} package (Cribari-Neto and Zeileis, 2010), the regression structure
#'     for \code{sigma}, say in terms of \code{s1, s2, ..., sq}, is specified as
#'     \code{y ~ x1 + x2 + ... + xp | s1 + s2 + ... + sq} using functionalities inherited
#'     from the \code{\link[Formula]{Formula}} package (Zeileis and Croissant, 2010).
#'
#' @return The \code{BCSreg} function returns an object of class \code{"BCSreg"},
#'      which consists of a list with the following components:
#'  \describe{
#'     \item{coefficients}{ a list containing the elements \code{"mu"} and
#'         \code{"sigma"} that consist of the estimates of the coefficients
#'          associated with the scale and relative dispersion regression structures, respectively.}
#'     \item{fitted.values}{ a vector with the fitted median responses. Not to be confused with the fitted values for \code{mu}.}
#'     \item{mu}{ a vector with the fitted scale parameters.}
#'     \item{sigma}{ a vector with the fitted relative dispersion parameters.}
#'     \item{lambda}{ the maximum likelihood estimate of the skewness parameter (\code{lambda}), or its fixed value
#'         specified in the \code{BCSreg} function.}
#'     \item{zeta}{ the specified value for the extra parameter of the corresponding BCS distribution, if applicable.}
#'     \item{family}{ the generating family of the fitted BCS distribution.}
#'     \item{link}{ a list with elements \code{"mu"} and \code{"sigma"} with the
#'         specified link functions for the \code{mu} and \code{sigma} regression structures, respectively.}
#'     \item{logLik}{ log-likelihood of the fitted model.}
#'     \item{vcov}{ asymptotic covariance matrix of the estimators. By default, the asymptotic
#'         covariance matrix is based on a analytical expression of the observed information matrix.
#'         It can be obtained numerically based on the Hessian matrix via \code{\link[stats]{optim}})
#'         if the argument \code{hessian = TRUE} is used in the \code{BCSreg} function.}
#'     \item{residuals}{ a vector of quantile residuals.}
#'     \item{pseudo.r.squared}{pseudo R-squared value.}
#'     \item{Upsilon.zeta}{an overall goodness-of-fit measure.}
#'     \item{v}{a vector with the \eqn{v(z)} values for all the observations.
#'         The \eqn{v} function is a weighting function that is involved in the
#'         parameter estimation process. It depends on the generating density function
#'         of the corresponding BCS distribution.}
#'     \item{nobs}{ number of observations.}
#'     \item{df.null}{ residual degrees of freedom in the null model (constant scale and relative dispersion).}
#'     \item{df.residual}{ residual degrees of freedom in the fitted model.}
#'     \item{control}{the control arguments passed to the optim call.}
#'     \item{start}{a vector with the starting values used in the iterative process.}
#'     \item{optim}{a list with the output from \code{\link[stats]{optim}}.}
#'     \item{converged}{logical indicating successful convergence of the iterative process.}
#'     \item{call}{the original function call.}
#'     \item{formula}{the formula used.}
#'     \item{terms}{a list with elements "\code{mu}", "\code{sigma}", and "\code{full}" containing
#'         the term objects for the respective models.}
#'     \item{levels}{a list with elements "\code{mu}", "\code{sigma}", and "\code{full}" containing
#'         the levels of the categorical regressors.}
#'     \item{contrasts}{a list with elements "\code{mu}" and "\code{sigma}"
#'         containing the contrasts corresponding to levels from the respective models.}
#'     \item{model}{the full model frame (if \code{y = TRUE}).}
#'     \item{y}{the response variable (if \code{y = TRUE}).}
#'     \item{x}{a list with elements "\code{mu}" and "\code{sigma}" with the model matrices from
#'         the \code{mu} and \code{sigma} submodels (if \code{x = TRUE}).}
#'    }
#'
#' @references Cribari-Neto F, Zeileis A (2010). Beta Regression in R. \emph{Journal of Statistical
#'     Software}, \bold{34}, 1-24
#'
#'     Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'     applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, \bold{101}, 321-344.
#'
#'     Medeiros, R. M. R., and Queiroz, F. F. (2025). Modeling positive continuous data:
#'     Box-Cox symmetric regression models and their extensions
#'
#'     Zeileis A, Croissant Y (2010). Extended Model Formulas in R: Multiple Parts and Multiple
#'     Responses. \emph{Journal of Statistical Software}, \bold{34}, 1-13.
#'
#'
#' @export
#'
#'
BCSreg <- function(formula, data, subset, na.action,
                   family = "NO", zeta,
                   link = "log", sigma.link = "log",
                   control = BCSreg.control(...), model = TRUE,
                   y = TRUE, x = FALSE, ...) {
  cl <- match.call()
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  oformula <- stats::as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if (length(formula)[2] < 2) {
    formula <- Formula::as.Formula(formula(formula), ~1)
    simple_formula <- TRUE
  } else {
    if (length(formula)[2] > 2) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts.", call. = FALSE) # RHS right-hand side
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- stats::terms(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1)
  mtS <- stats::delete.response(stats::terms(formula, data = data, rhs = 2))
  Y <- stats::model.response(mf, "numeric")
  X <- stats::model.matrix(mtX, mf)
  S <- stats::model.matrix(mtS, mf)
  if (length(Y) < 1) {
    stop("Empty model", call. = FALSE)
  }
  if (!(min(Y) > 0)) {
    stop("Invalid dependent variable, all observations must be positive.", call. = FALSE)
  }
  n <- length(Y)
  family <- match.arg(family, c("NO", "HP", "LOI", "LOII", "PE", "SN", "SL", "ST"))
  if (family == "SL" | family == "ST" | family == "SN" | family == "HP" | family == "PE") {
    if (missing(zeta)) {
      stop("For the family of distributions specified by the user an extra parameter is required.", call. = FALSE)
    } else {
      if (zeta <= 0) stop("Invalid extra parameter; zeta must be positive.", call. = FALSE)
    }
  } else {
    zeta <- NULL
  }
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
  val <- BCSreg.fit(
    X = X, y = Y, S = S, zeta = zeta, family = family,
    link = link, sigma.link = sigma.link, control = control
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
  class(val) <- "BCSreg"
  val
}
