#' The Box-Cox Symmetric (BCS) Distributions
#'
#' Density, distribution function, quantile function, and random generation
#'     for the class of the Box-Cox symmetric (BCS) distributions.
#'
#' @name bcs
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `n` is a vector, its length is used as the number of
#'     required observations.
#' @param mu vector of strictly positive scale parameters.
#' @param sigma vector of strictly positive relative dispersion parameters.
#' @param lambda vector of real-valued skewness parameters. If \code{lambda = 0}, the BCS
#'     distribution reduces to the corresponding log-symmetric distribution with parameters
#'     \code{mu}, \code{sigma}, and \code{zeta}.
#' @param zeta strictly positive extra parameter. It must be specified with only one value.
#' @param family a character that specifies the generating family of the BCS distribution.
#'     Available options are: \code{"NO"}, \code{"ST"}, \code{"LOI"}, \code{"LOII"},
#'     \code{"PE"}, \code{"SN"}, \code{"HP"}, and \code{"SL"}, corresponding to the normal,
#'     Student-\emph{t}, type I logistic, type II logistic, power exponential, sinh-normal,
#'     hyperbolic, and slash distributions, respectively.
#' @param log,log.p logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#'     Default is \code{FALSE}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'    \eqn{P(X \le x)} otherwise, \eqn{P(X > x)}.
#'
#' @details The class of the BCS distributions was introduced by Ferrari and
#'     Fumes (2017). It consists of a broad class of probability models for
#'     positive continuous data, which includes flexible distributions with
#'     different levels of skewness and tail-heaviness.
#'
#'
#'
#'     The BCS class includes, as special cases, the Box-Cox \emph{t} (Rigby and Stasinopoulos, 2006),
#'     Box-Cox normal (or Box-Cox Cole-Green; Cole and Green, 1992), Box-Cox power exponential
#'     (Rigby and Stasinopoulos, 2004) distributions, as well as the log-symmetric
#'     distributions (Vanegas and Paula, 2016).
#'
#'     The currently available BCS distributions in the \code{BCSreg} package are listed below:
#'
#' \tabular{llc}{
#'  \bold{Distribution}  \tab \bold{Family abbreviation} \tab \bold{Number of parameters}\cr
#'  Box-Cox Hyperbolic  \tab \code{"HP"}      \tab  4  \cr
#'  Box-Cox Type I Logistic  \tab \code{"LOI"}      \tab  3  \cr
#'  Box-Cox Type II Logistic  \tab \code{"LOII"}      \tab  3  \cr
#'  Box-Cox Normal  \tab \code{"NO"}      \tab  3  \cr
#'  Box-Cox Power Exponential  \tab \code{"PE"}      \tab  4  \cr
#'  Box-Cox Sinh-Normal  \tab \code{"SN"}      \tab  4  \cr
#'  Box-Cox Slash  \tab \code{"SL"}      \tab  4  \cr
#'  Box-Cox \emph{t}  \tab \code{"ST"}      \tab  4  \cr
#'  }
#'
#' @return
#' \code{dbcs} returns the density function, \code{pbcs} gives the cumulative distribution function,
#' \code{qbcs} provides the quantile function, and \code{rbcs} generates random variables.
#'
#' @references Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical
#'     properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#' @references
#'  Cole, T., and Green, P.J. (1992). Smoothing reference centile curves: the LMS
#'      method and penalized likelihood. \emph{Statistics in medicine}, 11, 1305-1319.
#'
#'  Ferrari, S. L., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.
#'
#'  Rigby, R. A., and Stasinopoulos, D. M. (2004). Smooth centile curves for skew
#'      and kurtotic data modelled using the Box-Cox power exponential
#'      distribution. \emph{Statistics in medicine}, 23, 3053-3076.
#'
#'  Rigby, R. A., and Stasinopoulos, D. M. (2006). Using the Box-Cox t
#'      distribution in GAMLSS to model skewness and kurtosis. \emph{Statistical Modelling}, 6, 209-229.
#'
#'  Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions:
#'      statistical properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#' @examples
#' # Density
#'
#' ## Right-skewed distributions
#' curve(dbcs(x, 3, 0.3, -1.5, family = "NO"), xlim = c(0, 7), ylim = c(0, 0.7), ylab = "Density")
#' curve(dbcs(x, 3, 0.3, -1.5, 4, family = "ST"), add = TRUE, col = 2)
#' curve(dbcs(x, 3, 0.3, -1.5, 5, family = "PE"), add = TRUE, col = 4)
#' legend("topright", legend = c("NO", "ST", "PE"), lty = 1, col = c(1, 2, 4))
#'
#' ## Truncated symmetric distributions (with support on (0, Inf))
#' curve(dbcs(x, 3, 0.3, 1, family = "NO"), xlim = c(0, 7), ylim = c(0, 0.7), ylab = "Density")
#' curve(dbcs(x, 3, 0.3, 1, 4, family = "ST"), add = TRUE, col = 2)
#' curve(dbcs(x, 3, 0.3, 1, 5, family = "PE"), add = TRUE, col = 4)
#' legend("topright", legend = c("NO", "ST", "PE"), lty = 1, col = c(1, 2, 4))
#'
#' ## Left-skewed distributions
#' curve(dbcs(x, 3, 0.3, 3, family = "NO"), xlim = c(0, 7), ylim = c(0, 0.7), ylab = "Density")
#' curve(dbcs(x, 3, 0.3, 3, 4, family = "ST"), add = TRUE, col = 2)
#' curve(dbcs(x, 3, 0.3, 3, 5, family = "PE"), add = TRUE, col = 4)
#' legend("topright", legend = c("NO", "ST", "PE"), lty = 1, col = c(1, 2, 4))
#'
#' # Random generation
#'
#' ## Parameter setting
#' mu <- 5 # scale parameter
#' sigma <- 0.2 # relative dispersion parameter
#' lambda <- -2 # skewness parameter
#' zeta <- 6 # extra parameter (if necessary)
#'
#' ## Generating family
#' family <- "NO"
#'
#' ## Visualization
#' x <- rbcs(10000, mu, sigma, lambda, zeta, family = family)
#'
#' hist(x, prob = TRUE, col = "white", main = "")
#' curve(dbcs(x, mu, sigma, lambda, zeta, family = family), col = "blue", add = TRUE)
#'
#' plot(ecdf(x), main = "")
#' curve(pbcs(x, mu, sigma, lambda, zeta, family = family), col = "blue", add = TRUE) #'
#'
#' @importFrom stats dnorm qnorm rnorm pnorm dlogis plogis qlogis runif dt pt qt
#'
#' @export
dbcs <- function(x, mu, sigma, lambda, zeta = NULL, family, log = FALSE) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(lambda)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  if (any(mu <= 0)) {
    warning(paste("mu must be positive ", "\n", ""))
  }
  if (any(sigma <= 0)) {
    warning(paste("sigma must be positive", "\n", ""))
  }
  if (any(zeta <= 0)) {
    warning(paste("zeta must be positive", "\n", ""))
  }

  l0 <- which(lambda == 0)
  l1 <- which(lambda != 0)

  z <- rep(NA, maxl)
  z[l0] <- log(x[l0] / mu[l0]) / sigma[l0]
  z[l1] <- ((x[l1] / mu[l1])^lambda[l1] - 1) / (sigma[l1] * lambda[l1])

  switch(family,
    NO = {
      r <- dnorm(z)
      R <- pnorm(1 / (sigma * abs(lambda)))
    },
    ST = {
      r <- dt(z, zeta)
      R <- pt(1 / (sigma * abs(lambda)), zeta)
    },
    LOI = {
      r <- dlogisI(z)
      R <- plogisI(1 / (sigma * abs(lambda)))
    },
    LOII = {
      r <- dlogis(z)
      R <- plogis(1 / (sigma * abs(lambda)))
    },
    SN = {
      r <- (2 * cosh(sqrt(z^2)) * exp(-2 * sinh(sqrt(z^2)) * sinh(sqrt(z^2)) / zeta^2) / (sqrt(2 * pi) * zeta))
      R <- pnorm((2 / zeta) * sinh(1 / (sigma * abs(lambda))))
    },
    PE = {
      r <- gamlss.dist::dPE(z, mu = 0, sigma = 1, nu = zeta)
      R <- gamlss.dist::pPE(1 / (sigma * abs(lambda)), mu = 0, sigma = 1, nu = zeta)
    },
    HP = {
      r <- GeneralizedHyperbolic::dhyperb(z, mu = 0, delta = 1, alpha = zeta, beta = 0)
      R <- GeneralizedHyperbolic::phyperb(1 / (sigma * abs(lambda)),
        mu = 0, delta = 1,
        alpha = zeta, beta = 0
      )
    },
    SL = {
      r <- dslash(z, zeta = zeta)
      R <- pslash(1 / (sigma * abs(lambda)), zeta = zeta)
    },
    stop(gettextf("%s family not recognised", sQuote(family)), domain = NA)
  )


  log.lik <- rep(-Inf, maxl)

  # NaN index
  log.lik[which(mu <= 0 | sigma <= 0 | zeta <= 0)] <- NaN

  # Positive density index
  id0 <- which(x > 0 & lambda == 0 & !is.nan(log.lik))
  id1 <- which(x > 0 & lambda != 0 & !is.nan(log.lik))

  log.lik[id0] <- -log(x[id0] * sigma[id0]) + log(r[id0])
  log.lik[id1] <- (lambda[id1] - 1) * log(x[id1]) - log(sigma[id1] * mu[id1]^lambda[id1]) + log(r[id1]) - log(R[id1])

  if (!log) log.lik <- exp(log.lik)
  if (d > 1L) matrix(log.lik, ncol = d) else log.lik
}

#' @rdname bcs
#' @export
pbcs <- function(q, mu, sigma, lambda, zeta = NULL, family, lower.tail = TRUE, log.p = FALSE) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(lambda)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  if (any(mu <= 0)) {
    warning(paste("mu must be positive ", "\n", ""))
  }
  if (any(sigma <= 0)) {
    warning(paste("sigma must be positive", "\n", ""))
  }
  if (any(zeta <= 0)) {
    warning(paste("zeta must be positive", "\n", ""))
  }

  l0 <- which(lambda == 0)
  l1 <- which(lambda != 0)

  z <- rep(NA, maxl)
  z[l0] <- log(q[l0] / mu[l0]) / sigma[l0]
  z[l1] <- ((q[l1] / mu[l1])^lambda[l1] - 1) / (sigma[l1] * lambda[l1])

  switch(family,
    NO = {
      Rz <- pnorm(z)
      Rpos <- pnorm(1 / (sigma * lambda))
      Rneg <- pnorm(-1 / (sigma * lambda))
    },
    ST = {
      Rz <- pt(z, zeta)
      Rpos <- pt(1 / (sigma * lambda), zeta)
      Rneg <- pt(-1 / (sigma * lambda), zeta)
    },
    LOI = {
      Rz <- plogisI(z)
      Rpos <- plogisI(1 / (sigma * lambda))
      Rneg <- plogisI(-1 / (sigma * lambda))
    },
    LOII = {
      Rz <- plogis(z)
      Rpos <- plogis(1 / (sigma * lambda))
      Rneg <- plogis(-1 / (sigma * lambda))
    },
    SN = {
      Rz <- pnorm((2 / zeta) * sinh(z))
      Rpos <- pnorm((2 / zeta) * sinh(1 / (sigma * lambda)))
      Rneg <- pnorm((2 / zeta) * sinh(-1 / (sigma * lambda)))
    },
    PE = {
      Rz <- gamlss.dist::pPE(z, mu = 0, sigma = 1, nu = zeta)
      Rpos <- gamlss.dist::pPE(1 / (sigma * lambda), mu = 0, sigma = 1, nu = zeta)
      Rneg <- gamlss.dist::pPE(-1 / (sigma * lambda), mu = 0, sigma = 1, nu = zeta)
    },
    HP = {
      Rz <- GeneralizedHyperbolic::phyperb(z, mu = 0, delta = 1, alpha = zeta, beta = 0)
      Rpos <- GeneralizedHyperbolic::phyperb(1 / (sigma * lambda), mu = 0, delta = 1, alpha = zeta, beta = 0)
      Rneg <- GeneralizedHyperbolic::phyperb(-1 / (sigma * lambda), mu = 0, delta = 1, alpha = zeta, beta = 0)
    },
    SL = {
      Rz <- pslash(z, zeta = zeta)
      Rpos <- pslash(1 / (sigma * lambda), zeta = zeta)
      Rneg <- pslash(-1 / (sigma * lambda), zeta = zeta)
    },
    stop(gettextf("%s family not recognised", sQuote(family)), domain = NA)
  )


  id0 <- which(q > 0 & mu > 0 & sigma > 0 & lambda == 0)
  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda < 0)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda > 0)

  cdf <- rep(NaN, length.out = maxl)
  cdf[id0] <- Rz[id0]
  cdf[id1] <- Rz[id1] / Rneg[id1]
  cdf[id2] <- (Rz[id2] - Rneg[id2]) / Rpos[id2]

  cdf[which(q <= 0 & mu > 0 & sigma > 0)] <- 0

  if (!lower.tail == TRUE) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Não mexi a partir daqui

#' @rdname bcs
#' @export
qbcs <- function(p, mu, sigma, lambda, zeta = NULL, family, lower.tail = TRUE, log.p = FALSE) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(length(p), length(mu), length(sigma), length(lambda))

  p <- rep_len(p, maxl)
  mu <- rep_len(mu, maxl)
  sigma <- rep_len(sigma, maxl)
  lambda <- rep_len(lambda, maxl)

  if (any(mu <= 0)) {
    warning(paste("mu must be positive ", "\n", ""))
  }
  if (any(sigma <= 0)) {
    warning(paste("sigma must be positive", "\n", ""))
  }
  if (any(zeta <= 0)) {
    warning(paste("zeta must be positive", "\n", ""))
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  l0 <- which(lambda <= 0)
  l1 <- which(lambda > 0)

  qtf <- z_p <- rep(NaN, maxl)

  switch(family,
    NO = {
      z_p[l0] <- qnorm(p[l0] * pnorm(1 / (sigma[l0] * abs(lambda[l0]))))
      z_p[l1] <- qnorm(1 - (1 - p[l1]) * pnorm(1 / (sigma[l1] * abs(lambda[l1]))))
    },
    ST = {
      z_p[l0] <- qt(p[l0] * pt(1 / (sigma[l0] * abs(lambda[l0])), zeta), zeta)
      z_p[l1] <- qt(1 - (1 - p[l1]) * pt(1 / (sigma * abs(lambda[l1])), zeta), zeta)
    },
    LOI = {
      z_p[l0] <- qlogisI(p[l0] * plogisI(1 / (sigma[l0] * abs(lambda[l0]))))
      z_p[l1] <- qlogisI(1 - (1 - p[l1]) * plogisI(1 / (sigma[l1] * abs(lambda[l1]))))
    },
    LOII = {
      z_p[l0] <- qlogis(p[l0] * plogis(1 / (sigma[l0] * abs(lambda[l0]))))
      z_p[l1] <- qlogis(1 - (1 - p[l1]) * plogis(1 / (sigma[l1] * abs(lambda[l1]))))
    },
    SN = {
      z_p[l0] <- asinh((zeta / 2) * qnorm(p[l0] * pnorm((2 / zeta) * sinh(1 / (sigma[l0] * abs(lambda[l0]))))))
      z_p[l1] <- asinh((zeta / 2) * qnorm(1 - (1 - p[l1]) * pnorm((2 / zeta) * sinh(1 / (sigma[l1] * abs(lambda[l1]))))))
    },
    PE = {
      z_p[l0] <- gamlss.dist::qPE(p[l0] * gamlss.dist::pPE(1 / (sigma[l0] * abs(lambda[l0])),
        mu = 0, sigma = 1, nu = zeta
      ), mu = 0, sigma = 1, nu = zeta)
      z_p[l1] <- gamlss.dist::qPE(1 - (1 - p[l1]) * gamlss.dist::pPE(1 / (sigma[l1] * abs(lambda[l1])),
        mu = 0, sigma = 1, nu = zeta
      ), mu = 0, sigma = 1, nu = zeta)
    },
    HP = {
      z_p[l0] <- GeneralizedHyperbolic::qhyperb(
        p[l0] * GeneralizedHyperbolic::phyperb(1 / (sigma[l0] * abs(lambda[l0])),
          mu = 0, delta = 1, alpha = zeta, beta = 0
        ),
        mu = 0, delta = 1, alpha = zeta, beta = 0
      )
      z_p[l1] <- GeneralizedHyperbolic::qhyperb(
        1 - (1 - p[l1]) * GeneralizedHyperbolic::phyperb(1 / (sigma[l1] * abs(lambda[l1])),
          mu = 0, delta = 1, alpha = zeta, beta = 0
        ),
        mu = 0, delta = 1, alpha = zeta, beta = 0
      )
    },
    SL = {
      z_p[l0] <- qslash(p[l0] * pslash(1 / (sigma[l0] * abs(lambda[l0])), zeta = zeta), zeta = zeta)
      z_p[l1] <- qslash(1 - (1 - p[l1]) * pslash(1 / (sigma[l1] * abs(lambda[l1])), zeta = zeta), zeta = zeta)
    },
    stop(gettextf("%s family not recognised", sQuote(family)), domain = NA)
  )


  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda != 0)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda == 0)
  id3 <- which(p == 0 & mu > 0 & sigma > 0)
  id4 <- which(p == 1 & mu > 0 & sigma > 0)

  qtf[id1] <- exp(log(mu[id1]) + (1 / lambda[id1]) * log1p(sigma[id1] * lambda[id1] * z_p[id1]))
  qtf[id2] <- exp(log(mu[id2]) + sigma[id2] * z_p[id2])
  qtf[id3] <- 0
  qtf[id4] <- Inf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

#' @rdname bcs
#' @export
rbcs <- function(n, mu, sigma, lambda, zeta = NULL, family) {
  if (any(mu <= 0)) {
    stop(paste("mu must be positive ", "\n", ""))
  }
  if (any(sigma <= 0)) {
    stop(paste("sigma must be positive", "\n", ""))
  }
  if (any(zeta <= 0)) {
    stop(paste("zeta must be positive ", "\n", ""))
  }
  if (any(n <= 0)) {
    stop(paste("n must be a positive integer", "\n", ""))
  }

  n <- ceiling(n)
  p <- runif(n)
  qbcs(p, mu = mu, sigma = sigma, lambda = lambda, zeta = zeta, family = family)
}



# #' @keywords internal
# pslash <- function(q, nu){
#   aux <-  function(x){
#     s_aux <- x^2/2
#     beta_aux <- nu + (1/2)
#     gama_aux <- zipfR::Igamma(beta_aux, s_aux)
#     r <- (nu/sqrt(2*pi))*(1/(s_aux^beta_aux))*zipfR::Igamma(beta_aux, s_aux)
#     return(r)
#   }
#
#   acumu_aux <- function(q) stats::integrate(aux, lower=-Inf, upper=q)$value
#   v.acumu_aux <- Vectorize(acumu_aux)
#   cdf <- v.acumu_aux(q)
#   cdf
# }
#
# #' @keywords internal
# qslash <- function(p, nu){
#
#   qtf <- function(input){
#     p <- input
#
#     if (!is.na(p)){
#       obj <- function(q){
#         pslash(q, nu) - p
#       }
#
#       nleqslv::nleqslv(stats::qnorm(p) / (0.5^(1/nu)), obj)$x
#     }else{
#       numeric(0)
#     }
#   }
#
#   q0 <- rep(0, length(p[p == 0.5]) )
#   q <- vector()
#   if (length(p[p != 0.5]) < 1){
#     q <- numeric(0)
#   } else{
#     q[p != 0.5 & p != 0 & p != 1] <-
#       as.numeric(apply(matrix(p[p != 0.5 & p != 0 & p != 1], ncol = 1),
#                        1, qtf))
#   }
#
#   q[p == 0] <- -Inf
#   q[p == 1] <- Inf
#   qtf <- c(q0, q)
#   index <- c(which(p == 0.5), which(p != 0.5))
#
#   qtf[sort(index, index.return = T)$ix]
# }


# The Slash distribution ----------------------------------------------------------------------

## Inverse gamma function
ig <- function(a, x) pmin(exp(lgamma(a) + stats::pgamma(x, a, scale = 1, log.p = TRUE)), .Machine$double.xmax)

## Probability density function
dslash <- function(x, mu = 0, sigma = 1, zeta, log = FALSE) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(length(x), length(mu), length(sigma), length(zeta))

  x <- rep_len(x, maxl)
  mu <- rep_len(mu, maxl)
  sigma <- rep_len(sigma, maxl)
  zeta <- rep_len(zeta, maxl)

  pmf <- rep_len(-Inf, maxl)

  # NaN index
  pmf[which(sigma <= 0 | zeta <= 0)] <- NaN

  u <- ((x - mu) / sigma)^2

  id1 <- which(u > 0 & !is.nan(pmf), arr.ind = TRUE)
  id2 <- which(u == 0 & !is.nan(pmf), arr.ind = TRUE)

  pmf[id1] <- log(ig(zeta[id1] + 0.5, u[id1] / 2)) + log(zeta[id1]) + zeta[id1] * log(2) -
    0.5 * log(pi) - (zeta[id1] + 0.5) * log(u[id1])
  pmf[id2] <- log(2 * zeta[id2]) - log(2 * zeta[id2] + 1) - 0.5 * log(2 * pi)

  if (!log) pmf <- exp(pmf)
  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Cumulative distribution function
pslash <- function(q, zeta, log.p = FALSE) {
  mu <- 0
  sigma <- 1

  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(length(q), length(mu), length(sigma), length(zeta))

  q <- rep_len(q, maxl)
  mu <- rep_len(mu, maxl)
  sigma <- rep_len(sigma, maxl)

  cdf <- rep_len(0, maxl)

  # NaN index
  cdf[which(sigma <= 0 | zeta <= 0)] <- NaN

  # Positive density index
  id1 <- which(is.finite(q) & q != mu & !is.nan(cdf))
  id2 <- which(q == mu & !is.nan(cdf))
  id3 <- which(q == -Inf)
  id4 <- which(q == Inf)

  # Constructing the Slash distribution
  W <- distr::AbscontDistribution(
    d = function(x) dslash(x, zeta = zeta),
    Symmetry = distr::SphericalSymmetry(0)
  )

  cdf[id1] <- distr::p(W)((q[id1] - mu[id1]) / sigma[id1])
  cdf[id2] <- 0.5
  cdf[id3] <- 0
  cdf[id4] <- 1


  if (log.p) cdf <- log(cdf)

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

# Quantile function
qslash <- function(p, zeta) {
  mu <- 0
  sigma <- 1

  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(length(p), length(mu), length(sigma), length(zeta))

  p <- rep_len(p, maxl)
  mu <- rep_len(mu, maxl)
  sigma <- rep_len(sigma, maxl)

  qtf <- rep_len(NA, maxl)

  # NaN index
  qtf[which(p < 0 | p > 1 | sigma <= 0 | zeta <= 0)] <- NaN

  # Positive density index
  id1 <- which(p != 0.5 & p > 0 & p < 1 & !is.nan(qtf), arr.ind = TRUE)
  id2 <- which(p == 0 & !is.nan(qtf), arr.ind = TRUE)
  id3 <- which(p == 0.5 & !is.nan(qtf), arr.ind = TRUE)
  id4 <- which(p == 1 & !is.nan(qtf), arr.ind = TRUE)

  # Constructing the Slash distribution
  W <- distr::AbscontDistribution(
    d = function(x) dslash(x, zeta = zeta),
    Symmetry = distr::SphericalSymmetry(0)
  )

  qtf[id1] <- distr::q(W)(p[id1])
  qtf[id2] <- -Inf
  qtf[id3] <- 0
  qtf[id4] <- Inf

  qtf <- mu + sigma * qtf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

## VERIFICAÇÃO VISUAL DE QUE A SLASH ESTÁ OK

## SLASH = mu + sigma * Z / sqrt(W), W ~ Beta(zeta, 1), Z ~ N(0, 1)

# # Random generation
# rslash <- function(n, mu = 0, sigma = 1, zeta) {
#   mu + sigma * stats::rnorm(n) / sqrt(stats::rbeta(n, zeta, 1))
# }
#
# x <- rslash(10000, zeta = 2)
#
# plot(density(x))
# curve(dslash(x, zeta = 2), add = TRUE, col = "blue")
#
# plot(ecdf(x))
# curve(pslash(x, zeta = 2), add = TRUE, col = "blue")
#
# plot(seq(0.001, 0.999, 0.001), qslash(seq(0.001, 0.999, 0.001), zeta = 2),
#      type = "l", xlab = "p", ylab = "Quantile")
# curve(qslash(x, zeta = 2), add = TRUE, col = "blue")


# The logistic type I distribution ------------------------------------------------------------

## Density function
dlogisI <- function(x, mu = 0, sigma = 1, log = FALSE) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(length(x), length(mu), length(sigma))

  x <- rep_len(x, maxl)
  mu <- rep_len(mu, maxl)
  sigma <- rep_len(sigma, maxl)

  pmf <- rep_len(-Inf, maxl)

  # NaN index
  pmf[which(sigma <= 0)] <- NaN

  u <- ((x - mu) / sigma)^2

  id <- which(u >= 0 & !is.nan(pmf), arr.ind = TRUE)

  const <- 1.484300029
  pmf[id] <- log(const) - u[id] - log((1 + exp(-u[id]))^2)

  if (!log) pmf <- exp(pmf)
  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

Wloi <- distr::AbscontDistribution(
  d = function(x) dlogisI(x),
  Symmetry = distr::SphericalSymmetry(0)
)

## Cumulative distribution function
plogisI <- function(q, log.p = FALSE) {
  mu <- 0
  sigma <- 1

  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(length(q), length(mu), length(sigma))

  q <- rep_len(q, maxl)
  mu <- rep_len(mu, maxl)
  sigma <- rep_len(sigma, maxl)

  cdf <- rep_len(0, maxl)

  # NaN index
  cdf[which(sigma <= 0)] <- NaN

  # Positive density index
  id1 <- which(is.finite(q) & !is.nan(cdf))
  id2 <- which(q == -Inf)
  id3 <- which(q == Inf)

  cdf[id1] <- distr::p(Wloi)((q[id1] - mu[id1]) / sigma[id1])
  cdf[id2] <- 0
  cdf[id3] <- 1

  if (log.p) cdf <- log(cdf)

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
qlogisI <- function(p) {
  mu <- 0
  sigma <- 1

  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(length(p), length(mu), length(sigma))

  p <- rep_len(p, maxl)
  mu <- rep_len(mu, maxl)
  sigma <- rep_len(sigma, maxl)

  qtf <- rep_len(NA, maxl)

  # NaN index
  qtf[which(p < 0 | p > 1 | sigma <= 0)] <- NaN

  # Positive density index
  id1 <- which(p != 0.5 & p > 0 & p < 1 & !is.nan(qtf), arr.ind = TRUE)
  id2 <- which(p == 0 & !is.nan(qtf), arr.ind = TRUE)
  id3 <- which(p == 1 & !is.nan(qtf), arr.ind = TRUE)

  qtf[id1] <- distr::q(Wloi)(p[id1])
  qtf[id2] <- -Inf
  qtf[id3] <- Inf

  qtf <- mu + sigma * qtf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# ## VERIFICAÇÃO VISUAL DA LOGISTICA TIPO i
# u <- runif(10000)
# x <- qlogisI(u)
#
# plot(density(x))
# curve(dlogisI(x), add = TRUE, col = "blue")
#
# plot(ecdf(x))
# curve(plogisI(x), add = TRUE, col = "blue")
#
# plot(seq(0.001, 0.999, 0.001), qlogisI(seq(0.001, 0.999, 0.001)),
#      type = "l", xlab = "p", ylab = "Quantile")
# curve(qlogisI(x), add = TRUE, col = "blue")
