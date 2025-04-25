#' The Zero-Adjusted Box-Cox Symmetric Distributions
#'
#' Density, distribution function, quantile function, and random generation
#'     for the class of the zero-adjusted Box-Cox symmetric (ZABCS) distributions.
#'
#' @name ZABCS
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `n` is a vector, its length is used as the number of
#'     required observations.
#' @param alpha vector of zero-adjusted parameters, with values on (0, 1).
#' @param mu vector of strictly positive scale parameters.
#' @param sigma vector of strictly positive relative dispersion parameters.
#' @param lambda vector of real-valued skewness parameters. If \code{lambda = 0}, the BCS
#'     distribution reduces to the corresponding log-symmetric distribution with parameters
#'     \code{mu} and \code{sigma} (and a possible extra parameter \code{zeta}).
#' @param zeta strictly positive extra parameter. It must be specified with only one value
#'     in cases where the BCS distribution has an extra parameter. See “Details” below.
#' @param family a character that specifies the symmetric generating family of the BCS distribution.
#'     Available options are: \code{"NO"} (default), \code{"ST"}, \code{"LOI"}, \code{"LOII"},
#'     \code{"PE"}, \code{"SN"}, \code{"HP"}, and \code{"SL"}, corresponding to the normal,
#'     Student-\emph{t}, type I logistic, type II logistic, power exponential, sinh-normal,
#'     hyperbolic, and slash distributions, respectively.
#' @param log,log.p logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#'     Default is \code{FALSE}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'    \eqn{P(X \le x)} otherwise, \eqn{P(X > x)}.
#'
#' @details The class of the ZABCS distributions was introduced by Medeiros and
#'     Queiroz (2025) as an extension of the Box-Cox symmetric (BCS) distributions
#'     (Ferrari and Fumes, 2017). The models consists of a broad class of probability
#'     distributions for positive continuous data which may include zeros.
#'
#'     Let \eqn{Y} be a positive continuous random variable with a ZABCS distribution
#'     with parameters \eqn{\alpha \in (0, 1)}, \eqn{\mu > 0}, \eqn{\sigma > 0},
#'     and \eqn{\lambda \in \mathbb{R}} and density generating function \eqn{r}.
#'     The probability density function of
#'     \eqn{Y} is given by
#'
#'     \eqn{
#'       f^{(0)}(y; \alpha, \mu, \sigma, \lambda)  = \left\{
#'         \begin{array}{rcl}
#'         \alpha, &  y=0,\\
#'         (1 - \alpha)f(y; \alpha, \mu, \sigma, \lambda), & y > 0,\\
#'         \end{array}
#'         \right.
#'     }
#'
#'    where
#'
#'     \eqn{
#'     f(y; \mu, \sigma, \lambda) = \left\{\begin{array}{ll}
#'       \dfrac{y^{\lambda-1}}{\mu^\lambda \sigma} \dfrac{r(z^2)}{R\left(\frac{1}{\sigma |\lambda|}\right)}, & \mbox{ if } \lambda \neq 0,\\
#'       \dfrac{1}{y\sigma} r(z^2), & \mbox{ if } \lambda = 0,
#'       \end{array}\right., \quad y > 0,
#'     }
#'
#'     with
#'
#'     \eqn{
#'     z = \left\{
#'     \begin{array}{ll}
#'     \dfrac{1}{\sigma \lambda} \left\{\left(\frac{y}{\mu}\right)^\lambda - 1 \right\}, & \mbox{ if } \lambda \neq 0, \\
#'     \dfrac{1}{\sigma} \log\left(\frac{y}{\mu}\right), & \mbox{ if } \lambda = 0,
#'     \end{array}
#'     \right.
#'     }
#'
#'     \eqn{r:[0,\infty) \longrightarrow [0, \infty)}
#'     satisfies \eqn{\int_0^\infty u^{-1/2}r(u)\textrm{d} u = 1}, and
#'     \eqn{R(x) = \int_{-\infty}^x r(u^2)\textrm{d} u, x \in \mathbb{R}}.
#'
#'     The function \eqn{r} is called density generating function, and it specifies the
#'     generating symmetric family of \eqn{Y} within the class of the ZABCS probability
#'     models. This function can also depend on extra parameters, such as the zero-adjusted
#'     Box-Cox \emph{t} and zero-adjusted Box-Cox power exponential distributions. We call
#'     these extra parameters \code{zeta}. The currently available ZABCS distributions in the
#'     \code{BCSreg} package are listed below:
#'     \tabular{llc}{
#'     \bold{Distribution}  \tab \bold{Family abbreviation} \tab \bold{N. of extra parameters}\cr
#'     Zero-adjusted Box-Cox Hyperbolic  \tab \code{"HP"}      \tab  1  \cr
#'     Zero-adjusted Box-Cox Type I Logistic  \tab \code{"LOI"}      \tab  0  \cr
#'     Zero-adjusted Box-Cox Type II Logistic  \tab \code{"LOII"}      \tab  0  \cr
#'     Zero-adjusted Box-Cox Normal  \tab \code{"NO"}      \tab  0  \cr
#'     Zero-adjusted Box-Cox Power Exponential  \tab \code{"PE"}      \tab  1  \cr
#'     Zero-adjusted Box-Cox Sinh-Normal  \tab \code{"SN"}      \tab  1  \cr
#'     Zero-adjusted Box-Cox Slash  \tab \code{"SL"}      \tab  1  \cr
#'     Zero-adjusted Box-Cox \emph{t}  \tab \code{"ST"}      \tab  1  \cr
#'     }
#'
#' @return
#' \code{dZABCS} returns the density function, \code{pZABCS} gives the cumulative distribution function,
#' \code{qZABCS} provides the quantile function, and \code{rZABCS} generates random variables.
#'
#' @references
#'  Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.
#'
#'  Medeiros, R. M. R., and Queiroz, F. F. (2025). Modeling positive continuous data:
#'      Box-Cox symmetric regression models and their extensions.
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @seealso \code{\link{BCS}} to access the density function, distribution
#'     function, quantile function, and a random number generator for the BCS
#'     distributions. \code{\link{BCSreg}} for estimating the parameters of a
#'     ZABCS regression model.
#'
#' @examples
#' # Probability density function
#'
#' ## Right-skewed distributions
#' curve(dZABCS(x, 0.4, 3, 0.3, -1.5, family = "NO"), from = 0.001, to = 7,
#'       xlim = c(0, 7), ylim = c(0, 0.5), ylab = "Density")
#' curve(dZABCS(x, 0.4, 3, 0.3, -1.5, 4, family = "ST"), add = TRUE, col = 2, from = 0.001)
#' curve(dZABCS(x, 0.4, 3, 0.3, -1.5, 5, family = "PE"), add = TRUE, col = 4, from = 0.001)
#' points(0, 0.4, type = "h", lty = 2)
#' points(0, 0.4, pch = 16, lty = 2)
#' legend("topright", legend = c("BCNO", "BCT", "BCPE"), lty = 1, col = c(1, 2, 4))
#'
#' ## Truncated symmetric distributions (with support on (0, Inf))
#' curve(dZABCS(x, 0.4, 3, 0.3, 1, family = "NO"), from = 0.001, to = 7,
#'       xlim = c(0, 7), ylim = c(0, 0.5), ylab = "Density")
#' curve(dZABCS(x, 0.4, 3, 0.3, 1, 4, family = "ST"), add = TRUE, col = 2, from = 0.001)
#' curve(dZABCS(x, 0.4, 3, 0.3, 1, 5, family = "PE"), add = TRUE, col = 4, from = 0.001)
#' points(0, 0.4, type = "h", lty = 2)
#' points(0, 0.4, pch = 16, lty = 2)
#' legend("topright", legend = c("BCNO", "BCT", "BCPE"), lty = 1, col = c(1, 2, 4))
#'
#' ## Left-skewed distributions
#' curve(dZABCS(x, 0.4, 3, 0.3, 3, family = "NO"), from = 0.001, to = 7,
#'       xlim = c(0, 7), ylim = c(0, 0.5), ylab = "Density")
#' curve(dZABCS(x, 0.4, 3, 0.3, 3, 4, family = "ST"), add = TRUE, col = 2, from = 0.001)
#' curve(dZABCS(x, 0.4, 3, 0.3, 3, 5, family = "PE"), add = TRUE, col = 4, from = 0.001)
#' points(0, 0.4, type = "h", lty = 2)
#' points(0, 0.4, pch = 16, lty = 2)
#' legend("topright", legend = c("BCNO", "BCT", "BCPE"), lty = 1, col = c(1, 2, 4))
#'
#'
#' # Random generation
#'
#' ## Parameter setting
#' alpha <- 0.2   # zero-adjustment parameter
#' mu <- 5        # scale parameter
#' sigma <- 0.2   # relative dispersion parameter
#' lambda <- -0.2 # skewness parameter
#'
#' ## Generating family
#' family <- "NO"
#'
#' ## Visualization
#' x <- rZABCS(10000, alpha, mu, sigma, lambda, family = family)
#'
#' hist(x, prob = TRUE, col = "white", main = "")
#'
#' points(0, mean(x == 0), type = "h", lty = 2)
#' points(0, mean(x == 0), pch = 16, lty = 2)
#' curve(dZABCS(x, alpha, mu, sigma, lambda, zeta, family = family), col = "blue", add = TRUE)
#'
#' plot(ecdf(x), main = "")
#' curve(pZABCS(x, alpha, mu, sigma, lambda, zeta, family = family), col = "blue", add = TRUE)
#' @importFrom stats dnorm qnorm rnorm pnorm dlogis plogis qlogis runif dt pt qt
#'
#' @export
dZABCS <- function(x, alpha, mu, sigma, lambda, zeta, family = "NO", log = FALSE) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(lambda)), length(alpha))

  x <- rep(x, length.out = maxl)
  alpha <- rep(alpha, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  if (any(alpha < 0 | alpha > 1)) {
    warning(paste("alpha must taken values on [0, 1] ", "\n", ""))
  }
  if (family %in% c("HP", "PE", "SN", "SL", "ST")) {
    if (missing(zeta))
      stop(gettextf("%s family has an extra parameter that must be specified via the 'zeta' argument",
                    sQuote(family)), domain = NA)

    if (zeta < 0)
      warning(paste("zeta must be positive", "\n", ""))
  } else {
    zeta <- NULL
  }

  log.lik <- rep(-Inf, maxl)

  # NaN index
  log.lik[which(alpha < 0 | alpha >  1)] <- NaN

  # Positive density index
  id0 <- which(x == 0 & !is.nan(log.lik))
  id1 <- which(x > 0 & !is.nan(log.lik))

  log.lik[id0] <- log(alpha[id0])
  log.lik[id1] <- log(1 - alpha[id1]) +
    dBCS(x = x[id1], mu = mu[id1], sigma = sigma[id1], lambda = lambda[id1],
         zeta = zeta, family = family, log = TRUE)

  if (!log) log.lik <- exp(log.lik)
  if (d > 1L) matrix(log.lik, ncol = d) else log.lik
}

#' @rdname ZABCS
#' @export
pZABCS <- function(q, alpha, mu, sigma, lambda, zeta, family = "NO", lower.tail = TRUE, log.p = FALSE) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(lambda)), length(alpha))

  q <- rep(q, length.out = maxl)
  alpha <- rep(alpha, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  if (any(alpha < 0 | alpha > 1)) {
    warning(paste("alpha must taken values on [0, 1] ", "\n", ""))
  }
  if (family %in% c("HP", "PE", "SN", "SL", "ST")) {
    if (missing(zeta))
      stop(gettextf("%s family has an extra parameter that must be specified via the 'zeta' argument",
                    sQuote(family)), domain = NA)

    if (zeta < 0)
      warning(paste("zeta must be positive", "\n", ""))
  } else {
    zeta <- NULL
  }

  id0 <- which(q < 0 & mu > 0 & sigma > 0 & alpha >= 0 & alpha <= 1)
  id1 <- which(q == 0 & mu > 0 & sigma > 0 & alpha >= 0 & alpha <= 1)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & alpha >= 0 & alpha <= 1)

  cdf <- rep(NaN, length.out = maxl)
  cdf[id0] <- 0
  cdf[id1] <- alpha[id1]
  cdf[id2] <- alpha[id2] +
    (1 - alpha[id2]) * pBCS(q = q[id2], mu = mu[id2], sigma = sigma[id2], lambda = lambda[id2],
         zeta = zeta, family = family)

  if (!lower.tail == TRUE) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

#' @rdname ZABCS
#' @export
qZABCS <- function(p, alpha, mu, sigma, lambda, zeta, family = "NO", lower.tail = TRUE, log.p = FALSE) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(length(p), length(mu), length(sigma), length(lambda), length(alpha))

  p <- rep_len(p, maxl)
  alpha <- rep_len(alpha, maxl)
  mu <- rep_len(mu, maxl)
  sigma <- rep_len(sigma, maxl)
  lambda <- rep_len(lambda, maxl)

  if (any(alpha < 0 | alpha > 1)) {
    warning(paste("alpha must taken values on [0, 1] ", "\n", ""))
  }
  if (family %in% c("HP", "PE", "SN", "SL", "ST")) {
    if (missing(zeta))
      stop(gettextf("%s family has an extra parameter that must be specified via the 'zeta' argument",
                    sQuote(family)), domain = NA)

    if (zeta < 0)
      warning(paste("zeta must be positive", "\n", ""))
  } else {
    zeta <- NULL
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  qtf <- rep(NaN, maxl)

  # Quantile function
  id0 <- which(p <= alpha)
  id1 <- which(p > alpha)

  qtf[id0] <- 0
  qtf[id1] <- qBCS(p = (p[id1] - alpha[id1]) / (1 - alpha[id1]),
                   mu = mu[id1], sigma = sigma[id1], lambda = lambda[id1],
                   zeta = zeta, family = family)

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

#' @rdname ZABCS
#' @export
rZABCS <- function(n, alpha, mu, sigma, lambda, zeta, family = "NO") {

  if (any(alpha < 0 | alpha > 1)) {
    stop(paste("alpha must taken values on [0, 1] ", "\n", ""))
  }

  n <- ceiling(n)
  p <- runif(n)
  qZABCS(p, alpha = alpha, mu = mu, sigma = sigma, lambda = lambda, zeta = zeta, family = family)
}
