% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_the_BCS_distributions.R
\name{bcs}
\alias{bcs}
\alias{dbcs}
\alias{pbcs}
\alias{qbcs}
\alias{rbcs}
\title{The Box-Cox Symmetric (BCS) Distributions}
\usage{
dbcs(x, mu, sigma, lambda, zeta = NULL, family, log = FALSE)

pbcs(
  q,
  mu,
  sigma,
  lambda,
  zeta = NULL,
  family,
  lower.tail = TRUE,
  log.p = FALSE
)

qbcs(
  p,
  mu,
  sigma,
  lambda,
  zeta = NULL,
  family,
  lower.tail = TRUE,
  log.p = FALSE
)

rbcs(n, mu, sigma, lambda, zeta = NULL, family)
}
\arguments{
\item{x, q}{vector of positive quantiles.}

\item{mu}{vector of strictly positive scale parameters.}

\item{sigma}{vector of strictly positive relative dispersion parameters.}

\item{lambda}{vector of real-valued skewness parameters. If \code{lambda = 0}, the BCS
distribution reduces to the corresponding log-symmetric distribution with parameters
\code{mu}, \code{sigma}, and \code{zeta}.}

\item{zeta}{strictly positive extra parameter. It must be specified with only one value.}

\item{family}{a character that specifies the generating family of the BCS distribution.
Available options are: \code{"NO"}, \code{"ST"}, \code{"LOI"}, \code{"LOII"},
\code{"PE"}, \code{"SN"}, \code{"HP"}, and \code{"SL"}, corresponding to the normal,
Student-\emph{t}, type I logistic, type II logistic, power exponential, sinh-normal,
hyperbolic, and slash distributions, respectively.}

\item{log, log.p}{logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
Default is \code{FALSE}.}

\item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
\eqn{P(X \le x)} otherwise, \eqn{P(X > x)}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If `n` is a vector, its length is used as the number of
required observations.}
}
\value{
\code{dbcs} returns the density function, \code{pbcs} gives the cumulative distribution function,
\code{qbcs} provides the quantile function, and \code{rbcs} generates random variables.
}
\description{
Density, distribution function, quantile function, and random generation
    for the class of the Box-Cox symmetric (BCS) distributions.
}
\details{
The class of the BCS distributions was introduced by Ferrari and
    Fumes (2017). It consists of a broad class of probability models for
    positive continuous data, which includes flexible distributions with
    different levels of skewness and tail-heaviness.



    The BCS class includes, as special cases, the Box-Cox \emph{t} (Rigby and Stasinopoulos, 2006),
    Box-Cox normal (or Box-Cox Cole-Green; Cole and Green, 1992), Box-Cox power exponential
    (Rigby and Stasinopoulos, 2004) distributions, as well as the log-symmetric
    distributions (Vanegas and Paula, 2016).

    The currently available BCS distributions in the \code{BCSreg} package are listed below:

\tabular{llc}{
 \bold{Distribution}  \tab \bold{Family abbreviation} \tab \bold{Number of parameters}\cr
 Box-Cox Hyperbolic  \tab \code{"HP"}      \tab  4  \cr
 Box-Cox Type I Logistic  \tab \code{"LOI"}      \tab  3  \cr
 Box-Cox Type II Logistic  \tab \code{"LOII"}      \tab  3  \cr
 Box-Cox Normal  \tab \code{"NO"}      \tab  3  \cr
 Box-Cox Power Exponential  \tab \code{"PE"}      \tab  4  \cr
 Box-Cox Sinh-Normal  \tab \code{"SN"}      \tab  4  \cr
 Box-Cox Slash  \tab \code{"SL"}      \tab  4  \cr
 Box-Cox \emph{t}  \tab \code{"ST"}      \tab  4  \cr
 }
}
\examples{
# Density

## Right-skewed distributions
curve(dbcs(x, 3, 0.3, -1.5, family = "NO"), xlim = c(0, 7), ylim = c(0, 0.7), ylab = "Density")
curve(dbcs(x, 3, 0.3, -1.5, 4, family = "ST"), add = TRUE, col = 2)
curve(dbcs(x, 3, 0.3, -1.5, 5, family = "PE"), add = TRUE, col = 4)
legend("topright", legend = c("NO", "ST", "PE"), lty = 1, col = c(1, 2, 4))

## Truncated symmetric distributions (with support on (0, Inf))
curve(dbcs(x, 3, 0.3, 1, family = "NO"), xlim = c(0, 7), ylim = c(0, 0.7), ylab = "Density")
curve(dbcs(x, 3, 0.3, 1, 4, family = "ST"), add = TRUE, col = 2)
curve(dbcs(x, 3, 0.3, 1, 5, family = "PE"), add = TRUE, col = 4)
legend("topright", legend = c("NO", "ST", "PE"), lty = 1, col = c(1, 2, 4))

## Left-skewed distributions
curve(dbcs(x, 3, 0.3, 3, family = "NO"), xlim = c(0, 7), ylim = c(0, 0.7), ylab = "Density")
curve(dbcs(x, 3, 0.3, 3, 4, family = "ST"), add = TRUE, col = 2)
curve(dbcs(x, 3, 0.3, 3, 5, family = "PE"), add = TRUE, col = 4)
legend("topright", legend = c("NO", "ST", "PE"), lty = 1, col = c(1, 2, 4))

# Random generation

## Parameter setting
mu <- 5 # scale parameter
sigma <- 0.2 # relative dispersion parameter
lambda <- -2 # skewness parameter
zeta <- 6 # extra parameter (if necessary)

## Generating family
family <- "NO"

## Visualization
x <- rbcs(10000, mu, sigma, lambda, zeta, family = family)

hist(x, prob = TRUE, col = "white", main = "")
curve(dbcs(x, mu, sigma, lambda, zeta, family = family), col = "blue", add = TRUE)

plot(ecdf(x), main = "")
curve(pbcs(x, mu, sigma, lambda, zeta, family = family), col = "blue", add = TRUE) #'

}
\references{
Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical
    properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.

Cole, T., and Green, P.J. (1992). Smoothing reference centile curves: the LMS
     method and penalized likelihood. \emph{Statistics in medicine}, 11, 1305-1319.

 Ferrari, S. L., and Fumes, G. (2017). Box-Cox symmetric distributions and
     applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.

 Rigby, R. A., and Stasinopoulos, D. M. (2004). Smooth centile curves for skew
     and kurtotic data modelled using the Box-Cox power exponential
     distribution. \emph{Statistics in medicine}, 23, 3053-3076.

 Rigby, R. A., and Stasinopoulos, D. M. (2006). Using the Box-Cox t
     distribution in GAMLSS to model skewness and kurtosis. \emph{Statistical Modelling}, 6, 209-229.

 Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions:
     statistical properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
}
