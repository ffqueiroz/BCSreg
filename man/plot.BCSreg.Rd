% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_methods.R
\name{plot.BCSreg}
\alias{plot.BCSreg}
\title{Diagnostic Plots for a Box-Cox Symmetric Regression Fit}
\usage{
\method{plot}{BCSreg}(
  x,
  which = 1:4,
  ask = prod(graphics::par("mfcol")) < length(which) && grDevices::dev.interactive(),
  pch = "+",
  las = 1,
  cex = 0.8,
  lwd = 2,
  ...
)
}
\arguments{
\item{x}{an object of class \code{"BCSreg"}.}

\item{which}{numeric; if a subset of the plots is required, specify a subset
of the numbers \code{1:7}.}

\item{ask}{logical; if \code{TRUE}, the user is asked before each plot.}

\item{pch, las, cex, lwd, ...}{graphical parameters (see \code{\link[graphics]{par}})}
}
\value{
\code{plot} method for \code{"\link{BCSreg}"} objects returns seven types
    of diagnostic plots.
}
\description{
This function provides plots for diagnostic analysis of a Box-Cox symmetric
    or a zero-adjusted regression fit.
}
\details{
The \code{plot} method for \code{\link{BCSreg}} objects provides seven types
    of diagnostic plots in the following order:
    \describe{
        \item{Residuals vs fitted values}{a plot of the residuals
            versus the fitted medians.}
        \item{Residuals vs observation indices.}{an index plot of the residuals
            versus the observation indices.}
        \item{Normal probability plot}{a normal probability plot of the residuals with a
            confidence region constructed according to Fox (2016) using the
            \code{\link[car]{qqPlot}} function.}
        \item{Case-weight perturbation}{An index plot of local influence based on the
            case-weight perturbation scheme.}
        \item{Density plot}{a graph that compares the empirical density of the residuals
            with the density of the standard normal distribution.}
        \item{Fitted vs observed values}{a dispersion diagram of the fitted values
            versus the observed values.}
        \item{Residuals vs v(z) function}{a dispersion diagram of the \eqn{v(z)} function
            versus the residuals. For some BCS models, the \eqn{v(z)} function
            may be interpreted as weights in the estimation process. If \code{family = "NO"},
            the \eqn{v(z)} function is constant.}
     }

     The \code{which} argument can be used to select a subset of the implemented plots.
     Default is \code{which = 1:4}. See \code{\link{residuals.BCSreg}} for details on
     the residuals.
}
\examples{
## Data set: fishery (for description, run ?fishery)
hist(fishery$cpue, xlab = "Catch per unit effort")
plot(cpue ~ tide_phase, fishery, pch = 16,
    xlab = "Tide phase", ylab = "Catch per unit effort")
plot(cpue ~ location, fishery, pch = 16,
    xlab = "Location", ylab = "Catch per unit effort")
plot(cpue ~ max_temp, fishery, pch = 16,
    xlab = "Maximum temperature", ylab = "Catch per unit effort")

## Fit a double Box-Cox normal regression model:
fit <- BCSreg(cpue ~ location + tide_phase |
                location + tide_phase + max_temp, fishery)

## Available plots:

### Residuals vs fitted values (fitted medians)
plot(fit, which = 1)

### Residuals vs observation indices
plot(fit, which = 2)

### Normal probability plot
plot(fit, which = 3)

### Local influence
plot(fit, which = 4)

### Density plot
plot(fit, which = 5)

### Fitted medians vs response
plot(fit, which = 6)

### v(z) function
plot(fit, which = 7)
}
\author{
Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>

Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
}
