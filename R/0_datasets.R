#' Fishing Data from Todos-os-Santos Bay
#'
#' This dataset contains information on fishing landings of the white ray
#' using the traditional "grozeira" method in Todos-os-Santos Bay, Brazil.
#' The data were collected from January 2012 to January 2013, covering 186 fishing
#' operations. The objective of the study was to analyze the relationship
#' between the catch per unit effort (CPUE) and various environmental factors.
#'
#' @format A data frame with 186 rows and 8 variables:
#'
#' \describe{
#'   \item{period}{fishing period, categorized as "Dry" or "Rainy".}
#'   \item{location}{fishing location, categorized as "Area 1", "Area 2", "Area 3", or "Area 4".}
#'   \item{wind_speed}{wind speed (m/s).}
#'   \item{tide_phase}{tide phase, categorized as "Quadrature" or "Spring tide".}
#'   \item{max_temp}{maximum temperature (°C).}
#'   \item{min_temp}{minimum temperature (°C).}
#'   \item{sunshine_duration}{duration of sunshine (hours).}
#'   \item{cpue}{represents the catch per unit effort, calculated as the
#'       productivity (in grams) divided by the product of the number of hooks and
#'       the soak time (in hours). .}
#' }
#'
#' @usage data(fishery)
#'
#' @details
#' The data set can be used to study the impact of environmental factors
#' on the catch efficiency of white ray fishing. It provides valuable information
#' for fisheries management and environmental studies.
#'
#' @source
#' Data collected from a field study on traditional ray fishing in Todos-os-Santos Bay, Brazil.
#'
#' @references
#'  Marion, C. (2015). \emph{Função da Baía de Todos os Santos, Bahia, no ciclo de
#'  vida de Dasyatis guttata (Chondrichthyes: Dasyatidae)}.
#'  Doctoral dissertation, University of São Paulo, Institute of Oceanography.
#'  Retrieved from \url{https://teses.usp.br/teses/disponiveis/21/21134/tde-22072015-154346/pt-br.php}
#'
#'  Paula, G. A., and Kumagaia, G. H. \emph{Relatório de análise estatística sobre
#'  o projeto: "Variabilidade espaço-temporal da captura da raia branca, Dasyattis
#'  guttata, na pesca artesanal da Baía de Todos os Santos, Bahia"}.
#'  São Paulo, IME-USP, 2014.
#'
#' @examples
#' data(fishery)
#' summary(fishery)
#' plot(cpue ~ tide_phase, fishery,
#'      xlab = "Tide phase", ylab = "CPUE (kg)",
#'      main = "Effect of tide phase on CPUE")
"fishery"
