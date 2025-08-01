#' Fishing Data from Todos-os-Santos Bay
#'
#' This dataset contains information on white ray landings using the traditional
#' "grozeira" fishing method in Baía de Todos os Santos, Brazil. Data were collected
#' from January 2012 to January 2013 by Marion (2015). The study aimed
#' to analyze the relationship between catch per unit effort (cpue) and various
#' environmental factors.
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
#' @usage data(raycatch)
#'
#' @details
#' In Brazil, marine fish are harvested through various fishing activities, with artisanal fishing
#' playing a particularly prominent role, especially in the state of Bahia, in the northeastern
#' region of the country. In this region, particularly within the Baía de Todos os Santos, most of
#' the marine production comes from small-scale operations. Among the species landed and sold,
#' rays constitute a significant portion of the catch. The data analyzed in this application were
#' collected by Marion (2015) over a 13-month period (January 2012 to January 2013), based on 231
#' fishing trips that employed the grozeira (a type of bottom longline) as gear. After eliminating
#' missing data, the data contains information on 186 fish landings. The objective is to identify
#' key factors influencing ray catch levels.  The response variable is the catch per unit effort
#' (\code{cpue}), defined as
#' \eqn{
#' \text{cpue} = \dfrac{\text{Productivity (g)}}{\text{Nº of hooks} \times \text{Immersion time (hours)}}.
#' }
#'
#'
#' @source
#' Data collected from a field study on traditional ray fishing in Baía de Todos os Santos, Brazil.
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
#' data(raycatch)
#' summary(raycatch)
#' plot(cpue ~ tide_phase, raycatch,
#'      xlab = "Tide phase", ylab = "CPUE (kg)",
#'      main = "Effect of tide phase on CPUE")
"raycatch"
