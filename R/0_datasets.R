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


#' Renewable Electricity Output for 186 Countries (2015)
#'
#' This dataset provides information on renewable electricity output and a set
#' of socioeconomic, environmental, and institutional variables for 186 countries
#' in the year 2015. The data were obtained from the World Bank
#' (\url{https://data.worldbank.org/}).
#'
#' @format A data frame with 186 rows and 9 variables:
#'
#' \describe{
#'   \item{\code{country_name}}{country name.}
#'   \item{\code{adj_sav_edu}}{adjusted savings: education expenditure (percent of GNI).}
#'   \item{\code{agri_land}}{natural logarithm of total agricultural land area (in sq. km).}
#'   \item{\code{elec_fossil}}{electricity production from fossil fuels (percent of total).
#'   Indicates the share of electricity generated from oil, gas, and coal, representing reliance on non-renewable energy.}
#'   \item{\code{forest_area}}{forest area (percent of land area). Reflects natural resource availability.}
#'   \item{\code{gov_effec}}{government effectiveness index (range: approximately -2.5 to 2.5). Measures public service quality, policy implementation, and government credibility.}
#'   \item{\code{pop_density}}{natural logarithm of population density (people per sq. km).}
#'   \item{\code{access_elec}}{access to electricity (percent of population).}
#'   \item{\code{renew_elec_output}}{renewable electricity output (in TWh).}
#' }
#'
#' @details
#'
#' The response variable is the renewable electricity output, measured in terawatt-hours (TWh),
#' representing the total electricity generated from renewable sources such as wind, solar
#' photovoltaic, solar thermal, hydro, marine, geothermal, solid biofuels, renewable municipal waste,
#' liquid biofuels, and biogas. Hydro pumped storage is excluded.
#'
#' @source \url{https://data.worldbank.org/}
#' @usage data(renewables2015)
#' @keywords datasets
#' @examples
#' data(renewables2015)
#' summary(renewables2015)
#' pairs(renewables2015[, -1], pch = 16, cex = 0.8, col = rgb(0,0,0, 0.5))
"renewables2015"
