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
#'       the soak time (in hours).}
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

#' Household Expenditures on Basic Education (POF 2017--2018, São Paulo)
#'
#' @description
#'
#' A dataset containing household-level information on expenditures on basic
#' education and associated socioeconomic characteristics, derived from the
#' 2017--2018 Brazilian Consumer Expenditure Survey (Pesquisa de Orçamentos
#' Familiares, POF), conducted by the Instituto Brasileiro de Geografia e
#' Estatística (IBGE).
#'
#' This dataset is used to illustrate regression modeling for zero-adjusted
#' data with a substantial proportion of zero observations. The response
#' variable corresponds to total household expenditure on basic education,
#' measured over the 12 months preceding the interview and assigned to the
#' household reference person. Expenditures include childcare, preschool,
#' regular primary and secondary education, youth and adult education, and
#' supplementary (equivalency) programs at the primary and secondary levels.
#'
#' Households reporting no expenditure on basic education are assigned a value
#' of zero. The final sample consists of 4,232 households residing in the state
#' of São Paulo, Brazil. Approximately 93\% of the observations correspond to
#' zero expenditure, indicating a highly zero-inflated distribution.
#'
#' @format
#' A data frame with 4,232 observations and 11 variables:
#'
#' \describe{
#'   \item{expense}{Total household expenditure on basic education (in BRL)
#'   over the 12 months preceding the interview. Values equal to zero indicate
#'   no reported expenditure.}
#'
#'   \item{residence_type}{Type of household residence, categorized as
#'   \code{"Urban"} or \code{"Rural"}.}
#'
#'   \item{age}{Age of the household reference person, in years.}
#'
#'   \item{sex}{Sex of the household reference person, coded as
#'   \code{"Male"} or \code{"Female"}.}
#'
#'   \item{race}{Self-reported race or ethnicity of the reference person,
#'   according to IBGE classification.}
#'
#'   \item{health_plan}{Indicator of whether the reference person is covered by
#'   a private health plan (\code{"Yes"} or \code{"No"}).}
#'
#'   \item{literacy}{Literacy status of the reference person
#'   (\code{"Yes"} or \code{"No"}).}
#'
#'   \item{years_schooling}{Number of completed years of formal education of
#'   the reference person.}
#'
#'   \item{education_level}{Highest educational level attained by the reference
#'   person.}
#'
#'   \item{income_pc}{Per capita disposable household income (in BRL),
#'   calculated as total disposable household income divided by the number of
#'   residents. Disposable income includes monetary and non-monetary earnings,
#'   net of direct taxes, social contributions, and other mandatory deductions.}
#'
#'   \item{sons}{Number of children living in the household, including children
#'   of the reference person and/or the spouse.}
#' }
#'
#' @details
#' In each household, the POF survey designates a reference person, typically
#' responsible for financial and administrative decisions. All individual-level
#' covariates refer to this reference person. Expenditure values are aggregated
#' at the household level but attributed to the reference person for modeling
#' purposes.
#'
#' Summary statistics indicate that the median education expenditure is zero,
#' while the mean expenditure is BRL 108.8, reflecting the presence of a small
#' number of households with substantially higher spending levels. The maximum
#' observed expenditure is BRL 48,000.
#'
#' @source
#' Instituto Brasileiro de Geografia e Estatística (IBGE).
#' Pesquisa de Orçamentos Familiares (POF) 2017--2018.
#' Available (in Portuguese) at:
#' \url{https://www.ibge.gov.br/estatisticas/sociais/populacao/24786-pesquisa-de-orcamentos-familiares-2.html}
#'
#'
#' @usage
#' data(pof_education_sp)
#'
#' @keywords datasets education expenditure zero-inflation Brazil
"education"

