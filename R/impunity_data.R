#' Impunity Index Dataset (2023)
#'
#' @description
#' This dataset contains information for 119 countries in the year 2023,
#' used to study the determinants of the Impunity Index, defined as one minus
#' the Rule of Law Index published by the World Justice Project (WJP).
#'
#' The Impunity Index captures the extent to which a country fails to uphold
#' the rule of law, reflecting weak government accountability, limited access
#' to justice, and other institutional deficiencies. Following Cribari-Neto
#' and Santos (2023), the dependent variable is modeled as a continuous
#' response in the interval (0, 1), making it suitable for the simplex regression
#' model with parametric mean links.
#'
#' @format A data frame with 119 observations and 11 variables:
#' \describe{
#'   \item{Country}{Country name}
#'   \item{Democracy}{Quality of Democracy Index (2020)}
#'   \item{EconomicFreedom}{Economic Freedom Index (2019)}
#'   \item{GDP}{Gross Domestic Product per capita (2018, in thousands of USD)}
#'   \item{Gini}{Gini coefficient (inequality measure, 2018)}
#'   \item{HDI}{Human Development Index (2019)}
#'   \item{HealthSpending}{Public health spending (\% of GDP, 2020)}
#'   \item{Impunity}{Impunity Index (1 - Rule of Law Index, WJP, 2023)}
#'   \item{Press}{Press Freedom Index (2023)}
#'   \item{dnordic}{Dummy variable equal to 1 for Nordic countries (Norway, Denmark, Sweden, Finland), 0 otherwise}
#'   \item{X}{Country identifier}}
#'
#' @source World Justice Project (2023), Freedom House (2023), UNDP (2019), World Bank (2018)
"impunity_dataset"
