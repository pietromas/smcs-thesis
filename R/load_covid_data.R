################################################################################
## load_covid_data.R -- downloads the Covid-19 Forecast Hub data used in
## Section 1.6 and saves them to data/.
##
## Run once before R/ch1_covid.R. The download takes several minutes and
## requires an internet connection. The resulting .rda files are not tracked in
## this repository.
##
## Usage:  Rscript R/load_covid_data.R
## Output: data/forecasts_death.rda, data/truth.rda
##
## First-time setup. The packages below are not on CRAN, so install them once
## by hand before running this script:
##
##   install.packages(c("remotes", "devtools", "dplyr", "tidyverse",
##                      "rvest", "here", "covidcast"))
##   remotes::install_github("reichlab/zoltr")
##   remotes::install_github("epiforecasts/scoringutils")
##   remotes::install_github("Chicago/RSocrata")
##   remotes::install_github("reichlab/covidData")
##   devtools::install_github("reichlab/covidHubUtils")
################################################################################

library(dplyr)
library(covidData)
library(covidHubUtils)

dir.create("data", showWarnings = FALSE)

## Download forecasts for the US.
forecasts_death <- load_forecasts(
  models = NULL,
  date_window_size = 0,
  locations = "US",
  hub = "US",
  types = c("quantile"),
  targets = paste(1, "wk ahead inc death"),
  source = "zoltar",
  verbose = FALSE,
  as_of = NULL
)

## Download true statistics from Johns Hopkins University.
truth <- load_truth(
  truth_source = "JHU",
  hub = "US",
  locations = "US",
  target_variable = c("inc death"),
  temporal_resolution = "weekly"
)

save(forecasts_death, file = "data/forecasts_death.rda")
save(truth, file = "data/truth.rda")

cat(sprintf("wrote data/forecasts_death.rda (%d rows) and data/truth.rda (%d rows)\n",
            nrow(forecasts_death), nrow(truth)))
