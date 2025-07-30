# heading ----
rm(list = ls())

Sys.setenv(LANG = "en")
options(scipen = 999)
message(
  paste0("Date: ", format(Sys.Date(), "%d %B %Y")), "\n",
  paste0("By: ", "Francesco Invernizzi"), "\n",
  paste0("Description: ", "code - master code to upload packages and set directories"), "\n",
  paste0("Version of R used: ", R.Version()$version.string)
)

# install load libraries ----
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse, # load multiple 'tidyverse' packages in a single step
  sf, # provides simple features access for R
  amt, #   manage and analyze animal movement data
  lubridate, # parse and manipulate dates
  circular, # for circular data
  knitr,
  kableExtra #
)

# set folder paths ----
getwd()
dir_data <- paste0(getwd(), "/data")

# print information on session (packages version, etc.) ----
# writeLines(
#   text = capture.output(
#     sessionInfo()
#   ),
#   con = paste0(
#     getwd(),
#     "/code/",
#     "00_log_session_information",
#     ".txt"
#   )
# )