# install load libraries ----
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse, # load multiple 'tidyverse' packages in a single step
  hardhat, # create high-quality R packages for modeling
  sf, # provides simple features access for R
  amt, #   manage and analyze animal movement data
  lubridate, # parse and manipulate dates
  circular, # for circular data
  knitr, # for dynamic report generation
  kableExtra, # for creating nice tables
  here # for file paths
)

# heading ----
rm(list = ls())
here::i_am("code/_master_code.R")

Sys.setenv(LANG = "en")
options(scipen = 999)
message(
  paste0("Date: ", format(Sys.Date(), "%d %B %Y")), "\n",
  paste0("By: ", "Francesco Invernizzi"), "\n",
  paste0(
    "Description: ",
    "code - master code to upload packages and set directories"
  ), "\n",
  paste0("Version of R used: ", R.Version()$version.string)
)

# set folder paths ----
# source("code/_functions.R")
dir_data <- here("data")
