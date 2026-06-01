library(dplyr)

# Una funzione che restituisce un dataframe di N righe
genera_dati <- function(n) {
    tibble(valore = rnorm(n), tipo = "Test")
}

dati_finti <- tibble(gruppo = c("A", "A", "B", "B", "B"))

dati_finti %>%
    reframe(
        suppressWarnings(as_tibble(sim.burst(n = n(), mod = mod1))),
        .by = gruppo
    )
