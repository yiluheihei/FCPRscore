here::i_am("script/07-risk_group-validation.R")

library(here)
library(dplyr)

load(here("output/data/risk.rda"))
source(here("R/surv.R"))


# cptac -------------------------------------------------------------------

load(here("output/data/cptac_hnsc.rda"))
plot_km(cpg_marker, cpg_marker_cox$coef,
        cptac3_dat$beta, cptac3_dat$clic,
        quantile(risk_score, 0.65))


# gse124052 ------------------------------------------------------------------

load(here("output/data/gse124052.rda"))
plot_km(cpg_marker, cpg_marker_cox$coef,
        gse124052_dat$beta, gse124052_dat$clic,
        quantile(risk_score, 0.65))


# gse75537 ----------------------------------------------------------------

load(here("output/data/gse75537.rda"))
plot_km(cpg_marker, cpg_marker_cox$coef,
        gse75537_beta, gse75537_meta,
        quantile(risk_score, 0.65))


# cptac + gse124052 + gse75537 -------------------------------------------------------
combined_cpg <- intersect(names(cptac3_dat$beta),
                          names(gse124052_dat$beta)) |>
    intersect(names(gse75537_beta))
combined_beta <- rbind(
    cptac3_dat$beta[combined_cpg],
    gse124052_dat$beta[combined_cpg],
    gse75537_beta[combined_cpg]
)
combined_clic <- rbind(
    cptac3_dat$clic[c("sample", "OS.time", "OS")],
    gse124052_dat$clic[c("sample", "OS.time", "OS")],
    gse75537_meta[c("sample", "OS.time", "OS")]
)
plot_km(cpg_marker, cpg_marker_cox$coef,
        combined_beta, combined_clic,
        quantile(risk_score, 0.65))
