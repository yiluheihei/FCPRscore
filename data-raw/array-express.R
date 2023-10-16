here::i_am("data-raw/array-express.R")

library(here)
library(ArrayExpress)
library(minfi)

source(here("R/idat-preprocess.R"))

emtab_10576 <- getAE("E-MTAB-10576", 
                     here("data/array_express/E-MTAB-10576/"))
emtab_10576 <- ae2bioc(mageFiles = emtab_10576)

options(warn = 0)
read.delim("data/array_express/E-MTAB-10576/E-MTAB-10576.sdrf.txt")

ematb_10576 <- minfi::read.metharray.exp(
    here("data/array_express/E-MTAB-10576/"),
    verbose = TRUE)

# 10576
emtab_10576_clic <- read.delim(
    here("data/array_express/E-MTAB-10576/E-MTAB-10576.sdrf.txt")) |> 
    dplyr::select(sample = Source.Name,
                  clinical = Factor.Value.clinical.information.) |> 
    dplyr::distinct()

emtab_10576_dat <- preprocess_dnam(
    here("data/array_express/E-MTAB-10576/"),
    emtab_10576_clic,
    cptac3_hnsc = FALSE
)

# 12431
emtab_12431_clic <- read.delim(
    here("data/array_express/E-MTAB-12431/E-MTAB-12431.sdrf.txt")) |> 
    dplyr::select(sample = Source.Name,
                  clinical = Factor.Value.clinical.information.) |> 
    dplyr::distinct()

emtab_12431_dat <- preprocess_dnam(
    here("data/array_express/E-MTAB-12431/"),
    emtab_12431_clic,
    cptac3_hnsc = FALSE
)

emtab_10577_clic <- 
    emtab_12431_clic <- read.delim(
        here("data/array_express/E-MTAB-10577/E-MTAB-10577.sdrf.txt")) |> 
    dplyr::select(sample = Source.Name,
                  clinical = Factor.Value.clinical.history.) |> 
    dplyr::distinct()
emtab_10577_dat <- preprocess_dnam(
    here("data/array_express/E-MTAB-10577/"),
    emtab_10577_clic,
    cptac3_hnsc = FALSE
)

emtab_10578_clic <- read.delim(
        here("data/array_express/E-MTAB-10578/E-MTAB-10578.sdrf.txt")) |> 
    dplyr::select(sample = Source.Name,
                  clinical = Factor.Value.clinical.history.) |> 
    dplyr::distinct()
emtab_10578_dat <- preprocess_dnam(
    here("data/array_express/E-MTAB-10578/"),
    emtab_10578_clic,
    cptac3_hnsc = FALSE
)

source(here("R/surv.R"))
load(here("output/data/risk_cpg_marker.rda"))

plot_surv(
    cpg_marker, cpg_marker_cox$coef,
    beta = emtab_10576_dat$beta, clic = emtab_10576_dat$clic,
    risk_score_cutoff = quantile(risk_score, 0.65),
    type = "group_bar"
)

plot_surv(
    cpg_marker, cpg_marker_cox$coef,
    beta = emtab_12431_dat$beta, clic = emtab_12431_dat$clic,
    risk_score_cutoff = quantile(risk_score, 0.65),
    type = "group_bar"
)

plot_surv(
    cpg_marker, cpg_marker_cox$coef,
    beta = emtab_10577_dat$beta, clic = emtab_10577_dat$clic,
    risk_score_cutoff = quantile(risk_score, 0.65),
    type = "group_bar"
)

plot_surv(
    cpg_marker, cpg_marker_cox$coef,
    beta = emtab_10578_dat$beta, clic = emtab_10578_dat$clic,
    risk_score_cutoff = quantile(risk_score, 0.65),
    type = "group_bar"
)



combined_cpg <- intersect(names(emtab_10576_dat$beta), 
                          names(emtab_10577_dat$beta)) |> 
    intersect(names(emtab_10578_dat$beta)) |> 
    intersect(names(emtab_12431_dat$beta))
combined_beta <- rbind(
    emtab_10576_dat$beta[combined_cpg],
    emtab_10577_dat$beta[combined_cpg],
    emtab_10578_dat$beta[combined_cpg],
    emtab_12431_dat$beta[combined_cpg]
)
combined_clic <- rbind(
    emtab_10576_dat$clic,
    emtab_10577_dat$clic,
    emtab_10578_dat$clic,
    emtab_12431_dat$clic
)
plot_surv(
    cpg_marker, cpg_marker_cox$coef,
    beta = combined_beta, clic = combined_clic,
    risk_score_cutoff = quantile(risk_score, 0.65),
    type = "group_bar"
)
