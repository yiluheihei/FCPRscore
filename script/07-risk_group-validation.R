here::i_am("script/07-risk_group-validation.R")

library(here)
library(dplyr)

load(here("output/data/risk.rda"))
source(here("R/surv.R"))

# cptac -------------------------------------------------------------------

# cptac_beta <- data.table::fread(
#     here("data/CPTAC/filter_CPTAC.csv")) |>
#     data.table::transpose(keep.names = "sample", make.names = "V1") |>
#     data.frame()
# 
# cptac_clic <- read.csv(here("data/CPTAC/filter_CPTAC_clin_1.csv"))
# names(cptac_clic)[1] <- "sample"
# plot_km(cpg_marker, cpg_marker_cox$coef,
#         cptac_beta, cptac_clic,
#         quantile(risk_score, 0.65))


load(here("output/data/cptac_hnsc.rda"))
plot_km(cpg_marker, cpg_marker_cox$coef, 
        cptac3_dat$beta, cptac3_dat$clic,
        quantile(risk_score, 0.65))


# gse124052 ------------------------------------------------------------------

##　from shensi
# gse124052_beta <- data.table::fread(here("data/GSE124052/filter_GSE124052.csv")) |>
#     data.table::transpose(keep.names = "sample", make.names = "V1") |>
#     data.frame()
# gse124052_clic <- read.csv(here("data/GSE124052/filter_GSE124052_clin_1.csv"))
# names(gse124052_clic)[1] <- "sample"
# plot_km(cpg_marker, cpg_marker_cox$coef,
#         gse124052_beta, gse124052_dat$clic,
#         quantile(risk_score, 0.65))

load(here("output/data/gse124052.rda"))
plot_km(cpg_marker, cpg_marker_cox$coef, 
        gse124052_dat$beta, gse124052_dat$clic,
        quantile(risk_score, 0.65))



# gse75537 ----------------------------------------------------------------
# data from gse75537 has been preprocessed using minfi

## data　from shensi
# gse75537_beta <- data.table::fread(here("data/GSE75537/filter_GSE75537.csv")) |>
#     data.table::transpose(keep.names = "sample", make.names = "V1") |>
#     data.frame()
# gse75537_clic <- read.csv(here("data/GSE75537/filter_GSE75537_clin_1.csv")) |>
#     dplyr::rename(sample = sampe)
# plot_km(cpg_marker, cpg_marker_cox$coef,
#         gse75537_beta, gse75537_clic,
#         quantile(risk_score, 0.65))

# data from geo
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


# gse41114 ----------------------------------------------------------------
# # 0.43

## data from shensi
# gse41114_beta <- data.table::fread(here("data/GSE41114/filter_GSE41114.csv")) |> 
#     data.table::transpose(keep.names = "sample", make.names = "V1") |> 
#     data.frame()
# gse41114_clic <- read.csv(here("data/GSE41114/filter_GSE41114_clin_1.csv")) |> 
#     dplyr::rename(sample = sampe)
# 
# plot_km(cpg_marker, cpg_marker_cox$coef, 
#         gse41114_beta, gse41114_clic,
#         quantile(risk_score, 0.5))


## data from geo
# gse41114 <- getGEO("GSE41114")[[1]]
# saveRDS(gse41114, here("output/data/gse41114_eset.rds"))

# gse41114 <- readRDS(here("output/data/gse41114_eset.rds"))
# gse41114_meta <- Biobase::pData(gse41114) |> 
#     dplyr::select(geo_accession,
#                   sample = `case:ch1`,
#                   sample_type = source_name_ch1) |> 
#     dplyr::filter(sample_type == "OSCC")
# ## clinical data supplementary table2 of paper PMID:23619168 
# gse41114_table_s2 <- readxl::read_xlsx(
#     here("data/GSE41114/supplementary_tables.xlsx"),
#     skip = 4,
#     sheet = 2) |> 
#     dplyr::select(sample = `Sample #`,
#            OS.time = `Overall Survival (years)`,
#            OS = `Overall Survival`) |> 
#     dplyr::mutate(sample = as.character(sample))
# gse41114_meta <- dplyr::inner_join(gse41114_meta, gse41114_table_s2)
# gse41114_meta$sample <- gse41114_meta$geo_accession
# gse41114_beta <- Biobase::exprs(gse41114)[, gse41114_meta$sample] |> 
#     t() |> as.data.frame()
# gse41114_beta$sample <- rownames(gse41114_beta)
load(here("output/data/gse41114.rda"))
plot_km(cpg_marker, cpg_marker_cox$coef, 
        gse41114_beta, gse41114_meta,
        quantile(risk_score, 0.35))



# gse52793----------------------------------------------------------------------

## from shensi
# gse52793_beta <- data.table::fread(here("data/GSE52793/filter_GSE52793.csv")) |>
#     data.table::transpose(keep.names = "sample", make.names = "V1") |>
#     data.frame()
# gse52793_clic <- read.csv(here("data/GSE52793/filter_GSE52793_clin_1.csv"))
# names(gse52793_clic)[1] <- "sample"
# 
# plot_km(cpg_marker, cpg_marker_cox$coef,
#         gse52793_beta, gse52793_clic,
#         quantile(risk_score, 0.65))

load(here("output/data/gse52793.rda"))
plot_km(cpg_marker, cpg_marker_cox$coef, 
        gse52793_beta, gse52793_meta,
        quantile(risk_score, 0.65))
