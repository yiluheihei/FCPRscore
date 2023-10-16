here::i_am("data-raw/geo-dat.R")
library(GEOquery)
library(here)


# gse124052 ---------------------------------------------------------------

source(here("R/idat-preprocess.R"))

gse124052_clic <- read.csv(here("data/GSE124052/filter_GSE124052_clin_1.csv"))
names(gse124052_clic)[1] <- "sample"

gse124052_dat <- preprocess_dnam(
    base_dir = here("data/GSE124052/GSE124052_RAW/"),
    clic = gse124052_clic
)
save(gse124052_dat, file = here("output/data/gse124052.rda"))


# gse75537 ----------------------------------------------------------------

gse75537 <- getGEO("GSE75537")
gse75537 <- gse75537[[1]]
gse75537_meta <- Biobase::pData(gse75537) |>
    dplyr::mutate(sample_type = ifelse(grepl("_T$", source_name_ch1),
                                       "tumor", "normal"),
                  OS.time = as.numeric(`overall survival (months):ch1`)/12,
                  OS = ifelse(`survival status:ch1` == "Dead", 1, 0)) |>
    dplyr::select(sample = geo_accession,
                  sample_type,
                  OS.time,
                  OS) |>
    dplyr::filter(sample_type == "tumor")
gse75537_beta <- Biobase::exprs(gse75537)[, gse75537_meta$sample]
gse75537_beta <- t(gse75537_beta) |> as.data.frame()
gse75537_beta$sample <- rownames(gse75537_beta)

save(gse75537_beta, gse75537_meta,
     file = here("output/data/gse75537.rda"))


