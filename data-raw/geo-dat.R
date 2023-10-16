here::i_am("data-raw/geo-dat.R")
library(GEOquery)
library(here)



# gse124052 ---------------------------------------------------------------

# zip_files <- list.files(
#     here("data/GSE124052/GSE124052_RAW/"),
#     full.names = TRUE
# )
# lapply(zip_files, R.utils::gunzip)


source(here("R/idat-preprocess.R"))

gse124052_clic <- read.csv(here("data/GSE124052/filter_GSE124052_clin_1.csv"))
names(gse124052_clic)[1] <- "sample"

gse124052_dat <- preprocess_dnam(
    base_dir = here("data/GSE124052/GSE124052_RAW/"),
    clic = gse124052_clic
)
save(gse124052_dat, file = here("output/data/gse124052.rda"))


# gse41114 ----------------------------------------------------------------

gse41114 <- getGEO("GSE4111i4")[[1]]
gse41114_meta <- Biobase::pData(gse41114) |> 
    dplyr::select(geo_accession,
                  sample = `case:ch1`,
                  sample_type = source_name_ch1) |> 
    dplyr::filter(sample_type == "OSCC")
## clinical data supplementary table2 of paper PMID:23619168 
gse41114_table_s2 <- readxl::read_xlsx(
    here("data/GSE41114/supplementary_tables.xlsx"),
    skip = 4,
    sheet = 2) |> 
    dplyr::select(sample = `Sample #`,
                  OS.time = `Overall Survival (years)`,
                  OS = `Overall Survival`) |> 
    dplyr::mutate(sample = as.character(sample))
gse41114_meta <- dplyr::inner_join(gse41114_meta, gse41114_table_s2)
gse41114_meta$sample <- gse41114_meta$geo_accession
gse41114_beta <- Biobase::exprs(gse41114)[, gse41114_meta$sample] |> 
    t() |> as.data.frame()
gse41114_beta$sample <- rownames(gse41114_beta)

save(gse41114_beta, gse41114_meta, 
     file = here("output/data/gse41114.rda"))


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


# gse52793 ----------------------------------------------------------------

gse52793 <- getGEO("GSE52793")[[1]]
gse52793_meta <- Biobase::pData(gse52793) |> 
    dplyr::mutate(OS.time = as.numeric(`survival_months:ch1`)/12,
                  OS = ifelse(`vital_status:ch1` == "Deceased", 1, 0)) |> 
    dplyr::select(sample = geo_accession,
                  OS.time,
                  OS)
gse52793_beta <- Biobase::exprs(gse52793)[, gse52793_meta$sample] |> 
    t() |> as.data.frame()
gse52793_beta$sample <- rownames(gse52793_beta)

save(gse52793_beta, gse52793_meta, 
     file = here("output/data/gse52793.rda"))


# gse87053 ----------------------------------------------------------------

gse87053 <- getGEO("GSE87053")
