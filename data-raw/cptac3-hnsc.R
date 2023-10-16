here::i_am("data-raw/cptac3-hnsc.R")

library(TCGAbiolinks)
library(here)

cptac3_clic <- read.csv(here("data/CPTAC/filter_CPTAC_clin_1.csv"))
names(cptac3_clic)[1] <- "sample"
cptac3_dnam <- GDCquery(project = "CPTAC-3",
         access = "open",
         data.category = "DNA Methylation",
         data.type = "Masked Intensities",
         data.format = "IDAT")

# filter hnsc
cptac3_hnsc_info <- dplyr::inner_join(
    cptac_clic,
    cptac3_dnam$results[[1]],
    by = c("sample" = "cases.submitter_id")) |> 
    dplyr::filter(grepl("Primary Tumor", sample_type)) |> 
    dplyr::mutate(sample_type = "tumor")

# join clinical information
cptac3_hnsc_meta <- dplyr::select(
    cptac3_hnsc_info,
    -file_name, -id, -submitter_id,-file_id,  -channel,
    -created_datetime, -md5sum) |> 
    dplyr::distinct()

barcode <- cptac3_hnsc_info$cases |> 
    strsplit(split = ";") |> 
    vapply(\(x)x[1], FUN.VALUE = character(1)) |> 
    unique()
#　query (barcode) for downloading
cptac3_hnsc <- GDCquery(project = "CPTAC-3",
                        access = "open",
                        data.category = "DNA Methylation",
                        data.type = "Masked Intensities",
                        data.format = "IDAT",
                        barcode = barcode)
GDCdownload(cptac3_hnsc, query = cptac3_hnsc, 
            directory = here("data/GDCdata/cptac3-hnsc/"))


# read and preprocess
source(here("R/idat-preprocess.R"))
cptac3_dat <- preprocess_dnam(
    here("data/GDCdata/cptac3-hnsc/CPTAC-3",
         "DNA_Methylation/Masked_Intensities/"),
    clic = cptac3_clic,
    cptac3_hnsc = TRUE,
    cptac3_hnsc_meta
)



save(cptac3_dat, file = here("output/data/cptac_hnsc.rda"))


