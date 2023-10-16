here::i_am("script/03-cpg-unicox.R")

library(here)
library(dplyr)

# download TCGA-HNSC clinical data
# library(TCGAbiolinks)
# hnsc_clic <- GDCquery_clinic("TCGA-HNSC", type = "clinical")
# saveRDS(hnsc_clic, here("output/data/hnsc_clic.rds"))

clic <- readxl::read_xlsx(here("data/clin.xlsx"))
hnsc_clic <- dplyr::filter(clic , type == "HNSC")
hnsc_clic <- hnsc_clic[-1]
hnsc_clic <- dplyr::rename(
    hnsc_clic, 
    barcode = bcr_patient_barcode,
    age = age_at_initial_pathologic_diagnosis) |>
    mutate(OS.time = OS.time/365)
saveRDS(hnsc_clic, here("output/data/hnsc_clic.rds"))

methy_sample_info <- readRDS(here("output/data/methy_sample_info.rds"))
beta <- readRDS(here("output/data/beta_normaltumor.rds"))
dmp_sig <- readRDS(here("output/data/dmp_sig.rds"))

dmp_beta <- beta[match(rownames(dmp_sig), rownames(beta)), ]
dmp_beta <- tibble::rownames_to_column(dmp_beta, "cpg")
dmp_beta <- data.table::transpose(
    dmp_beta,
    keep.names = "sample.submitter_id",
    make.names = "cpg"
)
dmp_beta$barcode <- TCGAutils::TCGAbarcode(dmp_beta$sample.submitter_id)
dmp_beta$Sample_Type <- methy_sample_info$Sample_Type[
    match(dmp_beta$sample.submitter_id, methy_sample_info$sample.submitter_id)
]

# select tumor samples
# dmp_beta_tumor <- dplyr::filter(dmp_beta, Sample_Type == "tumor")
# dmp_beta_surv <- left_join(hnsc_clic, dmp_beta_tumor, by = c("barcode" = "barcode"))

# shensi method: 使用`match`匹配barcode只是return the
# first position of the element, 所以可能选择了normal的样本, 应
# 选择tumor sample
# dmp_beta_tumor <- dmp_beta[match(hnsc_clic$barcode, dmp_beta$barcode), ]
dmp_beta_tumor <- dplyr::filter(dmp_beta, Sample_Type == "tumor")
dmp_beta_surv <- left_join(hnsc_clic, dmp_beta_tumor, 
                            by = c("barcode" = "barcode"))
save(dmp_beta_tumor, dmp_beta_surv,
        file = here("output/data/dmp_beta_tumor.rda"))

source(here("R/unicox.R"))

# about 15 minutes
cpg_unicox_res <- lapply(rownames(dmp_sig), 
                         run_uni_cox, df = dmp_beta_surv) |> 
    dplyr::bind_rows()
cpg_unicox_sig <- dplyr::filter(
    cpg_unicox_res,
    ph_pval < 0.01,
    ph_hypothesis_p > 0.05
)
save(cpg_unicox_res, cpg_unicox_sig,
     file = here("output/data/cpg_unicox.rda"))
