# RNA-seq中htseq参考基因组是用的hg38?

here::i_am("script/04-residual.R")

library(here)
library(dplyr)
library(MethReg)

load(here("output/data/cpg_unicox.rda"))
load(here("output/data/dmp_beta_tumor.rda"))
beta <- readRDS(here("output/data/beta_normaltumor.rds"))
methy_sample_info <- readRDS(here("output/data/methy_sample_info.rds"))

unicox_beta <- beta[match(cpg_unicox_sig$feature, rownames(beta)), ]
gene_fpkm <- data.table::fread(
    here("data/TCGA-HNSC.htseq_fpkm.tsv")
)

# samples with methylation and gene expression data
samples <- intersect(names(unicox_beta), names(gene_fpkm))
# keep the tumor samples: 500 samples
sample_type <- methy_sample_info$Sample_Type[
    match(samples, methy_sample_info$sample.submitter_id)
]
samples <- samples[sample_type == "tumor"]
# samples have both age and gender: 499 samples
hnsc_clic <- readRDS(here("output/data/hnsc_clic.rds"))
sample_barcode <- TCGAutils::TCGAbarcode(samples)
samples <- data.frame(barcode = sample_barcode, 
                      sample_id = samples)
sample_info <- left_join(samples, hnsc_clic, by = c("barcode" = "barcode"))
sample_info <- select(sample_info,
                      barcode, sample_id,
                      age, gender)
sample_info <- na.omit(sample_info)

unicox_beta <- select(unicox_beta, all_of(sample_info$sample_id))
unicox_gene_fpkm <- select(gene_fpkm, 
                           all_of(c("Ensembl_ID", sample_info$sample_id))) |> 
    tibble::column_to_rownames("Ensembl_ID")

rownames(sample_info) <- sample_info$sample_id
sample_meta <- select(sample_info, age, gender)
gene_residual <- get_residuals(
    data = log2(unicox_gene_fpkm + 1),
    metadata.samples = sample_meta , 
    cores = 5
)

# log trans
unicox_mval <- unicox_beta / (1 - unicox_beta)
cpg_residual <- get_residuals(
    data = unicox_mval,
    sample_meta
)

save(cpg_residual, gene_residual, 
     file = here("output/data/residual.rda"))
