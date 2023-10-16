
here::i_am("script/02-differential-analysis.R")
library(here)
library(dplyr)
library(limma)


# dmp ---------------------------------------------------------------------

source(here("R/dmpFinder2.R"))

beta <- readRDS(here("output/data/beta.rds"))
beta <- as.data.frame(na.omit(beta))

# keep tumor and normal samples
methy_sample_info <- readRDS(here("output/data/methy_sample_info.rds"))
methy_sample_info <- dplyr::filter(
    methy_sample_info, 
    Sample_Type %in% c("tumor", "normal")
)
beta <- beta[methy_sample_info$sample.submitter_id]
#names(beta) <- methy_sample_info$sample.submitter_id
saveRDS(beta, here("output/data/beta_normaltumor.rds"))

m_val <- log2(beta/(1 - beta))

dmp <- dmpFinder2(as.matrix(m_val), 
                  pheno = methy_sample_info$Sample_Type, 
                  type = "categorical" , 
                  shrinkVar = TRUE)

dmp_sig <- dplyr::filter(dmp, qval < 0.05, abs(diff_beta) > 0.2)
dmp_sig <- mutate(dmp_sig,
                  sig = ifelse(diff_beta > 0.2, "up", "down"))
saveRDS(dmp_sig, here("output/data/dmp_sig.rds"))


# deg ---------------------------------------------------------------------
library(DESeq2)

gene_count <- data.table::fread(here("data/TCGA-HNSC.htseq_counts.tsv")) |> 
    tibble::column_to_rownames("Ensembl_ID")
gene_count <- 2 ^ gene_count - 1

gene_sample_info <- dplyr::filter(
    methy_sample_info,
    sample.submitter_id %in% names(gene_count)
    
)
gene_count <- gene_count[gene_sample_info$sample.submitter_id]

dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(round(gene_count)), 
    colData =  gene_sample_info, design = ~ Sample_Type
)
#　remove gene with count lower than 4 ?
# dds <- dds[rowSums(counts(dds)) > 4, ] 
dds <- DESeq(dds, fitType = 'mean', 
             minReplicatesForReplace = 7, 
             parallel = FALSE)
deg <- results(dds, contrast = c("Sample_Type" , "tumor" , "normal"))
deg <- data.frame(deg)
saveRDS(deg, here("output/data/deg.rds"))
# padj ?
deg_sig <- dplyr::filter(
    deg,
    abs(log2FoldChange) > 0, 
    padj < 0.05) |> 
    mutate(sig = ifelse(log2FoldChange > 0, "up", "down")) |> 
    tibble::rownames_to_column("gene_id")
