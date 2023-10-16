here::i_am("script/05-methreg-lr.R")

library(here)
library(MethReg)
library(dplyr)


# methreg -----------------------------------------------------------------

load(here("output/data/residual.rda"))
rownames(gene_residual) <- strsplit(rownames(gene_residual), "\\.") |> 
    vapply(\(x) x[1], FUN.VALUE = character(1))

cpg_se <- make_dnam_se(
    dnam = cpg_residual,
    genome = "hg19",
    arrayType = "450k",
    betaToM = FALSE,
    verbose = TRUE
)
gene_se <- make_exp_se(
    gene_residual,
    genome = "hg19",
    verbose = TRUE
)

hg19_450k <- MethReg::get_met_probes_info(genome = "hg19", arrayType = "450k")
region_450k <- MethReg::make_names_from_granges(hg19_450k)

## create triplet --------------------------------------------------

library(JASPAR2022)
library(TFBSTools)
# library(motifmatchr)
# library(BSgenome.Hsapiens.UCSC.hg19)


# 1. gene promoter overlap
triplet_promoter <- create_triplet_distance_based(
    cpg_se,
    genome = "hg19",
    target.method = "genes.promoter.overlap",
    target.promoter.upstream.dist.tss = 1500,
    target.promoter.downstream.dist.tss = 1500,
    motif.search.window.size = 500,
    motif.search.p.cutoff = 1e-08,
    cores = 5
)


# # 2. window
triplet_distal_window <- create_triplet_distance_based(
    cpg_se,
    genome = "hg19",
    target.method = "window",
    target.window.size = 500 * 10^3,
    target.rm.promoter.regions.from.distal.linking = TRUE,
    motif.search.window.size = 500,
    motif.search.p.cutoff = 1e-08,
    cores = 5
)


# 3. nearby genes
triplet_distal_nearby_genes <- create_triplet_distance_based(
    cpg_se,
    genome = "hg19",
    target.method = "nearby.genes",
    target.num.flanking.genes = 5,
    target.window.size = 500 * 10^3,
    target.rm.promoter.regions.from.distal.linking = TRUE,
    motif.search.window.size = 500,
    motif.search.p.cutoff = 1e-08,
    cores = 5
)


# 4. regulon
triplet_regulon <- create_triplet_regulon_based(
    region = cpg_se,
    genome = "hg19",
    motif.search.window.size = 500,
    tf.target = dorothea::dorothea_hs,
    max.distance.region.target = 10^6
)

# extract probeID from regionID
triplet_regulon$probeID <- names(hg19_450k)[
    match(triplet_regulon$regionID, region_450k)
]

# merge
triplet_all <- list(triplet_promoter, triplet_distal_window,
                    triplet_distal_nearby_genes, triplet_regulon)
triplet_nms <- lapply(triplet_all, names) |> 
    Reduce("intersect", x = _)
triplet_all <- lapply(triplet_all, \(x) x[triplet_nms]) |> 
    bind_rows()
triplet_all_unique <- distinct(triplet_all, probeID, target)


## TF encichment score ------------------
# Using dorothea and viper
regulon_dorothea <- dorothea::dorothea_hs
rnaseq_tf_es <- get_tf_ES(
    exp = na.omit(SummarizedExperiment::assay(gene_se)),
    regulons = regulon_dorothea
)


## fit linear model: target gene, TF, DNAm -----------------------

interaction_promoter <- interaction_model(
    triplet = triplet_promoter,
    dnam = cpg_se,
    exp = gene_se,
    # dnam.group.threshold = 0.1,
    sig.threshold = 0.05,
    fdr = TRUE,
    stage.wise.analysis = TRUE,
    filter.correlated.tf.exp.dnam = FALSE,
    filter.triplet.by.sig.term = TRUE,
    cores = 6, 
    tf.activity.es = rnaseq_tf_es
)
interaction_distal_window <- interaction_model(
    triplet = triplet_distal_window, 
    dnam = cpg_se,
    exp = gene_se,
    # dnam.group.threshold = 0.1,
    sig.threshold = 0.05,
    fdr = TRUE,
    stage.wise.analysis = TRUE,
    filter.correlated.tf.exp.dnam = FALSE,
    filter.triplet.by.sig.term = TRUE ,
    cores = 5, 
    tf.activity.es = rnaseq_tf_es
)
interaction_distal_nearby_gene <- interaction_model(
    triplet = triplet_distal_nearby_genes, 
    dnam = cpg_se,
    exp = gene_se,
    # dnam.group.threshold = 0.1,
    sig.threshold = 0.05,
    fdr = TRUE,
    stage.wise.analysis = TRUE,
    filter.correlated.tf.exp.dnam = FALSE,
    filter.triplet.by.sig.term = TRUE ,
    cores = 5, 
    tf.activity.es = rnaseq_tf_es
)
interaction_regulon <- interaction_model(
    triplet = triplet_regulon, 
    dnam = cpg_se,
    exp = gene_se,
    # dnam.group.threshold = 0.1,
    sig.threshold = 0.05,
    fdr = TRUE,
    stage.wise.analysis = TRUE,
    filter.correlated.tf.exp.dnam = FALSE,
    filter.triplet.by.sig.term = TRUE ,
    cores = 5, 
    tf.activity.es = rnaseq_tf_es
)

interaction_all <- list(
    interaction_promoter, interaction_distal_window,
    interaction_distal_nearby_gene, interaction_regulon
)
interaction_nms <- lapply(interaction_all,names) |> 
    Reduce("intersect", x = _)
interaction_all <- lapply(
    interaction_all,
    \(x) x[interaction_nms]) |> 
    bind_rows()
interaction_sig <- dplyr::filter(
    interaction_all,
    `RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.05
)
save(interaction_promoter,
     interaction_distal_window,
     interaction_distal_nearby_gene,
     interaction_regulon,
     interaction_sig, 
     interaction_all,
     file = here("output/data/methreg_interaction.rda"))


# stage and trans  --------------------------------------------------------

hnsc_clic <- readRDS(here("output/data/hnsc_clic.rds"))
hnsc_clic_updated <- read.csv(
    here("data/clinical_HNSC.csv"),
    row.names = 1
)
hnsc_clic <- full_join(hnsc_clic, hnsc_clic_updated, 
                       by = c("barcode" = "submitter_id"))

## stage 
hnsc_clic_stage <- filter(
    hnsc_clic,
    clinical_stage != "[Not Available]") |> 
    mutate(stage_group = ifelse(clinical_stage %in% c("Stage I", "Stage II"),
                                "early", "advanced"),
           .after = clinical_stage)
## trans
hnsc_clic_trans <- filter(
    hnsc_clic,
    ajcc_clinical_n != "NX") |> 
    mutate(trans_group = ifelse(ajcc_clinical_n  == "N0", 
                                "nontrans", "trans"),
           .after = ajcc_clinical_n)
# 506 samples
hnsc_clic_stage_trans <- inner_join(
    hnsc_clic_stage, 
    select(hnsc_clic_trans, barcode, trans_group),
    by = c("barcode" = "barcode")
)


# linear regression: direct interaction of cpg with gene using triplet --------

# gene_residual_rlm <- SummarizedExperiment::assay(gene_se)
# gene_residual_rlm <- na.omit(gene_residual_rlm)
# gene_shared <- intersect(unique(triplet_all_unique$target), 
#                          rownames(gene_residual_rlm))
# gene_residual_rlm <- gene_residual_rlm[
#     match(gene_shared, rownames(gene_residual_rlm)), 
# ]
# 
# cpg_shared <- intersect(unique(triplet_all_unique$probeID),
#                         rownames(cpg_residual))
# cpg_residual_rlm <- cpg_residual[
#     match(cpg_shared, rownames(cpg_residual)), 
# ]
# triplet_all_unique_rlm <- dplyr::filter(
#     triplet_all_unique,
#     probeID %in% cpg_shared,
#     target %in% gene_shared) |> 
#     t() |> as.data.frame()
#     
# 
# source(here("R/rlm-gene-cpg.R"))
# cpg_gene_rlm_res <- auxfunction(
#     gene_residual_rlm,
#     cpg_residual_rlm,
#     gene_cpg = triplet_all_unique_rlm$V2,
#     clic = hnsc_clic_stage_trans
# )
# 
# cpg_gene_rlm_res <- lapply(
#     triplet_all_unique_rlm, 
#     \(x) auxfunction(
#         gene_residual = gene_residual_rlm,
#         cpg_residual = cpg_residual_rlm,
#         gene_cpg = x,
#         clic = hnsc_clic_stage_trans
#     )
# )
# cpg_gene_rlm_res <- bind_rows(cpg_gene_rlm_res)
# 
# cpg_gene_rlm_res$probe_id <- unlist(triplet_all_unique_rlm[1, ])
# cpg_gene_rlm_res$gene_id <- unlist(triplet_all_unique_rlm[2, ])
# rownames(cpg_gene_rlm_res) <- NULL
# saveRDS(cpg_gene_rlm_res, here("output/data/cpg_gene_rlm_res.rds"))
# 
# cpg_gene_rlm_res_sig <- dplyr::filter(
#     cpg_gene_rlm_res,
#     RLM_met_residual_pvalue <= 0.05
# )
# save(cpg_gene_rlm_res_sig, 
#      here("output/data/cpg_gene_rlm_res_sig.rds"))


# liner regression: use cpg gene pairs rather than triplet ---------

## 1. get linked genes of candidate cpgs: cpg - gene pairs
promoter_cpg_gene_pair <- MethReg::get_region_target_gene(
    cpg_se,
    method = "genes.promoter.overlap",
    genome = "hg19",
    promoter.upstream.dist.tss = 1500,
    promoter.downstream.dist.tss = 1500
)
promoter_cpg_gene_pair$probeID <- names(hg19_450k)[
    match(promoter_cpg_gene_pair$regionID, region_450k)
]

nearby_cpg_gene_pair <- MethReg::get_region_target_gene(
    cpg_se,
    method = "nearby.genes",
    num.flanking.genes = 5,
    rm.promoter.regions.from.distal.linking = TRUE,
    genome = "hg19"
)
nearby_cpg_gene_pair$probeID <- names(hg19_450k)[
    match(nearby_cpg_gene_pair$regionID, region_450k)
]

window_cpg_gene_pair <- MethReg::get_region_target_gene(
    cpg_se,
    method = "window",
    genome = "hg19",
    window.size = 500 * 10^3,
    rm.promoter.regions.from.distal.linking = TRUE
)
window_cpg_gene_pair$probeID <- names(hg19_450k)[
    match(window_cpg_gene_pair$regionID, region_450k)
]

cpg_gene_pair_all <- list(
    dplyr::select(promoter_cpg_gene_pair, probeID, target),
    dplyr::select(window_cpg_gene_pair, probeID, target),
    dplyr::select(nearby_cpg_gene_pair, probeID, target)) |> 
    do.call(rbind, args = _) |> 
    dplyr::distinct()

save(cpg_gene_pair_all,
     promoter_cpg_gene_pair,
     window_cpg_gene_pair,
     nearby_cpg_gene_pair,
     file = here("output/data/cpg_gene_pair.rda"))

gene_residual_rlm <- SummarizedExperiment::assay(gene_se)
gene_residual_rlm <- na.omit(gene_residual_rlm)
gene_shared <- intersect(unique(cpg_gene_pair_all$target), 
                         rownames(gene_residual_rlm))
gene_residual_rlm <- gene_residual_rlm[
    match(gene_shared, rownames(gene_residual_rlm)), 
]

# cpg_residual_rlm <- cpg_residual
cpg_shared <- intersect(unique(cpg_gene_pair_all$probeID),
                        rownames(cpg_residual))
cpg_residual_rlm <- cpg_residual[
    match(cpg_shared, rownames(cpg_residual)), 
]

cpg_gene_pair_all_rlm <- dplyr::filter(
    cpg_gene_pair_all,
    probeID %in% cpg_shared,
    target %in% gene_shared) |> 
    t() |> as.data.frame()

source(here("R/rlm-gene-cpg.R"))
cpg_gene_pair_rlm_res <- lapply(
    cpg_gene_pair_all_rlm, 
    \(x) auxfunction(
        gene_residual = gene_residual_rlm,
        cpg_residual = cpg_residual_rlm,
        gene_cpg = x,
        clic = hnsc_clic_stage_trans
    )) |> 
    bind_rows()
cpg_gene_pair_rlm_res$RLM_met_residual_padj <- 
    p.adjust(cpg_gene_pair_rlm_res$RLM_met_residual_pvalue, method = "fdr")
cpg_gene_pair_rlm_res$probe_id <- unlist(cpg_gene_pair_all_rlm[1, ])
cpg_gene_pair_rlm_res$gene_id <- unlist(cpg_gene_pair_all_rlm[2, ])
cpg_gene_pair_rlm_sig <- dplyr::filter(
    cpg_gene_pair_rlm_res,
    RLM_met_residual_padj < 0.01
)
save(cpg_gene_pair_rlm_res, 
     cpg_gene_pair_rlm_sig,
     file = here("output/data/cpg_gene_pair_rlm.rda"))


# best_marker <- c("cg02409878", "cg23867673", "cg01995815")
