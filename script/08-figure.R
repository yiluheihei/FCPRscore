library(survival)
library(survminer)
library(GSEABase)
library(ggpubr)
library(ggsci)
library(ggstatsplot)
library(GSVA)
library(GSEABase)
library(enrichplot)
library(here)
library(dplyr)
library(patchwork)

base_dir <- "shensi-figure"


clin_filter <-  readRDS("shensi-figure/clin_filter.rds")
source(here("shensi-figure/function.R"))

# Figure 2 -------------------------------------------------------------------
load(here("shensi-figure/risk_cpg_marker.rda"))
load(here("shensi-figure/TCGA_km.rda"))
load(here("shensi-figure/cptac_hnsc.rda"))
load(here("shensi-figure/gse124052.rda"))
load(here("shensi-figure/gse75537.rda"))
load(here("shensi-figure/test_RS.rda"))
load(here("shensi-figure/test_combined.rda"))

risk_score_cutoff <- quantile(risk_score, 0.65)

km_args <- data.frame(
  beta = I(list(mutate(beta_surv, sample = rownames(beta_surv)),
           combined_beta,
           cptac3_dat$beta,
           gse75537_beta,
           gse124052_dat$beta)),
  clic = I(list(dplyr::rename(hnsc_clic, sample = barcode),
           combined_clic,
           cptac3_dat$clic,
           gse75537_meta,
           gse124052_dat$clic)),
  plot_title = c("TCGA-HNSCC", "Validation dataset", "CPTAC-HNSCC",
                 "GSE75537", "GSE124052")
)
p1_survs <- purrr::pmap(
  km_args, 
  plot_km, 
  cpg_marker = cpg_marker,
  cpg_marker_coef = cpg_marker_cox$coef,
  risk_score_cutoff = risk_score_cutoff)

auc_args <- list(
  group_text = list(RS_combined, RS_cptac, RS_gse75537, RS_gse124052),
  plot_title = c("Validation dataset", "CPTAC-HNSCC",
                 "GSE75537", "GSE124052"))
# auc_args$clin_risk <- lapply(
#   auc_args$clin_risk,
#   \(x) if("risk_group" %in% names(x)) {
#     rename(x, group = risk_group) }
#   else {
#     x
#   }
# )
p1_auc1 <- plot_AUC(clin_filter, plot_title = "TCGA-HNSCC")
p1_aucs <- purrr::pmap(auc_args, plot_AUC_test, 
                       quan = quantile(clin_filter$OS.time,0.75))

p1_survs <- lapply(p1_survs, \(x) x + theme(legend.position = "none"))
# p1 <- c(p1_survs, list(guide_area()), list(p1_auc1), p1_aucs)


p2 <- ((p1_survs[[1]] + theme(legend.position = c(0.8, 0.9))) +
       p1_survs[[2]] + p1_survs[[3]] + p1_survs[[4]] + plot_layout(nrow = 1)) /
  (p1_survs[[5]] + p1_auc1 + p1_aucs[[1]] + plot_spacer() +
     plot_layout(nrow = 1)) /
  (p1_aucs[[2]] + p1_aucs[[3]] + p1_aucs[[4]] + plot_spacer() +
     plot_layout(nrow = 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
ggsave(filename = "output/figure/Figure2.tiff", p2,
       dpi = 300, width = 12, height = 8.6)

ggsave(filename = "output/figure/Figure2.pdf", p2,
       width = 12, height = 8.6)

# Figure 3 -------------------------------------------------------------------
load(here("shensi-figure/external_model.rda"))
load(here("shensi-figure/model_list.rda"))


p3a <- AUC_compare_cpg(ls_cpgmodel) +
  theme(legend.text = element_text(size = 8))

p3b <- AUC_compare_gene(ls_genemodel) +
  theme(legend.text = element_text(size = 8))

cindex_data <- readRDS(here("shensi-figure/cindex-compared model.rds"))
cindex_data$model <- paste0(
  c("FCPRscore", "ChenD", "GhafarpourV","XuY", "ShenS", 
    "BasuB", "HanY", "ZhuQ", "ShiC", "LiuB", "ChenY"),
  c("", rep("_DNAm", 4), rep("_RNA", 6)))
feature_n <- gsub("(.*)(\\(n = \\d+\\))", "\\2", cindex_data$legend)
cindex_data$legend <- paste(cindex_data$model, "\n", feature_n)
cindex_data$label <- round(cindex_data$Cindex, 3)
cindex_data$type <- ifelse(grepl("DNAm", cindex_data$model), 
                           "DNAm-based models", "RNA-based models")
cindex_data$type[cindex_data$model == "FCPRscore"] <- "FCPRscore"
cindex_data <- cindex_data %>%
  mutate(model = factor(model, levels = cindex_data$model)) %>%
  mutate(legend = factor(legend, levels = cindex_data$legend))
cindex_data <- dplyr::filter(cindex_data, model != "XuY_DNAm")
p3c <- ggplot(cindex_data, aes(x = legend, y = Cindex)) + geom_col(aes(fill=legend)) +
  geom_text(aes(label = label), y = 0.4, color = "white") +
  # facet_wrap(vars(type), scales = "free") +
  scale_fill_aaas() +
  labs(title= NULL,x = "", fill = NULL) +
  scale_y_continuous(limits = c(0, 0.68), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.6))

p3 <- (p3a | p3b) / (p3c + plot_spacer() + plot_layout(nrow = 1, widths = c(8, 2))) + 
  plot_layout(nrow = 2, heights = c(6, 8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(file = "output/figure/Figure3.tiff", p3,
       dpi = 300, width = 10, height = 10)
ggsave(file = "output/figure/Figure3.pdf", p3,
       width = 10, height = 10)


# Figure 4 -------------------------------------------------------------------
result1 <- readRDS(here("shensi-figure/forest_result.rds"))
cindex_data1 <- read.csv("shensi-figure/cindex_clinical-info.csv",header = TRUE,row.names = 1)
load(here("shensi-figure/risk_cpg_marker.rda"))

library(forestplot)
p4a <- forest_plot(result1)
ggsave("output/figure/f4a.pdf", p4a, width = 9, height = 9.6)

p4b <- nomo_coxyear(clin_filter)
ggsave("output/figure/f4b.pdf", p4c, width = 7, height = 9.6)

nomogram_models <- c("Combined", "FCPRscore", "HPV_status", "Prior_treatment", 
                     "Clinical_stage")
p4c <- cindex_data1 %>%
  mutate(model = factor(model, levels = cindex_data1$model),
         legend = factor(nomogram_models, levels = nomogram_models),
         label = round(Cindex, digits = 3)) %>%
  ggplot(aes(x = legend, y = Cindex)) + geom_col(aes(fill = legend)) +
  geom_text(aes(label = label), y = 0.4, color = "white") +
  scale_fill_aaas() +
  labs(title= NULL,x = "", fill = NULL) +
  scale_y_continuous(limits = c(0, 0.68), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.6))

p4d <- calibration(clin_filter)

p4e <- plot_nomogram_auc(clin_filter,0.75)

p4cde <- wrap_elements(full = p4c) + 
  wrap_elements(full = p4d) + 
  wrap_elements(full = p4e) + 
  plot_layout(nrow = 1, widths = c(1,3.3, 1))
ggsave("output/figure/f4cde.pdf", p4cde,
       width = 20, height = 4)


# Figure 5 -------------------------------------------------------------------
# HRD_clin <- readRDS(here("shensi-figure/HRD_clin.rds"))
# estimate_clin <- readRDS(here("shensi-figure/estimate_clin.rds"))
# aneuploidy_clin <- readRDS(here("shensi-figure/aneuploidy_clin.rds"))
load(here("shensi-figure/cell_type.rda"))

library(scales) # trans_new
p5_cell_prop <- cell_proportion(merge.data) +
  theme(legend.position = "top")

## tcga hnscc molecular and immune features -----------------------------

hnscc_clic <- readxl::read_xlsx(here("data/tcga-hnscc-feature/immune-landscape_immune-feature.xlsx")) |> 
  dplyr::select(patient = `TCGA Participant Barcode`, 
                # study = `TCGA Study`, 
                leukocyte_fraction = `Leukocyte Fraction`,
                lymphocyte_infiltration_score = `Lymphocyte Infiltration Signature Score`,
                stromal_fraction = `Stromal Fraction`,
                intratumor_heterogeneity = `Intratumor Heterogeneity`,
                proliferation = Proliferation,
                SNV_neoantigens = `SNV Neoantigens`,
                indel_neoantigens = `Indel Neoantigens`,
                aneuploidy_score = `Aneuploidy Score`,
                hrd_score = `Homologous Recombination Defects`) |> 
  dplyr::filter(patient %in% names(risk_score))


# absolute purity 
absolute_purity <- readr::read_tsv(here("data/tcga-hnscc-feature/TCGA_absolute_purity.txt")) |> 
  select(sample = array, absolute_purity = purity) |> 
  mutate(patient = gsub("^(.*)-\\d+\\w+$", "\\1", sample), .before = sample)

# tmb and msisensor
cbioportal_hnscc_clic <- readr::read_tsv(here("data/tcga-hnscc-feature/hnsc_tcga_tmb_msi.tsv")) |> 
  select(patient = `Patient ID`,
         sample = `Sample ID`,
         cbioportal_aneuploidy_score = `Aneuploidy Score`,
         msi_score = `MSIsensor Score`,
         TMB = `TMB (nonsynonymous)`)
# hnscc_clic <- left_join(hnscc_clic, cbioportal_hnscc_clic)

# TIDE and ESTIMATE scores
tide <- readr::read_tsv(here("output/tide_res.txt"))
names(tide)[1] <- "sample"
tide <- dplyr::select(tide, sample, tide_score = TIDE, tide_msi_score = `MSI Score`,
                      tide_dysfunction = Dysfunction,
                      tide_exclusion = Exclusion) |> 
  dplyr::mutate(sample = gsub('^(.*-\\d+)\\w+$', "\\1", sample),
                patient = gsub("^(.*)-\\d+\\w+$", "\\1", sample), .before = sample)

# The ESTIMATE (Estimation of STromal and Immune cells in MAlignant Tumors using
# Expression data) method was developed with the aim to estimate the fraction of
# tumor cells in a sample by using gene expression instead of copy number data.
# The fundamental assumption of this method is that the tumor microenvironment
# is a very rich and dynamic ecosystem, in which immune infiltrating cells and
# stroma play a major role. The ESTIMATE score is defined as the combination
# (i.e. sum) of immune and stromal scores and can be thought of as a "non-tumor
# score". Consequently, a high ESTIMATE enrichment gives a low tumor purity
# score and vice versa.
estimate <- readr::read_tsv(here("output/estimate_res.tsv")) |> 
  rename(sample = sample_id, 
         stromal_score = StromalScore,
         immune_score = ImmuneScore,
         estimate_score = ESTIMATEScore,
         estimate_purity = TumorPurity) |> 
  mutate(patient = gsub("^(.*)-\\d+\\w+$", "\\1", sample), .before = sample)

# stemness
dnass <- data.table::fread(
  here("data/tcga-hnscc-feature/StemnessScores_DNAmeth_20170210.tsv.gz")) |> 
    data.table::transpose(
      keep.names = "sample",
      make.names = "sample") |> 
  select(sample, dna_stemness = DNAss) |> 
  mutate(patient = gsub("^(.*)-\\d+\\w+$", "\\1", sample), .before = sample)

rnass <- data.table::fread(
  here("data/tcga-hnscc-feature/StemnessScores_RNAexp_20170127.2.tsv.gz")) |> 
    data.table::transpose(
      keep.names = "sample",
      make.names = "sample") |> 
  select(sample, rna_stemness = RNAss) |> 
  mutate(patient = gsub("^(.*)-\\d+\\w+$", "\\1", sample), .before = sample)

hnscc_clic <- Reduce(
  left_join,
  list(hnscc_clic, cbioportal_hnscc_clic, tide, 
       estimate, dnass, rnass))
hnscc_clic$risk_score <- risk_score[match(hnscc_clic$patient, names(risk_score))]
hnscc_clic$risk_group <- ifelse(hnscc_clic$risk_score > risk_score_cutoff, 
                                "High-risk", "Low-risk")
immune_features <- setdiff(
  names(hnscc_clic), 
  c("patient", "sample", "risk_score", "risk_group", "msi_score",
    "cbioportal_aneuploidy_score", "tide_msi_score", "proliferation"))
immune_features_main <- c(
  "leukocyte_fraction", "lymphocyte_infiltration_score", "stromal_fraction", "intratumor_heterogeneity", 
  "aneuploidy_score", "hrd_score", "TMB",  "tide_score", "immune_score", 
  "estimate_score", "estimate_purity", "dna_stemness")
immune_features_supp <- setdiff(immune_features, immune_features_main)
immune_features_ylab <- strsplit(immune_features, "_") |> 
  vapply(\(x) {
    if (x[1] %in% c("rna", "dna", "estimate", "absolute", "tide", "hrd")) {
      x[1] <- toupper(x[1])
    }
    substr(x[1], 1, 1) <- toupper(substr(x[1], 1, 1))
    paste(x, collapse = " ")
  }, FUN.VALUE = character(1))
names(immune_features_ylab) <- immune_features
immune_features_main_ylab <- immune_features_ylab[immune_features_main]

hnscc_clic <- mutate(hnscc_clic, across(all_of(immune_features), as.numeric))
hnscc_clic <- select(hnscc_clic, all_of(immune_features), risk_group)

p_immune_features_main <- lapply(
  immune_features_main, 
  ggbetweenstats2, 
  df = hnscc_clic)
p_immune_features_main <- mapply(
  \(x, y) {
    if (y == "TMB") {
      x <- x + scale_y_sqrt()
    }
    x + labs(y = y)
  },
  p_immune_features_main,
  immune_features_main_ylab,
  SIMPLIFY = FALSE)
names(p_immune_features_main) <- immune_features_main

## immune feature supp -----------------
immune_features_supp_ylab <- immune_features_ylab[immune_features_supp]
p_immune_features_supp <- lapply(
  immune_features_supp, 
  ggbetweenstats2, 
  df = hnscc_clic)
p_immune_features_supp <- mapply(
  \(x, y) {
    if (grepl("neoantigen", y)) {
      x <- x + scale_y_sqrt()
    }
    x + labs(y = y)
  },
  p_immune_features_supp,
  immune_features_supp_ylab,
  SIMPLIFY = FALSE)
names(p_immune_features_supp) <- immune_features_supp

ps_stromal_stemness <- p_immune_features_supp$stromal_score +
  p_immune_features_supp$rna_stemness + 
  plot_layout(nrow = 1) +
  plot_annotation(tag_level = "A") &
  theme(plot.tag = element_text(face = "bold"))

saveRDS(ps_stromal_stemness, "shensi-figure/fs/ps_stromal_stemness.rds")

ggsave(here("output/figure/fs_stromal_stemness.tiff"), 
       ps_stromal_stemness,
       width = 5, height = 3)
ps_ici <- Reduce("+", p_immune_features_supp[1:4]) +
  plot_layout(nrow = 1) +
  plot_annotation(tag_level = "A") &
  theme(plot.tag = element_text(face = "bold"))
ggsave(here("output/figure/fs_ici.tiff"), 
       ps_ici,
       width = 10, height = 4)
  

# Increased ITH associates with worse clinical outcomes or lower efficacy of
# immunomodulator (IM) therapy in a number of cancer types (McGranahan et al.,
# 2016; Morris et al., 2016).
p5_immune_feature <- p_immune_features_main[c(
  "leukocyte_fraction", "lymphocyte_infiltration_score", "immune_score", 
  "estimate_score", "estimate_purity", "stromal_fraction",
  "intratumor_heterogeneity", "aneuploidy_score", "hrd_score",
  "dna_stemness")
]


# geneset_list <- list(geneset_EMT,geneset_TGF,geneset_WNT)
ssgsea_data <- readRDS(here("shensi-figure/ssgeea_data.rds"))
load(here("shensi-figure/geneset_ssgeea.rda"))

enriched_gs <- data.frame(gs = c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_MYC_TARGETS_V1",
  "KEGG_GLUTATHIONE_METABOLISM",
  "HALLMARK_MTORC1_SIGNALING",
  "KEGG_WNT_SIGNALING_PATHWAY",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY"),
  type = rep(c("High-risk", "Low-risk"), each = 5))
enriched_gs_gene <- dplyr::filter(geneset_gsea, term %in% enriched_gs$gs)
enriched_gs_gene$term <- droplevels(enriched_gs_gene$term)
enriched_gs_ls <- split(
  enriched_gs_gene$gene,
  enriched_gs_gene$term
)
ps_gs_gsea <- plot_gs_gsea(ssgsea_data, enriched_gs_ls, clin_filter)
saveRDS(ps_gs_gsea, "shensi-figure/fs/ps_10_ssgsea.rds")
ggsave(here("output/figure/ps_10gs_ssgsea.pdf"), ps_gs_gsea, width = 10, height = 8)



# p5 <- (p5_cell_prop + p5_is1 + p5_is2 + plot_layout(nrow = 1, widths = c(5, 1, 1)))/
#   (Reduce("+", p5_immune_feature) + plot_layout(nrow = 2)) /
#   ((p5_three_pathway + theme(legend.position = "none")) + 
#      (p5_gsea_low) + (p5_gsea_high) + plot_layout(nrow = 1)) +
#   plot_annotation() &
#   theme(plot.tag = element_text(face = "bold"))

# p5 <- (p5_cell_prop + (p5_three_pathway + theme(legend.position = "none")) + 
#          plot_layout(nrow = 1, widths = c(5, 1))) /
#   (Reduce("+", p5_immune_feature) + plot_layout(nrow = 2)) /
#   (ggplotify::as.ggplot(p5_gsea_low) | ggplotify::as.ggplot(p5_gsea_high)) +
#   plot_layout(nrow = 3, heights = c(0.8, 1.5, 1)) +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(face = "bold"))
p5 <- (p5_cell_prop + theme(legend.position = "right")) /
  (Reduce("+", p5_immune_feature) + plot_layout(nrow = 2)) +
  plot_layout(nrow = 2, heights = c(0.8, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("output/figure/f5.pdf", p5, width = 18, height = 14)
ggsave("output/figure/f5.tiff", p5, dpi = 300, width = 18, height = 14)

ggsave(here("output/figure/fs_three_pathway.tiff"), p5_three_pathway,
       width = 5, height = 4)

  # (ggplotify::as.ggplot(p5_gsea_low) | ggplotify::as.ggplot(p5_gsea_high)) +
  # plot_layout(nrow = 3, heights = c(0.8, 1.5, 1)) +
  # plot_annotation(tag_levels = "A") &
  # theme(plot.tag = element_text(face = "bold"))

# Figure 6 immune subtype -------------------------------------------

load(here("shensi-figure/gsea_dt.rda"))
library(clusterProfiler)

geneset_gsea <- readRDS(here("shensi-figure/geneset_gsea.rds"))
geneset_gsea <- filter(
  geneset_gsea,
  grepl("(HALLMARK)|(KEGG)", term)) |> 
  mutate(term = droplevels(term))

p6_gsea_low <- plot_gsea_low(geneList_low,geneset_gsea)

high_gsea_res <- readRDS(here("shensi-figure/gsea_high-risk.rds"))
p6_gsea_high <- print(enrichplot::gseaplot2(
  high_gsea_res,
  c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_MYC_TARGETS_V1",
    "KEGG_GLUTATHIONE_METABOLISM",
    "HALLMARK_MTORC1_SIGNALING",
    "KEGG_WNT_SIGNALING_PATHWAY"), 
  pvalue_table = FALSE, 
  title = paste("Enriched in High-risk group")))

chisq_test_IS <- readRDS(here("shensi-figure/chisq_test_IS.rds"))
chisq_test_IE <- readRDS(here("shensi-figure/chisq_test_IE.rds"))
p6_is2 <- ggbarstats(
  chisq_test_IS, Subtype, Group,
  type = "nonparametric",
  bf.message = FALSE,
  # ggtheme = ggplot2::theme_classic(),
  proportion.test = FALSE,
  legend.title = "") + 
  scale_fill_aaas() + 
  scale_x_discrete(breaks = c("high", "low"), 
                   labels = c("High-risk", "Low-risk")) +
  labs(x = NULL)
pval_is2 <- extract_stats(p6_is2)$subtitle_data$p.value |> 
  formatC(format = "e", digits = 2)
p6_is2 <- p6_is2 + labs(
  subtitle = as.expression(bquote(italic(P) ~"=" ~ .(pval_is2)))) +
  theme(plot.subtitle = element_text(hjust = 0.5))

p6_is1 <- ggbarstats(
  chisq_test_IE[[1]], 
  Subtype, Group, 
  type = "nonparametric",
  bf.message = FALSE,
  legend.title = "",
  proportion.test = FALSE) + 
  scale_fill_aaas() +
  scale_x_discrete(breaks = c("high", "low"), 
                   labels = c("High-risk", "Low-risk")) +
  labs(x = NULL)
pval_is1 <- extract_stats(p6_is1)$subtitle_data$p.value |> 
  formatC(format = "e", digits = 2)
p6_is1 <- p6_is1 + labs(
  subtitle = as.expression(bquote(italic(P) ~"=" ~ .(pval_is1)))) +
  theme(plot.subtitle = element_text(hjust = 0.5))
# p6 <- p6_is1 + p6_is2 + plot_layout(nrow = 1) +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(face = "bold"))
# p6 <- (ggplotify::as.ggplot(p6_gsea_low) + p6_is2 + plot_layout(nrow = 1, widths = c(3, 1))) /
#   (ggplotify::as.ggplot(p6_gsea_high) + p6_is1 + plot_layout(nrow = 1, widths = c(3, 1))) +
#   plot_annotation(tag_levels = "A") & 
#   theme(plot.tag = element_text(face = "bold"))
p6 <- p6_is2 + p6_is1 + 
  ggplotify::as.ggplot(p6_gsea_low) + ggplotify::as.ggplot(p6_gsea_high) +  
  plot_layout(nrow = 2, byrow = FALSE, widths = c(1, 3)) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold"))
ggsave("output/figure/f6.pdf", p6, width = 14, height = 12)
ggsave("output/figure/f6.tiff", p6, width = 14, height = 12, dpi = 300)

# Figure 7 -------------------------------------------------------------------
load(here("shensi-figure/ICG29_exp.rda"))
# TIDE_filter <- readRDS(here("shensi-figure/TIDE_res.rds"))
# MSIsensor_score <- readRDS(here("shensi-figure/MSI.rds"))

igc_p <- read.csv(here("shensi-figure/IGC_deseq2.csv"), 
                  row.names = 1)
igcs <- igc_p$SYMBOL
igc_count <- filter(gene_ICG_all, Symbol %in% igcs)

p7_igc <- plot_icg_box(igc_count, igc_p)

p_immune_features_main$TMB +
  p_immune_features_main$tide_score +
  p7_igc + plot_layout(nrow = 1, widths = c(1, 1, 4))

p7_stages <- lapply(
  paste("Stage", c("I", "III", "IV")),
  plot_km_stage,
  group_text_RS = clin_filter
)

p7_chemo <- lapply(
  c("high", "low"),
  plot_km_chemo,
  group_text_RS = clin_filter
)
p7_surv <- c(p7_stages, p7_chemo)

p7 <- (p_immune_features_main$TMB +
  p_immune_features_main$tide_score +
  p7_igc + p7_surv[[1]] +
  plot_layout(nrow = 1, widths = c(1, 1, 4, 2))) /
  (Reduce("+", p7_surv[2:length(p7_surv)]) + plot_layout(nrow = 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(here("output/figure/f7.pdf"), p7, width = 16.8, height = 8)
ggsave(here("output/figure/f7.tiff"), p7, width = 16.8, height = 8)


