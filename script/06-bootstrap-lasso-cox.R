here::i_am("script/06-bootstrap-lasso-cox.R")

library(here)
library(dplyr)
library(glmnet)

# interaction_sig <- readRDS(here("output/data/methreg_interaction_sig.rds"))
# cpg_gene_rlm_sig <- readRDS(here("output/data/cpg_gene_rlm_res_sig.rds"))
load(here("output/data/methreg_interaction.rda"))
load(here("output/data/cpg_gene_pair_rlm.rda"))
# beta_residual <- readRDS(here("output/data/cpg_residual.rds"))
load(here("output/data/dmp_beta_tumor.rda"))
hnsc_clic <- readRDS(here("output/data/hnsc_clic.rds"))

cpg_candidate <- unique(c(interaction_sig$probeID,
                        cpg_gene_pair_rlm_sig$probe_id))

# cpg_interaction_candidate <- dplyr::select(
#     interaction_sig,
#     probe_id = probeID,
#     gene_id = target
# )
# cpg_rlm_candidate <- dplyr::select(
#     cpg_gene_pair_rlm_sig,
#     probe_id, gene_id
# )
# cpg_candidate_all <- rbind(cpg_interaction_candidate,
#                            cpg_rlm_candidate) |> 
#     dplyr::distinct()
# 
# deg <- readRDS(here("output/data/deg.rds")) 
# deg$gene_id <- strsplit(rownames(deg), "\\.") |> 
#     vapply(\(x)x[1], FUN.VALUE = character(1))
# deg_sig <- dplyr::filter(deg, 
#                          abs(log2FoldChange) > 1,
#                          padj < 0.1)
# cpg_candidate_all <- dplyr::filter(cpg_candidate_all, gene_id %in% deg_sig$gene_id)
# cpg_candidate <- unique(cpg_candidate_all$probe_id)
# 
# cpg_candidate <- c(cpg_candidate, 
#                    "cg02409878", "cg23867673", "cg01995815")

dmp_beta_surv <- inner_join(hnsc_clic, dmp_beta_tumor, 
                            by = c("barcode" = "barcode")) |> 
    tibble::column_to_rownames("barcode")
beta_surv <- dplyr::select(
    dmp_beta_surv,
    all_of(c(cpg_candidate, "OS.time", "OS"))
)
beta_surv <- na.omit(beta_surv)

cv_model <- glmnet::cv.glmnet(
    beta_surv[1:(ncol(beta_surv) - 2)] |> as.matrix(), 
    survival::Surv(beta_surv$OS.time, beta_surv$OS),
    family = "cox", 
    alpha = 1,
    nfolds = 10
)

boot_lamda <- function(df, indices) {
    d <- df[indices, ]
    x <- as.matrix(d[, 1:(ncol(d) - 2)])
    y <- survival::Surv(d$OS.time, d$OS)
    cv_model <- glmnet::cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
    
    return(cv_model$lambda.1se)
}


RNGkind("L'Ecuyer-CMRG")
# set.seed(123) # 1000, 995, 991, 980
set.seed(123)  
# seed: 1234, 2023
# # 1911 candidate cpgs: 100 minutes
# # 1140: 70 minutes
# # 963 : 70 minutes
# # 655：61 minutes
# # 658: 65 minutes
# # 678

system.time(boot_res <- boot::boot(
    data = beta_surv, statistic = boot_lamda, 
    R = 1000, ncpus = 12, parallel = "snow"
))

lambda_1se <- boot_res[["t"]]
cpg_coef <-coef(cv_model, s = c(lambda_1se[,1]))
cpg_coef <- as.matrix(cpg_coef)

cpg_nonzero_n <- apply(cpg_coef, 1, \(x) sum(x > 0)) |> 
    data.frame(n = _) |> 
    tibble::rownames_to_column("probe_id") |> 
    dplyr::arrange(desc(n))

# cpg_n <- apply(cpg_coef, 1, \(x) sum(x > 0))
# cpg_n <- data.frame(cpg = names(cpg_n), n = cpg_n)
# cpg_n <- cpg_n[order(cpg_n$n, decreasing = TRUE), ]


# multi cox ---------------------------------------------------------------

# nonzero > 980, "cg02409878", "cg23867673", "cg01995815" (cg01984743)
cpg_marker <- cpg_nonzero_n$probe_id[1:3]
perform_multicox <- function(cpg, df) {
    fml <- as.formula(
        paste0('survival::Surv(OS.time, OS) ~ ',
               paste(rev(cpg), collapse = " + "))
    )
    f <- survival::coxph(fml, data = df) |> summary()
    f_coef <- f$coefficients
    
    res <- data.frame(cpg = rownames(f_coef),
                      coef = f_coef[, "coef"],
                      hr = f_coef[, "exp(coef)"],
                      p = f_coef[, "Pr(>|z|)"],
                      lower_ci = f$conf.int[, "lower .95"],
                      upper_ci = f$conf.int[, "upper .95"])
    
    res
}

cpg_marker_cox <- perform_multicox(cpg_marker, beta_surv)
cpg_marker_cox <- cpg_marker_cox[cpg_marker, ]


# km survival -------------------------------------------------------------

# 527 samples have OS data
# na_idx <- attr(beta_surv, "na.action")
# dmp_beta_surv$barcode[-na_idx]
source(here("R/surv.R"))
risk_score <- vapply(
    beta_surv[cpg_marker] |> t() |> as.data.frame(),
    \(x) sum(x * cpg_marker_cox$coef),
    FUN.VALUE = numeric(1)
)
risk_score_cutoff <- quantile(risk_score, 0.5)

# cutoff determined using survminer::surv_cutpoint
risk_score_df <- data.frame(risk_score)
risk_score_df$barcode <- names(risk_score)
risk_score_df <- inner_join(risk_score_df, hnsc_clic)
risk_score_cp <- survminer::surv_cutpoint(
    risk_score_df,
    "OS.time", "OS",
    "risk_score",
    minprop = 0.3
)
risk_score_cutpoint <- risk_score_cp$cutpoint$cutpoint

cps <- lapply(
    seq(0.2, 0.8, 0.01),
    \(x) {
        cp <- quantile(risk_score_df$risk_score, x)
        risk_score_df$group <- ifelse(risk_score_df$risk_score >= cp,
                                      "High", "Low")
        
        fit <- survminer::surv_fit(survival::Surv(OS.time, OS) ~ group, data = risk_score_df)
        
        survminer::surv_pvalue(fit)
    }) |> bind_rows()
cps$quantile <- seq(0.2, 0.8, 0.01)
# lowest p value
arrange(cps, pval)[1, ]


plot_km(
    cpg_marker,
    cpg_marker_cox$coef,
    mutate(beta_surv, sample = rownames(beta_surv)),
    dplyr::rename(hnsc_clic, sample = barcode),
    quantile(risk_score, 0.65)
)


save(cpg_marker, cpg_marker_cox, cpg_nonzero_n,
     risk_score, risk_score_cutoff, risk_score_cutpoint,
     file = here("output/data/risk_cpg_marker.rda"))


# https://www.proteinatlas.org/
# cpg调控的基因有些与预后相关(如cg02409878:OSR2)，有些无关(cg01995815:CLSTN2)，
# 也正说明了cpg位点作为marker的重要性;
# 
# 虽然所有4个位点都是危险因素，但是调控的位点有的是favorable (如cg02409878:OSR2)，
# 有的是unfavorable (cg01984743: STC2), 说明单纯的通过基因及启动子区域位点负相关
# 可能忽略一些重要的预后相关的因素。
# 
# cg23867673:CDH23 Prognostic marker in urothelial cancer (favorable)
