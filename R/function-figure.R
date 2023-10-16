plot_km <- function(cpg_marker,
                    cpg_marker_coef,
                    beta,
                    clic,
                    risk_score_cutoff,
                    plot_title = "") {
  stopifnot("sample" %in% names(clic))
  stopifnot("sample" %in% names(beta))

  risk_score <- vapply(
    beta[cpg_marker] |> t() |> as.data.frame(),
    \(x) sum(x * cpg_marker_coef),
    FUN.VALUE = numeric(1)) |>
    data.frame(risk_score= _) |>
    mutate(risk_group = ifelse(risk_score  > risk_score_cutoff, "High", "Low"))
  risk_score$sample <- beta$sample
  df <- inner_join(clic, risk_score, by = c("sample" = "sample"))

  mul_cox_model_lasso<- as.formula(paste0('Surv(OS.time, OS)~',"risk_group"))
  mul_cox_lasso <- coxph(mul_cox_model_lasso,data = df )
  cox_lasso <- summary(mul_cox_lasso)
  mul_HR_lasso <- round(cox_lasso$conf.int[2],2)
  mul_5_lasso <- round(cox_lasso$conf.int[3],2)
  mul_95_lasso <- round(cox_lasso$conf.int[4],2)
  hr_label <- paste(
    paste("HR:", mul_HR_lasso),
    paste("(","95%CI, ",
          round(cox_lasso$conf.int[3],2),"-",
          round(cox_lasso$conf.int[4],2),")",sep = "")
  )

  fit <- survminer::surv_fit(survival::Surv(OS.time, OS) ~ risk_group, data = df)
  pval_info <- surv_pvalue(fit)
  pval <- pval_info$pval.txt
  p <- survminer::ggsurvplot(fit, pval = FALSE,pval.coord = c(0, 0.2),
                             pval.size = 4,xlab = "Overall survival (years)",
                             legend.labs = c("High-risk","Low-risk"),
                             legend.title = "") # risk_table = TRUE
  p <- p$plot +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_aaas() +
    labs(y = "Overall survival probability", color = NULL, title = plot_title) +
    annotate("text", x = 0.08, y = 0.04, size = 3,
             label = paste0(pval, "\n", hr_label),
             hjust = 0) +
    theme(legend.position = c(0.8, 0.9))

  print(p)
}


plot_AUC <- function(clin_risk,quan = 0.75,
                     color = ggsci::pal_aaas()(2)[1],
                     plot_title = ""){
  library(pROC)
  # ROC计算
  df <- clin_risk %>%
    dplyr::select(risk_score,OS,OS.time,group)
  df <- na.omit(df)
  df <- filter(df,df$OS.time <= quantile(df$OS.time,quan))

  rocobj <- roc(df$OS, df$risk_score,smooth = F)
  # 计算临界点/阈值
  cutOffPoint <- coords(rocobj, "best")
  cutOffPointText <- paste0(round(cutOffPoint[1],3),"(",round(cutOffPoint[2],3),",",round(cutOffPoint[3],3),")")

  # 计算AUC值
  auc <- sprintf("%0.2f", auc(rocobj)[1])#round(auc(rocobj)[1],2)
  # AUC的置信区间
  auc_low<-ci(rocobj,of="auc")[1]
  auc_high<-ci(rocobj,of="auc")[3]

  # 计算置信区间
  ciobj <- ci.se(rocobj,specificities=seq(0, 1, 0.01))
  data_ci<-ciobj[1:101,1:3]
  data_ci<-as.data.frame(data_ci)
  x=as.numeric(rownames(data_ci))
  data_ci<-data.frame(x,data_ci)

  # 绘图
  auc1 <- substr(as.character(auc) , 1 , 6)
  label <- paste("AUC = " , auc1 , sep = "")
  p <- ggroc(rocobj,
             color=color,
             size = 1,
             legacy.axes = TRUE)+
    theme_classic()+
    # theme(title = element_text(size=22,face = "bold"))+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
                 colour='grey',
                 linetype = 'dotdash',
                 linewidth = 0.4) +
    geom_ribbon(data = data_ci,
                aes(x= 1 - x, ymin = X2.5.,ymax = X97.5.),
                fill = 'lightblue',
                alpha=0.5)+
    annotate("text", x=0.7, y=0.10, label=label, size = 3) +
    labs(x = "False Positive Rate",
       y = "True Positive Rate",
       title = plot_title)

  print(p)
}


plot_AUC_test <- function(group_text,
                          quan = quantile(clin_filter$OS.time,0.75),
                          plot_title = "",
                          color = ggsci::pal_aaas()(2)[1]){
  library(pROC)
  # ROC计算
  df <- group_text %>%
    dplyr::select(OS,OS.time,risk_group,risk_score)
  df <- na.omit(df)
  df <- filter(df,df$OS.time <= quan)

  rocobj <- roc(df$OS, df$risk_score,smooth = F)
  # 计算临界点/阈值
  cutOffPoint <- coords(rocobj, "best")
  cutOffPointText <- paste0(round(cutOffPoint[1],3),"(",round(cutOffPoint[2],3),",",round(cutOffPoint[3],3),")")

  # 计算AUC值
  auc<-sprintf("%0.2f",auc(rocobj)[1])
  # AUC的置信区间
  auc_low<-ci(rocobj,of="auc")[1]
  auc_high<-ci(rocobj,of="auc")[3]

  # 计算置信区间
  ciobj <- ci.se(rocobj,specificities=seq(0, 1, 0.01))
  data_ci<-ciobj[1:101,1:3]
  data_ci<-as.data.frame(data_ci)
  x=as.numeric(rownames(data_ci))
  data_ci<-data.frame(x,data_ci)

  # 绘图
  auc1 <- substr(as.character(auc) , 1 , 6)
  label <- paste("AUC = " , auc1 , sep = "")

  p <- ggroc(rocobj,
             color=color,
             size = 1,
             legacy.axes = TRUE)+
    theme_classic()+
    # theme(title = element_text(size=22,face = "bold"))+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
                 colour='grey',
                 linetype = 'dotdash',
                 linewidth = 0.4) +
    geom_ribbon(data = data_ci,
                aes(x= 1 - x, ymin = X2.5.,ymax = X97.5.),
                fill = 'lightblue',
                alpha=0.5)+
    annotate("text", x=0.7, y=0.10, label=label, size = 3) +
    labs(x = "False Positive Rate",
         y = "True Positive Rate",
         title = plot_title)
  print(p)
}


AUC_compare_cpg <- function(ls_cpgmodel){
  library(pROC)
  # ROC计算
  rocobj1 <- roc(ls_cpgmodel[[1]]$OS, ls_cpgmodel[[1]]$risk_score,smooth = F)
  rocobj2 <- roc(ls_cpgmodel[[2]]$OS, ls_cpgmodel[[2]]$risk_score,smooth = F)
  rocobj3 <- roc(ls_cpgmodel[[3]]$OS, ls_cpgmodel[[3]]$risk_score,smooth = F)
  rocobj4 <- roc(ls_cpgmodel[[4]]$OS, ls_cpgmodel[[4]]$risk_score,smooth = F)
  rocobj5 <- roc(ls_cpgmodel[[5]]$OS, ls_cpgmodel[[5]]$risk_score,smooth = F)
  # 计算AUC值
  auc1 <- sprintf("%0.2f",auc(rocobj1)[1])
  auc2 <- sprintf("%0.2f",auc(rocobj2)[1])
  auc3 <- sprintf("%0.2f",auc(rocobj3)[1])
  auc4 <- sprintf("%0.2f",auc(rocobj4)[1])
  auc5 <- sprintf("%0.2f",auc(rocobj5)[1])
  # AUC的置信区间
  auc1_low<-ci(rocobj1,of="auc")[1]
  auc1_high<-ci(rocobj1,of="auc")[3]

  auc2_low<-ci(rocobj2,of="auc")[1]
  auc2_high<-ci(rocobj2,of="auc")[3]

  auc3_low<-ci(rocobj3,of="auc")[1]
  auc3_high<-ci(rocobj3,of="auc")[3]

  auc4_low<-ci(rocobj4,of="auc")[1]
  auc4_high<-ci(rocobj4,of="auc")[3]

  auc5_low<-ci(rocobj5,of="auc")[1]
  auc5_high<-ci(rocobj5,of="auc")[3]
  # 绘图
  auc1_1 <- substr(as.character(auc1) , 1 , 6)
  auc2_1 <- substr(as.character(auc2) , 1 , 6)
  auc3_1 <- substr(as.character(auc3) , 1 , 6)
  auc5_1 <- substr(as.character(auc5) , 1 , 6)
  aucs <- paste0("AUC = ", c(auc1_1, auc2_1, auc3_1, auc5_1))

  model_nms <- c("FCPRscore", "ChenD", "GhafarpourV", "ShenS")
  legend_label <- paste0(
    model_nms, c("", rep("_DNAm", 4)),
    " (n = ",
    c(4, 5, 14, 7),
    ", ",
    aucs,
    ")"
  )

  p <- ggroc(list(FCPRscore = rocobj1,
                  Ch_DNAm = rocobj2,
                  Gh_DNAm = rocobj3,
                  # Xu_DNAm = rocobj4 ,
                  Sh_DNAm = rocobj5),
             size = 0.8,
             legacy.axes = TRUE )+
    ggtitle("FCPRscore vs. DNAm-based models") +
    theme_classic()+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour='grey', linetype = 'dotdash') +
    theme(legend.position = c(0.7,0.3),legend.title = element_blank(),
          legend.background = element_rect(fill = alpha("white", 0))) +
    scale_color_aaas(alpha = 0.7, labels = legend_label) +
    labs(x = "False Positive Rate",
         y = "True Positive Rate")
}


AUC_compare_gene <- function(ls_genemodel){
  library(pROC)
  # ROC计算
  rocojb_ls <- lapply(ls_genemodel, \(x) roc(x$OS, x$risk_score, smooth = FALSE))
  auc_vals <- vapply(rocojb_ls,
                     \(x) sprintf("%0.2f", auc(x)[1]) |>  substr(1, 6),
                     FUN.VALUE = character(1))
  model_nms <- c("FCPRscore", "BasuB", "HanY", "ZhuQ", "ShiC", "LiuB", "ChenY")
  aucs <-  paste0("AUC = ", auc_vals)
  legend_label <- paste0(
    model_nms, c("", rep("_RNA", 6)), " (n = ",
    c(4, 8, 10, 5, 7, 8, 3),
    ", ",
    aucs,
    ")"
  )
  # legend_label <- paste0(model_nms, " (", aucs, ")")

  p <- ggroc(rocojb_ls,
             size=0.8,legacy.axes = TRUE )+
    ggtitle("FCPRscore vs. RNA-based models") +
    theme_classic()+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour='grey', linetype = 'dotdash') +
    theme(legend.position = c(0.7, 0.3),legend.title = element_blank(),
          legend.background = element_rect(fill = alpha("white", 0))) +
    scale_color_aaas(alpha = 0.7, labels = legend_label) +
    labs(x = "False Positive Rate",
         y = "True Positive Rate")
}


source(here("shensi-figure/draw-forestplot.R"))
forest_plot <- function(result1){
  p <- forestplot(result1[,c(1, 7,8)],
             mean=result1[,4],
             lower=result1[,5],
             upper=result1[,6],
             zero=1,
             boxsize=0.6,
             # labeltext =   ,
             # graph.pos= "right" ,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2,columns=c(1:4)),
                             "2" = gpar(lty=2,columns=c(1:4)),
                             "3" = gpar(lty=1,lwd=2,columns=c(1:4)),
                             "4" = gpar(lty=1,lwd=2,columns=c(1:4)),
                             "7" = gpar(lty=1,lwd=2,columns=c(1:4)),
                             "10" = gpar(lty=1,lwd=2,columns=c(1:4)),
                             "15" = gpar(lty=1,lwd=2,columns=c(1:4)),
                             "18" = gpar(lty=1,lwd=2,columns=c(1:4)),
                             "21" = gpar(lty=1,lwd=2,columns=c(1:4)),
                             "25"= gpar(lwd=2,lty=1,columns=c(1:4)),
                             "30"= gpar(lwd=2,lty=1,columns=c(1:4)),
                             "35"= gpar(lwd=2,lty=1,columns=c(1:4)) ,
                             "40"= gpar(lwd=2,lty=1,columns=c(1:4)),
                             "43"= gpar(lwd=2,lty=1,columns=c(1:4))),
             graphwidth = unit(.25,"npc"),
             xticks=c(-2,1,3,7,12,21) ,
             txt_gp=fpTxtGp(
               label=gpar(cex=1),
               ticks=gpar(cex=1),
               xlab=gpar(cex=1.5),
               title=gpar(cex=2)),
             lwd.zero=1,
             lwd.ci=1.5,
             lwd.xaxis=2,
             lty.ci=1.5,
             ci.vertices =T,
             ci.vertices.height=0.2,
             clip=c(0.1,8),
             ineheight=unit(8, 'mm'),
             line.margin=unit(8, 'mm'),
             colgap=unit(6, 'mm'),
             fn.ci_norm="fpDrawDiamondCI",
             mar = unit(rep(c(1, 0), times = 2), "mm"),
             # title="Independent Predict Factor",
             col=fpColors(box ='#021eaa',
                          lines ='#021eaa',
                          zero = "black"))


  draw_fp(p) |>
    cowplot::ggdraw()
}


nomo_coxyear <- function(clin_filter){
  library(rms)
  d2 <- clin_filter %>%
    dplyr::select(age,gender,pack_years_smoked,alcohol_history,HPV_status,risk_score,OS,OS.time)
  d2 <- na.omit(d2)
  d2$alcohol_history[which(d2$alcohol_history == "Not Reported")] <- "No"
  d2$pack_years_smoked[which(is.na(d2$pack_years_smoked) == TRUE)] <- 0
  d2$OS <- factor(d2$OS,
                  levels = c(0,1),
                  labels = c("Live", "Dead"))
  dd=datadist(d2)
  options(datadist = dd)
  f2 <- psm(Surv(OS.time,OS) ~ risk_score + age + gender + pack_years_smoked + alcohol_history + HPV_status ,# + prior_treatment
            data =  d2, dist='lognormal')
  surv <- Survival(f2) # 构建生存概率函数

  nom <- nomogram(f2, fun=list(function(x) surv(1, x),
                               function(x) surv(2, x),
                               function(x) surv(3, x)),
                  funlabel=c("1-year Survival probability",
                             "2-year Survival probability",
                             "3-year Survival probability"),
                  nint = 6)
  # save base r plot as ggplot

  plot(nom, xfrac=.4, cex.axis = 0.8)
  p <- recordPlot()
  cowplot::ggdraw(p)
}


calibration <- function(clin_filter){
  library(rms)
  d2 <- clin_filter %>%
    dplyr::select(age,gender,pack_years_smoked,alcohol_history,HPV_status,risk_score,OS,OS.time)
  d2 <- na.omit(d2)
  d2$alcohol_history[which(d2$alcohol_history == "Not Reported")] <- "No"
  d2$pack_years_smoked[which(is.na(d2$pack_years_smoked) == TRUE)] <- 0
  d2$OS.time <- d2$OS.time * 365

  dd=datadist(d2)
  options(datadist = dd)
  options(na.action="na.delete")

  f1 <- cph(formula = Surv(OS.time,OS) ~ risk_score + age + gender + pack_years_smoked + alcohol_history + HPV_status ,
            data =  d2, x=T,y=T,surv = T,na.action=na.delete,time.inc = 365)
  cal1 <- calibrate(f1, cmethod="KM", method="boot",u=365,m=50,B=1000)

  f2 <- cph(formula = Surv(OS.time,OS) ~ risk_score + age + gender + pack_years_smoked + alcohol_history + HPV_status ,
            data =  d2, x=T,y=T,surv = T,na.action=na.delete,time.inc = 730)
  cal2 <- calibrate(f2, cmethod="KM", method="boot",u=730,m=50,B=1000)

  f3 <- cph(formula = Surv(OS.time,OS) ~ risk_score + age + gender + pack_years_smoked + alcohol_history + HPV_status ,
            data =  d2, x=T,y=T,surv = T,na.action=na.delete,time.inc = 1095)
  cal3 <- calibrate(f3, cmethod="KM", method="boot",u=1095,m=50,B=1000)

  calibrate2df <- function(cali, year) {
    df <- data.frame(pred = cali[, "mean.predicted"],
                     cal = cali[, "KM"],
                     cal_corrected = cali[, "KM.corrected"],
                     se = cali[, "std.err"])
    df$time <- year

    ciupper <- function(surv, d) ifelse(surv==0, 0, pmin(1, surv*exp(d)))
    cilower <- function(surv, d) ifelse(surv==0, 0, surv*exp(-d))

    df$upper <- ciupper(df$cal, 1.959964*df$se)
    df$lower <- cilower(df$cal, 1.959964*df$se)

    df
  }

  df1 <- calibrate2df(cal1, "1-year")
  df2 <- calibrate2df(cal2, "2-year")
  df3 <- calibrate2df(cal3, "3-year")
  df <- list(df1, df2, df3) |>
    dplyr::bind_rows()

  p <- ggplot(df, aes(pred)) +
    geom_errorbar(aes(ymin = lower, ymax =  upper)) +
    facet_wrap(vars(time), scales = "free") +
    geom_point(aes(y = cal)) +
    geom_point(aes(y = cal_corrected), shape = 4) +
    geom_line(aes(pred, cal), col = ggsci::pal_aaas()(1), linewidth = 0.6) +
    geom_abline(slope = 1, color = "grey", linetype = "dotdash", linewidth = 0.6) +
    theme_bw() +
    labs(x = "Predicted probability", y = "Actua probability")

  print(p)
}


plot_nomogram_auc <- function(clin_risk, quan = 0.75){
  library(pROC)
  # ROC计算
  df <- clin_risk
  df <- filter(df,df$OS.time <= quantile(df$OS.time,quan))

  df$HPV_status[which(is.na(df$HPV_status) == TRUE)] <- "HNSC_HPV-"
  df$clinical_stage[which(is.na(df$clinical_stage) == TRUE)] <- "Stage I"

  df$HPV_status <- ifelse(df$HPV_status == "HNSC_HPV+" , 1 ,0)
  df$prior_treatment <- ifelse(df$prior_treatment == "Yes" ,1 ,0)
  df$alcohol_history <- ifelse(df$alcohol_history == "Yes" , 1 , 0)
  df$smoked_history <- ifelse(df$smoked_history == "high (>42py)" , 1 ,
                              ifelse(df$smoked_history == "median (<42py)" , 0 ,-1))
  df$clinical_stage <- ifelse(df$clinical_stage == "Stage I",-1 ,
                              ifelse(df$clinical_stage == "Stage II" ,0 ,
                                     ifelse(df$clinical_stage == "Stage III" ,1 , 2)))

  combine_all <- lrm(OS ~ risk_score + HPV_status + prior_treatment +
                       alcohol_history + smoked_history + clinical_stage,data = df)
  df$combine_all <- combine_all[["linear.predictors"]]

  rocobj_ls <- lapply(
    c("combine_all", "risk_score", "HPV_status",
      "prior_treatment", "clinical_stage"),
    \(x) roc(df$OS, df[[x]], smooth = FALSE)
  )

  auc_vals <- vapply(rocobj_ls,
                    \(x) sprintf("%0.2f", auc(x)[1]) |>  substr(1, 6),
                     FUN.VALUE = character(1))
  model_nms <- c("Combined", "FCPRscore", "HPV_status", "Prior_treatment",
                 "Clinical_stage")
  aucs <-  paste0("AUC = ", auc_vals)
  legend_label <- paste0(
    model_nms, " (",
    aucs,
    ")"
  )

  p <- ggroc(rocobj_ls,
             size=0.8,legacy.axes = TRUE )+
    # ggtitle("FCPRscore vs. RNA-based models") +
    theme_classic()+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour='grey', linetype = 'dotdash') +
    theme(legend.position = c(0.7, 0.3),legend.title = element_blank(),
          legend.background = element_rect(fill = alpha("white", 0))) +
    scale_color_aaas(alpha = 0.7, labels = legend_label) +
    labs(x = "False Positive Rate",
         y = "True Positive Rate")

  print(p)
}

compare_feature <- function(df, feature) {
  df <- dplyr::select(df, one_of(c("risk_group", feature)))
  df <- na.omit(df)

  p <- ggplot(df, aes(risk_group, .data[[feature]])) +
    geom_boxplot(aes(fill = risk_group), show.legend = FALSE, outlier.shape = NA) +
    stat_compare_means(comparisons = list(c("High-risk", "Low-risk")),
                       tip.length = 0,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                        symbols = c("***", "**", "*", ".", ""))) +
    scale_fill_aaas() +
    labs(x = NULL, y = feature) +
    theme_classic()

  print(p)
}

ggbetweenstats2 <- function(df, y, ...) {
  p <- ggstatsplot::ggbetweenstats(
    df, risk_group, !! y,
    bf.message = FALSE,
    centrality.plotting = FALSE,
    results.subtitle = FALSE,
    pairwise.comparisons = FALSE,
    ggtheme = theme_ggstatsplot(),
    ...
  )
  # pval <- ggstatsplot::extract_stats(p)$subtitle_data$p.value
  p + geom_signif(
    comparisons = list(c("High-risk", "Low-risk")),
    map_signif_level = function(p) symnum(
      p,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", "ns")) |>
      as.character()) +
    labs(x = NULL) +
    scale_color_aaas()
}


cell_proportion <- function(merge.data){
  p <- ggplot(merge.data, aes(Cell, Expression, fill = Group)) +
      geom_boxplot(outlier.shape = NA) +
      stat_compare_means(aes(group = Group),
                         method="wilcox.test",
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                          symbols = c("***", "**", "*", ".", "")),
                         label = "p.signif") +
      scale_y_continuous(trans = "sqrt") +
      scale_fill_aaas(breaks = c("high", "low"),
                      labels = c("High-risk", "Low-risk")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      labs(x = NULL, y = "Cell proportion", fill = NULL)

  print(p)
}


plot_gsea_low <- function(gene_ls ,gs, z = 0.5){
  List = gene_ls[order(gene_ls$log2FoldChange,decreasing = TRUE),]
  List = List[,4]
  names(List) <- as.character(gene_ls[,2])
  gseares <- GSEA(List,
                  TERM2GENE  = gs,
                  verbose = F,
                  pvalueCutoff = z)
  enriched_gs <- dplyr::filter(
    gseares@result,
    NES > 1,
    pvalue < 0.05,
    p.adjust < 0.25
  )
  gs_ids <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
              "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
              "HALLMARK_INTERFERON_ALPHA_RESPONSE",
              "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
              "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY"
              )
  p <- enrichplot::gseaplot2(gseares, gs_ids,
                 pvalue_table = FALSE,
                 title = paste("Enriched in Low-risk group"))
  print(p)
}

plot_gs_gsea <- function(expr, geneset_list,clin_filter) {
  if (!inherits(expr, "matrix")) expr <- as.matrix(expr)
  ssgsea_res <- gsva(expr, geneset_list, method = "ssgsea", min.sz = 1) |>
    t() |>
    as.data.frame()

  ssgsea_res$sample <- rownames(ssgsea_res)
  ssgsea_res$patient <- gsub("^(.*)-\\d+\\w+$", "\\1", ssgsea_res$sample)
  group_ssgsea <- inner_join(
    ssgsea_res, clin_filter,
    by = c("patient" = "barcode"))
  gs_nms <- names(geneset_list)
  tidy_group_ssgsea <- tidyr::pivot_longer(
    group_ssgsea,
    cols = all_of(gs_nms),
    names_to = "geneset",
    values_to = "enrichment_score"
  )
  tidy_group_ssgsea$type <- ifelse(
    tidy_group_ssgsea$geneset %in% c(
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
      "HALLMARK_MYC_TARGETS_V1",
      "KEGG_GLUTATHIONE_METABOLISM",
      "HALLMARK_MTORC1_SIGNALING",
      "KEGG_WNT_SIGNALING_PATHWAY"),
    "Enriched in High-risk group",
    "Enriched in Low-risk group")
  tidy_group_ssgsea$geneset <- gsub("(HALLMARK_)|(KEGG_)", "",
                                    tidy_group_ssgsea$geneset)
  p <- ggplot(tidy_group_ssgsea, aes(geneset, enrichment_score, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(vars(type), scales = "free") +
    stat_compare_means(aes(group = group),
                       method="wilcox.test",
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                        symbols = c("***", "**", "*", ".", "")),
                       label = "p.signif") +
    scale_fill_aaas(breaks = c("high", "low"),
                    labels = c("High-risk", "Low-risk")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = NULL, y = "Enrichment score", fill =  NULL)

  p
}

plot_icg_box <- function(igc_count, igc_p) {
  igc_p$label <- symnum(igc_p$padj,
                        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols = c("***", "**", "*", ".", NA))
  igc_p$group1 <- "High-risk"
  igc_p$group2 <- "Low-risk"
  igc_p$y.position <- 8
  igc_p <- arrange(igc_p, padj)
  igc_p$.y. <- "Expression"
  igc_p$Symbol <- igc_p$SYMBOL
  igc_p$xmin <- seq(0.8, nrow(igc_p) -1 + 0.8, 1)
  igc_p$xmax <- seq(1.2, nrow(igc_p) + 0.2, 1)
  igc_p$groups <- rep(list(c("High-risk", "Low-risk")), times = nrow(igc_p))
  igc_p$x <- seq_len(nrow(igc_p))
  igcs <- igc_p$SYMBOL

  igc_count$Symbol <- factor(igc_count$Symbol, levels = igcs)

  p <- ggplot(igc_count) +
    geom_boxplot(aes_string(x = "Symbol", y = "Expression", fill = "Group"),
                 outlier.shape = NA)  +
    stat_pvalue_manual(data = igc_p,
                       x = "Symbol",
                       label = "label",
                       tip.length = 0,
                       position = position_dodge(0.8)) +
    scale_fill_aaas(breaks = c("high", "low"),
                    labels = c("High-risk", "Low-risk")) +
    labs(x = NULL, y = "log2(FPKM+1", fill = NULL) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.9))

  print(p)
}

plot_km_stage <- function(group_text_RS ,
                          y = c("Stage I","Stage II","Stage III","Stage IV")){
  dt <- group_text_RS
  dt_OS <- dt %>%
    dplyr::select(OS,OS.time,group , clinical_stage)
  dt_OS <- na.omit(dt_OS)
  dt_OS$clinical_stagefilter <- ifelse(
    dt_OS$clinical_stage == "Stage I" , "Stage I" ,
    ifelse(dt_OS$clinical_stage == "Stage II" ,"Stage II" ,
           ifelse(dt_OS$clinical_stage == "Stage III" ,"Stage III" , "Stage IV")))
  dt_OS_filter <- filter(dt_OS , dt_OS$clinical_stagefilter == y)
  fit <- survminer::surv_fit(Surv(OS.time, OS) ~ group, data = dt_OS_filter)
  p <- survminer::ggsurvplot(
    fit, pval.coord = c(0.1, 0.2),pval = TRUE,
    pval.size = 4,xlab = "Overall survival (years)",
    break.time.by = 0.5,xlim=c(0,8),
    legend.labs = c("High-risk","Low-risk"),title = y,
    legend.title = "")
  p <- p$plot +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_aaas() +
    labs(y = "Overall survival probability", color = NULL, title = y) +
    theme(legend.position = c(0.8, 0.9))

  print(p)
}


plot_km_chemo <- function(group_text_RS , x = c("high","low")) {
  plot_title <- ifelse(x == "high", "High-risk", "Low-risk")

  dt <- group_text_RS
  dt_OS <- dt %>%
    dplyr::select(OS,OS.time,group , clinical_stage,chemical_treatment)
  dt_OS <- na.omit(dt_OS)
  dt_OS$chemical_treatment <- ifelse(dt_OS$chemical_treatment == "yes" , "yes" , "no")
  dt_OS$clinical_stagefilter <- ifelse(
    dt_OS$clinical_stage == "Stage I" , "Stage I" ,
    ifelse(dt_OS$clinical_stage == "Stage II" ,"Stage II" ,
           ifelse(dt_OS$clinical_stage == "Stage III" ,"Stage III" , "Stage IV")))
  dt_OS_filter <- filter(
    dt_OS ,
    dt_OS$clinical_stagefilter == "Stage III" |
      dt_OS$clinical_stagefilter == "Stage IV" )
  dt_OS_filter <- filter(dt_OS_filter , dt_OS_filter$group == x)
  fit <- survminer::surv_fit(Surv(OS.time, OS) ~ chemical_treatment, data = dt_OS_filter)
  p <- survminer::ggsurvplot(
    fit, pval.coord = c(0.1, 0.2),pval = TRUE,
    pval.size = 4,
    break.time.by = 0.5,xlim=c(0,8),
    legend.labs = c("Chemotherapy","No chemotherapy"))

  p <- p$plot +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_aaas() +
    labs(x = "Overall survival (years)", y = "Overall survival probability",
         color = NULL, title = plot_title) +
    theme(legend.position = c(0.8, 0.9))
  p$plot <- p$plot + scale_color_aaas()+scale_fill_aaas()
  p
}



