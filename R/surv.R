plot_km <- function(cpg_marker, 
                    cpg_marker_coef,
                    beta, 
                    clic, 
                    risk_score_cutoff) {
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
    
    fit <- survminer::surv_fit(survival::Surv(OS.time, OS) ~ risk_group, data = df)
    p <- survminer::ggsurvplot(fit, pval = TRUE,risk.table = TRUE)
    
    list(data = df, plot = p)
}
