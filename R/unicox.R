run_uni_cox <- function(feature, df) {
    formula <- as.formula(paste0('survival::Surv(OS.time, OS) ~ ', feature))
    fit <- survival::coxph(formula, data = df)
    fit_summary <- summary(fit)
    
    # Proportional Hazards Assumption: time independent
    ph_hypothesis_p <- survival::cox.zph(fit)$table[1, 3]
    
    ci5 <- round(fit_summary$conf.int[, 3], 2)
    ci95 <-round(fit_summary$conf.int[, 4], 2)
    res <- data.frame(feature = feature,
                      beta = fit_summary$coefficients[, 1],
                      hazard_ratio = fit_summary$coefficients[, 2],
                      ph_pval = fit_summary$coefficients[, 5],
                      ph_ci5 = ci5,
                      ph_ci95 = ci95,
                      ph_hypothesis_p = ph_hypothesis_p)
    
    res
}
