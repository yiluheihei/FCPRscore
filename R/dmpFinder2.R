dmpFinder2 <- function(dat, pheno, type = c("categorical", "continuous"),
                      qCutoff = 1, shrinkVar = FALSE) {
    
    # Check inputs
    type <- match.arg(type)
    if (is(dat, "MethylSet")) {
        .isMatrixBackedOrStop(dat, "dmpFinder")
        M <- getM(dat)
    } else if (is(dat, "DelayedMatrix")) {
        stop("This function does not yet support DelayedArray objects")
    } else {
        stopifnot(is.numeric(dat))
        M <- dat
        if (is.vector(M)) M <- matrix(M, nrow = 1)
    }
    n <- length(pheno)
    if (n != ncol(M)) stop("length of pheno does not equal number of samples")
    
    m2beta <- function(m) {
        return(2^m / (2^m + 1))
    }
    
    if (type == "categorical") {
        pheno <- factor(as.character(pheno))
        design <- model.matrix(~pheno)
        fit <- lmFit(M, design)
        coef1 <- fit$coefficients[, 1]
        if (shrinkVar) {
            fit <- contrasts.fit(fit, contrasts(pheno))
            fit <- eBayes(fit)
            tab <- data.frame(
                coef1 = coef1,
                diff_coef = fit$coefficients[, 1],
                f = fit[["F"]],
                pval = fit[["F.p.value"]])
        } else {
            fit1 <- lmFit(M)
            RSS1 <- rowSums2((M - fitted(fit1))^2)
            RSS <- rowSums2((M - fitted(fit))^2)
            df1 <- length(levels(pheno)) - 1
            df2 <- n - length(levels(pheno))
            Fstat <- ((RSS1 - RSS)/df1)/(RSS/df2)
            if (df2 > 1e+06) {
                F.p.value <- pchisq(df1 * Fstat, df1, lower.tail = FALSE)
            }
            else {
                F.p.value <- pf(Fstat, df1, df2, lower.tail = FALSE)
            }
            tab <- data.frame(
                coef1 = fit$coefficients[, 1],
                diff_coef = fit$coefficients[, 2],
                f = Fstat,
                pval = F.p.value)
        }
        
        tab$coef2 <- tab$coef1 + tab$diff_coef
        m1 <- tab$coef1
        m2 <- tab$coef2
        diff_beta <- m2beta(m2) - m2beta(m1)
        tab$diff_beta <- diff_beta
    }
    else if (type == "continuous") {
        design <- model.matrix(~pheno)
        fit <- lmFit(M, design)
        if (shrinkVar) {
            fit <- eBayes(fit)
            sigma <- sqrt(fit$s2.post)
            df <- fit$df.prior + fit$df.residual
        } else {
            sigma <- fit$sigma
            df <- fit$df.residual
        }
        coef <- fit$coefficients
        if (is.vector(coef)) coef <- matrix(coef, ncol = 2)
        stdev <- fit$stdev.unscaled
        if (is.vector(stdev)) stdev <- matrix(stdev, ncol = 2)
        t <- coef[, 2] / (stdev[, 2] * sigma)
        pval <- 2 * pt(abs(t), df = df, lower.tail = FALSE)
        tab <- data.frame(
            intercept = coef[, 1],
            beta = coef[, 2],
            t = t,
            pval = pval)
    }
    p0 <- siggenes::pi0.est(tab$pval[!is.na(tab$pval)])$p0
    tab$qval <- siggenes::qvalue.cal(tab$pval, p0)
    if (qCutoff < 1) tab <- tab[tab$qval <= qCutoff,]
    o <- order(tab$pval)
    tab <- tab[o, ]
    ## tab <- cbind(tab, annotate(rownames(tab)))
    if (nrow(tab) == 0) {
        message(sprintf("No significant DMPs at FDR cutoff of %s.", qCutoff))
    }
    tab
}
