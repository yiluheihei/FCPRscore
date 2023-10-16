auxfunction <- function(gene_residual, cpg_residual, gene_cpg, clic){
    probe_id <-gene_cpg[1]
    gene_id <-gene_cpg[2]
    
    gene_residual <- gene_residual[match(gene_id , rownames(gene_residual), 
                                         nomatch = 0), , drop = FALSE]
    
    cpg_residual <- cpg_residual[match(probe_id, rownames(cpg_residual),
                                       nomatch = 0), , drop = FALSE]

    
    residual_barcode <- TCGAutils::TCGAbarcode(colnames(gene_residual))
    colnames(gene_residual) <- residual_barcode
    colnames(cpg_residual) <- residual_barcode
    clic <- clic[match(residual_barcode, clic$barcode, nomatch = 0), ]
    gene_residual <- gene_residual[, match(clic$barcode, 
                                           colnames(gene_residual))]
    cpg_residual <- cpg_residual[, match(clic$barcode, 
                                         colnames(cpg_residual))]
    
    data <- data.frame(
        met_residual = cpg_residual %>% as.numeric(),
        rna_residual = gene_residual %>% as.numeric(),
        stage_group = clic$stage_group,
        trans_group = clic$trans_group
    )
    
    reg <- MASS::rlm(
        formula = rna_residual ~ met_residual + stage_group + trans_group, 
        data = data,
        psi = MASS::psi.bisquare,
        maxit = 100) |> 
        summary() |> coef() |> data.frame()
    
    degree_freedom <- nrow(data) - 3
    reg$pval <- 2 * (1 - pt(abs(reg$t.value), df = degree_freedom))
    
    quant_pval <- reg[-1, 4,drop = FALSE] |> 
        t() |>
        as.data.frame()
    colnames(quant_pval) <- paste0("RLM_",colnames(quant_pval),"_pvalue")
    
    quant_estimate <- reg[-1, 1, drop = FALSE] |> 
        t() |>
        as.data.frame()
    colnames(quant_estimate) <- paste0("RLM_",colnames(quant_estimate),"_estimate")
    
    return(cbind(quant_pval, quant_estimate))
}  
