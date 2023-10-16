preprocess_dnam <- function(base_dir, 
                            clic, 
                            cptac3_hnsc = TRUE,
                            cptac3_hnsc_meta = NULL) {
    stopifnot(dir.exists(base_dir))
    
    if (is.null(cptac3_hnsc_meta)) {
        if (cptac3_hnsc) {
            warning("set `cptac3_hnsc` as `FALSE` for other dataset")
            cptac3_hnsc = FALSE
        }
    }
    
    rgset <- minfi::read.metharray.exp(base = base_dir)
    
    # preprocess and normalize
    probe_pval <- minfi::detectionP(rgset)
    keep <-  apply(probe_pval, 1 , \(x) all(x < 0.01))
    rgset <- rgset[keep, ]
    # 过滤掉所有探针p值的均值大于0.05的样本
    rgset <- rgset[, colMeans(probe_pval) < 0.05]
    # 预处理，背景降噪和归一化，注意，此处可以根据情况，替换成另外的算法
    gr_set_normed <- minfi::preprocessFunnorm(rgset)
    # 获取注释文件
    methy_set <- minfi::preprocessRaw(rgset)
    gm_set <- minfi::mapToGenome(methy_set)
    
    annotation <- minfi::getAnnotation(gr_set_normed)
    snps <- minfi::getSnpInfo(gr_set_normed)
    
    # 探针过滤，去除在性染色体上的探针
    sex_probe <- rownames(annotation)[annotation$chr %in% c("chrX", "chrY")]
    keep <- !(featureNames(gr_set_normed) %in% sex_probe)
    gr_set_normed <- gr_set_normed[keep, ]
    
    # 去除覆盖了SNP位点的探针，如果感觉过滤掉的探针太多，可以适当上调SNP频率, 
    # 将maf的值变大 
    gr_set_normed <- dropLociWithSnps(gr_set_normed, snps = c("SBE", "CpG"), maf = 0)
    beta  <- getBeta(gr_set_normed)
    
    if (cptac3_hnsc) {
        colnames(beta) <- cptac3_hnsc_meta$sample[
            match(colnames(beta), 
                  paste0(cptac3_hnsc_meta$analysis_submitter_id, "_noid"))
        ]
        clic <- clic[
            match(colnames(beta), clic$sample), 
        ]
    } else {
        colnames(beta) <- vapply(
            strsplit(colnames(beta), "_"),
            \(x) x[1],
            FUN.VALUE = character(1)
        )
        beta <- beta[, clic$sample]
    }
    
    beta <- data.frame(beta, check.names = FALSE) |> 
        tibble::rownames_to_column("cpg") |> 
        data.table::transpose(make.names = "cpg", keep.names = "sample") |> 
        data.frame()
    
    return(list(beta = beta, clic = clic))
}
