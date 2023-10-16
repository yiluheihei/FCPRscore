here::i_am("script/01-methy-preprocess.R")

library(here)
library(dplyr)

# 450K: there are 12 arrays on a single physical “slide” 
# (organized in a 6 by 2 grid). Slides are organized into “plates” containing 
# at most 8 slides (96 arrays).

# methylation meta data -------------------

hnsc_info <- data.table::fread(here("data/tcga-hnsc-info.tsv"))

methy_files <- list.files(here("data/GDCdata/")) 
methy_files <- methy_files[grep(".idat" ,methy_files , fixed = TRUE)]

methy_meta <- data.frame(file = methy_files)
methy_ls <- strsplit(methy_files, "_")
methy_meta$Sentrix_ID <- vapply(methy_ls, "[", FUN.VALUE = character(1), 1)
methy_meta$Sentrix_Position <- vapply(methy_ls, "[", FUN.VALUE = character(1), 2)

methy_meta <- inner_join(methy_meta, hnsc_info, by = c("file" = "file_name"))
methy_meta$basename <- paste(
    methy_meta$Sentrix_ID , 
    methy_meta$Sentrix_Position , 
    sep = "_"
)

# pd_3 <- filter(pd_2 , pd_2$Sample_Type == "tumor" )
# pd_4 <- pd_3[!duplicated(pd_3$basename) , ]
# 
# pd_5 <- pd_2[!duplicated(pd_2$basename) , ]


# methylation 数据预处理 ------------------------------------------------------

library(minfi)

# 根据SampleSheet.csv 文件，读取所有样本的 .idat 文件
sample_info <- read.metharray.sheet(here("data/GDCdata"), 
                                    pattern="sample_sheet.csv")
#targets$Basename <- paste(myDir , targets$Sample_ID)
# basename <- sample_info[!duplicated(sample_info$Basename),]
rg_set <- read.metharray.exp(targets = sample_info)
pd <- pData(rgSet) 

# 计算探针的p值，过滤掉在任何一个样本中p值大于0.01的探针
probe_pval <- detectionP(rg_set)
keep <-  apply(probe_pval, 1 , \(x) all(x < 0.01))
rg_set <- rg_set[keep, ]

# 过滤掉所有探针p值的均值大于0.05的样本
rg_set <- rg_set[, colMeans(probe_pval) < 0.05]

# 预处理，背景降噪和归一化，注意，此处可以根据情况，替换成另外的算法
gr_set_normed <- preprocessFunnorm(rg_set)

# 获取注释文件
methy_set <- preprocessRaw(rg_set)
gm_set <- mapToGenome(methy_set)

annotation <- getAnnotation(gr_set_normed)
snps <- getSnpInfo(gr_set_normed)

# 探针过滤，去除在性染色体上的探针
sex_probe <- rownames(annotation)[annotation$chr %in% c("chrX", "chrY")]
keep <- !(featureNames(gr_set_normed) %in% sex_probe)
gr_set_normed <- gr_set_normed[keep, ]

# 去除覆盖了SNP位点的探针，如果感觉过滤掉的探针太多，可以适当上调SNP频率, 
# 将maf的值变大 
gr_set_normed <- dropLociWithSnps(gr_set_normed, snps = c("SBE", "CpG"), maf = 0)
beta  <- getBeta(gr_set_normed)
M <- getM(gr_set_normed) #, type = "beta", betaThreshold = 0.001

saveRDS(gr_set_normed , here("output/data/gr_set_normed.rds"))
saveRDS(M, here("output/data/m.rds"))
saveRDS(probe_pval, here("output/data/probe_pval.rds"))


# pca ------------------

library(TCGAbiolinks)
hnsc_biospecimen <- GDCquery_clinic("TCGA-HNSC", type = "Biospecimen")
sample_info <- inner_join(sample_info, hnsc_biospecimen,
                          by = c("sample.submitter_id" = "submitter_id"))
sample_info$barcode <- TCGAutils::TCGAbarcode(sample_info$sample.submitter_id)

saveRDS(sample_info, here("output/data/methy_sample_info.rds"))

# 没有ffpe 样本
# is_ffpe: A boolean value that denotes whether tissue samples used in the 
# analysis were formalin-fixed paraffin-embedded (FFPE).
sum(sample_info$is_ffpe)

beta <- as.data.frame(na.omit(beta))
beta <- beta[sample_info$basename]
names(beta) <- sample_info$sample.submitter_id
saveRDS(beta, here("output/data/beta.rds"))

beta_batch <- data.table::transpose(
    beta_batch,
    keep.names = "sample"
)
beta_batch$sample <- NULL
# use the 10000 most variably cpgs
cpg_sd_idx <- vapply(beta_batch, sd, FUN.VALUE = numeric(1)) |> 
    order(decreasing = TRUE)
beta_10k <- beta_batch[cpg_sd_idx[1:10000]]
beta_10k_dist <- dist(beta_10k)

pca_res <- prcomp(beta_10k, scale. = TRUE)

library(ggfortify)

# tumor vs. normal 显著差异
autoplot(pca_res,  data = sample_info, colour = "sample_type")
# p = 0.001
vegan::adonis2(beta_10k_dist ~ sample_type, sample_info)$`Pr(>F)`[1]

# array, 不需要去除batch
autoplot(pca_res,  data = sample_info, colour = "Array")
# p = 0.511
vegan::adonis2(beta_10k_dist ~ Array, sample_info)$`Pr(>F)`[1]

# slide, p = 0.001 ? remove batch?
autoplot(pca_res,  data = sample_info, colour = "Slide")
vegan::adonis2(beta_10k_dist ~ Slide, sample_info)
vegan::mrpp(beta_10k_dist, sample_info$Slide)
