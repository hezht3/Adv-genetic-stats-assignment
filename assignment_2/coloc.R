require(coloc)
require(tidyverse)

setwd("D:/OneDrive - Johns Hopkins/Course/140.686.01 - Advanced Methods for Statistical Genetics and Genomics/assignment/homework_2/files")

# load dataset

gwas_sum <- read.table("ISG15_protein_summary.txt", header = TRUE)
eqtl_adipo <- read.table("ISG15_ciseQTL_Adipose_summary.txt", header = TRUE)
eqtl_blood <- read.table("ISG15_ciseQTL_Blood_summary.txt", header = TRUE)

gwas_sum <- gwas_sum %>% 
    mutate(type = "quant") %>% 
    rename(pvalues = pval) %>%
    separate(col = SNP,
             into = c("Chr", "position", "allele_1", "allele_2"),
             sep = ":")
eqtl_adipo <- eqtl_adipo %>% 
    mutate(type = "quant") %>% 
    rename(pvalues = pval) %>%
    separate(col = SNP,
             into = c("Chr", "position", "allele_1", "allele_2"),
             sep = ":")
eqtl_blood <- eqtl_blood %>% 
    mutate(type = "quant") %>% 
    rename(pvalues = pval) %>%
    separate(col = SNP,
             into = c("Chr", "position", "allele_1", "allele_2"),
             sep = ":")

# check dataset

check_dataset(gwas_sum, warn.minp = 1e-10)
check_dataset(eqtl_adipo, warn.minp = 1e-10)
check_dataset(eqtl_blood, warn.minp = 1e-10)

# colocalisation analysis - GWAS summary statistics and adipose tissues eQTL

my.res <- coloc.abf(dataset1 = eqtl_adipo,
                    dataset2 = gwas_sum)
print(my.res)

subset(my.res$results, SNP.PP.H4 > 0.01)   # likely causal variant: SNP.PP.H4 > 0.01

o <- order(my.res$results$SNP.PP.H4, decreasing = TRUE)
cs <- cumsum(my.res$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
my.res$results[o, ][1:w, ]$snp   # 95% credible set

## sensitivity analysis

sensitivity(my.res, rule = "H4 > 0.5")

# colocalisation analysis - GWAS summary statistics and whole blood eQTL

my.res2 <- coloc.abf(dataset1 = eqtl_blood,
                     dataset2 = gwas_sum)
print(my.res2)

subset(my.res2$results, SNP.PP.H4 > 0.01)   # likely causal variant: SNP.PP.H4 > 0.01

o <- order(my.res2$results$SNP.PP.H4, decreasing = TRUE)
cs <- cumsum(my.res2$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
my.res2$results[o, ][1:w, ]$snp   # 95% credible set

## sensitivity analysis

sensitivity(my.res2, rule = "H4 > 0.5")
