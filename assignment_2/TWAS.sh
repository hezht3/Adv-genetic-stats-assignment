# ------------------------------------- 1. TWAS analysis of Breast Cancer -------------------------------------

# ------------------------------------- Input: GWAS summary statistics -------------------------------------
cd "C:\Users\Johnathan He"
python ./ldsc/munge_sumstats.py \
--sumstats ./fusion_twas-master/BC17/data_BC17_1 \
--out ./fusion_twas-master/BC17/BC17 \
tar xjf ./fusion_twas-master/BC17/BC17.sumstats.gz

# ------------------------------------- Performing the expression imputation -------------------------------------
cd "C:\Users\Johnathan He\fusion_twas-master"
export PATH="C:\Program Files\R\R-4.1.2\bin:$PATH"

# CHR 1
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 1 \
--out ./OUT/BC17.1.dat   # performing the expression imputation

cat ./OUT/BC17.1.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.1.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.1.top \
--out ./OUT/BC17.1.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 1 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 2
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 2 \
--out ./OUT/BC17.2.dat   # performing the expression imputation

cat ./OUT/BC17.2.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.2.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.2.top \
--out ./OUT/BC17.2.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 2 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 3
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 3 \
--out ./OUT/BC17.3.dat   # performing the expression imputation

cat ./OUT/BC17.3.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.3.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.3.top \
--out ./OUT/BC17.3.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 3 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 4
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 4 \
--out ./OUT/BC17.4.dat   # performing the expression imputation

cat ./OUT/BC17.4.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.4.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.4.top \
--out ./OUT/BC17.4.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 4 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 5
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 5 \
--out ./OUT/BC17.5.dat   # performing the expression imputation

cat ./OUT/BC17.5.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.5.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.5.top \
--out ./OUT/BC17.5.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 5 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 6
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 6 \
--out ./OUT/BC17.6.dat   # performing the expression imputation
                         # Results in the MHC are written to ./OUT/BC17.6.dat.MHC , evaluate with caution due to complex LD structure

cat ./OUT/BC17.6.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.6.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.6.top \
--out ./OUT/BC17.6.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 6 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 7
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 7 \
--out ./OUT/BC17.7.dat   # performing the expression imputation

cat ./OUT/BC17.7.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.7.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.7.top \
--out ./OUT/BC17.7.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 7 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 8
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 8 \
--out ./OUT/BC17.8.dat   # performing the expression imputation

cat ./OUT/BC17.8.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.8.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.8.top \
--out ./OUT/BC17.8.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 8 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 9
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 9 \
--out ./OUT/BC17.9.dat   # performing the expression imputation

cat ./OUT/BC17.9.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.9.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.9.top \
--out ./OUT/BC17.9.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 9 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 10
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 10 \
--out ./OUT/BC17.10.dat   # performing the expression imputation

cat ./OUT/BC17.10.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.10.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.10.top \
--out ./OUT/BC17.10.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 10 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 11
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 11 \
--out ./OUT/BC17.11.dat   # performing the expression imputation

cat ./OUT/BC17.11.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.11.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.11.top \
--out ./OUT/BC17.11.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 11 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 12
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 12 \
--out ./OUT/BC17.12.dat   # performing the expression imputation

cat ./OUT/BC17.12.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.12.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

# CHR 13
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 13 \
--out ./OUT/BC17.13.dat   # performing the expression imputation

cat ./OUT/BC17.13.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.13.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

# CHR 14
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 14 \
--out ./OUT/BC17.14.dat   # performing the expression imputation

cat ./OUT/BC17.14.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.14.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.14.top \
--out ./OUT/BC17.14.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 14 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 15
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 15 \
--out ./OUT/BC17.15.dat   # performing the expression imputation

cat ./OUT/BC17.15.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.15.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.15.top \
--out ./OUT/BC17.15.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 15 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 16
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 16 \
--out ./OUT/BC17.16.dat   # performing the expression imputation

cat ./OUT/BC17.16.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.16.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.16.top \
--out ./OUT/BC17.16.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 16 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 17
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 17 \
--out ./OUT/BC17.17.dat   # performing the expression imputation

cat ./OUT/BC17.17.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.17.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.17.top \
--out ./OUT/BC17.17.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 17 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 18
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 18 \
--out ./OUT/BC17.18.dat   # performing the expression imputation

cat ./OUT/BC17.18.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.18.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.18.top \
--out ./OUT/BC17.18.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 18 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 19
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 19 \
--out ./OUT/BC17.19.dat   # performing the expression imputation

cat ./OUT/BC17.19.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.19.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.19.top \
--out ./OUT/BC17.19.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 19 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 20
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 20 \
--out ./OUT/BC17.20.dat   # performing the expression imputation

cat ./OUT/BC17.20.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.20.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.20.top \
--out ./OUT/BC17.20.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 20 \
--plot --locus_win 100000   # joint/conditional tests and plots

# CHR 21
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 21 \
--out ./OUT/BC17.21.dat   # performing the expression imputation

cat ./OUT/BC17.21.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.21.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

# CHR 22
Rscript FUSION.assoc_test.R \
--sumstats ./BC17/BC17.sumstats \
--weights ./TCGA-BRCA.TUMOR/TCGA-BRCA.TUMOR.pos \
--weights_dir ./TCGA-BRCA.TUMOR/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 22 \
--out ./OUT/BC17.22.dat   # performing the expression imputation

cat ./OUT/BC17.22.dat | awk 'NR == 1 || $NF < 0.05/4463' > ./OUT/BC17.22.top   # output: gene-disease association, adjusting for 4463 genes in TCGA reference

Rscript FUSION.post_process.R \
--sumstats ./BC17/BC17.sumstats \
--input ./OUT/BC17.22.top \
--out ./OUT/BC17.22.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 22 \
--plot --locus_win 100000   # joint/conditional tests and plots