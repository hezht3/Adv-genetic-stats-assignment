# TWAS FUSION tutorial

# Installation

## Download and unpack the FUSION software package from github
cd "C:\Users\Johnathan He"
wget https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip -o master.zip
cd fusion_twas-master

## Download and unpack the (1000 Genomes) LD reference data
wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar xjvf LDREF.tar.bz2

## Download and unpack the plink2R library (by Gad Abraham)
wget https://github.com/gabraham/plink2R/archive/master.zip
unzip -o master.zip

## Launch R and install required libraries
export PATH="C:\Program Files\R\R-4.1.2\bin:$PATH"
R
install.packages(c('optparse','RColorBrewer'))
install.packages(c("Rcpp", "RcppEigen"))
devtools::install_local('D:/OneDrive - Johns Hopkins/Course/140.686.01 - Advanced Methods for Statistical Genetics and Genomics/assignment/homework_2/plink2R/plink2R/')
quit()

# Typical analysis and output
wget https://data.broadinstitute.org/alkesgroup/FUSION/SUM/PGC2.SCZ.sumstats

## Download and prepare the GWAS and GTEx whole blood data
mkdir WEIGHTS
cd WEIGHTS
wget https://data.broadinstitute.org/alkesgroup/FUSION/WGT/GTEx.Whole_Blood.tar.bz2
tar xjf GTEx.Whole_Blood.tar.bz2

## Performing the expression imputation
cd ..
Rscript FUSION.assoc_test.R \
--sumstats PGC2.SCZ.sumstats \
--weights ./WEIGHTS/GTEx.Whole_Blood.pos \
--weights_dir ./WEIGHTS/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 22 \
--out PGC2.SCZ.22.dat

## Joint/conditional tests and plots
cat PGC2.SCZ.22.dat | awk 'NR == 1 || $NF < 0.05/2058' > PGC2.SCZ.22.top

Rscript FUSION.post_process.R \
--sumstats PGC2.SCZ.sumstats \
--input PGC2.SCZ.22.top \
--out PGC2.SCZ.22.top.analysis \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 22 \
--plot --locus_win 100000