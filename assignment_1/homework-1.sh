# Set working directory
cd "D:\OneDrive - Johns Hopkins\Course\140.686.01 - Advanced Methods for Statistical Genetics and Genomics\assignment\homework_1"

# Python 2 environment
python ./ldsc/ldsc.py -h   # test ldsc

# Check Summary Statistics
head ./data/data_BC17.txt    # brest cancer
head ./data/data_BMI18.txt   # body mass index
head ./data/data_SCZ.txt     # schizophrenia

# Reformatting Summary Statistics
# 1. A unique identifier -> rsid
# 2. Allele 1 (effect allele) -> A1
# 3. Allele 2 (non-effect allele) -> A2
# 4. Sample size -> N
# 5. A P-value -> pval
# 6. A signed summary statistic -> z
python ./ldsc/munge_sumstats.py -h   # check munge_sumstats.py options
python ./ldsc/munge_sumstats.py --sumstats ./data/data_BC17.txt --out ./clean/bc       # brest cancer
python ./ldsc/munge_sumstats.py --sumstats ./data/data_BMI18.txt --ignore beta --out ./clean/bmi     # body mass index
python ./ldsc/munge_sumstats.py --sumstats ./data/data_SCZ.txt --ignore beta --out ./clean/scz     # schizophrenia

# Estimating Heritability
python ./ldsc/ldsc.py --h2 ./clean/bc.sumstats.gz --ref-ld-chr ./data/eur_w_ld_chr/ --w-ld-chr ./data/eur_w_ld_chr/ --out ./clean/bc_h2
python ./ldsc/ldsc.py --h2 ./clean/bmi.sumstats.gz --ref-ld-chr ./data/eur_w_ld_chr/ --w-ld-chr ./data/eur_w_ld_chr/ --out ./clean/bmi_h2
python ./ldsc/ldsc.py --h2 ./clean/scz.sumstats.gz --ref-ld-chr ./data/eur_w_ld_chr/ --w-ld-chr ./data/eur_w_ld_chr/ --out ./clean/scz_h2

# Conversion to Liability Scale
## Breast cancer sample prevalence: 106571/(106571 + 95762) = 0.5267109
python ./ldsc/ldsc.py --h2 ./clean/bc.sumstats.gz --ref-ld-chr ./data/eur_w_ld_chr/ --w-ld-chr ./data/eur_w_ld_chr/ --out ./clean/bc_h2_lia_eur_fem --samp-prev 0.5267109 --pop-prev 0.0008     # Europe, female
## Schizophrenia sample prevalence: 35476/(35476 + 46839) = 0.4309786
python ./ldsc/ldsc.py --h2 ./clean/scz.sumstats.gz --ref-ld-chr ./data/eur_w_ld_chr/ --w-ld-chr ./data/eur_w_ld_chr/ --out ./clean/scz_h2_lia --samp-prev 0.4309786 --pop-prev 0.0034     # Europe