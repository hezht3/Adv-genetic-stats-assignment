require(coloc)

# ---------------------------- format datasets ---------------------------- #

# A minimum coloc dataset

data(coloc_test_data)
attach(coloc_test_data)
minimum_data = D1[c("beta","varbeta","snp","position","type","sdY")]
str(minimum_data)

## check dataset

check_dataset(minimum_data, warn.minp = 1e-10)

plot_dataset(minimum_data)

## for case control trait

minimum_ccdata = D1[c("beta","varbeta","snp","position")]
minimum_ccdata$type = "cc"
str(minimum_ccdata)

check_dataset(minimum_ccdata)

## no sdY

nosdY_data = D1[c("beta","varbeta","snp","position","type","N","MAF")]
str(nosdY_data)

check_dataset(nosdY_data)

## no beta or var

nobeta_data = D1[c("MAF","snp","position","type","sdY","N")]
nobeta_data$pvalues = pnorm(-abs(D1$beta/sqrt(D1$varbeta)))*2
str(nobeta_data)

check_dataset(nobeta_data)

nobeta_ccdata = D1[c("MAF","snp","position","N")]
nobeta_ccdata$pvalues = pnorm(-abs(D1$beta/sqrt(D1$varbeta)))*2
nobeta_ccdata$type = "cc"
nobeta_ccdata$s = 0.5
str(nobeta_ccdata)

check_dataset(nobeta_ccdata)

# LD for coloc.susie

str(D1$LD)

# NB allele order

## check_alignment()

# ---------------------------- single variant assumption ---------------------------- #

## Fine mapping under a single causal variant assumption

plot_dataset(D1)

my.res <- finemap.abf(dataset = D1)
my.res[21:30,]

tail(my.res,3)

## (Approximate) Bayes Factor colocalisation analyses

my.res <- coloc.abf(dataset1 = D1,
                    dataset2 = D2)
print(my.res)

subset(my.res$results, SNP.PP.H4 > 0.01)   # extract more likely causal variant

o <- order(my.res$results$SNP.PP.H4, decreasing = TRUE)
cs <- cumsum(my.res$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
my.res$results[o,][1:w, ]$snp   # extract 95% credible set

# ---------------------------- sensitivity to prior values ---------------------------- #

## coloc explorer app: https://chr1swallace.shinyapps.io/coloc-priors/

# Sensitivity analysis

my.res <- coloc.abf(dataset1 = D1,
                    dataset2 = D2,
                    p12 = 1e-6)
my.res

sensitivity(my.res, rule = "H4 > 0.5")

sensitivity(my.res, rule = "H4 > 3*H3 & H0 < 0.1")

# ---------------------------- using SuSiE to relax the single causal variant assumption ---------------------------- #

# Multiple causal variants, using SuSiE to separate the signals

par(mfrow = c(2, 1))
plot_dataset(D3, main = "Dataset D3")
plot_dataset(D4, main = "Dataset D4")

## Standard coloc

my.res <- coloc.abf(dataset1 = D3, dataset2 = D4)
class(my.res)
my.res

sensitivity(my.res, "H4 > 0.9")

## Check LD matrix of data set

check_dataset(D3, req = "LD")
check_dataset(D4, req = "LD")

S3 = runsusie(D3)
summary(S3)
S4 = runsusie(D4)
summary(S4)

if(requireNamespace("susieR", quietly = TRUE)) {
    susie.res = coloc.susie(S3, S4)
    print(susie.res$summary)
}

### Sensitivity

if(requireNamespace("susieR", quietly = TRUE)) {
    sensitivity(susie.res, "H4 > 0.9", row = 1, dataset1 = D3, dataset2 = D4)
    sensitivity(susie.res, "H4 > 0.9", row = 2, dataset1 = D3, dataset2 = D4)
}
