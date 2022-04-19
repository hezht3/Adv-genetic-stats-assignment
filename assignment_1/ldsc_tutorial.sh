cd "C:\Users\Johnathan He\ldsc"
cd tutorial

# 1. Compute the genetic correlation between schizophrenia and bipolar disorder

# Download Data
# Need to install `wget`
wget https://4119f3fb-a-f9436c1e-s-sites.googlegroups.com/a/broadinstitute.org/pgc-summer-school-2015/lecture-materials/pgc.cross.bip.zip
wget https://4119f3fb-a-f9436c1e-s-sites.googlegroups.com/a/broadinstitute.org/pgc-summer-school-2015/lecture-materials/pgc.cross.scz.zip
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2

# Munge Data
# Need to run using Git Bash
tar -jxvf eur_w_ld_chr.tar.bz2
unzip -o pgc.cross.bip.zip
unzip -o pgc.cross.scz.zip
bunzip2 w_hm3.snplist.bz2

# In python2 environment
python munge_sumstats.py --sumstats "tutorial/pgc.cross.SCZ17.2013-05.txt" --N 17115 --out scz --merge-alleles "tutorial/w_hm3.snplist"
python munge_sumstats.py --sumstats "tutorial/pgc.cross.BIP11.2013-05.txt" --N 11810 --out bip --merge-alleles "tutorial/w_hm3.snplist"
# LD Score Regression
python ldsc.py --rg scz.sumstats.gz,bip.sumstats.gz --ref-ld-chr "tutorial/eur_w_ld_chr/" --w-ld-chr "tutorial/eur_w_ld_chr/" --out scz_bip
less scz_bip.log