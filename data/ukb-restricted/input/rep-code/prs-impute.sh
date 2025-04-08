#!/bin/bash
## Senan Hogan-Hennessy
## 27 March 2025.
# Impute the parents' (mean) genotype, and implied E[ Ed PGI | self, siblings].
# Young et al. (2022) (https://www.nature.com/articles/s41588-022-01085-0).
# https://snipar.readthedocs.io/en/latest/tutorial.html

# Note: this setup relies on files constructed in data-extract.py,
#       so that python3 code must be run first, before this works.
#       Also Snipar only works on python3.9, so must run in a 3.9 venv.


################################################################################
## IMport necessary files.

apt-get install plink2
# Snipar requires kinship + agesex file, computed in data-extract.py
#                     + chr_i BGEN files computed in prs-pipeline.sh.

# Get the pre-computed BGEN files, limited to Ed PGI relevant SNPs.
for i in {1..22}
do
    echo $i
    # Download BGEN
    dx download data-input/chr_$i.bgen
    # Get sample files.
    dx download "Bulk/Imputation/UKB imputation from genotype/ukb22828_c${i}_b0_v3.sample"
    mv "ukb22828_c${i}_b0_v3.sample" "chr_${i}.sample"
    # Convert BGEN -> BED (using BGEN filetype has a bug in Snipar package).
    plink2                            \
        --bgen chr_$i.bgen  ref-first \
        --sample chr_$i.sample        \
        --make-bed                    \
        --out chr_$i
done

# Get the Ed PGI GWAS weights
dx download data-input/EA4_additive_p1e-5_clumped.csv
# Get the computed Kinship file.
dx download data-input/kinship-adjusted.dat
# Get the computed age-sex file.
dx download data-input/phenotype-agesex.dat


################################################################################
## Install Snipar package

# Dependencies -> older version of python for this package.
versionNeeded="3.9.16"
apt-get install libssl-dev openssl python3-virtualenv libffi-dev cargo
wget https://www.python.org/ftp/python/$versionNeeded/Python-$versionNeeded.tgz
tar xzvf Python-$versionNeeded.tgz
cd Python-$versionNeeded
./configure --enable-optimizations
make
make install
cd ..
# Set the environment to use the older python version + rust manager, update pip
virtualenv --python="python$versionNeeded" "env-$versionNeeded"
source env-$versionNeeded/bin/activate
# Install the Snipar package, in the 3.9 venv.
git clone "https://github.com/AlexTISYoung/snipar"
pip install ./snipar
# Test the package has installed correctly.
python3 -m unittest snipar.tests


################################################################################
## Impute the missing parental genotypes, by Snipar package.

# See Snipar guide for the workflow.
# https://snipar.readthedocs.io/en/latest/guide.html#workflow

# Step 1. Infer IBD between siblings.
ibd.py                            \
    --bed chr_@                   \
    --chr_range 1-22              \
    --king kinship-adjusted.dat   \
    --agesex phenotype-agesex.dat \
    --threads 4

# Step 2. Impute parental genotypes.
impute.py                             \
    --bed       chr_@                 \
    --chr_range 1-22                  \
    --ibd       ibd_chr_@.ibd         \
    --king      kinship-adjusted.dat  \
    --agesex    phenotype-agesex.dat  \
    --out       parent-imputed-geno-@ \
    --threads   4

#! Command to list the bug for.
impute.py --bgen chr_@ --chr_range 1-22 --ibd ibd_chr_@.ibd --king kinship-adjusted.dat --agesex phenotype-agesex.dat --out parent-imputed-geno-@
#! Note:  32GB of RAM, to avoid getting memory killed.

https://github.com/AlexTISYoung/snipar/issues/new
Bug in pre_prepocessdata.py for .bgen input files.

Hi, hope you are well today.

I am listing a bug, related to line 600 in XXXYY.py.

When I run the imput.py script, I am greeted with the following error.

CODE HERE.

In addition, if I use plink2 to convert the 22 chromosome bgen files to bim, and uses those files as input for impute.py, then everything runs smoothly.

# Convert BGEN -> BED, then run impute.py --bed .
for i in {1..22}
do
    echo $i
    plink2                            \
        --bgen chr_$i.bgen  ref-first \
        --sample chr_$i.sample        \
        --make-bed                    \
        --out chr_$i
done
impute.py --bed chr_@ --chr_range 1-22 --ibd ibd_chr_@.ibd --king kinship-adjusted.dat --agesex phenotype-agesex.dat --out parent-imputed-geno-@


# Step 3. Calculate the PGI, with parental values too.
pgs.py imputed-ed-pgi                          \
    --bed       chr_@                          \
    --imp       parent-imputed-geno-@          \
    --chr_range 1-22                           \
    --weights   EA4_additive_p1e-5_clumped.csv \
    --SNP       "rsid"                         \
    --beta_col  "Beta"                         \
    --A1        "effect_allele"                \
    --A2        "noneffect_allele"             \
    --sep       ","                            \
    --fit_sib                                  \
    --parsum                                   \
    --threads   4

# Save the parental imputed Ed PGI file.
dx upload imputed-ed-pgi.pgs.txt --path data-clean/imputed-ed-pgi.pgs.txt
