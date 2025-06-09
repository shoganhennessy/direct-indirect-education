#!/bin/bash
## Senan Hogan-Hennessy
## 28 February 2025.
# This is an adjusted version of the PRS pipeline, provided by
# Collister Liu Clifton (2022, https://doi.org/10.3389/fgene.2022.818574) 
# https://github.com/2cjenn/PRS_Pipeline/blob/master/demo.sh

# Log the output of this file
# dx download data-input/prs-pipeline.sh
# bash -x prs-pipeline.sh &> prs-pipeline.log 
# bash prs-pipeline.sh 2>&1 | tee prs-pipeline.log

# Input: GWAS summary stats, from SSGAC (made into a csv).
#gwasResultsCsv=$1
#gwasResultsCsv="EA4_additive_p1e-5_clumped.csv"
#gwasResultsCsv="tan_2024_educational_attainment.csv"
gwasResultsCsv="EA4_additive_excl_23andMe.csv"
#gwasResultsCsv="ADHD1_multi_p5e-8_sumstats.csv"

dx download data-input/$gwasResultsCsv
head $gwasResultsCsv

# Name of the output file.
#gwasResultsOutput=$2
#gwasResultsOutput="raw-ed-pgi-okbay-2022"
gwasResultsOutput="raw-ed-pgi-okbay-exclude-2022"
#gwasResultsOutput="raw-ed-pgi-tan-2024"
#gwasResultsOutput="raw-adhd-pgi"


# Adjust GWAS files to format needed.
awk -F, '{ if (NR>1) { print $1 }}' $gwasResultsCsv > rsidlist.txt
awk -F, '{ if (NR>1) { print sprintf("%02d", $2)":"$3"-"$3 }}' $gwasResultsCsv > chrposlist.txt

## Get the BGEN software working on a UKB RAP terminal.
# Install Bgenix.
wget http://code.enkre.net/bgen/tarball/release/bgen.tgz
tar -xvzf bgen.tgz
cd bgen.tgz
# compile it (1) adapt python -> python3, (2) fix a bug in the Bgen C+ code.
apt-get install python-is-python3
sed -i -e "s/ios::streampos/streampos/g" src/View.cpp
./waf configure
./waf build
# test the Bgen compilation.
./build/test/unit/test_bgen
./build/apps/bgenix -g example/example.16bits.bgen -list
cd ..

# Extract the SNP-relevant regions from the 22 autosome data files with bgenix.
cmd=""
for i in {1..22}
do
    echo $i
    # Extract the SNP relevant sections from Bulk imputation BGEN
    bgen.tgz/build/apps/bgenix -g \
        "/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${i}_b0_v3.bgen" \
        -incl-rsids rsidlist.txt \
        -incl-range chrposlist.txt > chr_${i}.bgen &
    # Concatenate names of completed files.
    cmd=$cmd"chr_${i}.bgen "
done

# Combine the .bgen files for each chromosome into one BGEN file
bgen.tgz/build/apps/cat-bgen -g $cmd -og initial_chr.bgen
# Write index file .bgen.bgi
bgen.tgz/build/apps/bgenix -g initial_chr.bgen -index -clobber


# Import the GWAS estimates into the sqlite database as a table called Betas.
apt-get install sqlite3
sqlite3 initial_chr.bgen.bgi "DROP TABLE IF EXISTS Betas;"
sqlite3 -separator "," initial_chr.bgen.bgi ".import $gwasResultsCsv Betas"

sqlite3 initial_chr.bgen.bgi "DROP TABLE IF EXISTS Joined;"
# And inner join it to the index table (Variants), making a new table (Joined)
# By joining on alleles as well as chromosome and position 
# we can ensure only the relevant alleles from any multi-allelic SNPs are retained
sqlite3 -header -csv initial_chr.bgen.bgi \
    "CREATE TABLE Joined AS 
        SELECT Variant.*, Betas.chr_name, Betas.Beta FROM Variant INNER JOIN Betas 
            ON Variant.chromosome = printf('%02d', Betas.chr_name) 
            AND Variant.position = Betas.chr_position 
            AND Variant.allele1 = Betas.noneffect_allele 
            AND Variant.allele2 = Betas.effect_allele 
        UNION 
        SELECT Variant.*, Betas.chr_name, -Betas.Beta FROM Variant INNER JOIN Betas 
            ON Variant.chromosome = printf('%02d', Betas.chr_name) 
            AND Variant.position = Betas.chr_position 
            AND Variant.allele1 = Betas.effect_allele AND 
            Variant.allele2 = Betas.noneffect_allele;"

# Filter the .bgen file to include only alleles specified in Betas for each SNP.
bgen.tgz/build/apps/bgenix -g initial_chr.bgen -table Joined > single_allelic.bgen
# Produce an index file for the new .bgen
bgen.tgz/build/apps/bgenix -g single_allelic.bgen -index


########################################################
#----------------------- SNP QC -----------------------#
########################################################

# Get plink2, and save available ram (in Mb)
apt-get install plink2
freememory=$(free --mega | awk '/^Mem:/{print $7}')

# Convert to plink2 format.
plink2 --bgen single_allelic.bgen ref-first \
    --hard-call-threshold 0.1 \
    --sample "/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c1_b0_v3.sample" \
    --set-all-var-ids @:#_\$r_\$a \
    --freq \
    --make-pgen \
    --memory $freememory \
    --out raw

# Identify ambiguous SNPs
awk '/^[^#]/ { if( $5>0.4 && $5<0.6 && ( ($3=="A" && $4=="T") || ($4=="T" && $3=="A") || ($3=="C" && $4=="G") || ($4=="G" && $3=="C") ) ) { print $0 }}' \
    raw.afreq > exclrsIDs_ambiguous.txt

# Concatenate the MAF+Info files into one tsv file.
for i in {1..22}
do
    awk -v chr=$i 'BEGIN {FS="\t"; OFS="\t"} { print chr,$0,chr":"$3"_"$4"_"$5 }' "/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${i}_b0_v3.mfi.txt"
done > ukb_mfi_all_v3.tsv

# Exclude the ambiguous SNPs.
plink2 --pfile raw \
    --exclude exclrsIDs_ambiguous.txt \
    --extract-col-cond ukb_mfi_all_v3.tsv 9 10 --extract-col-cond-min 0.4 \
    --maf 0.005 \
    --write-snplist \
    --make-pgen \
    --out snpQC

# Write the samples.
plink2 --pfile raw \
    --extract snpQC.snplist \
    --write-samples \
    --out sampleQC


########################################################
#----------------- Calculate the PRS ------------------#
########################################################

# Put a file of the effect scores.
sqlite3 -separator " " -list initial_chr.bgen.bgi \
    "SELECT chr_name || ':' || position || '_' ||allele1 || '_' || allele2, allele2, Beta FROM Joined;" \
    > score.txt

# Calculate the weighted sum for each individual.
plink2 --pfile raw                              \
    --extract  snpQC.snplist                    \
    --keep     sampleQC.id                      \
    --score    score.txt     no-mean-imputation \
    --out      ed-pgi-score

# Save the calculated PRS.
dx upload ed-pgi-score.sscore --path data-clean/$gwasResultsOutput.tsv
