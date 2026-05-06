#!/bin/bash
set -eo pipefail

#Set an ID for each SNP so they can be filtered by name
bcftools annotate --set-id +"%CHROM:%POS" "$prefix".snp.filtered.pass.vcf > "$prefix".snp.filtered.nodot.vcf

#Subset to a specific species (keep only samples of specific species/taxa you want to do downstream analyses with)
vcftools --vcf "$prefix".snp.filtered.nodot.vcf --keep "$prefix"keep.txt  --recode --recode-INFO-all --out "$prefix".vcf

#Include only sites with Minor Allele Frequency greater than or equal to 0.05, with 20% missing data or less
vcftools --vcf "$prefix".vcf --maf 0.05 --max-missing 0.8 --recode --recode-INFO-all --out "$prefix"filtered_maf5_missing80.vcf

#Filter out linked genes to minimize linkage disequilibrium
plink --vcf "$prefix"filtered_maf5_missing80.vcf --indep-pairwise 50 5 0.2 --allow-extra-chr --recode --out "$prefix"
plink --extract "$prefix".prune.in --out "$prefix"_pruned --file "$prefix" --make-bed --allow-extra-chr --recode

#Generate eigenvalues and loadings for 20 PCA axes
plink --pca 20 --file "$prefix"_pruned --allow-extra-chr --out "$prefix"_pruned

#Generate basic stats (heterozygosity, inbreeding coefficient, allele frequencies)
plink --freq --het 'small-sample' --ibc --file "$prefix"_pruned --allow-extra-chr -out "$prefix"_pruned

#Calculate Weir and Cockerham FST 
plink2 --bfile "$prefix"_ss_pruned  --pheno "$prefix"_ss_fstpops.txt  --fst PHENO1 method=wc  --out "$prefix"_ss_fst --allow-extra-chr

#File conversions 
#If you have plink files (.bed/.bim/.fam) and you want to convert it to vcf:
plink --bfile parviflora_ss_pruned --recode vcf --allow-extra-chr --out parviflora_ss_pruned

#Then convert to a file for R
vcftools --vcf parviflora_ss_pruned.vcf --012 --out parviflora_ss_snp

#Cleaning
cut -f2- parviflora_ss_snp.012 | sed 's/-1/NA/g' >snp.temp
tr -d '\t' <parviflora_ss_snp.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
paste <(echo "ID" | cat - parviflora_ss_snp.012.indv) <(echo "" | cat header - snp.temp) > parviflora_ss_snp.forR
rm header snp.temp