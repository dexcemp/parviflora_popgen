reference=H230-3.ref.fasta
prefix=parviflora

# Make gvcf list (must end in .list or the next step will throw an error)
ls *-g.vcf > vcf.list

# Combine and Jointly call GVCFs
gatk CombineGVCFs -R $reference --variant vcf.list --output "$prefix".cohort.g.vcf
gatk GenotypeGVCFs -R $reference -V "$prefix".cohort.g.vcf -O "$prefix".cohort.unfiltered.vcf

# Keep only SNPs passing a hard filter
#outputs vcf containing only SNPS (removes indels)
time gatk SelectVariants -V "$prefix".cohort.unfiltered.vcf -R $reference -select-type-to-include SNP -O "$prefix".SNPall.vcf
#hardfiltering: flags variants as PASS or hardfilter (does not remove variants, only flags)
time gatk VariantFiltration -R $reference -V "$prefix".SNPall.vcf --filter-name "hardfilter" -O "$prefix".snp.filtered.vcf --filter-expression "QD < 5.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 & ReadPosRankSum < -8.0"
#keep SNPs that only passed the hardfilter
gatk SelectVariants -V "$prefix".snp.filtered.vcf --exclude-filtered -O "$prefix".snp.filtered.pass.vcf 