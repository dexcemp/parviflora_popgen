# parviflora_popgen
Pipeline for genetic structure and environmental association analyses (RDA) for the Heuchera parviflora group.

* 1_gatk_variantcall.sh: Script to perform variant calling in GATK4.
* 2_gatk_combinefilter.sh: Script to remove indels and perform hard filtering in GATK4.
* 3_plink_filtering&stats.sh: Script to subset SNP dataset and filter out sites with minor allele frequency < 0.05 and missing data > 20% using vcftools; filter out linked genes, generate loadings for PCA, generate basic popgen stats, and calculate FST using PLINK.
* 4_fastSTRUCTURE.sh: Script to run fastSTRUCTURE, finding maximum likelihood K, and visualization with DISTRUCT.
* 5_genetic_structure_visualization.R: Script to plot PCA and FST heatmaps in R.
* 6_RDA.R: Script to extract climate data, perform correlation and redundancy analyses, and generate RDA plots.
* 7_partialRDA_variancepartitioning.R: Script to extract dbMEMs and run partial redundancy analyses.
