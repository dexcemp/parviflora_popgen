#PCA
#Read eigenvec file
pca <- read.table("parviflora_ss_pruned.eigenvec", header = FALSE)

#Give columns names
colnames(pca) <- c("FID","IID", paste0("PC",1:(ncol(pca)-2)))

#Add population data
meta <- read.table("parviflora_ss_popdata.txt", header=TRUE)
pca2 <- merge(pca, meta, by="IID")

#Add variance explained in axis labels
eig <- scan("parviflora_ss_pruned.eigenval")
var_exp <- eig/sum(eig)*100

#Plot in ggplot2
library(ggplot2)
ggplot(pca2, aes(PC1, PC2, fill=Pop)) +
  geom_point(shape = 21, size=3, color="black") +
  scale_fill_manual(values = c(
    "puberula" = "#00946a",
    "missouriensisN" = "#facc14",
    "missouriensisS" = "#a98600",
    "parvifloraN" = "#305cde",
    "parvifloraS" = "#6fd0ee",
    "saurensis" = "#ec5800")) +
  xlab(paste0("PC1 (", round(var_exp[1],2),"%)")) +
  ylab(paste0("PC2 (", round(var_exp[2],2),"%)")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))


#FST Heat map
#Load the PLINK .summary file
fst <- read.table("parviflora_ss_fst.fst.summary", header=TRUE, stringsAsFactors = FALSE)

#convert to an FST matrix
library(reshape2)
fst_matrix <- acast(fst, POP1 ~ POP2, value.var="WC_FST")

#fill the symmetric part
fst_matrix[lower.tri(fst_matrix)] <- t(fst_matrix)[lower.tri(fst_matrix)]
diag(fst_matrix) <- 0

write.csv(fst_matrix, "parviflora_ss_newfst.csv") 

#plot heatmap
library(pheatmap)
pheatmap(fst_matrix,
         display_numbers = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white", "yellow", "red"))(256))