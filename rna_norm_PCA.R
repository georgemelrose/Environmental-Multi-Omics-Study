install.packages("ggplot2")
install.packages("stringr")
install.packages("dplyr")
library(ggplot2)
library(stringr)
library(dplyr)

#Removing the gene names, to just look at gene names#
rna_norm_numeric <- select(rna_norm_counts, -X)


#PCA on normalised RNA-seq data#
normrna_pca <- prcomp(rna_norm_numeric, center = TRUE, scale = TRUE)
pca_summary <- summary(normrna_pca)
print(pca_summary)

pca_df <- data.frame (normrna_pca$x, make = stringr::word(rownames(rna_norm_numeric), 1))
print(normrna_pca$x)

ggplot(pca_df, aes(x = PC1, y = PC2, col = make)) +
  geom_point(size = 3) +
  labs(x = "PC1 81.93%", y= "PC2 12.37%", title = "PCA for normalised RNA data" ) + theme(legend.position = "bottom")

