#heatmap to show variants cluster by species

#load libraries
library(pheatmap)
library(dendextend)

#start with data frame containing cleaned up variant table with combined replicates
df <- as.data.frame(variant_table_cleaned_combined)

#NAs to 0s for clustering purposes
df[is.na(df)] <- 0

#get just frequency data into a matrix, leave out inoculums (the ferret would cluster with P3 for example, and other things could be thrown off)
mat <- as.matrix(df[,5:17])
row.names(mat) <- df$variant

#scale data to a distribution with mean 0 and SD as 1
mat_scale <- scale(mat)

#default heatmap
pdf(file="animal_heatmap.pdf")
pheatmap(mat, fontsize_row = 5)
dev.off()

#heatmap where you only cluster based on animals (columns) is the same 
pheatmap(mat, cluster_rows = FALSE, fontsize_row = 5)

#default heatmap with scaling
pdf(file="animal_heatmap_scaled.pdf")
pheatmap(mat_scale, fontsize_row = 5)
dev.off()

#default clustering
t.mat <- t(mat)

clust <- hclust(dist(t.mat), method="complete")

pdf(file="animal_dendrogram.pdf")
as.dendrogram(clust) %>%
  plot(horiz=TRUE)
dev.off()


