#heatmap to show variants do not cluster by species

#load libraries
library(pheatmap)
library(dendextend)
library(openxlsx)

#start with data frame containing cleaned up variant table with combined replicates
df <- read.csv("variant_summary_processed.csv")

#NAs to 0s for clustering purposes
df[is.na(df)] <- 0

#get just frequency data into a matrix, leave out inoculums (the ferret would cluster with P3 for example, and other things could be thrown off)
#mat <- as.matrix(df[,11:41])
mat <- as.matrix(df[,9:24])
row.names(mat) <- df$variant

#default heatmap with all the variants on there (and all the replicates if you use that data)

pdf(file="heatmap.pdf")
pheatmap(mat, fontsize_row = 5)
dev.off()

#scale data to a distribution with mean 0 and SD as 1
mat_scale <- scale(mat)

#default heatmap with scaling (looks pretty much the same)
pdf(file="heatmap_scaled.pdf")
pheatmap(mat_scale, fontsize_row = 5)
dev.off()

#default clustering, also doesn't cluster by species, but arrangement does seem a little different
t.mat <- t(mat)
clust <- hclust(dist(t.mat), method="complete")

pdf(file="clustered_dendrogram.pdf")
as.dendrogram(clust) %>%
  plot(horiz=TRUE)
dev.off()


