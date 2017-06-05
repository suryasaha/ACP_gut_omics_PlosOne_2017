#using all replicates Min1-4, Pos1-4

#clean
rm(list = ls())

source("https://bioconductor.org/biocLite.R")
library(DESeq2)


counts <- read.delim("PE_SE_data_matrix.txt",header = TRUE, row.names = 1)
meta_data <- read.delim("PE_SE_meta_data_matrix.txt",header = TRUE, row.names = 1)

# round off decimals
counts  <- data.frame(row.names = rownames(counts), mapply(function(id) round(id), counts))

head(counts)

# DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_data,
  design = ~ Condition)
dds$Condition
# [1] pos pos pos pos pos pos pos pos min min min min min min min min pos pos pos pos pos pos pos pos min min
# [27] min min min min min min
# Levels: min pos

# Collapsing all technical replicates
as.data.frame( colData( dds )[ ,c("Sample","Run","Condition","Lib") ] )
ddsCollapsed <- collapseReplicates(dds, groupby = dds$Sample)
as.data.frame( colData( ddsCollapsed )[ ,c("Sample","Run", "Condition","Lib") ] )
as.data.frame( colData( ddsCollapsed )[ ,c("Sample", "Condition") ] )

# write to disk for edgeR processing
write.csv(counts(ddsCollapsed), file = "All_ddsCollapsed.csv")

# collapseRep puts the Min samples first
head(counts(ddsCollapsed))
# min1 min2 min3 min4 pos1 pos2 pos3 pos4
# gi|662182959|ref|XM_008472436.1|  219  240  224  245  154  230  295  311
# gi|662182961|ref|XM_008482233.1|   23   30   30   28   20   33   30   44

ddsCollapsed$Condition
# [1] min min min min pos pos pos pos
# Levels: min pos

library("edgeR")
# removing low counts, is 3 but should be 4 as we have 4 replicates
dim(ddsCollapsed[rowSums(cpm(counts(ddsCollapsed)) >1) >=3,])
# [1] 15385     8
# 15385/21986 genes have > 1 read/million in 3 or more libs
ddsCollapsed <- ddsCollapsed[rowSums(cpm(counts(ddsCollapsed)) >1) >=3,]
# write to disk for edgeR processing
write.csv(counts(ddsCollapsed), file = "Allnolow_ddsCollapsed.csv")

########################################################################################
#### Pos(experimental) vs Min (control). Min is reference ####
########################################################################################

ddsCollapsedPvM <- ddsCollapsed

#  if you never tell the DESeq2 functions  which  level  you  want  to  compare  against  (e.g.  which  level  represents  the  control  group), the comparisons will be based on the alphabetical order of the levels. Min will be ref here but still  specifying the reference level explicitly
ddsCollapsedPvM$Condition <- relevel(ddsCollapsedPvM$Condition, ref="min")

library(DESeq2)

boxplot(counts(ddsCollapsedPvM), las = 2, main = "Spread of counts per library after collapsing PE and SE")

library("corrplot")
corrplot(cor(counts(ddsCollapsedPvM)), method="square", tl.col="black", addgrid.col="black", is.corr=FALSE, main="Collapsed counts")

#DESeq
ddsCollapsedPvM <- DESeq(ddsCollapsedPvM)

# results with FDR alpha = 0.1
ddsCollapsedPvMResults <- results(ddsCollapsedPvM)
metadata(ddsCollapsedPvMResults)$alpha
# [1] 0.1

# results with FDR alpha = 0.01 or 1%
ddsCollapsedPvMResultsFDR0.01 <- results(ddsCollapsedPvM, alpha = 0.01)
metadata(ddsCollapsedPvMResultsFDR0.01)$alpha
# [1] 0.01

summary(ddsCollapsedPvMResultsFDR0.01)
# out of 15385 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)     : 672, 4.4% 
# LFC < 0 (down)   : 574, 3.7% 
# outliers [1]     : 10, 0.065% 
# low counts [2]   : 597, 3.9% 
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# sort by adjusted pvalue
ddsCollapsedPvMResultsFDR0.01SortedPadj <- ddsCollapsedPvMResultsFDR0.01[order(ddsCollapsedPvMResultsFDR0.01$padj),]

# padj < 0.1
sum(ddsCollapsedPvMResultsFDR0.01SortedPadj$padj < 0.1, na.rm = TRUE)
# 2191
# 1411
# [1] 771

# padj < 0.05
sum(ddsCollapsedPvMResultsFDR0.01SortedPadj$padj < 0.05, na.rm = TRUE)
# 1832
# 1135
# [1] 600


# Not screening for log2FoldChange as down regulated genes have -ive values. padj < 0.05
dim(subset(ddsCollapsedPvMResultsFDR0.01SortedPadj, ddsCollapsedPvMResultsFDR0.01SortedPadj$padj < 0.05))
# [1] 1832    6
# [1] 1135    6

# Output
write.csv(ddsCollapsedPvMResultsFDR0.01SortedPadj,quote = FALSE, file = "ddsCollapsedPvMResultsFDR1percSortedPadj_nolow-POSvsMIN.csv")

write.csv(subset(ddsCollapsedPvMResultsFDR0.01SortedPadj, ddsCollapsedPvMResultsFDR0.01SortedPadj$log2FoldChange > 0.5 & ddsCollapsedPvMResultsFDR0.01SortedPadj$padj < 0.1),quote = FALSE, file = "DESeq2_FDR1perc_LogFC0.5_nolow-POSvsMIN.csv")

write.csv(subset(ddsCollapsedPvMResultsFDR0.01SortedPadj, ddsCollapsedPvMResultsFDR0.01SortedPadj$padj < 0.05),quote = FALSE, file = "DESeq2_FDR1perc_nolow-POSvsMIN.csv")

###################
##### VIZ #########
###################

# shows the log2 fold changes attributable to a given variable over the mean of normalized counts. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
# The shrinkage is greater for the log2 fold change estimates from genes with low counts and high dispersion, as can be seen by the narrowing of spread of leftmost points
# need to use prefix to call correct method
DESeq2::plotMA(ddsCollapsedPvMResultsFDR0.01SortedPadj, main="DESeq2 All nolow collapsed", ylim=c(-2,2))
ddsCollapsedPvMResultsFDR0.01SortedPadj

# dispersion estimate plot shows the gene-wise estimates (black), the fitted values (red), and the final maximum a posteriori estimates used in testing (blue). blue circles are genes with high dispersion that are outliers
plotDispEsts(ddsCollapsedPvM, ylim = c(1e-6, 1e1), main="DESeq2 All nolow")

# gene with lowest $padj 
# counts of reads for a gene across the groups. plotCounts normalizes counts by sequencing depth and adds a pseudocountof 1/2 to allow for log scale plotting, looks great
plotCounts(ddsCollapsedPvM, gene=which.min(ddsCollapsedPvMResultsFDR0.01$padj), intgroup="Condition", main="DESeq2 trimmed\n gene with lowest adj pvalue")
# gene with max $padj, looks ok
plotCounts(ddsCollapsedPvM, gene=which.max(ddsCollapsedPvMResultsFDR0.01$padj), intgroup="Condition", main="DESeq2 trimmed\n gene with highest adj pvalue")

# PCA
rldCollapsed <- rlog(ddsCollapsedPvM)
plotPCA(rldCollapsed, intgroup = c("Condition"))
plotPCA(rldCollapsed, intgroup = c("Sample"))
# This looks great. Now the PC2 is only 8% and all Min samples group together. Pos1 is still different but all Pos samples are spread out.

# Heatmap showing the Euclidean distances between the samples as calculated from the regularized log transformation.
sampleDists <- dist(t(assay(rldCollapsed)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rldCollapsed$Condition, rldCollapsed$Sample, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library("pheatmap")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

