#clean
rm(list = ls())

source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library("edgeR")

# Reading in collapsed counts from DESeq2 so that PE and SE data is merged by sample!!
counts <- read.delim("Allnolow_ddsCollapsed.csv", sep = ",", header = TRUE, row.names = 1)
colSums(counts)
# min1    min2    min3    min4    pos1    pos2    pos3    pos4 
# 2993728 3612630 3296422 3630525 2665335 3063938 3556758 3945265 
colSums(counts)/1e06
# min1     min2     min3     min4     pos1     pos2     pos3     pos4 
# 2.993728 3.612630 3.296422 3.630525 2.665335 3.063938 3.556758 3.945265 

################################
##### Removing low counts ######
################################

dim(counts[rowSums(cpm(counts) > 1) >= 3,])
# [1] 15385     6
# 15385/21986 genes have > 1 read/million in 3 or more libs
counts_nolow <- counts[rowSums(cpm(counts) > 1) >= 3,]
colSums(counts_nolow)
# min1    min2    min3    min4    pos1    pos2    pos3    pos4 
# 2988747 3606323 3291364 3623719 2660719 3058679 3550961 3939584 
colSums(counts_nolow)/1e06
# min1     min2     min3     min4     pos1     pos2     pos3     pos4 
# 2.988747 3.606323 3.291364 3.623719 2.660719 3.058679 3.550961 3.939584 

library("corrplot")
corrplot(cor(counts_nolow), method="square", tl.col="black", addgrid.col="black", is.corr=FALSE, main="Collapsed counts (no low counts)")


groups <- c("min", "min", "min", "min", "pos", "pos", "pos", "pos")

barplot(colSums(counts_nolow),names.arg=groups, xlab="Library name", ylab="Read count",col="yellow", main="Nof reads mapped per sample after collapsing")

dge <- DGEList(counts = counts_nolow, group = groups)
names(dge)
head(dge$counts) # original count matrix
levels((dge$samples$group))
# [1] "min" "pos"
dge$samples # contains a summary of your samples
# group lib.size norm.factors
# min1   min  2988747            1
# min2   min  3606323            1
# min3   min  3291364            1
# min4   min  3623719            1
# pos1   pos  2660719            1
# pos2   pos  3058679            1
# pos3   pos  3550961            1
# pos4   pos  3939584            1

dge <- calcNormFactors(dge) # normalize libs to prevent over expressed genes from blanking out rest
dge$samples
# group lib.size norm.factors
# min1   min  2988747    1.0148539
# min2   min  3606323    1.0144485
# min3   min  3291364    0.9907325
# min4   min  3623719    1.0622451
# pos1   pos  2660719    0.9646374
# pos2   pos  3058679    1.0039838
# pos3   pos  3550961    1.0011110
# pos4   pos  3939584    0.9519459

dge$samples$lib.size * dge$samples$norm.factors # effective library sizes
# [1] 3033142 3658429 3260861 3849278 2566629 3070864 3554906 3750271

design <- model.matrix(~groups, data = dge$samples)
design
# (Intercept) groupspos
# min1           1         0
# min2           1         0
# min3           1         0
# min4           1         0
# pos1           1         1
# pos2           1         1
# pos3           1         1
# pos4           1         1
# attr(,"assign")
# [1] 0 1
# attr(,"contrasts")
# attr(,"contrasts")$groups
# [1] "contr.treatment"

dge <- estimateDisp(dge, design = design) #Estimate Common, Trended and Tagwise Negative Binomial dispersions by weighted likelihood empirical Bayes

##########################################
#### VIZ Pre-DE analysis #####
##########################################

barplot(dge$counts,las=2,main="Normalized counts/library", names.arg = as.character(dge$samples$group))

#http://cgrlucb.wikispaces.com/edgeR+spring2013
library(RColorBrewer)
colors <- brewer.pal(9, "Set1")
#If we want the normalized pseudo-counts, useful for instance for cluster analysis,
scale = dge$samples$lib.size * dge$samples$norm.factors
normCounts = round(t(t(dge$counts)/scale)*mean(scale))
boxplot(log2(normCounts+1), las=2, col=colors[dge$samples$group],
        main="Normalized counts", names.arg = as.character(dge$samples$group))

# PCA
plotMDS(dge,labels=groups,col=c("darkblue","darkgreen")[factor(groups)]) #diff colors for each group, logFC

plotMeanVar(dge,show.tagwise.vars=TRUE,NBline=TRUE, main="edgeR all nolow collapsed") #each dot represents the estimated mean and variance for each gene, with binned variances as well as the trended common dispersion overlaid. Explore the mean-variance relationship for DGE data

plotBCV(dge) # the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million). higher counts  => lower variation

##############################################################
#### Pos(experimental) vs Min (control). Min is reference ####
##############################################################

# dispersion can be "common", "trended", "tagwise" or "auto". Default behavior ("auto" is to use most complex dispersions found in data object. table from exactTest() does not contain p-values adjusted for multiple testing
dge_et_PvM <- exactTest(dge, pair = c("min","pos")) # find DE genes in pos vs min
# adjusts the raw p-values using the False Discovery Rate (FDR) correction, and returns the top 10 differentially expressed genes
topTags(dge_et_PvM, sort.by = "logFC")
# Comparison of groups:  pos-min 
# logFC   logCPM       PValue          FDR
# gi|662190810|ref|XM_008469966.1| 11.188808 5.467616 3.743889e-31 1.371422e-28
# gi|662194745|ref|XM_008472112.1|  9.013972 3.343713 1.931989e-07 5.650885e-06
# gi|662209799|ref|XM_008480292.1| -8.691748 3.064760 3.166431e-08 1.107172e-06
# gi|662209681|ref|XR_541691.1|    -8.274905 2.702527 1.371027e-13 1.165373e-11
# gi|662224600|ref|XM_008488429.1| -8.242495 2.645404 4.842098e-07 1.313329e-05
# gi|662222418|ref|XM_008487230.1| -8.194321 2.645114 1.764837e-09 7.870151e-08
# gi|662218885|ref|XM_008485272.1| -7.673639 2.159983 2.357102e-27 6.593457e-25
# gi|662217236|ref|XM_008484370.1|  7.413035 1.886825 1.991258e-23 4.084735e-21
# gi|662220329|ref|XM_008486073.1|  7.035372 4.988612 6.045834e-97 2.325379e-93
# gi|662208136|ref|XR_541655.1|     6.979616 1.532070 4.257322e-06 8.899307e-05

dge_tt_PvM <- topTags(dge_et_PvM, sort.by = "logFC", n = nrow(dge_et_PvM))

dge_nc_PvM <- cpm(dge,normalized.lib.sizes=TRUE) #Computes counts per million (CPM)
head(dge_nc_PvM[rownames(dge_tt_PvM$table),order(groups)],5) # depth-adjusted reads per million for some of the top 5 differentially expressed genes

# DE gene counts for log fold change cutoffs
dim(subset(dge_tt_PvM$table, dge_tt_PvM$table$logFC > 1))
# [1] 888   4

# All DE genes, no need for log2FC filter as genes will be +ive for up and -ive for down regulated.
dim(subset(dge_tt_PvM$table, dge_tt_PvM$table$FDR < 0.01))
# [1] 1549    4

# output
write.csv(subset(dge_tt_PvM$table, dge_tt_PvM$table$logFC > 0.5 & dge_tt_PvM$table$FDR < 0.01), "edgeR_FDR1perc_LogFC0.5_nolow-POSvsMIN.csv")
write.csv(subset(dge_tt_PvM$table, dge_tt_PvM$table$FDR < 0.01), "edgeR_FDR1perc_nolow-POSvsMIN.csv")

##########################################
#### VIZ Min as reference #####
##########################################

de_genes <- rownames( dge_tt_PvM )[ dge_tt_PvM$table$PValue <= 0.01 ]
hist( dge_tt_PvM[de_genes[1:100],] , breaks=10 , xlab="Log Concentration" , col="red" , xlim=c(-18,-6) , ylim=c(0,0.4) , freq=FALSE , main="Poisson: Top 100" )# DOES NOT WORK

# MA  plot  that  shows  the  relationship  between  concentration  and  fold-change  across  the  genes.   The differentially expressed genes are colored red and the non-differentially expressed are colored black.  The orange dots represent genes in which the counts were zero in all samples of one of the groups

# first 100 DE
plotSmear(dge,de.tags=rownames(dge_tt_PvM[de_genes[1:100],]),main="log-Fold Change vs log-Conc FDR < .01")

# all DE genes 
plotSmear(dge,de.tags=rownames(dge_tt_PvM[de_genes,]),main="log-Fold Change vs log-Conc FDR < .01")


top <- topTags(dge_et_PvM, n=nrow(dge_et_PvM$counts))$table
de <- rownames(top[top$FDR<0.01,])
#volcano plot
plot(top$logFC, -log10(top$PValue), pch=20, cex=.5, ylab="-log10(p-value)",
     xlab="logFC", main="DE genes\nPos vs Min, FDR 0.01", col=as.numeric(rownames(top) %in% de)+1)
abline(v=c(-2, 2), col=colors[2])


# gene to gene heatmap
library(gplots)
heatmap.2(log(normCounts[de[1:400],]+1), ColSideColor=colors[groups],
          labCol=groups, main="Top 400 genes DE Pos vs Min\nFDR < 0.01", labRow=NA)
# sample to sample heatmap
heatmap.2(log(normCounts[de,]+1), ColSideColor=colors[groups],
          labCol=groups, main="All genes DE Pos vs Min\nFDR < 0.01", labRow=NA)

