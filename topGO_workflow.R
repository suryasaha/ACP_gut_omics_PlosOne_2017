rm(list=ls())
library(topGO)

acpgeneID2GO <- readMappings(file = "formatted.out.protein.fa.tsv.GO")
str(head(acpgeneID2GO))
# List of 6
# $ XP_008467331.1: chr "GO:0008270"
# $ XP_008467333.1: chr "GO:0016021"
# $ XP_008467334.1: chr [1:5] "GO:0004713" "GO:0004871" "GO:0005925" "GO:0006468" ...
# $ XP_008467335.1: chr [1:3] "GO:0005215" "GO:0006810" "GO:0016020"
# $ XP_008467336.1: chr "GO:0016020"
# $ XP_008467338.1: chr [1:2] "GO:0005509" "GO:0005515"

acpgeneNames <- names(acpgeneID2GO)
summary(acpgeneNames)
# Length     Class      Mode 
# 11637 character character 
head(acpgeneNames)
# [1] "XP_008467331.1" "XP_008467333.1" "XP_008467334.1" "XP_008467335.1"
# [5] "XP_008467336.1" "XP_008467338.1"

DEgenes <- read.delim(file = "DE_All_DESeq2_edgeR.names", header = FALSE)
as.character(DEgenes$V1)
# [1] "XP_008467337.1" "XP_008467367.1" "XP_008467374.1" "XP_008467392.1"
# [5] "XP_008467411.1" "XP_008467414.1" "XP_008467418.1" "XP_008467420.1"
# [9] "XP_008467426.1" "XP_008467428.1" "XP_008467430.1" "XP_008467431.1"
DEgeneListFactor <- factor(as.integer(acpgeneNames %in% as.character(DEgenes$V1)))
# [1] 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1
# [39] 0 0 0 0 0 0 0 1 1 0 0 0 0 1 0 1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0

#assigning names to factor
names(DEgeneListFactor) <- acpgeneNames
str(DEgeneListFactor)
# Factor w/ 2 levels "0","1": 1 2 1 1 1 1 1 1 1 1 ...
# - attr(*, "names")= chr [1:11637] "XP_008467331.1" "XP_008467333.1" "XP_008467334.1" "XP_008467335.1" ...

### creating object with molecular function ontology GO hierarchy
acpGOdata <- new("topGOdata", ontology = "MF", allGenes = DEgeneListFactor, annot = annFUN.gene2GO, gene2GO = acpgeneID2GO)
acpGOdata
# ------------------------- topGOdata object -------------------------
#   
#   Description:
#   -   
#   
#   Ontology:
#   -  MF 
# 
# 11637 available genes (all genes from the array):
#   - symbol:  XP_008467331.1 XP_008467333.1 XP_008467334.1 XP_008467335.1 XP_008467336.1  ...
# - 1152  significant genes. 
# 
# 10188 feasible genes (genes that can be used in the analysis):
#   - symbol:  XP_008467331.1 XP_008467334.1 XP_008467335.1 XP_008467338.1 XP_008467339.1  ...
# - 1015  significant genes. 
# 
# GO graph (nodes with at least  1  genes):
#   - a graph with directed edges
# - number of nodes = 1223 
# - number of edges = 1562 
# 
# ------------------------- topGOdata object -------------------------

# One important point to notice is that not all the genes that are provided by geneList, the initial gene universe, can be annotated to the GO. This can be seen by comparing the number of all available genes, the genes present in geneList, with the number of feasible genes.

numGenes(acpGOdata)
# [1] 10188

length(sigGenes(acpGOdata))
# [1] 1015
numSigGenes(acpGOdata)

graph(acpGOdata)
# A graphNEL graph with directed edges
# Number of Nodes = 1223 
# Number of Edges = 1562 


### creating object with molecular function ontology GO hierarchy and min node size to avoid artifacts
acpGOnode5MFdata <- new("topGOdata", description = "MF GO and min node size 5", ontology = "MF", allGenes = DEgeneListFactor, annot = annFUN.gene2GO, gene2GO = acpgeneID2GO, nodeSize = 5)
acpGOnode5MFdata
# ------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  MF GO and min node size 5 
# 
# Ontology:
#   -  MF 
# 
# 11637 available genes (all genes from the array):
#   - symbol:  XP_008467331.1 XP_008467333.1 XP_008467334.1 XP_008467335.1 XP_008467336.1  ...
# - 1152  significant genes. 
# 
# 10188 feasible genes (genes that can be used in the analysis):
#   - symbol:  XP_008467331.1 XP_008467334.1 XP_008467335.1 XP_008467338.1 XP_008467339.1  ...
# - 1015  significant genes. 
# 
# GO graph (nodes with at least  5  genes):
#   - a graph with directed edges
# - number of nodes = 606 
# - number of edges = 785 
# 
# ------------------------- topGOdata object -------------------------

numGenes(acpGOnode5MFdata)
# [1] 10188

length(sigGenes(acpGOnode5MFdata))
# [1] 1015

graph(acpGOnode5MFdata)
# A graphNEL graph with directed edges
# Number of Nodes = 606 
# Number of Edges = 785 

### creating object with molecular function ontology GO hierarchy and min node size to avoid artifacts
acpGOnode5BPdata <- new("topGOdata", description = "BP GO and min node size 5", ontology = "BP", allGenes = DEgeneListFactor, annot = annFUN.gene2GO, gene2GO = acpgeneID2GO, nodeSize = 5)
acpGOnode5BPdata
# ------------------------- topGOdata object -------------------------
#   
#   Description:
#   -  BP GO and min node size 5 
# 
# Ontology:
#   -  BP 
# 
# 11637 available genes (all genes from the array):
#   - symbol:  XP_008467331.1 XP_008467333.1 XP_008467334.1 XP_008467335.1 XP_008467336.1  ...
# - 1152  significant genes. 
# 
# 6658 feasible genes (genes that can be used in the analysis):
#   - symbol:  XP_008467334.1 XP_008467335.1 XP_008467339.1 XP_008467346.1 XP_008467352.1  ...
# - 769  significant genes. 
# 
# GO graph (nodes with at least  5  genes):
#   - a graph with directed edges
# - number of nodes = 1151 
# - number of edges = 2456 
# 
# ------------------------- topGOdata object -------------------------
numGenes(acpGOnode5BPdata)
# [1] 6658

graph(acpGOnode5BPdata)
# A graphNEL graph with directed edges
# Number of Nodes = 1151 
# Number of Edges = 2456

# nof sig genes
str(sigGenes(acpGOnode5MFdata)) # 1015
str(sigGenes(acpGOnode5BPdata)) # 769

# GO terms used
length(usedGO(acpGOnode5MFdata))# 606
length(usedGO(acpGOnode5BPdata))# 1151

#select 10 random GO terms, count the number of genes with those terms and obtain their names
sel.terms <- sample(usedGO(acpGOnode5MFdata), 10)
num.ann.genes <- countGenesInTerm(acpGOnode5MFdata, sel.terms) ## the number of genes with selected GO terms
num.ann.genes
ann.genes <- genesInTerm(acpGOnode5MFdata, sel.terms) ## get the annotations
head(ann.genes)

termStat(acpGOnode5MFdata)
# Annotated Significant Expected
# GO:0000030         7           1     0.70
# GO:0000049        12           2     1.20
# GO:0000062        12           1     1.20
termStat(acpGOnode5BPdata)
# Annotated Significant Expected
# GO:0000041         9           3     1.04
# GO:0000042         5           0     0.58
# GO:0000070        35           1     4.04

# one first defines a test statistic for the chosen algorithm, meaning that an instance of object specific for the algorithm is constructed in which only the test statistic must be specified, and then calls a generic function (interface) to run the algorithm.

whichTests()
# [1] "fisher"     "ks"         "t"          "globaltest" "sum"        "ks.ties"   
whichAlgorithms()
# [1] "classic"     "elim"        "weight"      "weight01"    "lea"         "parentchild"


fishertest.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_acpGOnode5BPdata <- getSigGroups(acpGOnode5BPdata, fishertest.stat) 
resultFisher_acpGOnode5BPdata
# Description: BP GO and min node size 5 
# Ontology: BP 
# 'classic' algorithm with the 'Fisher test' test
# 1151 GO terms scored: 54 terms with p < 0.01
# Annotation data:
#   Annotated genes: 6658 
# Significant genes: 769 
# Min. no. of genes annotated to a GO: 5 
# Nontrivial nodes: 730 

resultFisher_acpGOnode5MFdata <- getSigGroups(acpGOnode5MFdata, fishertest.stat) 
resultFisher_acpGOnode5MFdata
# Description: MF GO and min node size 5 
# Ontology: MF 
# 'classic' algorithm with the 'Fisher test' test
# 606 GO terms scored: 46 terms with p < 0.01
# Annotation data:
#   Annotated genes: 10188 
# Significant genes: 1015 
# Min. no. of genes annotated to a GO: 5 
# Nontrivial nodes: 438 

weightcounttest.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test")
resultweightFisher_acpGOnode5BPdata <- getSigGroups(acpGOnode5BPdata, weightcounttest.stat) 
resultweightFisher_acpGOnode5BPdata
# Description: BP GO and min node size 5 
# Ontology: BP 
# 'weight' algorithm with the 'Fisher test : ratio' test
# 1151 GO terms scored: 12 terms with p < 0.01
# Annotation data:
#   Annotated genes: 6658 
# Significant genes: 769 
# Min. no. of genes annotated to a GO: 5 
# Nontrivial nodes: 730 

resultweightFisher_acpGOnode5MFdata <- getSigGroups(acpGOnode5MFdata, weightcounttest.stat) 
resultweightFisher_acpGOnode5MFdata
# Description: MF GO and min node size 5 
# Ontology: MF 
# 'weight' algorithm with the 'Fisher test : ratio' test
# 606 GO terms scored: 17 terms with p < 0.01
# Annotation data:
#   Annotated genes: 10188 
# Significant genes: 1015 
# Min. no. of genes annotated to a GO: 5 
# Nontrivial nodes: 438 

elimcounttest.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
resultelimFisher_acpGOnode5BPdata <- getSigGroups(acpGOnode5BPdata, elimcounttest.stat)
resultelimFisher_acpGOnode5BPdata
# Description: BP GO and min node size 5 
# Ontology: BP 
# 'elim' algorithm with the 'Fisher test : 0.01' test
# 1151 GO terms scored: 13 terms with p < 0.01
# Annotation data:
#   Annotated genes: 6658 
# Significant genes: 769 
# Min. no. of genes annotated to a GO: 5 
# Nontrivial nodes: 730 

resultelimFisher_acpGOnode5MFdata <- getSigGroups(acpGOnode5MFdata, elimcounttest.stat) 
resultelimFisher_acpGOnode5MFdata
# Description: MF GO and min node size 5 
# Ontology: MF 
# 'elim' algorithm with the 'Fisher test : 0.01' test
# 606 GO terms scored: 24 terms with p < 0.01
# Annotation data:
#   Annotated genes: 10188 
# Significant genes: 1015 
# Min. no. of genes annotated to a GO: 5 
# Nontrivial nodes: 438

############################################
##
## Using resultFisher_acpGOnode5BPdata and resultFisher_acpGOnode5MFdata
##
############################################

head(score(resultFisher_acpGOnode5BPdata))
# GO:0000041 GO:0000070 GO:0000122 GO:0000154 GO:0000184 GO:0000226 
# 0.0757147  0.9865310  0.4530642  0.4587429  0.1459316  0.9537745
head(score(resultFisher_acpGOnode5MFdata))
# GO:0000030 GO:0000049 GO:0000062 GO:0000166 GO:0000287 GO:0000988 
# 0.52042327 0.33930236 0.71636624 0.98608333 0.04497358 0.99423055

hist(score(resultFisher_acpGOnode5BPdata), xlab = "BP p-values")
hist(score(resultFisher_acpGOnode5MFdata), xlab = "MF p-values")
