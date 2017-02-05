# R-commands used in BrÃ¥te et al. 2015 Regulatory RNA at the root of animals: dynamic expression of developmental lncRNAs in the calcisponge Sycon ciliatum 

# 1. Correlation test of lincRNAs and neighboring genes
# 2. DESeq2 analysis (PCA plots and Differential expression tests)
# 3. WGCNA analysis


### 1. Correlation test of lincRNAs and neighboring genes ###

#Import count matrices of lincRNA and nearest neighbor. The neighbor should be on the line directly underneath the RNA. 
d_right = read.table(file="lincRNAs.right.neighbor.counts.csv", header=T, sep="\t")
d_left = read.table(file="lincRNAs.left.neighbor.counts.csv", header=T, sep="\t")

#remove gene names
d_right$id = NULL
d_left$id = NULL

#convert to matrix
d_right = as.matrix(d_right)
d_left = as.matrix(d_left)

#run correlation test
cor_right = sapply(2*(1:(nrow(d_right)/2)), function(pair) unname(cor.test(d_right[pair-1,], d_right[pair,], method="spearman")$estimate))
cor_left = sapply(2*(1:(nrow(d_left)/2)), function(pair) unname(cor.test(d_left[pair-1,], d_left[pair,], method="spearman")$estimate))

#get p-values
Pvalue_right = sapply(2*(1:(nrow(d_right)/2)), function(pair) unname(cor.test(d_right[pair-1,], d_right[pair,], method="spearman")$p.value))
Pvalue_left = sapply(2*(1:(nrow(d_left)/2)), function(pair) unname(cor.test(d_left[pair-1,], d_left[pair,], method="spearman")$p.value))

#adjust p-values for multiple comparisons
Adjusted_right = p.adjust(Pvalue_right, "BH")
Adjusted_left = p.adjust(Pvalue_left, "BH")
#######################


### 2. DESeq2 analysis ###
library(DESeq2)
#Same commands run for coding and lincRNAs. Only lincRNAs shown here.

#Import the lincRNA matrix, with the correlated lincRNAs (rho > 0.6) removed. 
lincRNAcountMatrix = read.csv("lincRNA_counts.csv", header=TRUE, sep="\t", check.names=FALSE, row.names=1)

#Setup metadata
samples = data.frame(row.names=c("05.apr", "11.apr", "28.apr", "11.mai", "13.mai", "16.mai", "18.mai", "23-05fert-1", "23-05fert-2", "26-05-early-cleav-1", "26-05early-cleav-2", "30-05late-cleav-1", "30-05late-cleav-2", "06-06preinv-1", "06-06preinv-2", "09-06preinv-1", "09-06preinv-2", "17-06ear-postinv-1", "17-06-ear-postinv-2", "R1D0B1", "R2D0B1", "R3D0B1"), 
                     condition=as.factor(c("Early vitellogenesis", "Early vitellogenesis", "Early vitellogenesis", "Late vitellogenesis", "Late vitellogenesis", "Late vitellogenesis", "Late vitellogenesis", "Fertilization", "Fertilization", "Early cleavage", "Early cleavage", "Late cleavage", "Late cleavage", "Early preinversion", "Early preinversion", "Late preinversion", "Late preinversion", "Early postinversion", "Early postinversion", "Non reproductive", "Non reproductive", "Non reproductive")))

#Create the DESeqDataSet
ddslincRNA = DESeqDataSetFromMatrix(countData=lincRNAcountMatrix, colData=samples, design= ~ condition)
colData(ddslincRNA)$condition=factor(colData(ddslincRNA)$condition, 
                              levels=c("Early vitellogenesis", "Late vitellogenesis", "Fertilization", 
                                       "Early cleavage", "Late cleavage", "Early preinversion", "Late preinversion", 
                                       "Early postinversion", "Late postinversion", "Non reproductive"))

#Relevel - set the Non reproductive samples as the reference condition
ddslincRNA$condition = relevel(ddslincRNA$condition, ref="Non reproductive")

#Run the DE analysis
ddslincRNA = estimateSizeFactors(ddslincRNA)
ddslincRNA = estimateDispersions(ddslincRNA)
ddslincRNA = nbinomWaldTest(ddslincRNA, maxit = 500)

#regularized log transformation
rld = rlog(ddslincRNA)

#PCA plot
print(plotPCA(rld, intgroup="condition"))


##Extract the results from all conditions vs Non reproductive from the Wald test##

#Early vitellogenesis vs Non reproductive
reslincRNAEarlyVitVsNonRepr = results(ddslincRNA, contrast=c("condition", "Early vitellogenesis", "Non reproductive"))
reslincRNAEarlyVitVsNonRepr = reslincRNAEarlyVitVsNonRepr[order(reslincRNAEarlyVitVsNonRepr$padj),]
reslincRNASigEarlyVitVsNonRepr = reslincRNAEarlyVitVsNonRepr[which(reslincRNAEarlyVitVsNonRepr$padj < 0.1), ]
reslincRNASigEarlyVitVsNonReprUpReg = reslincRNASigEarlyVitVsNonRepr[which(reslincRNASigEarlyVitVsNonRepr$log2FoldChange > 0), ]

#Late vitellogenesis vs Non reproductive
reslincRNALateVitVsNonRepr = results(ddslincRNA, contrast=c("condition", "Late vitellogenesis", "Non reproductive"))
reslincRNALateVitVsNonRepr = reslincRNALateVitVsNonRepr[order(reslincRNALateVitVsNonRepr$padj),]
reslincRNASigLateVitVsNonRepr = reslincRNALateVitVsNonRepr[which(reslincRNALateVitVsNonRepr$padj < 0.1), ]
reslincRNASigLateVitVsNonReprUpReg = reslincRNASigLateVitVsNonRepr[which(reslincRNASigLateVitVsNonRepr$log2FoldChange > 0), ]

#Fertilization vs Non reproductive
reslincRNAFertVsNonRepr = results(ddslincRNA, contrast=c("condition", "Fertilization", "Non reproductive"))
reslincRNAFertVsNonRepr = reslincRNAFertVsNonRepr[order(reslincRNAFertVsNonRepr$padj),]
reslincRNASigFertVsNonRepr = reslincRNAFertVsNonRepr[which(reslincRNAFertVsNonRepr$padj < 0.1), ]
reslincRNASigFertVsNonReprUpReg = reslincRNASigFertVsNonRepr[which(reslincRNASigFertVsNonRepr$log2FoldChange > 0), ]

#Early cleavage vs Non reproductive
reslincRNAEarlyClVsNonRepr = results(ddslincRNA, contrast=c("condition", "Early cleavage", "Non reproductive"))
reslincRNAEarlyClVsNonRepr = reslincRNAEarlyClVsNonRepr[order(reslincRNAEarlyClVsNonRepr$padj),]
reslincRNASigEarlyClVsNonRepr = reslincRNAEarlyClVsNonRepr[which(reslincRNAEarlyClVsNonRepr$padj < 0.1), ]
reslincRNASigEarlyClVsNonReprUpReg = reslincRNASigEarlyClVsNonRepr[which(reslincRNASigEarlyClVsNonRepr$log2FoldChange > 0), ]

#Late cleavage vs Non reproductive
reslincRNALateClVsNonRepr = results(ddslincRNA, contrast=c("condition", "Late cleavage", "Non reproductive"))
reslincRNALateClVsNonRepr = reslincRNALateClVsNonRepr[order(reslincRNALateClVsNonRepr$padj),]
reslincRNASigLateClVsNonRepr = reslincRNALateClVsNonRepr[which(reslincRNALateClVsNonRepr$padj < 0.1), ]
reslincRNASigLateClVsNonReprUpReg = reslincRNASigLateClVsNonRepr[which(reslincRNASigLateClVsNonRepr$log2FoldChange > 0), ]

#Early preinversion vs Non reproductive
reslincRNAEarlyPreVsNonRepr = results(ddslincRNA, contrast=c("condition", "Early preinversion", "Non reproductive"))
reslincRNAEarlyPreVsNonRepr = reslincRNAEarlyPreVsNonRepr[order(reslincRNAEarlyPreVsNonRepr$padj),]
reslincRNASigEarlyPreVsNonRepr = reslincRNAEarlyPreVsNonRepr[which(reslincRNAEarlyPreVsNonRepr$padj < 0.1), ]
reslincRNASigEarlyPreVsNonReprUpReg = reslincRNASigEarlyPreVsNonRepr[which(reslincRNASigEarlyPreVsNonRepr$log2FoldChange > 0), ]

#Late preinversion vs Non reproductive
reslincRNALatePreVsNonRepr = results(ddslincRNA, contrast=c("condition", "Late preinversion", "Non reproductive"))
reslincRNALatePreVsNonRepr = reslincRNALatePreVsNonRepr[order(reslincRNALatePreVsNonRepr$padj),]
reslincRNASigLatePreVsNonRepr = reslincRNALatePreVsNonRepr[which(reslincRNALatePreVsNonRepr$padj < 0.1), ]
reslincRNASigLatePreVsNonReprUpReg = reslincRNASigLatePreVsNonRepr[which(reslincRNASigLatePreVsNonRepr$log2FoldChange > 0), ]

#Early postinversion vs Non reproductive
reslincRNAEarlyPostVsNonRepr = results(ddslincRNA, contrast=c("condition", "Early postinversion", "Non reproductive"))
reslincRNAEarlyPostVsNonRepr = reslincRNAEarlyPostVsNonRepr[order(reslincRNAEarlyPostVsNonRepr$padj),]
reslincRNASigEarlyPostVsNonRepr = reslincRNAEarlyPostVsNonRepr[which(reslincRNAEarlyPostVsNonRepr$padj < 0.1), ]
reslincRNASigEarlyPostVsNonReprUpReg = reslincRNASigEarlyPostVsNonRepr[which(reslincRNASigEarlyPostVsNonRepr$log2FoldChange > 0), ]

#Late postinversion vs Non reproductive
reslincRNALatePostVsNonRepr = results(ddslincRNA, contrast=c("condition", "Late postinversion", "Non reproductive"))
reslincRNALatePostVsNonRepr = reslincRNALatePostVsNonRepr[order(reslincRNALatePostVsNonRepr$padj),]
reslincRNASigLatePostVsNonRepr = reslincRNALatePostVsNonRepr[which(reslincRNALatePostVsNonRepr$padj < 0.1), ]
reslincRNASigLatePostVsNonReprUpReg = reslincRNASigLatePostVsNonRepr[which(reslincRNASigLatePostVsNonRepr$log2FoldChange > 0), ]

#######################


### 3. WGCNA analysis ###
library(WGCNA)


options(stringsAsFactors = FALSE)

#Import normalized gene count matrix of developmental coding genes and lincRNAs, count values averaged between replicates

countMatrix = read.csv("count_matrix.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)

#Filter out lowly expressed genes, i.e. keep genes expressed more than 5 in at least 3 samples
use = rowSums(countMatrix > 5) >=3
countMatrixFiltered = countMatrix[ use, ]

#log2 transform (x+1 to avoid zeroes)
log2countMatrixFiltered = log2(countMatrixFiltered+1)

#Extract the genes with the highest variance (upper 75 percentile). 
library(genefilter)
rv = rowVars(log2countMatrixFiltered)
q75 = quantile(rowVars(log2countMatrixFiltered), .75)
q75log2countMatrixFiltered = log2countMatrixFiltered[rv > q75, ]

# Flip the count matrix
Flippedq75log2countMatrixFiltered = as.data.frame(t(q75log2countMatrixFiltered))
gsg = goodSamplesGenes(Flippedq75log2countMatrixFiltered, verbose = 3)
gsg$allOK

# Define the adjancency matric
ADJ1 = abs(cor(Flippedq75log2countMatrixFiltered, use="p"))^18 # Using the Pearson correlation

#Turn adjacency into a measure of dissimilarity
dissADJ=1-ADJ1

#Use of topological overlap to define dissimilarity
dissTOM = TOMdist(ADJ1)
collectGarbage()

## Using topological overlap
#Calculate the dendrogram
hierTOM = hclust(as.dist(dissTOM), method="average")

#Define modules
colorDynamicTOM = labels2colors(cutreeDynamic(hierTOM, method="tree"))

#plot dendrogram and modules
sizeGrWindow(10,5)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE,
                    main = "Gene dendrogram, dynamicTOM, power18")

