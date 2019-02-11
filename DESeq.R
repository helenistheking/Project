# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2", version = "3.8")

#Import DESeq2 library in R
library("DESeq2")

#Load gene(/transcript) count matrix and labels

countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
colData <- read.csv("PHENO_DATA.csv", sep=",", row.names=1)
##Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData))

#should be TRUE

countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))


#Create a DESeqDataSet from count matrix and labels
#technically a subclass of summarizedexperiment class
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, design = ~ cell_type + treatment)

#Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "treatment_DMSO_vs_CHX")
chx_dmsoresults

#Sort by adjusted p-value and display
(resOrdered <- res[order(res$padj), ])

#log fold change shrinkage is useful for visualization and ranking of genes
#use apeglm method for size shrinkage
resLFC <- lfcShrink(dds, coef="treatment_DMSO_vs_CHX", type="apeglm")
resLFC

#could use parallelizaition- if i had lots of samples, but I dont
#independent hypothesis weighting 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("IHW", version = "3.8")

library("IHW")
resIHW <- results(dds, filterFun= ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult


#plotMA shows the log2 fold changes atributated to a given varable over the mean normalized counts for all the samples
plotMA(resLFC, ylim=c(-2,2))

#large LFC due to being solo samples, AND no repeat :(

#identify the row number of indiviuals resulting indices 
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

write.csv(as.data.frame(resLFC), file="condition_CHXvsDMSO.csv")
