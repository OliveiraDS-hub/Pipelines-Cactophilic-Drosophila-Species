library(DESeq2)
library(ggplot2)
####### CONTROL RNA-SEQ

counttable = read.csv("raw_counts_genes.csv", header = TRUE, row.names = 1, sep = ",")
head(counttable)
dim(counttable)

colData = data.frame(row.names=colnames(counttable))
colData$condition = factor(c(rep("dbuz",3), rep("dkoep",3)))                     
colData

dataset <- DESeqDataSetFromMatrix(countData = counttable, 
                                  colData = colData, 
                                  design = ~condition)
dds = DESeq(dataset)

plotDispEsts(dds)
resultsNames(dds)

#Get normalized counts
#norm_counts = counts(dds, normalized=TRUE)
#write.csv(norm_counts, "buzzatii_head_normcounts.csv", quote = FALSE)

dds$condition
resultsNames(dds)
## PCA ##

vsd <- vst(nsub= 150, dds, blind=TRUE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
head(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=5) +
  ggtitle("PCA buzzatii species - HEAD - Genes") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


##### POPULATIONS CONSTRASTS #####

## dbuz vs dkoep ##
dbuz_vs_dkoep <- results(dds, parallel=T,  contrast=c("condition","dbuz","dkoep"))

ix = which.max(dbuz_vs_dkoep$log2FoldChange)
res <- dbuz_vs_dkoep[order((dbuz_vs_dkoep$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c(rep("blue",3), rep("red",3),rep("blue",3), rep("red",3), rep("blue",3)), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dbuz_vs_dkoep, "DESEQout_dbuzvsdkoep.csv", quote = FALSE)

upreg = subset(dbuz_vs_dkoep, padj<0.05 & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dbuzvsdkoep_UP.csv", quote = FALSE)

downreg = subset(dbuz_vs_dkoep, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dbuzvsdkoep_DOWN.csv", quote = FALSE)

