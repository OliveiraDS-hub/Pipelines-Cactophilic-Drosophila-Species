library(DESeq2)
library(ggplot2)
####### CONTROL RNA-SEQ

counttable = read.csv("raw_counts_genes_TEs.csv", header = TRUE, row.names = 1, sep = ",")
head(counttable)

colData = data.frame(row.names=colnames(counttable))
colData$condition = factor( c( rep("dari_head",3), rep("dmoj01_head",3), rep("dmoj20_head",3), rep("dmoj22_head",3), rep("dmoj26_head",3)))                     

dataset <- DESeqDataSetFromMatrix(countData = counttable, 
                                  colData = colData, 
                                  design = ~condition)
dds = DESeq(dataset)

dds = DESeq(ddsHTSeq, parallel = T)
plotDispEsts(dds)
resultsNames(dds)

#Get normalized counts
#norm_counts = counts(dds, normalized=TRUE)
#write.csv(norm_counts, "mulleri_head_normcounts.csv", quote = FALSE)

dds$condition
resultsNames(dds)
## PCA ##

vsd <- vst(dds, nsub=200, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
head(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=5) +
  ggtitle("PCA mulleri species - TEs") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


##### POPULATIONS CONSTRASTS #####
dds$condition
## dmoj01 vs dmoj20 ##

dmoj_01vsdmoj_20 <- results(dds, parallel=T,  contrast=c("condition","dmoj01_head","dmoj20_head"))
#dmoj_01vsdmoj_20 <- lfcShrink(dds, parallel=T,type="apeglm", contrast=c("condition","dmoj01_head","dmoj20_head"))

ix = which.max(dmoj_01vsdmoj_20$log2FoldChange)
res <- dmoj_01vsdmoj_20[order((dmoj_01vsdmoj_20$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c("blue","blue","blue","blue","red","red","red","red"), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dmoj_01vsdmoj_20, "DESEQout_dmoj_01vsdmoj_20.csv")

claudia = subset(dmoj_01vsdmoj_20, padj<0.05 & log2FoldChange > | 1 & log2FoldChange < -1)
dim(claudia)
upreg = subset(dmoj_01vsdmoj_20, padj<0.05 & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_01vsdmoj_20_UP.csv")

downreg = subset(dmoj_01vsdmoj_20, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_01vsdmoj_20_DOWN.csv")



## dmoj01 vs dmoj22 ##

dmoj_01vsdmoj_22 <- results(dds, parallel=T,  contrast=c("condition","dmoj01_head","dmoj22_head"))
ix = which.max(dmoj_01vsdmoj_22$log2FoldChange)
res <- dmoj_01vsdmoj_22[order((dmoj_01vsdmoj_22$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c("blue","blue","blue","blue","red","red","red","red"), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dmoj_01vsdmoj_22, "DESEQout_dmoj_01vsdmoj_22.csv")

upreg = subset(dmoj_01vsdmoj_22, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_01vsdmoj_22_UP.csv")

downreg = subset(dmoj_01vsdmoj_22, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_01vsdmoj_22_DOWN.csv")





## dmoj01 vs dmoj26 ##

dmoj_01vsdmoj_26 <- results(dds, parallel=T,  contrast=c("condition","dmoj01_head","dmoj26_head"))
ix = which.max(dmoj_01vsdmoj_26$log2FoldChange)
res <- dmoj_01vsdmoj_26[order((dmoj_01vsdmoj_26$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c("blue","blue","blue","blue","red","red","red","red"), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dmoj_01vsdmoj_26, "DESEQout_dmoj_01vsdmoj_26.csv")

upreg = subset(dmoj_01vsdmoj_26, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_01vsdmoj_26_UP.csv")

downreg = subset(dmoj_01vsdmoj_26, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_01vsdmoj_26_DOWN.csv")




##### POPULATIONS CONSTRASTS #####


## dmoj20 vs dmoj01 ##

dmoj_20vsdmoj_01 <- results(dds, parallel=T,  contrast=c("condition","dmoj20","dmoj01"))
ix = which.max(dmoj_20vsdmoj_01$log2FoldChange)
res <- dmoj_20vsdmoj_01[order((dmoj_20vsdmoj_01$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c(rep("blue",3), rep("red",3),rep("blue",3), rep("red",3), rep("blue",3)), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dmoj_20vsdmoj_01, "DESEQout_dmoj_20vsdmoj_01.csv")

upreg = subset(dmoj_20vsdmoj_01, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_20vsdmoj_01_UP.csv")

downreg = subset(dmoj_20vsdmoj_01, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_20vsdmoj_01_DOWN.csv")


## dmoj20 vs dmoj22 ##

dmoj_20vsdmoj_22 <- results(dds, parallel=T,  contrast=c("condition","dmoj20","dmoj22"))
ix = which.max(dmoj_20vsdmoj_22$log2FoldChange)
res <- dmoj_20vsdmoj_22[order((dmoj_20vsdmoj_22$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c(rep("blue",3), rep("red",3),rep("blue",3), rep("red",3), rep("blue",3)), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dmoj_20vsdmoj_22, "DESEQout_dmoj_20vsdmoj_22.csv")

upreg = subset(dmoj_20vsdmoj_22, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_20vsdmoj_22_UP.csv")

downreg = subset(dmoj_20vsdmoj_22, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_20vsdmoj_22_DOWN.csv")




## dmoj20 vs dmoj26 ##

dmoj_20vsdmoj_26 <- results(dds, parallel=T,  contrast=c("condition","dmoj20","dmoj26"))
ix = which.max(dmoj_20vsdmoj_26$log2FoldChange)
res <- dmoj_20vsdmoj_26[order((dmoj_20vsdmoj_26$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c(rep("blue",3), rep("red",3),rep("blue",3), rep("red",3), rep("blue",3)), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dmoj_20vsdmoj_26, "DESEQout_dmoj_20vsdmoj_26.csv")

upreg = subset(dmoj_20vsdmoj_26, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_20vsdmoj_26_UP.csv")

downreg = subset(dmoj_20vsdmoj_26, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_20vsdmoj_26_DOWN.csv")


## dmoj20 vs dari ##

dmoj20vsdari <- results(dds, parallel=T,  contrast=c("condition","dmoj20","dari"))

# Save result
write.csv(dmoj20vsdari, "DESEQout_dmoj_20vsdari.csv", quote = FALSE)

upreg = subset(dmoj20vsdari, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_20vsdari_UP.csv", quote = FALSE)

downreg = subset(dmoj20vsdari, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_20vsdari_DOWN.csv", quote = FALSE)










##### POPULATIONS CONSTRASTS #####

## dmoj22 vs dmoj01 ##

dmoj_22vsdmoj_01 <- results(dds, parallel=T,  contrast=c("condition","dmoj22","dmoj01"))
ix = which.max(dmoj_22vsdmoj_01$log2FoldChange)
res <- dmoj_22vsdmoj_01[order((dmoj_22vsdmoj_01$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c(rep("blue",3), rep("red",3),rep("blue",3), rep("red",3), rep("blue",3)), las=2, main=rownames(dds)[ ix  ]  )



# Save result
write.csv(dmoj_22vsdmoj_01, "DESEQout_dmoj_22vsdmoj_01.csv")

upreg = subset(dmoj_22vsdmoj_01, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_22vsdmoj_01_UP.csv")

downreg = subset(dmoj_22vsdmoj_01, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_22vsdmoj_01_DOWN.csv")





## dmoj22 vs dmoj20 ##

dmoj_22vsdmoj_20 <- results(dds, parallel=T,  contrast=c("condition","dmoj22","dmoj20"))
ix = which.max(dmoj_22vsdmoj_20$log2FoldChange)
res <- dmoj_22vsdmoj_20[order((dmoj_22vsdmoj_20$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c(rep("blue",3), rep("red",3),rep("blue",3), rep("red",3), rep("blue",3)), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dmoj_22vsdmoj_20, "DESEQout_dmoj_22vsdmoj_20.csv")

upreg = subset(dmoj_22vsdmoj_20, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_22vsdmoj_20_UP.csv")

downreg = subset(dmoj_22vsdmoj_20, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_22vsdmoj_20_DOWN.csv")







## dmoj22 vs dmoj26 ##

dmoj_22vsdmoj_26 <- results(dds, parallel=T,  contrast=c("condition","dmoj22","dmoj26"))
ix = which.max(dmoj_22vsdmoj_26$log2FoldChange)
res <- dmoj_22vsdmoj_26[order((dmoj_22vsdmoj_26$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c(rep("blue",3), rep("red",3),rep("blue",3), rep("red",3), rep("blue",3)), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dmoj_22vsdmoj_26, "DESEQout_dmoj_22vsdmoj_26.csv")

upreg = subset(dmoj_22vsdmoj_26, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_22vsdmoj_26_UP.csv")

downreg = subset(dmoj_22vsdmoj_26, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_22vsdmoj_26_DOWN.csv")

## dmoj22 vs dari ##

dmoj22vsdari <- results(dds, parallel=T,  contrast=c("condition","dmoj22","dari"))

# Save result
write.csv(dmoj22vsdari, "DESEQout_dmoj_22vsdari.csv", quote = FALSE)

upreg = subset(dmoj22vsdari, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_22vsdari_UP.csv", quote = FALSE)

downreg = subset(dmoj22vsdari, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_22vsdari_DOWN.csv", quote = FALSE)






## dmoj26 vs dmoj01 ##

dmoj_26vsdmoj_01 <- results(dds, parallel=T,  contrast=c("condition","dmoj26","dmoj01"))
ix = which.max(dmoj_26vsdmoj_01$log2FoldChange)
res <- dmoj_26vsdmoj_01[order((dmoj_26vsdmoj_01$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c(rep("blue",3), rep("red",3),rep("blue",3), rep("red",3), rep("blue",3)), las=2, main=rownames(dds)[ ix  ]  )


# Save result
write.csv(dmoj_26vsdmoj_01, "DESEQout_dmoj_26vsdmoj_01.csv", quote = FALSE)

upreg = subset(dmoj_26vsdmoj_01, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_26vsdmoj_01_UP.csv", quote = FALSE)

downreg = subset(dmoj_26vsdmoj_01, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_26vsdmoj_01_DOWN.csv", quote = FALSE)





## dmoj26 vs dmoj20 ##

dmoj_26vsdmoj_20 <- results(dds, parallel=T,  contrast=c("condition","dmoj26","dmoj20"))

# Save result
write.csv(dmoj_26vsdmoj_20, "DESEQout_dmoj_26vsdmoj_20.csv", quote = FALSE)

upreg = subset(dmoj_26vsdmoj_20, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_26vsdmoj_20_UP.csv", quote = FALSE)

downreg = subset(dmoj_26vsdmoj_20, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_26vsdmoj_20_DOWN.csv", quote = FALSE)




## dmoj26 vs dmoj22 ##

dmoj_26vsdmoj_22 <- results(dds, parallel=T,  contrast=c("condition","dmoj26","dmoj22"))

# Save result
write.csv(dmoj_26vsdmoj_22, "DESEQout_dmoj_26vsdmoj_22.csv", quote = FALSE)

upreg = subset(dmoj_26vsdmoj_22, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_26vsdmoj_22_UP.csv", quote = FALSE)

downreg = subset(dmoj_26vsdmoj_22, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_26vsdmoj_22_DOWN.csv", quote = FALSE)

## dmoj26 vs dari ##

dmoj26vsdari <- results(dds, parallel=T,  contrast=c("condition","dmoj26","dari"))

# Save result
write.csv(dmoj26vsdari, "DESEQout_dmoj_26vsdari.csv", quote = FALSE)

upreg = subset(dmoj26vsdari, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_26vsdari_UP.csv", quote = FALSE)

downreg = subset(dmoj26vsdari, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_26vsdari_DOWN.csv", quote = FALSE)








## dmoj01 vs dmoj20 ##

dmoj_01vsdmoj_20 <- results(dds, parallel=T,  contrast=c("condition","dmoj01","dmoj20"))

# Save result
write.csv(dmoj_01vsdmoj_20, "DESEQout_dmoj_01vsdmoj_20.csv", quote = FALSE)

upreg = subset(dmoj_01vsdmoj_20, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_01vsdmoj_20_UP.csv", quote = FALSE)

downreg = subset(dmoj_01vsdmoj_20, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_01vsdmoj_20_DOWN.csv", quote = FALSE)









## dmoj01 vs dmoj22 ##

dmoj_01vsdmoj_22 <- results(dds, parallel=T,  contrast=c("condition","dmoj01","dmoj22"))

# Save result
write.csv(dmoj_01vsdmoj_22, "DESEQout_dmoj_01vsdmoj_22.csv", quote = FALSE)

upreg = subset(dmoj_01vsdmoj_22, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_01vsdmoj_22_UP.csv", quote = FALSE)

downreg = subset(dmoj_01vsdmoj_22, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_01vsdmoj_22_DOWN.csv", quote = FALSE)




## dmoj01 vs dmoj26 ##

dmoj_01vsdmoj_26 <- results(dds, parallel=T,  contrast=c("condition","dmoj01","dmoj26"))

# Save result
write.csv(dmoj_01vsdmoj_26, "DESEQout_dmoj_01vsdmoj_26.csv", quote = FALSE)

upreg = subset(dmoj_01vsdmoj_26, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_01vsdmoj_26_UP.csv", quote = FALSE)

downreg = subset(dmoj_01vsdmoj_26, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_01vsdmoj_26_DOWN.csv", quote = FALSE)


## dmoj01 vs dari ##

dmoj_01vsdari <- results(dds, parallel=T,  contrast=c("condition","dmoj01","dari"))

# Save result
write.csv(dmoj_01vsdari, "DESEQout_dmoj_01vsdari.csv", quote = FALSE)

upreg = subset(dmoj_01vsdari, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dmoj_01vsdari_UP.csv", quote = FALSE)

downreg = subset(dmoj_01vsdari, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dmoj_01vsdari_DOWN.csv", quote = FALSE)












## dari vs dmoj01 ##

darivsdmoj_01 <- results(dds, parallel=T,  contrast=c("condition","dari","dmoj01"))
ix = which.min(darivsdmoj_01$log2FoldChange)
res <- darivsdmoj_01[order((darivsdmoj_01$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c("blue","blue","blue","blue","red","red","red","red"), las=2, main=rownames(dds)[ ix  ]  )

# Save result
write.csv(darivsdmoj_01, "DESEQout_dari_vsdmoj_01.csv", quote = FALSE)

upreg = subset(darivsdmoj_01, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dari_vsdmoj_01UP.csv", quote = FALSE)

downreg = subset(darivsdmoj_01, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dari_vsdmoj_01DOWN.csv", quote = FALSE)








## dari vs dmoj01 ##

darivsdmoj_01 <- results(dds, parallel=T,  contrast=c("condition","dari","dmoj01"))
ix = which.min(darivsdmoj_01$log2FoldChange)
res <- darivsdmoj_01[order((darivsdmoj_01$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c("blue","blue","blue","blue","red","red","red","red"), las=2, main=rownames(dds)[ ix  ]  )

# Save result
write.csv(darivsdmoj_01, "DESEQout_dari_vsdmoj_01.csv", quote = FALSE)

upreg = subset(darivsdmoj_01, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dari_vsdmoj_01UP.csv", quote = FALSE)

downreg = subset(darivsdmoj_01, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dari_vsdmoj_01DOWN.csv", quote = FALSE)






## dari vs dmoj20 ##

darivsdmoj_20 <- results(dds, parallel=T,  contrast=c("condition","dari","dmoj20"))
ix = which.max(darivsdmoj_20$log2FoldChange)
res <- darivsdmoj_20[order((darivsdmoj_20$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c("blue","blue","blue","blue","red","red","red","red"), las=2, main=rownames(dds)[ ix  ]  )

# Save result
write.csv(darivsdmoj_20, "DESEQout_dari_vsdmoj_20.csv", quote = FALSE)

upreg = subset(darivsdmoj_20, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dari_vsdmoj_20UP.csv", quote = FALSE)

downreg = subset(darivsdmoj_20, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dari_vsdmoj_20DOWN.csv", quote = FALSE)







## dari vs dmoj22 ##

darivsdmoj_22 <- results(dds, parallel=T,  contrast=c("condition","dari","dmoj22"))
ix = which.max(darivsdmoj_22$log2FoldChange)
res <- darivsdmoj_22[order((darivsdmoj_22$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c("blue","blue","blue","blue","red","red","red","red"), las=2, main=rownames(dds)[ ix  ]  )

# Save result
write.csv(darivsdmoj_22, "DESEQout_dari_vsdmoj_22.csv", quote = FALSE)

upreg = subset(darivsdmoj_22, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dari_vsdmoj_22UP.csv", quote = FALSE)

downreg = subset(darivsdmoj_22, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dari_vsdmoj_22DOWN.csv", quote = FALSE)







## dari vs dmoj26 ##

darivsdmoj_26 <- results(dds, parallel=T,  contrast=c("condition","dari","dmoj26"))
ix = which.max(darivsdmoj_26$log2FoldChange)
res <- darivsdmoj_26[order((darivsdmoj_26$log2FoldChange),decreasing = TRUE),]
head(res)
barplot(assay(dds)[ix,], cex.names = 0.80, col = c("blue","blue","blue","blue","red","red","red","red"), las=2, main=rownames(dds)[ ix  ]  )

# Save result
write.csv(darivsdmoj_26, "DESEQout_dari_vsdmoj_26.csv", quote = FALSE)

upreg = subset(darivsdmoj_26, padj<0.05  & log2FoldChange > 1)
dim(upreg)
mean(upreg$log2FoldChange)
write.csv(upreg, "dari_vsdmoj_26UP.csv", quote = FALSE)

downreg = subset(darivsdmoj_26, padj<0.05  & log2FoldChange < -1)
dim(downreg)
mean(downreg$log2FoldChange)
write.csv(downreg, "dari_vsdmoj_26DOWN.csv", quote = FALSE)

