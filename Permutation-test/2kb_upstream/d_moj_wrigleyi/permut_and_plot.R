#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library("regioneR")))
suppressWarnings(suppressMessages(library("BiocManager")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("ggplot2")))

args <- commandArgs(trailingOnly = TRUE)
all_genes <- args[1] #### All 2kb upstream/downstream regions
TEs_bed <- args[2] #### All TEs in the genome
species_name <- args[3] #### species suffix

### This script must have HLAU genes in the following structure:
### dmoj01_or_2kbUP.bed
### {species_name}_{hs_name}_2kbUP.bed
### hs_name are items from the list_HS_genes

list_HS_genes <- list("obp","or", "gustatory","ionotropic", "p450", "ugt", "gst", "esterase", "abc")

permut_test <- function (gene_bedf, TEs_bedf, spc_name) {
  all_TEs <- read.table(paste0(TEs_bedf), sep = "\t", header = FALSE)
  all_TEsGR <- toGRanges(all_TEs)
  
  all_2kb_regions <- read.table(gene_bedf, sep = "\t", header = FALSE)
  head(all_2kb_regions)
  columns = c("hs_gene", "pvalue", "zscore" ,"direction") 
  df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(df) = columns
  
  for (hs_name in list_HS_genes) {
    HS_gene <- read.table(paste0(species_name, "_", hs_name, "_2kbUP.bed"), sep = "\t", header = FALSE)
    HS_geneGR <- toGRanges(HS_gene)
    head(HS_gene)
    all_minus_HSgenes <- all_2kb_regions[!all_2kb_regions$V4 %in% HS_gene$V4, ]
    head(all_minus_HSgenes)
    all_2kb_regionsGR <- toGRanges(all_minus_HSgenes)
    
    numOverlaps(HS_geneGR, all_TEsGR, count.once = TRUE)
    
    pdf(paste0("perm_", hs_name, ".pdf"))
    pt = permTest(A = HS_geneGR, ntimes = 2000, randomize.function = resampleRegions, universe = all_2kb_regionsGR,
                  evaluate.function = numOverlaps, count.once = TRUE, B = all_TEsGR, verbose = FALSE, alternative = "auto")
    plot(pt)
    dev.off()
    
    sum_res = summary(pt)
    res_pvalue = format(round(sum_res$pvalue, 4), nsmall = 4)
    zscore = format(round(sum_res$zscore, 4), nsmall = 4)
    paste("ZSCORE=", zscore)
    direction = sum_res$test
    df[nrow(df) + 1,] <- list(hs_name, res_pvalue, zscore, direction)
  }
  write.table(df, file = "perm-test_res.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
}


permut_test(all_genes, TEs_bed, species_name)



#### Make z-score plot ###

###Import dataframe and fix gene names
data <- read.csv("perm-test_res.tsv", sep = "\t")
data$hs_gene[data$hs_gene == "gustatory"] <- "GR"
data$hs_gene[data$hs_gene == "ionotropic"] <- "IR"
data$hs_gene[data$hs_gene == "esterase"] <- "EST"
data$hs_gene[data$hs_gene == "p450"] <- "CYP"
data$hs_gene <- toupper(data$hs_gene)

### Set the order of hs_gene in the opposite order
data$hs_gene <- factor(data$hs_gene, levels = rev(data$hs_gene))

### Create the plot
plot <- ggplot(data, aes(x = zscore, y = hs_gene)) +
  geom_point(aes(color = ifelse(pvalue < 0.05, "red", "grey")), size = 4, shape = 21, stroke = 1, fill = ifelse(data$pvalue < 0.05, "red", "grey"), color = "black") +
  scale_color_identity() +
  theme_minimal() +
  scale_x_continuous(limits = c(-4, 4)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  labs(x = "Z-score", y = "", title = paste0(species_name))


ggsave(paste0(species_name,"_enrichment_permut.pdf"), plot = plot, width = 3, height = 5)
