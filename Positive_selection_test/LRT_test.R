args <- commandArgs(TRUE)
fname <- args[1]
#paste0("output_",fname,)
lnl_data = read.csv(paste0("output_",fname,"/results_selection.tsv"), sep ="\t", header = TRUE)
options(digits=9)

LRT_test <- function (output_name) {
  for (row in 1:nrow(lnl_data)) {
    gene_id <- as.character(lnl_data[row, 1])
    lnl_null <- lnl_data[row,5]
    lnl_alt <- lnl_data[row,10]
    subtraction_lnl = print(lnl_alt - lnl_null, 6)
    delta <- 2*(subtraction_lnl)
    {
      pvalue = pchisq(delta, df = 1, lower.tail=F)
      res <- data.frame(gene_id, pvalue)
      write.table(res, file = paste0(output_name, "_result_LRT.tsv"), sep = "\t",
                  append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }
  }
}

LRT_test(fname)
