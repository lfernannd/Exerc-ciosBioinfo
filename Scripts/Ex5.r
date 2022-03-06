library("affy")
library("GEOquery")

# getting file
# setwd("Files/")
# getGEOSuppFiles("GSE45276")
# untar("GSE45276_RAW.tar", exdir = "data")
# cels = list.files("data/", pattern = "CEL")
# sapply(paste("data", cels, sep = "/"), gunzip)

cels = list.files("data/", pattern = "CEL")

setwd("GSE45276//data//")

full_data = ReadAffy(filenames=cels)
full_data_rma_norm = rma(full_data)
data_rma = exprs(full_data_rma_norm)

write.table(data_rma, file = "data_rma.txt", sep = "\t")

library("arrayQualityMetrics")
arrayQualityMetrics(expressionset = full_data_rma_norm,
                    outdir = "Report_for_data_rma",
                    force = TRUE)
# no outliers detected

table <- read.table(file = "data_rma.txt", sep = "\t", header = TRUE)
means <- rowMeans(table)
var <- apply(table,1, var)
vm_matrix <- cbind(table, means, var)

selected_index <- which((vm_matrix$means > 6) & (vm_matrix$var > 0.2))
selected_probes <- vm_matrix[selected_index,] 

n_probes <- nrow(selected_probes)
# 4128 probes