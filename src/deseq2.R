#!/usr/bin/env Rscript

library(data.table)
library(tximport)
library(DESeq2)
library(VennDiagram)

#get list of quant files
file_list <-list.files("output/salmon", 
                       full.names=TRUE, 
                       recursive = TRUE, 
                       pattern = "quant.sf")
names(file_list) <- gsub(".*/", "", dirname(file_list))


#generate txtogene
quant_results <- fread(file_list[[1]])
txtogene <- quant_results[,.(
    TXNAME = Name, 
    GENEID = gsub("t[[:digit:]]+$", "", Name))]

#import quant files
txi <- tximport(file_list, type = "salmon", tx2gene = data.frame(txtogene))

#filter lowly expressed genes
mean_filter <- rowMeans(txi$counts) > 5
max_filter <- rowMax(txi$counts) > 10


#generate coldata
col_dt<-fread("data/sample_names.csv")
col_dt[grep("^TollA", sample), trt:="TollA"]
col_dt[grep("^h2o", sample), trt:="h2o"]
col_dt[grepl("^dpp", sample) | grepl("^gbb", sample), trt:="dpp"]
col_dt[, record_id:= factor(record_id, levels=colnames(txi$counts))]
setorder(col_dt, record_id)
col_data <- data.frame(col_dt, row.names = "record_id")

#generate deseq
dds<- DESeqDataSetFromTximport(txi, col_data, ~trt)
dds <- dds[max_filter | mean_filter]
dds <-DESeq(dds)

#extract results
alpha <- 0.05
threshold <- 1
res_TollA <- results(dds, 
                     contrast = c("trt", "TollA", "h2o"), 
                     lfcThreshold = threshold, alpha= alpha)

res_dpp <- results(dds, 
                     contrast = c("trt", "dpp", "h2o"), 
                     lfcThreshold = threshold, alpha= alpha)

# write results
fwrite(data.table(data.frame(res_TollA), keep.rownames = TRUE),
       "output/deseq2/res_tolla.csv")
fwrite(data.table(data.frame(res_dpp), keep.rownames = TRUE),
       "output/deseq2/res_dpp.csv")

# draw venn diagram
sig_dpp <- rownames(subset(res_dpp, padj < alpha))
sig_tolla <- rownames(subset(res_TollA, padj < alpha))
Set1 <- RColorBrewer::brewer.pal(2, "Set1")
vd <- venn.diagram(x = list(
    "DPP" = sig_dpp,
    "TollA" = sig_tolla),
    filename = NULL,
    fill = Set1[c(1:2)],
    lty = "solid",
    lwd = 1,
    cex = 1,
    cat.cex = 2,
    fontfamily = 'Sans',
    cat.fontfamily = 'Sans',
    alpha = 0.5,
    margin = 0
)
grid.newpage()
grid.draw(vd)

# ZOMG which gene
sig_dpp[sig_dpp %in% sig_tolla]

# how many DE genes?
subset(res_TollA, padj < alpha)
subset(res_dpp, padj < alpha)

res_TollA[order(res_TollA$padj),]
plotCounts(dds, "Nasvi2EG004500", intgroup = "trt")
plotCounts(dds, "Nasvi2EG001135", intgroup = "trt")

plotCounts(dds, "Nasvi2EG006445", intgroup = "trt")
