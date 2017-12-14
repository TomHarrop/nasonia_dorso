#!/usr/bin/env Rscript

library(data.table)
library(VennDiagram)
library(compare)

#read file (from previous analysis)
#from Shannon's computer
res_dpp <- read.csv('R/res_dpp.csv', header=TRUE,sep = ',')
res_TollA <- read.csv('R/res_TollA.csv', header=TRUE,sep = ',')
res_pers <- read.csv('R/pers_data.csv', header=TRUE, sep = ',')
with_names <- read.csv('R/with_nams.csv', header = TRUE, sep = ',')
alpha <- 0.05
threshold <- 1

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
TollA_sig <- subset(res_TollA, padj < alpha)
TollA_names <- as.character(TollA_sig$rn)

dpp_sig <- subset(res_dpp, padj < alpha)
dpp_names <- as.character(dpp_sig$rn)

#pull significant genes from pers dataset
pers_sig <- subset(res_pers, significant == "yes")
pers_names <- as.character(pers_sig$gene)

res_TollA[order(res_TollA$padj),]

TollA_and_pers <- TollA_names[TollA_names %in% pers_names]
dpp_and_pers <- dpp_names[dpp_names %in% pers_names]

TollA_genes <- with_names[with_names$ï..gene %in% TollA_and_pers, "Best.ortholog.description"]
dpp_genes <- with_names[with_names$ï..gene %in% dpp_and_pers, "Best.ortholog.description"]
