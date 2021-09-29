## Dependencies

library(Seurat)
library(cowplot)
library(dplyr)
library(monocle)
library(patchwork)
library(ggplot2)
library(edgeR)
library(ggrepel)
library(ggpubr)
library(enrichplot)
library(msigdbr)
library(pheatmap)
library(stringr)
library(reshape2)
library(DESeq2)
library(clusterProfiler)
library(org.Mmu.eg.db)
library(org.Hs.eg.db)
library(xlsx)
library(RColorBrewer)
library(BuenColors)
library(systemPipeR)





## load the fire

combined_cca_40_300 <- readRDS("20210928cornea1.RDS")
combined_cca_40_300 <- subset(combined_cca_40_300,subset = var != "Old")
combined_cca_40_300 <- subset(combined_cca_40_300,subset = var != "DM")





## Fig. 1a

DimPlot(combined_cca_40_300,label = T) +
  scale_color_manual(values=jdb_color_map(c("B","GMP","MPP","mono","CD4","Ery","CD8","GMP-B"))) 





## Fig. 1b

DefaultAssay(combined_cca_40_300) <- "RNA"
combined_cca_40_300 <- ScaleData(combined_cca_40_300)
combined_cca_40_300 <- NormalizeData(combined_cca_40_300,assay = "RNA")
mk <- FindAllMarkers(combined_cca_40_300,assay = "RNA",only.pos = T)
ave_exp <- AverageExpression(combined_cca_40_300,assays = "RNA")[[1]]

gene_list <- c()
for (i in unique(mk$cluster)) {
  temp <- mk[mk$cluster==i,]
  temp <- temp[order(temp$avg_logFC,decreasing = T),]
  gene_list <- c(gene_list,temp[1:30,]$gene)
}

df <- data.frame()
for (i in gene_list) {
  df <- rbind(df,
              ave_exp[rownames(ave_exp)==i,])
}
pheatmap(df,
         scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         #col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         show_rownames = T,
         treeheight_row = F,
         treeheight_col = 100,
         border = F) 





## Fig. 1d

immune_gene_list <- read.table("immuneGeneList.txt",sep = "\t",header = 1)

gene_list <- c()
for (i in unique(mk$cluster)) {
  temp <- mk[mk$cluster==i,]
  temp <- temp[order(temp$avg_logFC,decreasing = T),]
  temp <- temp[temp$gene %in% immune_gene_list$Symbol,]
  gene_list <- c(gene_list,temp$gene)
}


df <- data.frame()
for (i in gene_list) {
  df <- rbind(df,
              ave_exp[rownames(ave_exp)==i,])
}

pheatmap(df[,7:1],
         scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         #col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         show_rownames = F,
         treeheight_row = F,
         treeheight_col = 100,
         fontsize_col = 13,
         border = F) 





## Fig. 1f
