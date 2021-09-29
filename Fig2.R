## Untitled 50
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




## load the file

combined_cca_40_300 <- readRDS("20210928cornea1.RDS")
combined_cca_40_300 <- subset(combined_cca_40_300,subset = var != "DM")
combined_cca_40_300 <- subset(combined_cca_40_300,subset = var != "Old")
combined_cca_40_300 <- subset(combined_cca_40_300, idents = c("LSC/LPC","TAC","PMC", "TDC"))







load("/home/xyh/LEC_analysis/20210531_monocle/all_cell_atlas/all__remove_cell_cycle.monocle_by_seurat.RData")
pdata <- pData(monocle_obj)

bc <- pdata[pdata$var == "YOUNG",] %>% rownames()
monocle_subset <- monocle_obj[,bc]
plot_cell_trajectory(monocle_subset, color_by = "var")
diff_test_res <- differentialGeneTest(monocle_subset[c("KRT12", "MKI67", "KRT15", "KRT3"),],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")





p <- plot_pseudotime_heatmap(monocle_subset[combined_cca_40_300@assays$integrated@var.features,],
                             num_clusters = 3,
                             cores = 1,
                             return_heatmap=T,
                             trend_formula = "~sm.ns(Pseudotime, df=2)",
                             show_rownames = F)


p$tree_row


clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering[,2] <- rownames(clustering)





go <- enrichGO(clustering$V2[clustering$Gene_Clusters=="1"],pvalueCutoff = 0.05,pAdjustMethod = "BH",OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
enrichplot::dotplot(go,showCategory=5) +
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        )

go <- enrichGO(clustering$V2[clustering$Gene_Clusters=="2"],pvalueCutoff = 0.05,pAdjustMethod = "BH",OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
enrichplot::dotplot(go,showCategory=5) +
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        )

go <- enrichGO(clustering$V2[clustering$Gene_Clusters=="3"],pvalueCutoff = 0.05,pAdjustMethod = "BH",OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
enrichplot::dotplot(go,showCategory=5) +
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        )


