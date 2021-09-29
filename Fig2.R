### LECC食蟹猴LSC分化 调控




options(stringsAsFactors = FALSE)
options(scipen = 100)


library(monocle)
library(SeuratWrappers)

load("/home/xyh/LEC_analysis/20200611_rep/20200209/LEC_ncbi_300_40.RData")
combined_cca_40_300 <- combined_cca_40


combined_cca_40_300 <-  RunUMAP(object = combined_cca_40_300, reduction = "pca", dims = 1:35)
combined_cca_40_300 <-  FindNeighbors(object = combined_cca_40_300, reduction = "pca", dims = 1:35)
combined_cca_40_300 <-  FindClusters(combined_cca_40_300, resolution = 0.1)

DimPlot(combined_cca_40_300,label = T)






combined_cca_40_300 <- RenameIdents(combined_cca_40_300, `0` = "CEC", `1` = " LEC", `2` = "CTAC", 
                                    `3` = "LSC", `4` = "CTAC", `5` = "ETAC",`6` = "Conjunctival",`7` = "Melanocyte", `8` = "IC")

combined_cca_40_300$celltype <- Idents(combined_cca_40_300)

combined_cca_40_300 <- subset(combined_cca_40_300, idents = c("LSC","ETAC","CTAC", "CEC"))

combined_cca_40_300$var[combined_cca_40_300$var=="CTRL"] <- "AGED"
combined_cca_40_300$var[combined_cca_40_300$var=="STIM"] <- "DM"

combined_cca_40_300.YOUNG <- subset(combined_cca_40_300,subset = var != "DM")
combined_cca_40_300.YOUNG <- subset(combined_cca_40_300.YOUNG,subset = var != "AGED")

combined_cca_40_300.AGED <- subset(combined_cca_40_300,subset = var != "DM")
combined_cca_40_300.AGED <- subset(combined_cca_40_300.AGED,subset = var != "YOUNG")

combined_cca_40_300.DM <- subset(combined_cca_40_300,subset = var != "YOUNG")
combined_cca_40_300.DM <- subset(combined_cca_40_300.DM,subset = var != "AGED")











# load("/home/xyh/LEC_analysis/20210531_monocle/CEC.test1.monocle_by_seurat.RData")
# load("/home/xyh/LEC_analysis/20210531_monocle/Dm/cellcycle.test1.monocle_by_seurat.RData")


load("/home/xyh/LEC_analysis/20210531_monocle/all_cell_atlas/all__remove_cell_cycle.monocle_by_seurat.RData")
pdata <- pData(monocle_obj)

bc <- pdata[pdata$var == "YOUNG",] %>% rownames()
monocle_subset <- monocle_obj[,bc]
plot_cell_trajectory(monocle_subset, color_by = "var")
diff_test_res <- differentialGeneTest(monocle_subset[c("KRT12", "MKI67", "KRT15", "KRT3"),],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")





p <- plot_pseudotime_heatmap(monocle_subset[combined_cca_40@assays$integrated@var.features,],
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



# enrich <- enricher(clustering$V2[clustering$Gene_Clusters=="3"], TERM2GENE = metabolism_df[,c(1,2)],TERM2NAME = metabolism_df[,c(1,3)],pvalueCutoff = 1,pAdjustMethod = "none",maxGSSize = 10000,minGSSize = 1)
# a <- data.frame(enrich)




go <- enrichGO(clustering$V2[clustering$Gene_Clusters=="3"],pvalueCutoff = 0.05,pAdjustMethod = "BH",OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
a <- data.frame(go)


enrichplot::dotplot(go,showCategory=5) +
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  
  
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.text.x = element_text(face = "plain",size = 10.5,angle=0,hjust = 0.5,vjust = 0.5,color="black"),
        axis.text.y = element_text(face = "plain",size = 10.5,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        axis.title.y = element_text(size = 12,color="black"),
        )



a
