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
combined_cca_40_300 <- subset(combined_cca_40_300,subset = var != "Young")
call <- data.frame(msigdbr())
call <- call[call$gs_cat %in% c("H","C2","C3","C5"),c(3,5,13)]





## Fig. 5a
DimPlot(combined_cca_40_300,group.by = "var",split.by = "var")





## Fig. 5b
DefaultAssay(combined_cca_40_300) <- "RNA"
combined_cca_40_300 <- ScaleData(combined_cca_40_300)
combined_cca_40_300 <- NormalizeData(combined_cca_40_300,assay = "RNA")

a <- combined_cca_40_300@meta.data
a <- plyr::count(a[,c(5,9)])
a$celltype <- factor(a$celltype,levels = c("LSC/LPC","TAC","PMC","TDC","ConjEC-1","ConjEC-2","MC"))
ggplot(a,aes(x = celltype, y = freq,fill = var)) + 
  geom_bar(position = "fill",stat = "identity",width = 0.8,color="black") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=45,hjust = 1,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x="",y="Proportion")
  
  
  
  
  
## Fig. 5c

combined_cca_40_300 <- CellCycleScoring(combined_cca_40_300,s.features = cc.genes.updated.2019$s.genes,g2m.features = cc.genes.updated.2019$g2m.genes)
combined_cca_40_300@meta.data$Phase <- factor(combined_cca_40_300@meta.data$Phase,levels = c("G1","S","G2M"))

df <- combined_cca_40_300@meta.data
df <- plyr::count(df[,c(5,9,12)])
df$celltype <- factor(df$celltype,levels = c("LSC/LPC","TAC","PMC","TDC","ConjEC-1","ConjEC-2","MC"))
ggplot(df,aes(x = var, y = freq,fill = Phase)) + 
  geom_bar(position = "fill",stat = "identity",width = 0.8,color="black") +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=45,hjust = 1,vjust = 1,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x="",y="Proportion") +
  facet_wrap(~celltype,nrow = 1) +
  scale_fill_manual(values=c("#4169E1","#DC143C","#48D1CC"))





## Fig. 5d

combined_cca_40_300.gd <- combined_cca_40_300
combined_cca_40_300.gd$celltype.stim <- paste(Idents(combined_cca_40_300.gd), combined_cca_40_300.gd$var, sep = "_")
combined_cca_40_300.gd$celltype <- Idents(combined_cca_40_300.gd)
Idents(combined_cca_40_300.gd) <- "celltype.stim"
DefaultAssay(combined_cca_40_300.gd) <- "RNA"
combined_cca_40_300.gd <- ScaleData(combined_cca_40_300.gd,assay = "RNA")
combined_cca_40_300.gd <- NormalizeData(combined_cca_40_300.gd,assay = "RNA")


LSCLPC.diff <- FindMarkers(combined_cca_40_300.gd, ident.1 = "LSC/LPC_DM", ident.2 = "LSC/LPC_Old", verbose = FALSE)
TAC.diff <- FindMarkers(combined_cca_40_300.gd, ident.1 = "TAC_DM", ident.2 = "TAC_Old", verbose = FALSE)
PMC.diff <- FindMarkers(combined_cca_40_300.gd, ident.1 = "PMC_DM", ident.2 = "PMC_Old", verbose = FALSE)
TDC.diff <- FindMarkers(combined_cca_40_300.gd, ident.1 = "TDC_DM", ident.2 = "TDC_Old", verbose = FALSE)





## Fig. 5d & Fig. 5e

LSCLPC.diff.all <- FindMarkers(combined_cca_40_300.gd, ident.1 = "LSC/LPC_DM", ident.2 = "LSC/LPC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)
TAC.diff.all <- FindMarkers(combined_cca_40_300.gd, ident.1 = "TAC_DM", ident.2 = "TAC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)
PMC.diff.all <- FindMarkers(combined_cca_40_300.gd, ident.1 = "PMC_DM", ident.2 = "PMC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)
TDC.diff.all <- FindMarkers(combined_cca_40_300.gd, ident.1 = "TDC_DM", ident.2 = "TDC_CTRL", verbose = FALSE, logfc.threshold = 0, min.pct = 0)

# LSC/LPC KEGG
trans <- bitr(rownames(LSCLPC.diff.all),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
trans <- merge(trans,data.frame(SYMBOL=rownames(LSCLPC.diff.all),
                                avg_logFC=LSCLPC.diff.all$avg_logFC))
glist <- trans$avg_logFC
names(glist) <- trans$ENTREZID
glist <- sort(glist,decreasing = T)
gsea.kegg.LSCLPC <- gseKEGG(glist,organism = "hsa",pvalueCutoff = 1,pAdjustMethod = "BH") 

# TAC KEGG
trans <- bitr(rownames(TAC.diff.all),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
trans <- merge(trans,data.frame(SYMBOL=rownames(TAC.diff.all),
                                avg_logFC=TAC.diff.all$avg_logFC))
glist <- trans$avg_logFC
names(glist) <- trans$ENTREZID
glist <- sort(glist,decreasing = T)
gsea.kegg.TAC <- gseKEGG(glist,organism = "hsa",pvalueCutoff = 1,pAdjustMethod = "BH") 

# PMC KEGG
trans <- bitr(rownames(PMC.diff.all),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
trans <- merge(trans,data.frame(SYMBOL=rownames(PMC.diff.all),
                                avg_logFC=PMC.diff.all$avg_logFC))
glist <- trans$avg_logFC
names(glist) <- trans$ENTREZID
glist <- sort(glist,decreasing = T)
gsea.kegg.PMC <- gseKEGG(glist,organism = "hsa",pvalueCutoff = 1,pAdjustMethod = "BH") 

# TDC KEGG
trans <- bitr(rownames(TDC.diff.all),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
trans <- merge(trans,data.frame(SYMBOL=rownames(TDC.diff.all),
                                avg_logFC=TDC.diff.all$avg_logFC))
glist <- trans$avg_logFC
names(glist) <- trans$ENTREZID
glist <- sort(glist,decreasing = T)
gsea.kegg.TDC <- gseKEGG(glist,organism = "hsa",pvalueCutoff = 1,pAdjustMethod = "BH") 

a <- as.data.frame(table(c(data.frame(gsea.kegg.LSCLPC)[data.frame(gsea.kegg.LSCLPC)$pvalue<=0.05 & data.frame(gsea.kegg.LSCLPC)$NES>0,]$Description,
                           data.frame(gsea.kegg.TAC)[data.frame(gsea.kegg.TAC)$pvalue<=0.05 & data.frame(gsea.kegg.TAC)$NES>0,]$Description,
                           data.frame(gsea.kegg.PMC)[data.frame(gsea.kegg.PMC)$pvalue<=0.05 & data.frame(gsea.kegg.PMC)$NES>0,]$Description,
                           data.frame(gsea.kegg.TDC)[data.frame(gsea.kegg.TDC)$pvalue<=0.05 & data.frame(gsea.kegg.TDC)$NES>0,]$Description
)))

vennPlot(overLapper(list("LSCLPC upregulated pathway"=data.frame(gsea.kegg.LSCLPC)[data.frame(gsea.kegg.LSCLPC)$pvalue<=0.05 & data.frame(gsea.kegg.LSCLPC)$NES>0,]$Description,
                         "TAC upregulated pathway"=data.frame(gsea.kegg.TAC)[data.frame(gsea.kegg.TAC)$pvalue<=0.05 & data.frame(gsea.kegg.TAC)$NES>0,]$Description,
                         "PMC upregulated pathway"=data.frame(gsea.kegg.PMC)[data.frame(gsea.kegg.PMC)$pvalue<=0.05 & data.frame(gsea.kegg.PMC)$NES>0,]$Description,
                         "TDC upregulated pathway"=data.frame(gsea.kegg.TDC)[data.frame(gsea.kegg.TDC)$pvalue<=0.05 & data.frame(gsea.kegg.TDC)$NES>0,]$Description), type="vennsets"
)
)

a <- as.data.frame(table(c(data.frame(gsea.kegg.LSCLPC)[data.frame(gsea.kegg.LSCLPC)$pvalue<=0.05 & data.frame(gsea.kegg.LSCLPC)$NES<0,]$Description,
                           data.frame(gsea.kegg.TAC)[data.frame(gsea.kegg.TAC)$pvalue<=0.05 & data.frame(gsea.kegg.TAC)$NES<0,]$Description,
                           data.frame(gsea.kegg.PMC)[data.frame(gsea.kegg.PMC)$pvalue<=0.05 & data.frame(gsea.kegg.PMC)$NES<0,]$Description,
                           data.frame(gsea.kegg.TDC)[data.frame(gsea.kegg.TDC)$pvalue<=0.05 & data.frame(gsea.kegg.TDC)$NES<0,]$Description
)))

vennPlot(overLapper(list("LSCLPC downregulated pathway"=data.frame(gsea.kegg.LSCLPC)[data.frame(gsea.kegg.LSCLPC)$pvalue<=0.05 & data.frame(gsea.kegg.LSCLPC)$NES<0,]$Description,
                         "TAC downregulated pathway"=data.frame(gsea.kegg.TAC)[data.frame(gsea.kegg.TAC)$pvalue<=0.05 & data.frame(gsea.kegg.TAC)$NES<0,]$Description,
                         "PMC downregulated pathway"=data.frame(gsea.kegg.PMC)[data.frame(gsea.kegg.PMC)$pvalue<=0.05 & data.frame(gsea.kegg.PMC)$NES<0,]$Description,
                         "TDC downregulated pathway"=data.frame(gsea.kegg.TDC)[data.frame(gsea.kegg.TDC)$pvalue<=0.05 & data.frame(gsea.kegg.TDC)$NES<0,]$Description), type="vennsets"
)
)





## Fig. 5f & Fig. 5g & Fig. 5h

combined_cca_40_300.gd <- AddModuleScore(combined_cca_40_300.gd, features = list(call[call$gs_name=="HALLMARK_OXIDATIVE_PHOSPHORYLATION",]$gene_symbol),name="OXIDATIVE_PHOSPHORYLATION",assay = "RNA") 
combined_cca_40_300.gd <- AddModuleScore(combined_cca_40_300.gd, features = list(call[call$gs_name=="HALLMARK_GLYCOLYSIS",]$gene_symbol),name="GLYCOLYSIS",assay = "RNA") 
combined_cca_40_300.gd <- AddModuleScore(combined_cca_40_300.gd, features = list(call[call$gs_name=="KEGG_CITRATE_CYCLE_TCA_CYCLE",]$gene_symbol),name="TCA",assay = "RNA") 

df <- combined_cca_40_300.gd@meta.data

ggboxplot(df, "celltype", "GLYCOLYSIS1", fill = "var" ) +
  stat_compare_means(aes(group=var), label = "..p.format..",method = "t.test")+
  labs(y="Glycolysis gene set score",x="")

ggboxplot(df, "celltype", "OXIDATIVE_PHOSPHORYLATION1", fill = "var" ) +
  stat_compare_means(aes(group=var), label = "..p.format..",method = "t.test")+
  labs(y="OXPHOS gene set score",x="")

ggboxplot(df, "celltype", "TCA1", fill = "var" ) +
  stat_compare_means(aes(group=var), label = "..p.format..",method = "t.test")+
  labs(y="TCA cycle gene set score",x="")





## Supplementary Fig. 5 a

hm_up <- rbind(data.frame(gene=LSCLPC.diff[LSCLPC.diff$p_val_adj<=0.05 & LSCLPC.diff$avg_logFC>=0.3,] %>% rownames(),celltype="LSCLPC"),
               data.frame(gene=TAC.diff[TAC.diff$p_val_adj<=0.05 & TAC.diff$avg_logFC>=0.3,] %>% rownames(),celltype="TAC"),
               data.frame(gene=PMC.diff[PMC.diff$p_val_adj<=0.05 & PMC.diff$avg_logFC>=0.3,] %>% rownames(),celltype="PMC"),
               data.frame(gene=TDC.diff[TDC.diff$p_val_adj<=0.05 & TDC.diff$avg_logFC>=0.3,] %>% rownames(),celltype="TDC")
)
hm_up_new <- data.frame(gene=unique(hm_up$gene))
hm_up_new[,2:5] <- 0
hm_up_new <- hm_up_new[,-1]
rownames(hm_up_new) <- unique(hm_up$gene)
colnames(hm_up_new) <- c("LSCLPC","TAC","PMC","TDC")
for (i in 1:nrow(hm_up)) {
  hm_up_new[rownames(hm_up_new)==hm_up$gene[i],colnames(hm_up_new)==hm_up$celltype[i]] <- 1
}
hm_up_new <-rbind(hm_up_new[rowSums(hm_up_new)==4,],
                  hm_up_new[rowSums(hm_up_new)==3,],
                  hm_up_new[rowSums(hm_up_new)==2,],
                  hm_up_new[rowSums(hm_up_new)==1,])
hm_up_new <- hm_up_new[,1:4]
pheatmap(hm_up_new,
         color = colorRampPalette(c("grey", "red"))(5),cluster_cols = F,cluster_rows = F,show_rownames = F)

hm_dn <- rbind(data.frame(gene=LSCLPC.diff[LSCLPC.diff$p_val_adj<=0.05 & LSCLPC.diff$avg_logFC<= -0.3,] %>% rownames(),celltype="LSCLPC"),
               data.frame(gene=TAC.diff[TAC.diff$p_val_adj<=0.05 & TAC.diff$avg_logFC<= -0.3,] %>% rownames(),celltype="TAC"),
               data.frame(gene=PMC.diff[PMC.diff$p_val_adj<=0.05 & PMC.diff$avg_logFC<= -0.3,] %>% rownames(),celltype="PMC"),
               data.frame(gene=TDC.diff[TDC.diff$p_val_adj<=0.05 & TDC.diff$avg_logFC<= -0.3,] %>% rownames(),celltype="TDC")
)
hm_dn_new <- data.frame(gene=unique(hm_dn$gene))
hm_dn_new[,2:5] <- 0
hm_dn_new <- hm_dn_new[,-1]
rownames(hm_dn_new) <- unique(hm_dn$gene)
colnames(hm_dn_new) <- c("LSCLPC","TAC","PMC","TDC")
for (i in 1:nrow(hm_dn)) {
  hm_dn_new[rownames(hm_dn_new)==hm_dn$gene[i],colnames(hm_dn_new)==hm_dn$celltype[i]] <- 1
}
hm_dn_new <-rbind(hm_dn_new[rowSums(hm_dn_new)==4,],
                  hm_dn_new[rowSums(hm_dn_new)==3,],
                  hm_dn_new[rowSums(hm_dn_new)==2,],
                  hm_dn_new[rowSums(hm_dn_new)==1,])
hm_dn_new <- hm_dn_new[,1:4]
pheatmap(hm_dn_new,
         color = colorRampPalette(c("grey", "blue"))(5),cluster_cols = F,cluster_rows = F,show_rownames = F)





## Supplementary Fig. 5 b

go_LSCLPC_up_BP <- enrichGO(LSCLPC.diff[LSCLPC.diff$p_val_adj<=0.05 & LSCLPC.diff$avg_logFC>=0.3,] %>% rownames(), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                            pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
go_TAC_up_BP <- enrichGO(TAC.diff[TAC.diff$p_val_adj<=0.05 & TAC.diff$avg_logFC>=0.3,] %>% rownames(), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
go_PMC_up_BP <- enrichGO(PMC.diff[PMC.diff$p_val_adj<=0.05 & PMC.diff$avg_logFC>=0.3,] %>% rownames(), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
go_TDC_up_BP <- enrichGO(TDC.diff[TDC.diff$p_val_adj<=0.05 & TDC.diff$avg_logFC>=0.3,] %>% rownames(), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)

df_up <- data.frame()
for (i in ls() %>% grep(pattern = "^go",value = T) %>% grep(pattern = "up",value = T)) {
  a <- paste0('temp <- data.frame(',i,')')
  eval(parse(text=a))
  print(a)
  df_up <- rbind(df_up,
                 data.frame(temp,celltype=(i %>% strsplit("_") %>% unlist())[2]))
}
df_up$celltype <- factor(df_up$celltype,levels = c("LSCLPC","TAC","PMC","TDC"))

df_up_new <- df_up[df_up$Description %in% c("epithelial cell apoptotic process","mRNA catabolic process","cellular response to oxidative stress","response to oxidative stress","response to reactive oxygen species"),]
df_up_new$bubblesize <- df_up_new$Count/max(df_up_new$Count)
df_up_new$log <- -log10(df_up_new$p.adjust)

ggplot(data=df_up_new, mapping=aes(x=celltype,y=Description,color=log))+
  geom_point(stat= "identity",aes(size=Count),alpha=0.7,show.legend = TRUE)+ 
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=90,hjust = 1,vjust = 0.5,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x=element_blank(),y=element_blank(),col="-log10 adjusted p-value") +
  scale_color_gradientn(colors = c(colorRampPalette(as.character(jdb_palette("brewer_heat")))(250))[126:250])





## Supplementary Fig. 5 c

go_LSCLPC_dn_BP <- enrichGO(LSCLPC.diff[LSCLPC.diff$p_val_adj<=0.05 & LSCLPC.diff$avg_logFC<=-0.3,] %>% rownames(), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                            pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
go_TAC_dn_BP <- enrichGO(TAC.diff[TAC.diff$p_val_adj<=0.05 & TAC.diff$avg_logFC<=-0.3,] %>% rownames(), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
go_PMC_dn_BP <- enrichGO(PMC.diff[PMC.diff$p_val_adj<=0.05 & PMC.diff$avg_logFC<=-0.3,] %>% rownames(), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
go_TDC_dn_BP <- enrichGO(TDC.diff[TDC.diff$p_val_adj<=0.05 & TDC.diff$avg_logFC<=-0.3,] %>% rownames(), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)

df_dn <- data.frame()
for (i in ls() %>% grep(pattern = "^go",value = T) %>% grep(pattern = "dn",value = T)) {
  a <- paste0('temp <- data.frame(',i,')')
  eval(parse(text=a))
  print(a)
  df_dn <- rbind(df_dn,
                 data.frame(temp,celltype=(i %>% strsplit("_") %>% unlist())[2]))
}
df_dn$celltype <- factor(df_dn$celltype,levels = c("LSCLPC","TAC","PMC","TDC","Conjunctival","Melanocyte"))

df_dn_new <- df_dn[df_dn$Description %in% c("cornification","epidermal cell differentiation","keratinocyte differentiation","positive regulation of response to wounding","positive regulation of wound healing"),]
df_dn_new$bubblesize <- df_dn_new$Count/max(df_dn_new$Count)
df_dn_new$log <- -log10(df_dn_new$p.adjust)

ggplot(data=df_dn_new, mapping=aes(x=celltype,y=Description,color=log))+
  geom_point(stat= "identity",aes(size=Count),alpha=0.7,show.legend = TRUE)+ 
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme_set(theme_bw()) +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
  theme(panel.border = element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(face = "plain",size = 13,angle=90,hjust = 1,vjust = 0.5,colour = "black"),
        axis.text.y = element_text(face = "plain",size = 12,colour = "black")) +
  labs(x=element_blank(),y=element_blank(),col="-log10 adjusted p-value") +
  scale_color_gradientn(colors = c(colorRampPalette(as.character(jdb_palette("brewer_marine")))(250))[126:250])
