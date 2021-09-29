## Untitled 16
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


ENSM2symbol <- read.table("/home/zjw/zjw/annotation/gencode.v30.annotation.sorted.all.gene.gtf",sep = "\t")
call <- data.frame(msigdbr())
call <- call[call$gs_cat %in% c("H","C2","C3","C5"),c(3,5,13)]


gene_id <- strsplit(ENSM2symbol[,9],split = "gene_id ")
gene_id <- lapply(gene_id, function(x) x[2])
gene_id <- unlist(gene_id)
gene_id <- strsplit(gene_id,split = ";")
gene_id <- lapply(gene_id, function(x) x[1])
gene_id <- unlist(gene_id)

gene_name <- strsplit(ENSM2symbol[,9],split = "gene_name ")
gene_name <- lapply(gene_name, function(x) x[2])
gene_name <- unlist(gene_name)
gene_name <- strsplit(gene_name,split = ";")
gene_name <- lapply(gene_name, function(x) x[1])
gene_name <- unlist(gene_name)

gene_type <- strsplit(ENSM2symbol[,9],split = "gene_type ")
gene_type <- lapply(gene_type, function(x) x[2])
gene_type <- unlist(gene_type)
gene_type <- strsplit(gene_type,split = ";")
gene_type <- lapply(gene_type, function(x) x[1])
gene_type <- unlist(gene_type)



ENSM2symbol <- data.frame(gene_id,gene_name,gene_type)


ENSM2symbol <- unique(ENSM2symbol)



count <- data.frame(shCtrl_rep1=read.table("/home/zjw/zjw/20210109cornea/Scr-1_L1_339X39_Aligned.count")[,2],
                    shCtrl_rep2=read.table("/home/zjw/zjw/20210109cornea/Scr-2_L1_345X45_Aligned.count")[,2],
                    shZNF281_rep1=read.table("/home/zjw/zjw/20210109cornea/V2-1_L1_337X37_Aligned.count")[,2],
                    shZNF281_rep2=read.table("/home/zjw/zjw/20210109cornea/V2-2_L1_343X43_Aligned.count")[,2]
                    
)


rownames(count) <- read.table("/home/zjw/zjw/20210109cornea/Scr-1_L1_339X39_Aligned.count")[,1]

count <- count[!rowSums(count)==0,]







condition <- factor(c(rep("shCtrl",2),rep("shZNF281",2)))

dds <- DESeqDataSetFromMatrix(count, DataFrame(condition), design= ~ condition )

dds <- DESeq(dds)

res <- DESeq2::results(dds)

mcols(res)

summary(res)

all_gene <- as.data.frame(res)




all_gene$gene <- rownames(all_gene)
for(i in 1:nrow(all_gene)) {
  if(all_gene[i,7] %in% ENSM2symbol$gene_id) all_gene[i,7] <- ENSM2symbol[ENSM2symbol$gene_id==all_gene[i,7],]$gene_name
}




all_gene$padj[is.na(all_gene$padj)] <- 1
all_gene <- na.omit(all_gene)

all_gene[abs(all_gene$log2FoldChange)>1 & all_gene$padj<0.05,]





et_c2 <- all_gene



et_c2 <- et_c2[et_c2$gene %in% ENSM2symbol$gene_name[ENSM2symbol$gene_type %in% c("protein_coding",
                                                                                  "lincRNA",
                                                                                  "IG_C_gene",
                                                                                  "IG_D_gene",
                                                                                  "IG_J_gene",
                                                                                  "IG_LV_gene",
                                                                                  "IG_V_gene",
                                                                                  "IG_V_pseudogene",
                                                                                  "IG_J_pseudogene",
                                                                                  "IG_C_pseudogene",
                                                                                  "TR_C_gene",
                                                                                  "TR_D_gene",
                                                                                  "TR_J_gene",
                                                                                  "TR_V_gene",
                                                                                  "TR_V_pseudogene",
                                                                                  "TR_J_pseudogene")],]



glist <- et_c2$log2FoldChange
names(glist) <- as.character(et_c2$gene)
glist <- sort(glist,decreasing = T)
glist <- glist[!duplicated(names(glist))]


glist <- glist[names(glist) %in% et_c2$gene[ et_c2$baseMean>6]]   

# glist <- glist[names(glist) %in% regulons$KDM5B_extended]
# 
# gsea <- GSEA(glist, TERM2GENE = data.frame(gs_name="SASP",gene_symbol=as.character(sasp$Gene),gs_description="Senescence-associated secretory phenotype")[,c(1,2)],TERM2NAME = data.frame(gs_name="SASP",gene_symbol=as.character(sasp$Gene),gs_description="Senescence-associated secretory phenotype")[,c(1,3)],pvalueCutoff = 1,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0,exponent = 1)
# gseaplot2(gsea,geneSetID = "SASP",pvalue_table = T)
# gsea.list <- data.frame(gsea)



gsea <- GSEA(glist, TERM2GENE = call[,c(1,2)],
             #TERM2NAME = call_new[,c(1,3)],
             pvalueCutoff = 1,pAdjustMethod = "BH",maxGSSize = 10000,minGSSize = 0,exponent = 1,nPerm = 5000)

gsea.list <- data.frame(gsea)



gseaplot2(gsea,geneSetID = "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",pvalue_table = F,title = "Glycolysis / Gluconeogenesis",subplots = c(1,2)) +
  gseaplot2(gsea,geneSetID = "KEGG_OXIDATIVE_PHOSPHORYLATION",pvalue_table = F,title = "Oxidative phosphorylation",subplots = c(1,2)) +
  gseaplot2(gsea,geneSetID = "GO_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION",pvalue_table = F,title = "Positive regulation of epithelial cell proliferation",subplots = c(1,2)) +
  gseaplot2(gsea,geneSetID = "GO_TISSUE_MIGRATION",pvalue_table = F,title = "Tissue migration",subplots = c(1,2)) 
