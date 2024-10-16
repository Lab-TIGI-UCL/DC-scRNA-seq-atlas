########### R script for Figure3. ###########
########### Created by Danwen Qian on 01-05-2024. ###########
########### Last modified by Danwen Qian on 01-05-2024. ###########

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(ggthemes)
library(ggbeeswarm)
library(scales)
library(forcats)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(ggsci)
library(rcna)
library(glue)

###Fifure3A
##MANA score analysis
combined_MANA<-readRDS("combined_MANA.rds")
skin<-subset(combined_MANA, cancer_type %in% c("BCC","SKCM","SCC"))
non_skin<-subset(combined_MANA, cancer_type %in% c("BRCA","LC","RCC"))

plot<-ggplot(data=non_skin,aes(y=CXCL13_MANA_score,x=response,color=response))+theme_classic()+geom_boxplot(outlier.shape = NA,width=0.9) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+scale_y_continuous(trans=pseudo_log_trans(base = 10))
ggsave("MANA/CXCL13_MANA_non_skin.pdf", plot, height=3, width=3, dpi=600)

plot<-ggplot(data=non_skin,aes(y=CXCL13_MANA_score,x=response,color=response))+theme_classic()+geom_boxplot(outlier.shape = NA,width=0.9) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+scale_y_continuous(trans=pseudo_log_trans(base = 10))+facet_wrap(~cancer_type,nrow = 1,scales = "free_y")
ggsave("MANA/CXCL13_MANA_non_skin_cancer_type.pdf", plot, height=3, width=7, dpi=600)

plot<-ggplot(data=skin,aes(y=TCF7_MANA_score,x=response,color=response))+theme_classic()+geom_boxplot(outlier.shape = NA,width=0.91) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+scale_y_continuous(trans=pseudo_log_trans(base = 10))
ggsave("MANA/TCF7_MANA_skin.pdf", plot, height=3, width=3, dpi=600)

plot<-ggplot(data=skin,aes(y=TCF7_MANA_score,x=response,color=response))+theme_classic()+geom_boxplot(outlier.shape = NA,width=0.9) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+scale_y_continuous(trans=pseudo_log_trans(base = 10))+facet_wrap(~cancer_type,nrow = 1,scales = "free_y")
ggsave("MANA/TCF7_MANA_skin_cancer_type.pdf", plot, height=3, width=7, dpi=600)

library(nlme)
non_skin$response[non_skin$response=="responder"]<-1
non_skin$response[non_skin$response=="non_responder"]<-0
non_skin$response<-as.numeric(non_skin$response)
skin$response[skin$response=="responder"]<-1
skin$response[skin$response=="non_responder"]<-0
skin$response<-as.numeric(skin$response)

effect_model<-lme(response~CXCL13_MANA_score,random=~1|sample,data=non_skin,method="ML")
null_model<-lme(response~1,random=~1|sample,data=non_skin,method="ML")
anova(null_model,effect_model)

effect_model<-lme(response~TCF7_MANA_score,random=~1|sample,data=skin,method="ML")
null_model<-lme(response~1,random=~1|sample,data=skin,method="ML")
anova(null_model,effect_model)

###Figure3C heatmap
patient<-read.csv("MANA_score.csv")

cat("add DC")
combined_DC_integrated<-readRDS("combined_DC_integrated_16dataset_v2.rds")
DC_proportion<-as.data.frame.array(prop.table(table(combined_DC_integrated$sample_ID,combined_DC_integrated$clusters_0.4),margin=1))
DC_proportion$sample_ID<-rownames(DC_proportion)
DC_tumor_type<-combined_DC_integrated@meta.data[,5:6] %>% distinct(sample_ID, .keep_all = TRUE)
DC_proportion<-merge(DC_proportion,DC_tumor_type,by="sample_ID")
cor_table<-merge(patient,DC_proportion,by="sample_ID")
rownames(cor_table)<-cor_table[,1]
cor_table<-cor_table[,-1]

library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(corrplot)

R<-data.frame()
P<-data.frame()
for(i in unique(cor_table$cancer_type)){
  cor<-cor_table[cor_table$cancer_type==i,]
  cor_data <- rcorr(as.matrix(cor[,-(17)]),type ="spearman")
  cor_data[["P"]][is.na(cor_data[["P"]])]<-"1"
  assign(paste0("plot", i),ComplexHeatmap::pheatmap(cor_data[["r"]][1:2,-(1:2)],name="correlation",cluster_row = FALSE, cluster_col = FALSE,row_names_side = "left", row_dend_side = "right",
                         cellwidth=20,cellheight=20,column_names_side = "top", column_dend_side = "bottom",col = colorRamp2(seq(-1,1,0.01), colors = rev(COL2('RdBu',201))),
                         display_numbers =matrix(ifelse(cor_data[["P"]][1:2,-(1:2)] <=0.05, "*", ""), nrow(cor_data[["P"]][1:2,-(1:2)]))))
  r<-as.data.frame(cor_data$r[1:2,])
  r$cancer_type<-i
  R<-rbind(R,r)
  p<-as.data.frame(cor_data$P[1:2,])
  p$cancer_type<-i
  P<-rbind(P,p)
}

CXCL13_two<-cor_table[cor_table$cancer_type %in% c("BRCA","RCC"),]
CXCL13_two<- rcorr(as.matrix(CXCL13_two[,-(17)]),type ="spearman")
CXCL13_two[["P"]][is.na(CXCL13_two[["P"]])]<-"1"
CXCL13<-cor_table[cor_table$cancer_type %in% c("BRCA","RCC","LC","CRC","OV"),]
CXCL13<- rcorr(as.matrix(CXCL13[,-(17)]),type ="spearman")
CXCL13[["P"]][is.na(CXCL13[["P"]])]<-"1"
TCF7<-cor_table[cor_table$cancer_type %in% c("BCC","SKCM"),]
TCF7<- rcorr(as.matrix(TCF7[,-(17)]),type ="spearman")
TCF7[["P"]][is.na(TCF7[["P"]])]<-"1"

RCXCL13<-R[c(1,3,5,9,11),]
rownames(RCXCL13)<-RCXCL13[,17]
RCXCL13<-rbind(RCXCL13,CXCL13_two[["r"]][1,],CXCL13[["r"]][1,])
rownames(RCXCL13)[6:7]<-c("BRCA+RCC","combined")
PCXCL13<-P[c(1,3,5,9,11),]
rownames(PCXCL13)<-PCXCL13[,17]
PCXCL13<-rbind(PCXCL13,CXCL13_two[["P"]][1,],CXCL13[["P"]][1,])
rownames(PCXCL13)[6:7]<-c("BRCA+RCC","combined")
plot1<-ComplexHeatmap::pheatmap(RCXCL13[,3:16],name="correlation",cluster_row = FALSE, cluster_col = FALSE,row_names_side = "left", row_dend_side = "right",
                                                  cellwidth=20,cellheight=20,column_names_side = "top", column_dend_side = "bottom",col = colorRamp2(seq(-1,1,0.01), colors = rev(COL2('RdBu',201))),
                                                  display_numbers =matrix(ifelse( PCXCL13[,3:16]<=0.05, "*", ""), nrow(PCXCL13[,3:16])),row_order = c("RCC","BRCA","LC","CRC","OV","BRCA+RCC","combined"))

RTCF7<-R[c(8,14),]
rownames(RTCF7)<-RTCF7[,17]
RTCF7<-rbind(RTCF7,TCF7[["r"]][2,])
rownames(RTCF7)[3]<-"combined"
PTCF7<-P[c(8,14),]
rownames(PTCF7)<-PTCF7[,17]
PTCF7<-rbind(PTCF7,TCF7[["P"]][2,])
rownames(PTCF7)[3]<-"combined"
plot2<-ComplexHeatmap::pheatmap(RTCF7[,3:16],name="correlation",cluster_row = FALSE, cluster_col = FALSE,row_names_side = "left", row_dend_side = "right",
                                cellwidth=20,cellheight=20,column_names_side = "top", column_dend_side = "bottom",col = colorRamp2(seq(-1,1,0.01), colors = rev(COL2('RdBu',201))),
                                display_numbers =matrix(ifelse( PTCF7[,3:16]<=0.05, "*", ""), nrow(PTCF7[,3:16])),row_order = c("BCC","SKCM","combined"))

pdf("MANA/MANA_Proxy.pdf",width = 8, height = 8)
draw(plot1%v%plot2)
dev.off()


###Figure3D-E
CXCL9<-cor_table[cor_table$cancer_type %in% c("LC","CRC"),]
plot<-ggplot(CXCL9, aes(x=CXCL13_MANA_score,y=`cDC2_CXCL9+`,col=cancer_type)) +geom_point(aes())+ geom_smooth(method = 'lm', se = F,linewidth = 0.5)+theme_classic()+scale_color_d3()+ stat_cor(data=CXCL9, method = "spearman")
ggsave("MANA/CXCL9_CXCL13.pdf", plot, height=4, width=6, dpi=600)

dys<-cor_table[cor_table$cancer_type %in% c("BRCA","CRC","OV"),]
plot<-ggplot(dys, aes(CXCL13_MANA_score,`cDC2_CD1C+_B`,col=cancer_type)) +geom_point(aes())+ geom_smooth(method = 'lm', se = F,linewidth = 0.5)+theme_classic()+scale_color_d3()+ stat_cor(data=dys, method = "spearman")
ggsave("MANA/dys_CXCL13_non_skin.pdf", plot, height=4, width=6, dpi=600)

###Figure3F
###CNA analysis
MANA<-read.csv("MANA_score.csv")
CNA_combined_DC<-subset(combined_DC_integrated,sample_ID %in% MANA$sample_ID)
metadata<-CNA_combined_DC@meta.data
metadata$cells<-rownames(metadata)
metadata<- merge(metadata,MANA,by="sample_ID",all.x=TRUE)
order<-rownames(CNA_combined_DC@meta.data)
metadata <- metadata %>% slice(match(order, cells))
rownames(metadata)<-metadata$cells
CNA_combined_DC@meta.data <-metadata
CNA_combined_DC<-FindNeighbors(CNA_combined_DC, reduction = "pca", dims = 1:30)

CNA_combined_DC_non_skin<-subset(CNA_combined_DC,cancer_type %in% c("BC","RCC","LC","CRC","OvC"))
obj <- association.Seurat(
  seurat_object = CNA_combined_DC_non_skin,
  test_var = 'CXCL13_MANA_score',
  samplem_key = 'sample_ID',
  graph_use ="integrated_nn",
  verbose = TRUE,
  batches = NULL, ## no batch variables to include
  covs = NULL ## no covariates to include
)
p<-FeaturePlot(obj, features = c('cna_ncorrs'))[[1]] +
  scale_color_gradient2_tableau() +
  labs(
    title = 'CXCL13_MANA in non-skin cancer', color = 'Correlation',
    subtitle = sprintf('global p=%0.3f', obj@reductions$cna@misc$p)
  )
ggsave("CNA_MANA_CXCL13_non_skin.pdf", p, height=5, width=5, dpi=600)

CNA_combined_DC_skin<-subset(CNA_combined_DC,cancer_type %in% c("BCC","Melanoma"))
obj <- association.Seurat(
  seurat_object = CNA_combined_DC_skin,
  test_var = 'TCF7_MANA_score',
  samplem_key = 'sample_ID',
  graph_use ="integrated_nn",
  verbose = TRUE,
  batches = NULL, ## no batch variables to include
  covs = NULL ## no covariates to include
)
p<-FeaturePlot(obj, features = c('cna_ncorrs'))[[1]] +
  scale_color_gradient2_tableau() +
  labs(
    title = 'TCF7_MANA in skin cancer', color = 'Correlation',
    subtitle = sprintf('global p=%0.3f', obj@reductions$cna@misc$p)
  )
ggsave("CNA_MANA_TCF7_skin.pdf", p, height=5, width=5, dpi=600)



###Figure3G
#######psudebulk DE https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
library(SingleCellExperiment)
library(Matrix.utils)
library(tidyverse)
library(magrittr)
library(DESeq2)
library(apeglm)

combined_TCELL_integrated<-subset(combined_TCELL_integrated,CD8A>0,slot="counts")
TCELL<-data.frame(table(combined_TCELL_integrated$sample_ID))
filter<-TCELL[TCELL$Freq>=10,]
combined_TCELL_integrated<-subset(combined_TCELL_integrated,sample_ID %in% filter$Var1)
non_skin_DC<-subset(combined_DC_integrated,cancer_type %in% c("BC","CRC","LC","OvC","RCC")&(sample_ID %in% rownames(cor_table)))
non_skin_TCELL<-subset(combined_TCELL_integrated,cancer_type %in% c("BC","CRC","LC","OvC","RCC")&(sample_ID %in% cor_table$sample_ID))
CXCL13_patient<-aggregate(non_skin_TCELL$CXCL13_score,by=list(sample_ID=non_skin_TCELL$sample_ID),mean)
colnames(CXCL13_patient)[2]<-"CXCL13_MANA_score"
CXCL13_patient$CXCL13_MANA_class <- ifelse(CXCL13_patient$CXCL13_MANA_score>= median(CXCL13_patient$CXCL13_MANA_score),"High","Low")
CXCL13_high_patient<-CXCL13_patient$sample_ID[CXCL13_patient$CXCL13_MANA_class=="High"]
non_skin_DC[["CXCL13_MANA_class"]] <- ifelse(non_skin_DC$sample_ID %in% CXCL13_high_patient,"High","Low")

counts <-non_skin_DC@assays$RNA@counts
metadata <- non_skin_DC@meta.data

sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)
groups <- colData(sce)[, c("clusters_0.4", "sample_ID")]

cluster_names <- levels(colData(sce)$clusters_0.4)
sample_names <- levels(colData(sce)$sample_ID)

aggr_counts <- aggregate.Matrix(t(counts(sce)),
                       groupings = groups, fun = "sum")

aggr_counts <- t(aggr_counts)

splitf <- str_match(colnames(aggr_counts), "(cDC2_CD1C\\+_A|cDC2_CD1C\\+_B|cDC2_C1Q\\+|pDC|cDC2_HSP\\+|CCR7\\+_DC|cDC1|cDC2_S100B\\+|cDC2_CD207\\+|cDC2_ISG15\\+|Proliferating_DC|T_cDC2_mixed|cDC2_CXCL9\\+|AS_DC)_")[,2]

counts_ls <- list()

for (i in 1:length(cluster_names)) {

  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(splitf == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
}

##Generating matching metadata at the sample-level
metadata <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::select(CXCL13_MANA_class, sample_ID)

metadata <- metadata[!duplicated(metadata), ]

t <- table(colData(sce)$sample_ID,
           colData(sce)$clusters_0.4)

metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
    ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
    
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$cluster_id <- names(counts_ls)[i]
    df$sample_ID  <- gsub("(cDC2_CD1C\\+_A|cDC2_CD1C\\+_B|cDC2_C1Q\\+|pDC|cDC2_HSP\\+|CCR7\\+_DC|cDC1|cDC2_S100B\\+|cDC2_CD207\\+|cDC2_ISG15\\+|Proliferating_DC|T_cDC2_mixed|cDC2_CXCL9\\+|AS_DC)_","",df$cluster_sample_id)
    
    
    ## Retrieve cell count information for this cluster from global cell count table
    idx <- which(colnames(t) == unique(df$cluster_id))
    cell_counts <- t[, idx]
    
    ## Remove samples with zero cell contributing to the cluster
    cell_counts <- cell_counts[cell_counts > 0]
    
    ## Match order of cell_counts and sample_ids
    sample_order <- match(df$sample_ID, names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    
    ## Append cell_counts to data frame
    df$cell_count <- cell_counts
    
    
    ## Join data frame (capturing metadata specific to cluster) to generic metadata
    df <- plyr::join(df, metadata,
                     by = intersect(names(df), names(metadata)))
    
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df) <- df$cluster_sample_id
    
    ## Store complete metadata for cluster i in list
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- unique(df$cluster_id)

}

##Creating a DESeq2 object
idx <- which(names(counts_ls) == "cDC2_CD1C+_B")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

dds <- DESeqDataSetFromMatrix(cluster_counts,
                              colData = cluster_metadata,
                              design = ~ CXCL13_MANA_class)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor, annotation = cluster_metadata[, c("CXCL13_MANA_class"), drop=F])

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)
# Generate results object
res <- results(dds,
               name = "CXCL13_MANA_class_Low_vs_High",
               alpha = 0.05)
res <- lfcShrink(dds,
                 coef = "CXCL13_MANA_class_Low_vs_High",
                 res=res,
                 type = "apeglm")

res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

write.csv(res_tbl,"MANA/pseudobulk_dys_DE.csv")

res_tbl<-res_tbl[!is.na(res_tbl$padj),]

res_tbl$Group<-"not-significant"
#res_tbl$Group[which((res_tbl$padj < 0.05) & (res_tbl$log2FoldChange > 1))]="up-regulated"
#res_tbl$Group[which((res_tbl$padj < 0.05) & (res_tbl$log2FoldChange < -1))]="down-regulated"
res_tbl$Group[which((res_tbl$padj < 0.05) & (abs(res_tbl$log2FoldChange) > 0.25))]="significant"
res_tbl$logP<--log10(res_tbl$padj)
res_tbl$Lable<-""
up.gene<-(res_tbl[((res_tbl$padj < 0.05) & (res_tbl$log2FoldChange > 0.25)),] %>% slice_max(n = 10, order_by = log2FoldChange))$gene
down.gene<-(res_tbl[((res_tbl$padj < 0.05) & (res_tbl$log2FoldChange < -0.25)),]%>% slice_min(n = 10, order_by = log2FoldChange))$gene
top10<-c(as.character(up.gene),as.character(down.gene))
res_tbl$Lable[match(top10,res_tbl$gene)]<-top10
p<-ggscatter(res_tbl,x="log2FoldChange", y="logP", color="Group",palette=c("#BBBBBB","#CC0000"),size=2, label=res_tbl$Lable, font.label=10,repel=T,
xlab="log2FoldChange",
ylab="-log10 (P_value_adjust)")+theme_classic()+
geom_hline(yintercept=1.30, linetype="dashed")+
geom_vline(xintercept=c(-0.25,0.25), linetype="dashed")+
theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
ggsave(paste0("MANA/psedubulk_dys_volcano.pdf"), p, height=4, width=6, dpi=600)
