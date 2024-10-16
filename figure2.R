########### R script for trajectory analysis. ###########
########### Created by Danwen Qian on 29-04-2024. ###########
########### Last modified by Danwen Qian on 29-04-2024. ###########
###monocle3

library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(ggsci)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)


combined_DC_integrated<-readRDS("combined_DC_integrated_16dataset_v2.rds")
Idents(object = combined_DC_integrated) <- "clusters_0.4"

##Figure2A Pseudotime plot
cds <- as.cell_data_set(combined_DC_integrated)
integrated.sub<-subset(combined_DC_integrated,clusters_0.4 %in% c("cDC2_CD1C+_A","cDC2_CD1C+_B","cDC2_C1Q+","cDC2_HSP+","cDC2_CD207+","cDC2_ISG15+","cDC2_CXCL9+","cDC2_S100B+","cDC1","CCR7+_DC","AS_DC"))
cds_sub <- as.cell_data_set(integrated.sub)
rowData(cds_sub)$gene_name <- rownames(cds_sub)
rowData(cds_sub)$gene_short_name <- rowData(cds_sub)$gene_name
cds_sub <- cluster_cells(cds_sub)

cds_sub <- learn_graph(cds_sub)
p<-plot_cells(cds_sub,
              color_cells_by = "clusters_0.4",
              label_groups_by_cluster=FALSE,
              label_cell_groups = TRUE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              trajectory_graph_segment_size = 0.25,
              trajectory_graph_color = "grey50",
              group_label_size = 4)+scale_color_d3("category20c")
ggsave("monocle3/trajectory_11subsets.pdf", p, height=5, width=5, dpi=600)
###here we set AS_DC as the root
cds_sub <- order_cells(cds_sub)
p<-plot_cells(cds_sub,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=3,
              label_principal_points = FALSE,trajectory_graph_segment_size= 0.5)
ggsave("monocle3/pseudotime_11subsets.pdf", p, height=5, width=5, dpi=600)

##Figure2B
cds_sub <- estimate_size_factors(cds_sub)
mature_genes<-c("CD274","CD200","CD40","CCR7")
mature_cds <- cds_sub[rowData(cds_sub)$gene_short_name %in% mature_genes,]
mature_cds <- order_cells(mature_cds)
p<-plot_genes_in_pseudotime(mature_cds,
                         color_cells_by="clusters_0.4",
                         min_expr=0.5,ncol=2)+scale_color_d3("category20c")
ggsave("monocle3/genes_in_pseudotime.pdf", p, height=4, width=6, dpi=600)

##Figure2C Branch analysis
##Analyzing branches in single-cell trajectories,run for direct and indirect
cds_subset <- choose_graph_segments(cds_sub, clear_cds = F, return_list = T)
colData(cds_sub)$chosen<-"Unchosen"
colData(cds_sub)$chosen[colnames(cds_sub) %in% cds_subset$cells]<-"Chosen"
p<-plot_cells(cds_sub,
              color_cells_by = "chosen",
              label_groups_by_cluster=FALSE,
              label_cell_groups = T,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              group_label_size = 5,trajectory_graph_segment_size= 0.5)+scale_color_manual(values=c("purple","grey"))
ggsave("monocle3/branch_directly.pdf", p, height=5, width=5, dpi=600)

##Figure2D-E
Branch1<-choose_graph_segments(cds_sub, clear_cds = F, return_list = T)
Branch2<-choose_graph_segments(cds_sub, clear_cds = F, return_list = T)

overlap<-Branch1$cells[Branch1$cells %in% Branch2$cells]
Branch1<-Branch1$cells[!(Branch1$cells %in% overlap)]
Branch2<-Branch2$cells[!(Branch2$cells %in% overlap)]

DC_branch<-subset(combined_DC_integrated, cells = c(Branch1,Branch2))
DC_branch$Branch<-"Direct_branch"
DC_branch$Branch[colnames(DC_branch)  %in% Branch2]<-"Indirect_branch"

DC_branch_dys<-subset(DC_branch,clusters_0.4=="cDC2_CD1C+_B")
Idents(object=DC_branch_dys)<-"Branch"
de.markers <- FindMarkers(DC_branch_dys, ident.1 = "Direct_branch", ident.2 = "Indirect_branch",logfc.threshold=0)
write.csv(de.markers,file="monocle3/Branch_dys_DE.csv")
de.markers$Group<-"not-significant"
de.markers$Group[which((de.markers$p_val_adj < 0.05) & (abs(de.markers$avg_log2FC) > 1))]="significant"
de.markers$logP<--log10(de.markers$p_val_adj)
de.markers<-de.markers[de.markers$logP!="Inf",]
de.markers$Lable<-""
up.gene<-rownames(de.markers[((de.markers$p_val_adj < 0.05) & de.markers$avg_log2FC > 1),] %>% slice_max(n = 5, order_by = avg_log2FC))
down.gene<-rownames(de.markers[((de.markers$p_val_adj < 0.05) & de.markers$avg_log2FC < -1),] %>% slice_min(n = 5, order_by = avg_log2FC))
top10<-c(as.character(up.gene),as.character(down.gene))
de.markers$Lable[match(top10,rownames(de.markers))]<-top10
p<-ggscatter(de.markers,x="avg_log2FC", y="logP", color="Group",palette=c("#BBBBBB","#CC0000"),size=2, label=de.markers$Lable, font.label=8,repel=T,
xlab="log2FoldChange",
ylab="-log10 (P_value_adjust)")+
geom_hline(yintercept=1.30, linetype="dashed")+
geom_vline(xintercept=c(-1,1), linetype="dashed")+
theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
ggsave("monocle3/Branch_dys_volcano.pdf", p, height=5, width=5, dpi=600)

DC_branch_CCR7<-subset(DC_branch,clusters_0.4=="CCR7+_DC")
Idents(object=DC_branch_CCR7)<-"Branch"
de.markers <- FindMarkers(DC_branch_CCR7, ident.1 = "Direct_branch", ident.2 = "Indirect_branch",logfc.threshold=0)
write.csv(de.markers,file="monocle3/Branch_CCR7_DE.csv")
de.markers$Group<-"not-significant"
de.markers$Group[which((de.markers$p_val_adj < 0.05) & (abs(de.markers$avg_log2FC) > 1))]="significant"
de.markers$logP<--log10(de.markers$p_val_adj)
de.markers<-de.markers[de.markers$logP!="Inf",]
de.markers$Lable<-""
up.gene<-rownames(de.markers[((de.markers$p_val_adj < 0.05) & de.markers$avg_log2FC > 1),] %>% slice_max(n = 5, order_by = avg_log2FC))
down.gene<-rownames(de.markers[((de.markers$p_val_adj < 0.05) & de.markers$avg_log2FC < -1),] %>% slice_min(n = 5, order_by = avg_log2FC))
top10<-c(as.character(up.gene),as.character(down.gene))
de.markers$Lable[match(top10,rownames(de.markers))]<-top10
p<-ggscatter(de.markers,x="avg_log2FC", y="logP", color="Group",palette=c("#BBBBBB","#CC0000"),size=2, label=de.markers$Lable, font.label=8,repel=T,
xlab="log2FoldChange",
ylab="-log10 (P_value_adjust)")+
geom_hline(yintercept=1.30, linetype="dashed")+
geom_vline(xintercept=c(-1,1), linetype="dashed")+
theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
ggsave("monocle3/Branch_CCR7_volcano.pdf", p, height=6, width=6, dpi=600)


##Figure2F-H Branch analysis
cds_subset <- cds_sub[,cds_subset$cells]

subset_pr_test_res <- graph_test(cds_subset,neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
cds_subset <- preprocess_cds(cds_subset)
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.01,random_seed=1)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)),
                                cell_group=colData(cds_subset)$clusters_0.4)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
write.csv(agg_mat,"monocle3/modules_directly.csv")

agg_mat <-agg_mat[,c(1:2,5:6)]
p<-pheatmap::pheatmap(agg_mat,
                      scale="column", clustering_method="ward.D2")
ggsave("monocle3/module_directly.pdf", p, height=10, width=3, dpi=600)

#annotation of module
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
m_t2g <- msigdbr(species = "Homo sapiens") %>% filter(gs_cat %in% c("C5","H"))
m_t2g<-m_t2g %>% dplyr::select(gs_name, gene_symbol)
results<-list()

for (i in 1:63){
    genelist<-gene_module_df[gene_module_df$module==i,]$id
    em<-enricher(genelist, TERM2GENE=m_t2g)
    pathway<-em@result %>% filter(p.adjust<0.05)
    results[[paste0("Module_",i)]]<-pathway$Description
}

n.obs <- sapply(results, length)
seq.max <- seq_len(max(n.obs))
results <- sapply(results, "[", i = seq.max)
write.csv(results,"monocle3/module_directly_annotation.csv",row.names = F)

CCR7_dir<-(gene_module_df %>% subset(module %in% c("21","52","17","54","20","2")))$id
em<-enricher(CCR7_dir, TERM2GENE=m_t2g)
p<-barplot(em,showCategory = 20)
ggsave("monocle3/CCR7_dir_pathway.pdf", p, height=8, width=7, dpi=600)

dys_dir<-(gene_module_df %>% subset(module %in% c("46","49","34","18","6","30")))$id
em<-enricher(dys_dir, TERM2GENE=m_t2g)
p<-barplot(em,showCategory = 20)
ggsave("monocle3/dys_dir_pathway.pdf", p, height=8, width=7, dpi=600)

###indirectly way
cds_subset <- choose_graph_segments(cds_sub, clear_cds = F, return_list = T)
cds_subset <- cds_sub[,cds_subset$cells]

subset_pr_test_res <- graph_test(cds_subset,neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
cds_subset <- preprocess_cds(cds_subset)
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.01,random_seed=1)
write.csv(gene_module_df,"monocle3/gene_module_indirectly.csv")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)),
                                cell_group=colData(cds_subset)$clusters_0.4)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
write.csv(agg_mat,"monocle3/modules_indirectly.csv")

agg_mat <-agg_mat[,c(1:2,4:5,7)]
p<-pheatmap::pheatmap(agg_mat,
                      scale="column", clustering_method="ward.D2")
ggsave("monocle3/module_indirectly.pdf", p, height=10, width=4, dpi=600)

#annotation of module
results<-list()

for (i in 1:50){
    genelist<-gene_module_df[gene_module_df$module==i,]$id
    em<-enricher(genelist, TERM2GENE=m_t2g)
    pathway<-em@result %>% filter(p.adjust<0.05)
    results[[paste0("Module_",i)]]<-pathway$Description
}

n.obs <- sapply(results, length)
seq.max <- seq_len(max(n.obs))
results <- sapply(results, "[", i = seq.max)
write.csv(results,"monocle3/module_indirectly_annotation.csv",row.names = F)

CCR7_indir<-(gene_module_df %>% subset(module %in% c("43","3","36","32","5","23","29")))$id
em<-enricher(CCR7_indir, TERM2GENE=m_t2g)
p<-barplot(em,showCategory = 20)
ggsave("monocle3/CCR7_indir_pathway.pdf", p, height=8, width=7, dpi=600)

dys_indir<-(gene_module_df %>% subset(module %in% c("22","24","4","20")))$id
em<-enricher(dys_indir, TERM2GENE=m_t2g)
p<-barplot(em,showCategory = 20)
ggsave("monocle3/dys_indir_pathway.pdf", p, height=8, width=7, dpi=600)













