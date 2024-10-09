########### R script for identifying DC subtypes. ###########
########### Created by Danwen Qian on 29-10-2021. ###########
########### Last modified by Danwen Qian on 17-07-2021. ###########

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(scDblFinder)
library(scater)
library(ggrepel)
library(ggsci)
library(RColorBrewer)
library(patchwork)



combined_DC_integrated<-readRDS("combined_DC_integrated_16dataset_v2.rds")

DefaultAssay(combined_DC_integrated) <- "RNA"
Idents(combined_DC_integrated)<-"clusters_0.4"

##Figure1b summary plot
p1 <- combined_DC_integrated@meta.data %>% mutate(cancer_type=factor(cancer_type,levels=c("ESCA","PAAD","OV","THCA","UCEC","BCC","SKCM","BRCA","CRC","LC","RCC"))) %>% ggplot(aes(y=cancer_type,fill= cancer_type))+geom_bar( col="black")+theme_classic()+ggtitle("Cells")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.title.x=element_blank())+theme(axis.title.y=element_text(size = 16),axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))+scale_fill_d3("category20c")+NoLegend() ##fill= "#FFFFCC"##

sample<-combined_DC_integrated@meta.data %>% distinct(sample_ID, .keep_all = TRUE)
p2 <- sample %>% mutate(cancer_type=factor(cancer_type,levels=c("ESCA","PAAD","OV","THCA","UCEC","BCC","SKCM","BRCA","CRC","LC","RCC"))) %>% ggplot(aes(y=cancer_type,fill= cancer_type))+geom_bar(col="black")+theme_classic()+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+ggtitle("Samples")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.title.x=element_blank())+scale_fill_d3("category20c")+NoLegend()+theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))

patient<-combined_DC_integrated@meta.data %>% distinct(patient_ID, .keep_all = TRUE)
p3 <- patient %>% mutate(cancer_type=factor(cancer_type,levels=c("ESCA","PAAD","OV","THCA","UCEC","BCC","SKCM","BRCA","CRC","LC","RCC"))) %>% ggplot(aes(y=cancer_type,fill= cancer_type))+geom_bar(col="black")+theme_classic()+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+ggtitle("Patients")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.title.x=element_blank())+scale_fill_d3("category20c")+NoLegend()+theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))

p<-p1+p2+p3+plot_annotation(caption ="Counts",theme = theme(plot.caption = element_text(hjust = 0.5,size = 14)))

ggsave("plots/summary_barplot.pdf", p, height=8, width=20, dpi=600)

#pie-chart
cell_number<-data.frame(table(combined_DC_integrated$cancer_type))
cell_number<-cell_number[order(cell_number$Freq),]
piepercent<- round(100*cell_number$Freq/sum(cell_number$Freq), 1)
labels<-paste0(cell_number$Var1," ",piepercent,"%")
myPalette <- pal_d3("category20c")(20)
png(file = "plots/pie_chart_cells.png",height=480, width=480)
pie(cell_number$Freq, labels=labels,col =myPalette,cex=1)
dev.off()

number<-data.frame(table(sample$cancer_type))
number<-number[order(number$Freq),]
piepercent<- round(100*number$Freq/sum(number$Freq), 1)
labels<-paste0(number$Var1," ",piepercent,"%")
png(file = "plots/pie_chart_samples.png",height=960, width=960)
pie(number$Freq, labels=labels,col =myPalette,cex=1)
dev.off()

number<-data.frame(table(patient$cancer_type))
number<-number[order(number$Freq),]
piepercent<- round(100*number$Freq/sum(number$Freq), 1)
labels<-paste0(number$Var1," ",piepercent,"%")
png(file = "plots/pie_chart_patients.png",height=960, width=960)
pie(number$Freq, labels=labels,col =myPalette,cex=1)
dev.off()

##Figure1c UMAP
umap.df <- FetchData(combined_DC_integrated, vars = c("UMAP_1", "UMAP_2", "clusters_0.4"))
umap.cols <- DiscretePalette(length(unique(umap.df$clusters_0.4)), palette = "polychrome")
ggrepel.df <- umap.df %>% group_by(clusters_0.4) %>% summarise_at(c("UMAP_1", "UMAP_2"), mean)
ggrepel.df <- rbind(ggrepel.df, umap.df %>% mutate(clusters_0.4 = ""))
p<-ggplot(umap.df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(fill = clusters_0.4), size = 1,  pch = 21)+
    geom_label_repel(umap.df %>% group_by(clusters_0.4) %>% summarise_at(c("UMAP_1", "UMAP_2"), mean), size = 2, mapping = aes(x = UMAP_1, y = UMAP_2, label = clusters_0.4)) +guides(colour = guide_legend(override.aes = list(size=2))) + theme_classic() + NoLegend() +theme(panel.border =  element_rect(colour = "black", fill = NA, size = 1),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(text = element_text(size = 15))+scale_fill_d3("category20c")
ggsave("plots/custom_UMAP.pdf", p, height=4, width=4, dpi=600)

##Figure1d FeaturePlot
p1<-FeaturePlot(combined_DC_integrated,features=c("CLEC9A","FCER1A","CCR7","LILRA4","MKI67","SIGLEC6"),ncol=3,cols = c("lightgrey","coral2"))
ggsave("plots/FeaturePlot.pdf", p1, height=10, width=15, dpi=600)

##Figure1e bubble plot
markers.to.plot <- c("CLEC10A","CEBPD","CD1C","FCER1A","HLA-DRB5","S100A9","S100A8","C1QA","C1QB","C1QC","HSPA6","BAG3","DNAJB1","CD1A","CD207","ISG15","IFI6","IFITM1","CXCL9","CXCL10","S100B","LTB","IL22RA2","CLEC9A","XCR1","C1orf54","CPNE3","LAMP3","CCR7","FSCN1","CD274","CD200","CD40","LILRA4","GZMB","CLEC4C","IL3RA","STMN1","MKI67","AXL","SIGLEC6","PPP1R14A","CD3D","CD3E")
p<-DotPlot(combined_DC_integrated, features = markers.to.plot, cols = c("yellow","blue"), dot.scale = 8) +
   RotatedAxis()
ggsave("plots/Bubble_heatmap.png", p, height=9, width=20, dpi=600)


##Figure1F gene set enrichment analysis
library(msigdbr)
library(fgsea)
library(tibble)
library(pheatmap)

m_df <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

markers <- FindAllMarkers(combined_DC_integrated, min.pct = 0, logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)

fgseaRes<-data.frame()
for(i in c("cDC2_CD1C+_A","cDC2_CD1C+_B","cDC2_C1Q+","cDC2_HSP+","cDC2_CD207+","cDC2_ISG15+","cDC2_CXCL9+","cDC2_S100B+","cDC1","CCR7+_DC","pDC","Proliferating_DC","AS_DC","T_cDC2_mixed")){
cluster.genes<- markers %>% dplyr::filter(cluster == i) %>%arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC)%>% distinct(gene, .keep_all = TRUE)
ranks<- deframe(cluster.genes)
fgseaRes1<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaRes1$cluster<-i
fgseaRes<-bind_rows(fgseaRes, fgseaRes1)
}

top<- fgseaRes %>%
filter(pval < 0.05)%>%
filter(padj < 0.1)%>%
group_by(cluster)

pathway<-as.data.frame(unique(top$pathway))
colnames(pathway)<-"pathway"


for(i in c("cDC2_CD1C+_A","cDC2_CD1C+_B","cDC2_C1Q+","cDC2_HSP+","cDC2_CD207+","cDC2_ISG15+","cDC2_CXCL9+","cDC2_S100B+","cDC1","CCR7+_DC","pDC","Proliferating_DC","AS_DC","T_cDC2_mixed")){
    cluster.pathway<- fgseaRes %>% dplyr::filter(cluster == i)%>%dplyr::select(pathway,NES)
    pathway<-merge(pathway, cluster.pathway,by ="pathway")
}

colnames(pathway)[-1]<-c("cDC2_CD1C+_A","cDC2_CD1C+_B","cDC2_C1Q+","cDC2_HSP+","cDC2_CD207+","cDC2_ISG15+","cDC2_CXCL9+","cDC2_S100B+","cDC1","CCR7+_DC","pDC","Proliferating_DC","AS_DC","T_cDC2_mixed")
pathway<-na.omit(pathway)
rownames(pathway)<-pathway$pathway
pathway<-pathway[,-1]

p <- pheatmap(pathway,
cluster_rows = T,
cluster_cols = T,
show_rownames = T,
show_colnames = T,
fontsize = 8)
ggsave("plots/KEGG_heatmap.pdf", p, height=12, width=8, dpi=600)

##Figure1G stacted barplot
plot<-dittoBarPlot(combined_DC_integrated, "clusters_0.4", group.by = "cancer_type",color.panel  = pal_d3("category20c")(14),var.labels.reorder=c(5,6,4,12,2,9,3,7,14,10,13,8,11,1))
ggsave("plots/barplot_cancers.png", plot, height=9, width=9, dpi=600)
P<-propeller(clusters = combined_DC_integrated$clusters_0.4, sample = combined_DC_integrated$sample_ID,
           group = combined_DC_integrated$cancer_type)


##FigureF tissue barplot
DC_proportion<-data.frame(table(combined_DC_integrated$clusters_0.4,combined_DC_integrated$tissue))
DC_proportion$freq<-DC_proportion$Freq/5750*100
DC_proportion$freq[15:28]<-DC_proportion$Freq[15:28]/28453*100
colnames(DC_proportion)[c(1:2,4)]<-c("clusters","tissue","Cell fraction")
plot<-DC_proportion%>%
ggplot(aes(x=clusters,y=`Cell fraction`,fill=tissue)) +
geom_bar(stat = 'identity',position = 'dodge') +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
theme(plot.title = element_text(hjust=0.5, face="bold"))+scale_fill_manual(values=c("steelblue","#E69F00"))
ggsave("plots/barplot_tissue.pdf", plot, height=8, width=10, dpi=600)
