########### R script for running Ligand-receptor analysis ###########
########### Created by Danwen Qian on 28-08-2022. ###########
########### Last modified by Danwen Qian on 17-08-2023. ###########

library(tidyverse)
library(magrittr)
library(liana)
library(CellChat)
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(patchwork)

data_Tcell<-readRDS("TCELL_DC_interation.rds")


##Liana
liana_results_Tcell <- liana_wrap(data_Tcell,method = c("natmi", "connectome", "logfc", "sca", "cellphonedb","call_cellchat"),resource = c("Consensus"))
liana_results_Tcell <- liana_results_Tcell %>%
  liana_aggregate()

write.csv(liana_results_Tcell,"Tcell_DC_liana_results.csv")

liana_results_Tcell$target<-factor(liana_results_Tcell$target,levels=c("AS_DC","CCR7+_DC","cDC1","cDC2_C1Q+","cDC2_CD1C+_A","cDC2_CD1C+_B","cDC2_CD207+","cDC2_CXCL9+","cDC2_HSP+","cDC2_ISG15+","cDC2_S100B+","pDC","STMN1+_DC","T_cDC2_Doublets","Tnaive","CD4/CD8_Tpro","CD4/CD8_Tstr","CD4_Th","CD4_Tcm","CD8_Teff","CD8_Tem","Tisg","CD4_Treg","CD8_GZMK_Tex","Ttex","NK_like_T" ))

##panel 4A
pdf("Receptor/chord_Tcells.pdf", width = 6, height = 6)
liana_results_Tcell %>%
     filter(aggregate_rank <= 0.05) %>%
     chord_freq(source_groups = c("CCR7+_DC","cDC2_CD1C+_B","cDC2_CXCL9+","cDC1"),
                target_groups = c("CD8_GZMK_Tex","Tnaive","Ttex","CD8_Tem","CD4_Th", "CD4_Treg","CD4_Tcm","CD8_Teff","Tisg"))
dev.off()



##cellchat panel 4B
data_Tcell_cellchat<-subset(data_Tcell,clusters %in% c("CCR7+_DC","cDC2_CD1C+_B","cDC2_CXCL9+","cDC1","CD8_GZMK_Tex","Tnaive","Ttex","CD8_Tem","CD4_Th", "CD4_Treg","CD4_Tcm","CD8_Teff","Tisg"))

cellchat <- createCellChat(object = data_Tcell_cellchat, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB)
cellchat@DB <-CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

options(future.globals.maxSize = 5 * 1024^3)
cellchat <- computeCommunProb(cellchat, type = "triMean")

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

saveRDS(cellchat,"cellchat_T_DC.rds")

cellchat@idents<-factor(cellchat@idents,levels=c("CCR7+_DC","cDC2_CD1C+_B","cDC2_CXCL9+","cDC1","CD8_GZMK_Tex","Tnaive","Ttex","CD8_Tem","CD4_Th", "CD4_Treg","CD4_Tcm","CD8_Teff","Tisg"))

for(i in 1:33) {
    pathways.show <- cellchat@netP[["pathways"]][i]

    pdf(paste0("Receptor/cellchat_DC_T_",pathways.show,".pdf"), width = 6, height = 6)

    netVisual_aggregate(cellchat, signaling =pathways.show, layout = "circle",sources.use = c(10:13),targets.use = c(1:9))
    dev.off()
}


##Branch analysis panel 4C-D
##the same process for dir and indir branch, then combined them
liana_results_Tcell_branch<-read.csv("liana_results_branch_T.csv")
liana_results_Tcell_branch$target<-factor(liana_results_Tcell_branch$target,levels=c("cDC1_dir","CCR7+_DC_dir","cDC2_CD207+_dir","cDC2_C1Q+_dir","cDC2_S100B+_dir","cDC2_CD1C+_A_dir","cDC2_HSP+_dir","cDC2_CXCL9+_dir", "cDC2_ISG15+_dir","cDC2_CD1C+_B_dir", "AS_DC_dir","cDC2_CD207+_indir","CCR7+_DC_indir", "cDC2_HSP+_indir","cDC2_CD1C+_A_indir","cDC2_S100B+_indir","cDC2_CXCL9+_indir","cDC2_ISG15+_indir","cDC2_C1Q+_indir", "cDC2_CD1C+_B_indir","AS_DC_indir","Tnaive","CD4/CD8_Tpro","CD4/CD8_Tstr","CD4_Th","CD4_Tcm","CD8_Teff","CD8_Tem","Tisg","CD4_Treg","CD8_GZMK_Tex","Ttex","NK_like_T"))

pdf("Receptor/chord_Tcells_branch.pdf", width = 6, height = 6)
liana_results_Tcell_branch %>%
     filter(aggregate_rank <= 0.05) %>%
     chord_freq(source_groups = c("CCR7+_DC_dir","cDC2_CD1C+_B_dir","CCR7+_DC_indir","cDC2_CD1C+_B_indir"),
                target_groups = c("CD8_GZMK_Tex","Tnaive","Ttex","CD8_Tem","CD4_Th", "CD4_Treg","CD4_Tcm","CD8_Teff","Tisg"))
dev.off()

for(i in c("CCR7+_DC","cDC2_CD1C+_B","cDC2_CXCL9+","cDC1")){
    p<-liana_results_Tcell_branch %>%
        filter(aggregate_rank <= 0.05) %>%
        liana_dotplot(source_groups = c(paste0(i,"_dir"),paste0(i,"_indir")),
                      target_groups = c("CD8_GZMK_Tex","Tnaive","Ttex","CD8_Tem","CD4_Th", "CD4_Treg","CD4_Tcm","CD8_Teff","Tisg"),
                      show_complex = F)+ theme(axis.text.x = element_text(size=16,hjust = 1,vjust = 0.5,angle = 90,color = "black"))+scale_colour_gradientn(colours = colorRampPalette(brewer.pal(9, "OrRd"))(100))
    ggsave(paste0("Receptor/",i,"_branch_Tcells.pdf"), p, height=18, width=14, dpi=600)
    
}






