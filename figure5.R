library(ggplot2)
library(forcats)
library(ggbeeswarm)
library(ggpubr)
library(scales)

##Panel 5C
data<-read.csv("/Users/qiandanwen/Downloads/mIF.csv")


p<-data %>% mutate( Area = fct_reorder(Area, dysfunction.cDC2)) %>% ggplot(aes(y = dysfunction.cDC2, x = Area)) + geom_col(fill = "#ff7f0e")+theme_classic()+NoLegend()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5,angle = 90,color = "black"))+ ylab("VCAN+_cDC2/cDC2")

ggsave("/Users/qiandanwen/Downloads/dysfunctional_DC_barplot.pdf", p, height=4, width=8, dpi=600)


##Panel 5D
cbPalette2 <- c("#377EB8","#E41A1C")


p<-ggplot(data=data,aes(x=group, y=CD8.CD3_ave,color=group))+theme_classic()+geom_boxplot(lwd=0.2,outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+
  theme(strip.background = element_rect(size = 0.4),axis.line = element_line(size = 0.2),legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+
  stat_compare_means(hjust = 0)+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+ylab("CD8+ T cells/mm2")+scale_y_continuous(trans=pseudo_log_trans(base = 10),breaks = c( 5, 10, 20, 50, 100,200,400,800,1200))
ggsave("/Users/qiandanwen/Downloads/Tcells_boxplot.pdf", p, height=4, width=4, dpi=600)

##Panel 5E
Tcells<-fread("Tcells.csv")
Tumor<-fread("Tumor.csv")
coordinates<- merge(Tcells, Tumor, by = "Metadata_sample_name", allow.cartesian = T, suffixes = c("1", "2"))
coordinates[, distance := ((Location_Center_X1 - Location_Center_X2)^2 + (Location_Center_Y1 - Location_Center_Y2)^2)^0.5]

coordinates <- coordinates[, .SD[distance == min(distance)], by = id2]

coordinates[, id1 := NULL]
coordinates[, id2 := NULL]
fwrite(file = "tumor-T_Cells-distances.csv", x = coordinates)

source("Plot_Functions.R")

distances <- fread("tumor-T_Cells-distances.csv")
samples <- fread("sample_metadata.csv")
distances <- merge(distances, samples, by.x = "Metadata_sample_name", by.y = "sample_name")
distances[, pair := paste0(spatial_analysis_cell_type2, "-", spatial_analysis_cell_type1)]

p<-histogram_density(distances, "distance", "Tumor-CD8_T_cells", "Distance", "Cells", category_col = "comparison",
                  geom = "density",color = "black", fill = "grey", alpha = 0.3, log_x = T, log_y = F,x_breaks = c(1, 10, 100, 1000))+
                  theme(strip.background = element_rect(size = 0.4),axis.ticks = element_line(size = 0.2),axis.line = element_line(size = 0.2))

ggsave("/Users/qiandanwen/Downloads/Tcells_distance.pdf", p, height=4, width=5, dpi=600)

library(lmerTest)
library(lme4)
distances$comparison[distances$comparison=="Low"]<-0
distances$comparison[distances$comparison=="High"]<-1

model <- lmer(distance ~ comparison + (1 | Metadata_sample_name), data = distances)
summary(model)




