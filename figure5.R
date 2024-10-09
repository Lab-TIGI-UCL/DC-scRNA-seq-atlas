########### R script for fig5. ###########
########### Created by Danwen Qian on 20-02-2024. ###########
########### Last modified by Danwen Qian on 20-02-2024. ###########


library(ggplot2)
library(dplyr)
library(ggsci)
library(forestplot)
library(ggbeeswarm)
library(ggpubr)
library(scales)
library(meta)



DC<-read.csv("CIBERSORTx_DC.csv",row.names = 1)
DC_col<-read.csv("CIBERSORTx_DC.csv",row.names = 1,header  = F)
colnames(DC)<-DC_col[1,]
DC$response[DC$response=="non_responder"]<-0
DC$response[DC$response=="responder"]<-1
DC$response<-as.numeric(DC$response)


###Forest plot,
results<-data.frame(matrix(ncol=6,nrow=14))
results[,1]<-colnames(DC)[1:14]
colnames(results)<-c("Var","OR_mean","OR_1","OR_2","Pvalue","OR")
for(i in 1:14){
    this_dat<-DC[,c(i,18,19)]
    names(this_dat)[1]<-"DC"
    fit_results<-glm(response ~ DC,
    data = this_dat,
    family = binomial(link = "logit"))
    fit.result<-summary(fit_results)
    results[i,2]<-exp(fit.result$coefficients[2,1])
    results[i,3:4]<-confint(fit_results)[2,]
    results[i,3]<-exp(results[i,3])
    results[i,4]<-exp(results[i,4])
    results[i,5]<-fit.result$coefficients[2,4]
    }
results$OR<-paste0(round(results$OR_mean,2),
                        " (",
                        round(results$OR_1,2),
                        "-",
                        round(results$OR_2,2),
                        ")")
results$Pvalue<-round(results$Pvalue,3)
results[2:15,]<-results[1:14,]
results[1,]<-c("Variable","","","","P-value","OR (95% CI)")

write.csv(results,"Cibersortx_DC/OR.csv",row.names = F)

pdf("Cibersortx_DC/forestplot.pdf",  height=6, width=6)
forestplot(labeltext=as.matrix(results[,c(1,6,5)]),
           mean=as.numeric(results$OR_mean) ,
           lower=as.numeric(results$OR_1),
           upper=as.numeric(results$OR_2),
           zero=1,
           boxsize=0.2,
           graph.pos=2,clip=c(0,100),
           lineheight = unit(7,'mm'),
           colgap=unit(2,'mm'),
           lwd.zero=1.5,
           lwd.ci=2,
           col=fpColors(box='#458B00',
                        summary='#8B008B',
                        lines = 'black',
                        zero = '#7AC5CD'),
           xlab="OR",
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab  = gpar(cex = 0.8),
                            cex = 0.9),
           lty.ci = "solid",
           title = "Forestplot",
           line.margin = 0.08)
dev.off()

DC<-read.csv("CIBERSORTx_DC.csv",row.names = 1)

cbPalette2 <- c("#377EB8","#E41A1C")
for(i in c(1:14)){
  this_data<-DC[,c(i,18,19)]
  y<-colnames(this_data)[1]

  p1<-ggplot(data=this_data,aes_string(x="response", y=y,color="response"))+theme_classic()+geom_boxplot(lwd=0.2,outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(strip.background = element_rect(size = 0.4),axis.line = element_line(size = 0.2),legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means(aes(label = sprintf("P = %.3f", as.numeric(..p.format..))),hjust = 0)+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+scale_y_continuous(trans=pseudo_log_trans(base = 10))+ylab(DC_col[1,i])
  ggsave(paste0("Cibersortx_DC/",y,".pdf"),plot = p1, height=4, width=4,dpi=600)
  
  p1<-ggplot(data=this_data,aes_string(x="response", y=y,color="response"))+theme_classic()+geom_boxplot(lwd=0.2,outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(strip.background = element_rect(size = 0.4),axis.line = element_line(size = 0.2),legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means(aes(label = sprintf("P = %.3f", as.numeric(..p.format..))),hjust = 0)+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+facet_wrap(~cancer_type,nrow = 2,scales = "free_y")+scale_y_continuous(trans=pseudo_log_trans(base = 10))+ylab(DC_col[1,i])
  ggsave(paste0("Cibersortx_DC/",y,"_by_cancer.pdf"),plot = p1,height=6, width=8,dpi=600)
  
}

###KM polt
library("survival")
library("survminer")

for(i in 1:14){
    this_dat<-na.omit(DC[,c(i,19,23,24)])
    this_dat$DC<-ifelse(this_dat[,1]>median(this_dat[,1]),"High","Low")
    fit <- survfit(Surv(OS_time, OS_status) ~ DC, data=this_dat)
    p <-ggsurvplot(fit, data = this_dat,
               palette=c("red","blue"),  #更改线的颜色
               legend.labs=c("High","Low"), #标签
               legend.title=colnames(this_dat)[1],
               ylab="Cumulative survival (percentage)",
               xlab = " Time (Months)",
               censor.shape = 124,
               censor.size = 2,
               conf.int = FALSE,
               pval = TRUE,
               pval.coord = c(0, 0.03)
    )
    pdf(paste0("Cibersortx_DC/",colnames(this_dat)[1],"_forestplot.pdf"),height=4, width=4)
    print(p, newpage = FALSE)
    dev.off()
    }

for(i in 1:14){
    p<-ggplot(DC, aes_string("CD8_T",colnames(DC)[i])) +geom_point(aes())+theme_classic()+scale_color_d3()+geom_smooth(method = 'lm', se = T)+
        coord_cartesian(ylim =c(0,max(DC[,i])))+stat_cor(data=DC, method = "spearman")+ylab(DC_col[1,i])
    ggsave(paste0("Cibersortx_DC/",colnames(DC)[i],"_lm.pdf"),p,height=4, width=4,dpi=600)
}

##MANA score
library(GSVA)

exp<-read.csv("/Users/qiandanwen/Documents/postdoc in UCL/CPI2500/Combined_TPM_matrix_all_studies_n1003/Combined_TPM_matrix_all_studies_n1003.csv",header = TRUE)
exp<-exp %>% distinct(gene_id, .keep_all = TRUE)
row.names(exp)<-exp[,1]
exp<-exp[,-1]
exp<-log2(exp+1)
gs<-read.csv("MANA_geneset.csv")
gs<-as.list(gs)
ssGSEA<-gsva(as.matrix(exp),gs,method = "ssgsea",kcdf='Gaussian',ssgsea.norm = TRUE)

for(i in 1:14){
    p<-ggplot(DC, aes_string("MANA_score",colnames(DC)[i])) +geom_point(aes())+theme_classic()+scale_color_d3()+geom_smooth(method = 'lm', se = T)+
        coord_cartesian(ylim =c(0,max(DC[,i])))+stat_cor(data=DC, method = "spearman")+ylab(DC_col[1,i])
    ggsave(paste0("Cibersortx_DC/",colnames(DC)[i],"_MANA.pdf"),p,height=4, width=4,dpi=600)
}


##Tx
Tx<-read.csv("CIBERSORTx_Tx.csv")
for(i in 3:16){
    p<-ggplot(Tx, aes_string("CD8_T",colnames(Tx)[i])) +geom_point(aes())+theme_classic()+geom_smooth(method = 'lm', se = T)+
    coord_cartesian(ylim =c(0,max(Tx[,i])))+stat_cor(data=Tx, method = "spearman")+ylab(DC_col[1,i-2])
    ggsave(paste0("Cibersortx_DC/",colnames(Tx)[i],"_Tx_CD8.pdf"),p,height=4, width=4,dpi=600)
    p<-ggplot(Tx, aes_string("Morisita_LymTumour",colnames(Tx)[i])) +geom_point(aes())+theme_classic()+geom_smooth(method = 'lm', se = T)+
    coord_cartesian(ylim =c(0,max(Tx[,i])))+stat_cor(data=Tx, method = "spearman")+ylab(DC_col[1,i-2])
    ggsave(paste0("Cibersortx_DC/",colnames(Tx)[i],"_Tx_Morisita.pdf"),p,height=4, width=4,dpi=600)
}

Tx_col<-read.csv("CIBERSORTx_Tx.csv",header  = F)
colnames(Tx)<-Tx_col[1,]
for(i in 3:16){
    this_dat<-na.omit(Tx[,c(i,24,25)])
    this_dat$DC<-ifelse(this_dat[,1]>median(this_dat[,1]),"High","Low")
    fit <- survfit(Surv(os_time, os_status) ~ DC, data=this_dat)
    p <-ggsurvplot(fit, data = this_dat,
               palette=c("red","blue"),  #更改线的颜色
               legend.labs=c("High","Low"), #标签
               legend.title=colnames(this_dat)[1],
               ylab="Cumulative survival (percentage)",
               xlab = " Time (Months)",
               censor.shape = 124,
               censor.size = 2,
               conf.int = FALSE,
               pval = TRUE,
               pval.coord = c(0, 0.03)
    )
    pdf(paste0("Cibersortx_DC/",colnames(this_dat)[1],"_Tx_KM.pdf"),height=4, width=4)
    print(p, newpage = FALSE)
    dev.off()
    
    }


for(i in 3:16){
    this_dat<-na.omit(Tx[,c(i,22,24,25)])
    this_dat$DC<-ifelse(this_dat[,1]>median(this_dat[,1]),"High","Low")
    fit <- survfit(Surv(os_time, os_status) ~ DC, data=this_dat)
    p <-ggsurvplot(fit, data = this_dat,
               palette=c("red","blue"),  #更改线的颜色
               legend.labs=c("High","Low"), #标签
               legend.title=colnames(this_dat)[1],
               ylab="Cumulative survival (percentage)",
               xlab = " Time (Months)",
               censor.shape = 124,
               censor.size = 2,
               conf.int = FALSE,
               pval = TRUE,
               pval.coord = c(0, 0.03)
    )
    pdf(paste0("Cibersortx_DC/",colnames(this_dat)[1],"_Tx_KM.pdf"),height=4, width=4)
    print(p, newpage = FALSE)
    dev.off()
    p1<-ggplot(data=this_dat,aes_string(x="DC", y="Morisita_LymTumour",color="DC"))+theme_classic()+geom_boxplot(lwd=0.2,outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(strip.background = element_rect(size = 0.4),axis.line = element_line(size = 0.2),legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means(aes(label = sprintf("P = %.3f", as.numeric(..p.format..))),hjust = 0)+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+scale_y_continuous(trans=pseudo_log_trans(base = 10))
    ggsave(paste0("Cibersortx_DC/",Tx_col[1,i],"_Tx_Morisita_box.pdf"),plot = p1, height=4, width=4,dpi=600)
    }

