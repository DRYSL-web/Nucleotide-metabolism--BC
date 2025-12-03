remove(list = ls()) ##清空当前环境
#跑GSVA用TPM数据
######---------GSvA富集分析--------####
#针对单样品，又称SSGSEA分析
library(ggplot2) #绘图使用
library(ComplexHeatmap) #绘图使用
library(clusterProfiler) #数据处理使用
library(GSVA) #GSVA使用
library(GSEABase) #数据处理使用
library(dplyr) #数据处理使用
library(data.table) #数据读取使用

#读取数据
TCGA_BRCA_Exp_tpm <-read.csv("TCGA_brca_Tpm.csv",header = T)
colnames(TCGA_BRCA_Exp_tpm)[1] <- "gene_name"

## 去除重复基因
TCGA_BRCA_Exp_tpm <- TCGA_BRCA_Exp_tpm[!duplicated(TCGA_BRCA_Exp_tpm$gene_name),]    
#行名为基因名
rownames(TCGA_BRCA_Exp_tpm) <- TCGA_BRCA_Exp_tpm$gene_name 
#去掉第一列样本名
gsva_data <- TCGA_BRCA_Exp_tpm[,-1] 


##---------1、GSVA富集分析-----------#
gene_setKEGG <- getGmt("KEGG_metabolism_nc.gmt")
gene_set<-gene_setKEGG
gsva_resultKEGG <- gsva(as.matrix(matrix_gsva_data), 
                    gene_set, 
                    method = "gsva",
                    min.sz=1,
                    max.sz=Inf,
                    kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                    verbose=T, 
                    parallel.sz = parallel::detectCores())#调用所有核
write.csv(gsva_resultKEGG, file = './BRCA_GSVA_KEGG.csv')


#BiocManager::install("limma")
library(limma)
group <- substring(colnames(TCGA_BRCA_Exp_tpm), 14, 15)
table(group)
groups <- factor(c(rep("tumor",1113),rep("normal",113))) 
design <- model.matrix(~ groups + 0)
rownames(design) <- colnames(gsva_resultKEGG)
compareE <- makeContrasts(groupstumor-groupsnormal, levels = design) #groupstumor-groupsnormal注意替换
fit <- lmFit(gsva_resultKEGG, design)
fit <- contrasts.fit(fit, compareE)
fit <- eBayes(fit)
#得到未筛选的limma两组间比较结果
fit_all <- topTable(fit, coef = 1, adjust = 'BH', n = Inf)
fit_all <- cbind(rownames(fit_all),fit_all)
colnames(fit_all) = c("geneset", "logFC", "avgExp", "t", "pval", "FDR", "B")

fit_all$group <- ifelse(fit_all$logFC > 0 & fit_all$pval < 0.05 ,"up" ,
                      ifelse(fit_all$logFC < 0 & fit_all$pval < 0.05 ,"down","noSig")
)
write.csv(fit_all, file = './BRCA_KEGG_DEG.csv',row.names = F)


#自定义阈值筛选部分结果
diff_result = fit_all[fit_all$pval < 0.05,]
diff_result$up_down = ifelse(diff_result$t<0,'Down','Up')

#绘制条形图
library(ggplot2)
pp = ggplot(diff_result, aes(reorder(geneset, t), t))  + geom_col(aes(fill=up_down)) + scale_fill_manual(values=c("#6CC570","#2A5078")) +  coord_flip() +
  labs(x="Gene set", y="t value of GSVA Score" ) +
  theme_minimal() +
  geom_text(data = diff_result[diff_result$up_down == 'Down',], aes(x= geneset, y=0, label = geneset), hjust = 0, size = 2.5)+
  geom_text(data = diff_result[diff_result$up_down == 'Up',], aes(x= geneset, y=0, label = geneset), hjust = 1, size = 2.5)+ theme_bw() +
  theme(axis.line.x = element_line(colour = "black"))+
  theme(axis.line = element_blank()) +
  theme(axis.text.y=element_blank()) +
  theme(panel.grid =element_blank(),panel.border = element_blank(),axis.text.x=element_text(size=6),legend.text=element_text(size=8),axis.ticks = element_blank())+ labs( fill = "Up_Down")
#存图
ggsave("./KEGG.barplot.pdf",
       height = 10,width = 8,plot = pp, limitsize = FALSE)

