#GSEA富集分析

load('TCGA_exp.rda')
data=TCGA
colnames(data) <- gsub("\\.", "-", colnames(data))
tcga.risk=read.delim('risk.train.txt',sep='\t',header = T)
tcga.risk<-tcga.risk[,c("id","Risk")]
row.names(tcga.risk)<-tcga.risk$id
# tcga.risk$id<-NULL
head(tcga.risk)

#确保样品名一致
samesample=intersect(colnames(data),rownames(tcga.risk))
data=data[,samesample]
tcga.risk=tcga.risk[samesample,]

#表达数据，要求行名是基因名，列名是样品名
exprSet=data
#分组数据
group_list=tcga.risk$Risk
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
colnames(deg)
logFC=0.5
P.Value = 0.05
type1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
type2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(type1,"down",ifelse(type2,"up","stable"))
table(deg$change)
allDiff=rownames(deg[deg$change=='up'|deg$change=='down',])
deg_se=deg[allDiff,]
deg=deg[order(abs(deg$logFC),decreasing = T),]

library(tidyverse)
library(clusterProfiler)
library (org.Hs.eg.db)
library(enrichplot)
library(GOplot)
library(ggpubr)
deg_se=rownames_to_column(deg_se,'SYMBOL')
deg_se=deg_se[,-8]

gene.df <- bitr(deg_se$SYMBOL,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df)


library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析

#导入差异分析结果
length(deg_se$SYMBOL)#共19712个基因

#添加entrez ID列：
##：symbol转entrez ID：
symbol <- deg_se$SYMBOL
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")

#准备genelist文件：
##需要的genelist格式：entrez ID+log2fc
genelist <- deg_se$logFC
names(genelist) <- deg_se$SYMBOL
head(genelist)

#genelist过滤(ID转换中丢失的部分基因)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)#过滤后剩18138个基因
head(genelist)

#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)

#GSEA_KEGG富集分析：报错！！！！
KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "hsa",#不同物种选择可官方查询：https://www.genome.jp/kegg/catalog/org_list.html
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#提取结果表：
KEGG_ges_result <- KEGG_ges@result

#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
KEGG_ges_result[1:3,1:10]

#GSEA_GO富集分析：
Go <- gseGO(geneList = genelist,
            ont = 'BP',
            OrgDb = 'org.Hs.eg.db',
            keyType = 'ENTREZID',
            pvalueCutoff = 0.05)

#提取结果表：
Go_result <- Go@result
#显示前10列（第十一列为富集到该基因集的核心基因列表，这里过长不展示）：
Go_result[1:3,1:10]

#GSEA可视化：
##山峦图：
library(ggridges)
library(ggplot2)
library(enrichplot)

p <- ridgeplot(Go,
               showCategory = 10,
               fill = "p.adjust",
               decreasing  = F)
p

#同时展示多个基因集：
p2 <- gseaplot2(Go,
                geneSetID = c('GO:0022409','GO:0019221','GO:0006953'),#或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = c("#003399", "#FFFF00", "#FF6600"),
                pvalue_table = F,
                ES_geom = "line")#"dot"将线转换为点
p2

hallmark_list <- getGmt('h.all.v7.4.symbols.gmt') 

dat <- as.matrix(exprSet)

geneset <- hallmark_list
#gsva_mat =tcga.gsva
gsva_mat <- gsva(expr=dat, 
                 gset.idx.list=geneset, 
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())#调用所有核

#### 进行limma差异处理 ####
##设定 实验组exp / 对照组ctr

design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(gsva_mat)
contrast.matrix <- makeContrasts(contrasts= 'low-high',  
                                 levels = design)

fit1 <- lmFit(gsva_mat,design)                 #拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)                         #修正

summary(decideTests(efit,lfc=1, p.value=0.05)) #统计查看差异结果
tempOutput <- topTable(efit, coef= 'low-high', n=Inf)
degs <- na.omit(tempOutput) 

#degs=pathwy_p
rownames(degs)=gsub('HALLMARK_','',rownames(degs))
#str_replace(rownames(degs), "HALLMARK_","")


#### 发散条形图绘制 ####
library(tidyverse)  
library(ggthemes)
library(ggprism)
p_cutoff=0.001
degs <- degs
Diff <- rbind(subset(degs,logFC>0)[1:20,], subset(degs,logFC<0)[1:20,]) #选择上下调前20通路     
dat_plot <- data.frame(id  = row.names(Diff),
                       p   = Diff$P.Value,
                       lgfc= Diff$logFC)
dat_plot$group <- ifelse(dat_plot$lgfc>0 ,1,-1)    # 将上调设为组1，下调设为组-1
dat_plot$lg_p <- -log10(dat_plot$p)*dat_plot$group # 将上调-log10p设置为正，下调-log10p设置为负
dat_plot$threshold <- factor(ifelse(abs(dat_plot$p) <= p_cutoff,
                                    ifelse(dat_plot$lgfc >0 ,'Up','Down'),'Not'),
                             levels=c('Up','Down','Not'))

dat_plot <- dat_plot %>% arrange(lg_p)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

## 设置不同标签数量
low1 <- dat_plot %>% filter(lg_p < log10(p_cutoff)) %>% nrow()
low0 <- dat_plot %>% filter(lg_p < 0) %>% nrow()
high0 <- dat_plot %>% filter(lg_p < -log10(p_cutoff)) %>% nrow()
high1 <- nrow(dat_plot)

p <- ggplot(data = dat_plot,aes(x = id, y = lg_p, 
                                fill = threshold)) +
  geom_col()+
  coord_flip() + 
  scale_fill_manual(values = c('Up'= '#D43E3C','Not'='#cccccc','Down'= '#4F94C4')) +
  geom_hline(yintercept = c(-log10(p_cutoff),log10(p_cutoff)),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('-log10(P.Value) of GSVA score, low-risk versus high-risk') + 
  guides(fill="none")+
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'black') + #黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 黑色标签

ggsave("Fig 7c GSVA_barplot_pvalue.pdf",p,width = 10,height  = 15)


#相关性分析
tcga.gsva.risk=cbind.data.frame(t(tcga.gsva[,tcga.risk$Samples]),
                                risksocre=tcga.risk$riskscorez)
colnames(tcga.gsva.risk)=gsub('HALLMARK_','',colnames(tcga.gsva.risk))
cor_matr = cor(tcga.gsva.risk)
library(corrplot)
pdf('results/Fig6d.pdf',height = 15,width = 15)
corrplot(cor_matr, type="upper",method = 'ellipse',
         order="hclust", tl.col="black", tl.srt=45)
dev.off()
#
key.pathway=rownames(tcga.gsva)[1:4]
ggsurvplotKM<-function(dat, title = 'Groups', 
                       lables = c(), col = mycolor,
                       risk.table = TRUE,
                       tables.height = 0.25) {
  # dat：数据框，行为样本，列为时间、状态以及分组
  # the color palette to be used. Allowed values include "hue" for the default hue color scale; "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red"); and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty". See details section for more information. Can be also a numeric vector of length(groups); in this case a basic color palette is created using the function palette.
  library(ggplot2)
  library(survminer)
  library(survival)
  library(ggthemes)
  library(ggplotify)
  library(ggsci)
  
  colnames(dat) <- c("OS.time", "OS", "Groups")
  # # 将时间转换成年份
  # if (max(dat[, 1]) > 365) {
  #   dat[, 1] <- dat[, 1] / 365
  # } else if (max(dat[, 1]) > 24) {
  #   dat[, 1] <- dat[, 1] / 12
  # } else {
  #   dat <- dat
  # }
  fit <- survfit(Surv(OS.time, OS) ~ Groups,data=dat)
  surv.fit <- ggsurvplot(fit, data = dat, palette = col,
                         pval = TRUE, 
                         pval.method = T,
                         pval.method.size = 4,
                         pval.method.coord = c(0, 0.15),
                         surv.median.line='hv',
                         linetype = 1, 
                         pval.coord=c(0, 0.05), 
                         pval.size = 4,
                         risk.table = risk.table,
                         risk.table.y.text = FALSE,
                         legend.title = title,
                         legend.labs = lables,
                         xlab = 'Time(years)',
                         ggtheme=theme_bw(),
                         tables.height = tables.height)
  # 将图形转换为 ggplot 对象
  if (risk.table) {
    surv.fit1 <- surv.fit$plot + 
      theme(legend.position=c(1,1), 
            legend.justification=c(1,1),
            plot.margin=unit(c(0.1, 0.15, 0, 0.15), "inches"),
            legend.background = element_rect(fill = NA, colour = NA)
            # ,
            # axis.text.x=element_blank(),
            # axis.title.x=element_blank()
      )
    
    surv.fit2 <- surv.fit$table + 
      theme(plot.title=element_blank(),
            plot.margin=unit(c(0, 0.15, 0, 0.15), "inches")) +
      ylab('')
    surv.fit <- ggpubr::ggarrange(surv.fit1,
                                  surv.fit2, 
                                  ncol = 1, 
                                  nrow = 2,
                                  heights = c(1 - tables.height, 
                                              tables.height),
                                  align = "hv")
  } else {
    surv.fit <- surv.fit$plot + 
      theme(legend.position=c(1,1), 
            legend.justification=c(1,1),
            # plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
            legend.background = element_rect(fill = NA, colour = NA))
  }
  return(surv.fit)
}
mycolor <- ggsci::pal_npg(palette = c("nrc"), alpha =1)(8)
tcga.risk.pathway=cbind.data.frame(t(tcga.gsva[key.pathway,tcga.risk$Samples]),
                                   OS=tcga.risk$OS,OS.time=tcga.risk$OS.time)
tcga.risk.pathway$HALLMARK_TNFA_SIGNALING_VIA_NFKB=ifelse(tcga.risk.pathway$HALLMARK_TNFA_SIGNALING_VIA_NFKB>median(tcga.risk.pathway$HALLMARK_TNFA_SIGNALING_VIA_NFKB),'High','Low')

tcga.risk.pathway$HALLMARK_HYPOXIA=ifelse(tcga.risk.pathway$HALLMARK_HYPOXIA>median(tcga.risk.pathway$HALLMARK_HYPOXIA),'High','Low')

tcga.risk.pathway$HALLMARK_CHOLESTEROL_HOMEOSTASIS=ifelse(tcga.risk.pathway$HALLMARK_CHOLESTEROL_HOMEOSTASIS>median(tcga.risk.pathway$HALLMARK_CHOLESTEROL_HOMEOSTASIS),'High','Low')


fig6e <- ggsurvplotKM(tcga.risk.pathway[, c("OS.time", "OS", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")],
                      lables = c('High', 'Low'),risk.table = F,
                      title = 'TNFA_SIGNALING_VIA_NFKB')
fig6e
ggsave('results/fig6e.pdf',fig6e,height = 7,width = 7)

