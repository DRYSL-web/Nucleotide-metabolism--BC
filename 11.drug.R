
dir.create('results')
library(ggpubr)
library(ggplot2)

tcga.risk.cli<-read.csv('TCGA_risk.csv',header = T,row.names = 1)
row.names(tcga.risk.cli) <- gsub("-01$", "", row.names(tcga.risk.cli))
tcga.risk.cli$Samples<-row.names(tcga.risk.cli)

#http://www.sxdyc.com/drugPrrophetic
#运行后的结果进行后续分析
tcga.drug.socre=read.delim('results/predicted.drugs.score.txt',sep='\t',header = T,check.names = F,row.names = 1)
tcga.drug.socre[1:4,1:4]
key.drug=colnames(tcga.drug.socre)
dir.create('results/drugs')

sameSample=intersect(row.names(tcga.drug.socre), row.names(tcga.risk.cli))
tcga.drug.socre=tcga.drug.socre[sameSample,,drop=F]
tcga.risk.cli=tcga.risk.cli[sameSample,,drop=F]

tcga.drug.socre1=cbind.data.frame(tcga.drug.socre[tcga.risk.cli$Samples,],
                                  riskscore=tcga.risk.cli$RS,
                                  risk=tcga.risk.cli$group)
for (i in key.drug){
  #相关性散点图
  p1=ggplot(tcga.drug.socre1, aes_string(x = 'riskscore', y = i))+
    geom_point(color = "#8ab6d6") +  # 改用图片中的蓝色
    geom_smooth(method = "lm", color = "#1C9C56", fill = "#D1E5E0") +  # 改为图片中的绿色和浅灰绿色
    stat_cor(method = "spearman")+
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  
  
  #组间差异图
  p2=ggplot(tcga.drug.socre1,aes_string(x='risk',y = i))+
    geom_jitter(aes_string(colour='risk'),size=1,position = position_jitter(seed=1))+
    geom_boxplot(fill='white',colour='grey60',alpha=0)+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),  # 移除整个面板边框
      axis.line.x = element_line(color = "black"),  # 添加底部x轴线
      axis.line.y = element_line(color = "black"),   # 添加左侧y轴线
      axis.line.x.top = element_blank(),      # 移除顶部x轴线
      axis.line.y.right = element_blank()    # 移除右侧y轴线
    ) +
    xlab('')+ylab(i)+labs(colour='risk')+
    scale_color_manual(values = c("#1C9C55","#2C6A99"))+
    stat_compare_means(comparisons=list(c('High','Low')),
                       method="wilcox.test",label = "p.signif")
  ggsave(paste0('results/drugs/',i,'_cor.pdf'),p1,height = 6,width = 6)
  ggsave(paste0('results/drugs/',i,'_boxplot.pdf'),p2,height = 6,width = 6)
  
}
write.csv(tcga.drug.socre1,"tcga.drug.socre1.csv")


#绘制单个药物的敏感性
library(ggplot2)
library(data.table)
library(aplot)
library(ggpubr)
library(dplyr)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 加载模拟数据(因为是用ggplot做，所以多个组都可以)
All <- read.csv("tcga.drug.socre1.csv",header = T,row.names = 1)
colnames(All)
#提取药物和分组
df<-All[,c("Cisplatin","risk")]
library(dplyr)
df <- df %>% 
  rename(AUC = Cisplatin, Subtype = risk) %>%  # 重命名列
  mutate(Drug = "Cisplatin") %>%               # 添加Drug列
  select(Drug, everything())                   # ▼ 关键操作：将Drug移到第一列

head(df)
# 计算p值
# 这里是多组，用kruskal.test函数，还可以改为aov
# 如果是两组，可以改为t.test，或wilcox.test
p.val <- t.test(AUC ~ Subtype,
                      data = df)
p.lab <- paste0("P",
                ifelse(p.val$p.value < 0.001, " < 0.001",
                       paste0(" = ",round(p.val$p.value, 3)))) 
p.lab
# 设置颜色
green <- "#C7EAB2"
cyan <- "#5FC1C2"
blue <- "#1B90BE"

# 绘制上半部分密度图
p_top <- ggplot(df, aes(x = AUC, color = Subtype, fill = Subtype)) +
  geom_density() +
  # 让箱子的所有位置都颜色统一，如例文所示
  scale_color_manual(values = c(alpha(green,0.7),alpha(cyan,0.7),alpha(blue,0.7))) + # 设置透明色
  scale_fill_manual(values = c(alpha(green,0.7),alpha(cyan,0.7),alpha(blue,0.7))) +
  theme_classic() + # 如果显示采用这一行
  
  # 这里提取输入文件的第一列药物名称，写入x轴标签
  xlab(paste0("Estimated AUC of ", unique(df$Drug))) + 
  # 第一列非必需，可以像下面这样直接写xlab
  #xlab("Estimated AUC of Cisplatin") +
  
  ylab(NULL) + 
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_blank(), # 原文不显示纵轴的密度
        #axis.text.y = element_text(size = 12,color = "black"), # 如果要显示采用这一行
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rug()
p_top
# 画box plot
p_bot <- ggplot(df, aes(Subtype, AUC, fill = Subtype)) + 
  geom_boxplot(aes(col = Subtype)) + 
  scale_fill_manual(values = c(green, cyan, blue)) + 
  scale_color_manual(values = c(green, cyan, blue)) + 
  xlab(NULL) + ylab("Estimated AUC") + 
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_blank(), # 原文不显示箱线图的x轴
        #axis.text.x = element_text(size = 12,color = "black"), # 如要显示箱线图x轴采用这一行
        axis.text.y = element_text(size = 11,color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  annotate(geom="text",
           x = 1.5,
           hjust = 1,
           y = max(df$AUC),
           size = 4, angle = 270, fontface = "bold",
           label = p.lab) +
  coord_flip() # 翻转图像

# 用白色标记箱子的基本统计量
dat <- ggplot_build(p_bot)$data[[1]]
p_bot <- p_bot + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p_bot
library(aplot)
# 使用aplot拼图，底部箱型图稍微小一些
p <- p_top %>% insert_bottom(p_bot, height = 0.4)
pdf(file = "boxdensity.pdf", width = 6,height = 3)
p
invisible(dev.off())










#cmap

#差异分析
limma_DEG=function(exp,group,ulab,dlab){
  library(limma)
  ind1=which(group==ulab)
  ind2=which(group==dlab)
  sml <- c(rep('A',length(ind1)),rep('B',length(ind2)))    # set group names
  eset=exp[,c(ind1,ind2)]
  fl <- as.factor(sml)
  
  design <- model.matrix(~fl+0)
  colnames(design) <- levels(fl)
  cont.matrix<-makeContrasts(contrasts='A-B',levels=design)
  #print(head(eset))
  fit<-lmFit (eset,design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  #print(sml)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(eset))
  return(tT)
}

tcga.diff<-limma_DEG(exp = tcga.dat[,tcga.risk.cli$Samples],
                     group = tcga.risk.cli$risk,
                     ulab = 'High',dlab = 'Low')
fc.fit=log2(1.5);p.fit=0.05
tcga.diff.fit=tcga.diff[which(abs(tcga.diff$logFC)>fc.fit & tcga.diff$adj.P.Val<p.fit),]
tcga.diff.fit$type=ifelse(tcga.diff.fit$logFC>0,'Up','Down')
table(tcga.diff.fit$type)
write.table(tcga.diff.fit,'results/tcga.diff.fit.txt',quote = F,row.names = T,sep='\t')
#https://clue.io/
