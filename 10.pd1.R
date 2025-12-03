
dir.create('results')
library(survival)
library(ggplot2)
library(ggpubr)
#读取评分文件
risk=read.table("risk.train.txt", header=T, sep="\t", check.names=F,row.names=1)
head(risk)
riskGenes <- c("TAGLN2", "PCMT1", "PTMA", "CXCL13", "TUBA3D", "DCTPP1") 

#读取临床数据
clinical <- read.table('GSE91061_time.txt', header = T)
clinical <- clinical[, c("OS.time", "OS","Response")]
table(clinical$Response)
clinical <- na.omit(clinical)


load("GSE91061_exp.rda")


intersect(rownames(clinical), colnames(exp))

#构建模型
exp_model_data <- cbind(clinical[, c("OS.time", "OS")],
                               t(exp[, rownames(clinical)]))

#提取riskGenes的表达谱
exp.genes <- intersect(riskGenes, colnames(exp_model_data))
exp.genes

#
table(clinical$Response)
clinical1 <- clinical[clinical$Response != 'NE', ]
colnames(clinical1)

exp_model_data <- exp_model_data[rownames(clinical1), ]
exp_model_data<-exp_model_data[,c("OS.time", "OS",exp.genes)]
exp.genes=gsub('-','__',exp.genes)
colnames(exp_model_data)=gsub('-','__',colnames(exp_model_data))
fmla.exp <- as.formula(paste0("Surv(OS.time, OS) ~"
                                     ,paste0(exp.genes,collapse = '+')))
library(survival)
cox.exp <- coxph(fmla.exp, data =as.data.frame(exp_model_data))
exp_lan <- coef(cox.exp)

risk.exp=as.numeric(exp_lan%*%as.matrix(t(exp_model_data[,names(exp_lan)])))


#KM曲线
dat_km<-data.frame(time=exp_model_data$OS.time/365,
                   status=exp_model_data$OS,
                   Risk=ifelse(risk.exp>=median(risk.exp),'High','Low'),
                   riskscore=risk.exp)
rownames(dat_km)=rownames(exp_model_data)
dat_km=merge(dat_km,clinical1,by=0)
head(dat_km)
write.csv(dat_km,"dat_km.csv")

dat_km<-read.csv("dat_km.csv",header = T,row.names = 1)
###绘制K-M曲线####
library(survival)
library(survminer)

rt=dat_km[,c("Row.names","OS.time","OS","Risk")]
row.names(rt)<-rt[,1]
rt$Row.names<-NULL
rt$OS.time<-rt$OS.time/30
#比较高低风险组生存差异，得到显著性p???
  diff=survdiff(Surv(OS.time, OS) ~Risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS.time, OS) ~ Risk, data = rt)
#绘制K-M曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(months)",
                     ylab="Overall survival",
                     break.time.by = 5,
                     palette=c('#CD3333','#1874CD'),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
#输出图形
pdf(file="surv.train.pdf", width = 6.5, height =5.5, onefile = FALSE)
print(surPlot)
dev.off()






#占比图
library(ggplot2)
library(scales)  # 用于百分比格式化
# 转换为长格式（若尚未完成）
results1
results1 <- data.frame(
  type = c("CR/PR", "PD/SD", "CR/PR", "PD/SD"),
  Risk = c("High", "High", "Low", "Low"),
  Percentage = c(0.07, 0.93, 0.32, 0.68)
)
# 转换为百分比数值（15% → 15）
results1$Percentage <- results1$Percentage * 100
# 设置因子顺序（控制堆叠顺序）
results1$type <- factor(results1$type, levels = c("PD/SD", "CR/PR"))
results1$Risk <- factor(results1$Risk, levels = c("Low", "High"))
library(ggplot2)
library(scales)  # 用于百分比格式化
# 临床响应配色方案
clinical_colors <- c("CR/PR" = "#3474b5",  # 蓝色（响应组）
                     "PD/SD" = "#ef766d")  # 红色（非响应组）

fig <- ggplot(results1, aes(x = Risk, y = Percentage, fill = type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5),
            color = "white", size = 5, fontface = "bold") +
  scale_fill_manual(values = clinical_colors) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(x = "", y = "Percentage of Patients") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# 导出高清PDF
fig
ggsave("exp_response.pdf", fig, width = 5, height = 6,device = cairo_pdf)  # 确保矢量图清晰度
dev.off()




#箱线图
head(dat_km)
dat_km$Response <- factor(dat_km$Response, levels = c("CR/PR","PD/SD" ))
fig10h=ggboxplot(dat_km, x='Response', y='riskscore', 
                 fill = "Response", color = "black",
                 palette = c("#3474b5","#ef766d"), 
                 ylab='riskscore',xlab='',
                 add = "boxplot")+ 
  stat_compare_means(comparisons=list(c('CR/PR','PD/SD'))
                     ,method="wilcox.test",label = "p.signif")+labs(fill='Response')
fig10h
ggsave('fig10h.pdf',fig10h,height = 6,width = 6)



dat_km$Response1 <- factor(dat_km$Response1, levels = c('CR','PR','SD','PD' ))
fig10i=ggboxplot(dat_km, x='Response1', y='riskscore', 
                 fill = "Response1", color = "black",
                 palette = c("#3474b5","#ef766d","#00a087","#3c5488"), 
                 ylab='riskscore',xlab='',
                 add = "boxplot")+ 
  stat_compare_means(comparisons=list(c('CR','PD'),c('CR','SD'),
                                      c('PR','PD'),c('PR','SD'))
                     ,method="wilcox.test",label = "p.signif")+labs(fill='Response')
fig10i
ggsave('fig10i.pdf',fig10i,height = 6,width = 6)


