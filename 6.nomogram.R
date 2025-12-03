library(readr)
library(stringr)
library(ggsci)
library(plyr)
library(ggplot2)
library(ggpubr)

#读取临床信息
clin<-read.csv("combin.csv",row.names = 1)
rownames(clin) <- sub("-01$", "", rownames(clin)) #去除行名后缀-01
clin<-clin[c("Age","Stage_T","N","M","Stage","Subtype_mRNA","OS","group")]
clin <- clin[order(clin$group,clin$Age,clin$Stage,clin$Subtype_mRNA,clin$Stage_T,clin$N,clin$M,clin$OS),]
clin=clin[clin$Stage_T != 'TX',]
clin=clin[clin$N != 'NX',]
clin=clin[clin$M != 'MX',]
#临床相关性分析，得到显著性标记
sigVec=c("group")
for(clinical in colnames(clin[, 1:(ncol(clin)-1)])){
  data=clin[c("group", clinical)]
  colnames(data)=c("group", "clinical")
  #data=data[(data[,"clinical"]!="unknow"),]
  tableStat=table(data)
  stat=chisq.test(tableStat)
  pvalue=stat$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(clinical, Sig))
  print(tableStat)
  #print(paste(clinical, pvalue, Sig, sep="\t"))
}
#sigVec=c(sigVec,"Risk")
# 原始 sigVec 示例: c("group", "OS**", "Age*", "Stage_T", "N**", "M", "Stage***", "Subtype_mRNA")
# 生成带星号的变量名映射
sig_labels <- setNames(
  c( "Age*","Stage_T", "N**","M", "Stage***", "Subtype_mRNA","OS**","group"),  # 对应原始列名中的星号
  c( "Age", "Stage_T","N", "M","Stage", "Subtype_mRNA","OS","group")           # 原始变量名（不带星号）
)


clin$Age <- factor(clin$Age)
clin$Stage <- factor(clin$Stage,levels = c('Stage I','Stage II','Stage III','Stage IV'))
clin$Subtype_mRNAtage <- factor(clin$Subtype_mRNA)
clin$Stage_T <- factor(clin$Stage_T)
clin$N <- factor(clin$N)
clin$M <- factor(clin$M)
clin$group <- factor(clin$group)
clin$OS <- factor(clin$OS)


Age <- c(pal_nejm(alpha = 0.9)(8)[3],'#CF4E27')
names(Age) <- levels(clin$Age)
table(clin$Subtype_mRNAtage)
Subtype_mRNAtage <- c('#E0864A','rosybrown',"#3b4992",'#b09c85','goldenrod')
names(Subtype_mRNAtage) <- levels(clin$Subtype_mRNAtage)
table(clin$Stage_T)
Stage_T<- c('#5f559b','#a20056','#808180','#1b1919')
names(Stage_T) <- levels(clin$Stage_T)
table(clin$N)
N<- c('#008b45','#631879','#008280','#bb0021')
names(N) <- levels(clin$N)
table(clin$M)
M<- c('#3c5488','#f39b7f')
names(M) <- levels(clin$M)
table(clin$Stage)
Stage <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(Stage) <- levels(clin$Stage)
table(clin$OS)
OS <- c('#00a087','#e64b35')
names(OS) <- levels(clin$OS)
table(clin$group)
group <- c('#0072b5','#bc3c29')
names(group) <- levels(clin$group)


order=rownames(clin)
table(clin$group)

score=read.table("risk.train.txt",header = TRUE,row.names = 1)
rownames(score) <- sub("-01A$", "", rownames(score))
exp<-score[,c(3:8)]
exp=t(exp)

##绘制临床特征与基因表达热图
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(scales)

# -------------------------------------------------------------------------
Top = HeatmapAnnotation(
  Age = clin$Age,
  Subtype_mRNAtage = clin$Subtype_mRNAtage,
  Stage = clin$Stage,
  OS = clin$OS,
  Stage_T = clin$Stage_T,
  N = clin$N,
  M = clin$M,
  group = clin$group,
  annotation_legend_param = list(
    labels_gp = gpar(fontsize = 10),
    border = TRUE,
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    ncol = 1
  ),
  border = TRUE,
  col = list(
    group = group,
    Age = Age,
    Stage_T = Stage_T,
    Subtype_mRNAtage = Subtype_mRNAtage,
    N = N,
    M = M,
    Stage = Stage,
    OS = OS
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  # 动态添加星号标签
  annotation_label = sig_labels
)

#Figure5d
Heatmap(
  exp[, order],
  name = 'Z-score',
  top_annotation = Top,
  cluster_rows = TRUE,
  col = colorRamp2(c(-2, 0, 2), c('#82cfe1', 'white', '#e64c4c')),
  color_space = "RGB",
  cluster_columns = FALSE,
  border = TRUE,
  row_order = NULL,
  row_names_side = 'left',
  column_order = NULL,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 9),
  column_split = clin$group[order],
  gap = unit(1, "mm"),
  column_title = NULL,
  column_title_gp = gpar(fontsize = 10),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 10),
    border = TRUE,
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  column_gap = unit(2, 'mm')
)
dev.off()


##绘制小提琴图ggviolion
library(ggpubr)
sameSample=intersect(row.names(clin), row.names(score))
clin=clin[sameSample,,drop=F]
score=score[sameSample,,drop=F]
score<-score[,c("RS")]
rt=cbind(clin,score)

rt1=rt[rt$M != 'MX',]
comparisons <- list(c('M0','M1'))
ggviolin(rt1, x = "M", y = "score", fill = "M", palette = c("npg"),
         add = "boxplot", add.params = list(fill = "white"), order = c("M0", "M1")) +
  stat_compare_means(comparisons = comparisons)

rt2=rt
rt2$Stage=ifelse(rt2$Stage %in% c('Stage I','Stage II'),'Stage I-II','Stage III-IV')
comparisons <- list(c('Stage I-II','Stage III-IV'))
pdf('Fig 5e.pdf',width = 6,height = 6)
ggviolin(rt2, x = "Stage", y = "score", fill = "Stage", palette = c("npg"),
         add = "boxplot", add.params = list(fill = "white"), order = c('Stage I-II','Stage III-IV')) +
  stat_compare_means(comparisons = comparisons)
dev.off()

#绘制Figure5h 预测某项临床特征的ROC诊断曲线
library(pROC)
rt1<-read.csv("combin.csv",row.names = 1)
table(rt1$N)
rt1$N_group <- ifelse(rt1$N == "N0", "N0", "N1")
roc <-  roc(N_group ~ RS, data = rt1)
auc(roc)##ROC的曲线下面积
ci(roc,of = 'auc')#置信区间
#绘制绘ROC图
pdf("Fig 5h.pdf",width = 5,height = 5)
plot(roc,
     print.auc=F,   #输出AUC值
     print.thres=F,   #输出cut-off值
     main = "",   #设置图形的标题
     col= "#2fa1dd",   #曲线颜色
     legacy.axes=TRUE,#使横轴从0到1，表示为1-特异度
     lwd=3,##线条粗细
     identity.col="grey",   #对角线颜色
     identity.lty=1,identity.lwd=2)
#text(0.6,0.3,'AUC:0.710 (0.647, 0.772)',adj = 0,col='#2fa1dd')
legend(0.5,0.5,paste0('AUC=0.710'),col = "#2fa1dd",lty = 1,lwd = 2,bty = 'n')


dev.off()
#绘制Figure5f
library(ggplot2)
library(dplyr)
clin<-read.csv("combin.csv",row.names = 1)
rownames(clin) <- sub("-01$", "", rownames(clin)) #去除行名后缀-01
dat<-clin[c("OS","Age","Stage","Stage_T","M","Subtype_mRNA","group")]
# dat=dat[dat$Grade != 'GX',]##还有5个GX的，去掉之后 重新按照median分高低组
table(dat$group)

dat=dat[,-1]
# 按Risk分成High和Low，计算各列数值。
gname <- "group"
vname <- setdiff(colnames(dat), gname)
pie.high <- pie.low <- list()
fisher.p <- c()
for (i in vname) {
  tmp <- table(dat[,gname], dat[,i])
  p <- format(fisher.test(tmp)$p.value,digits = 2)
  names(p) <- i
  fisher.p <- c(fisher.p, p)
  
  pie.dat <- 
    tmp %>% as.data.frame() %>% group_by(Var1) %>% mutate(Pct = Freq/sum(Freq)) %>% as.data.frame()
  
  # 表格内的两行对应Risk的两类：Risk high和Risk low
  pie.high[[i]] <- pie.dat[which(pie.dat$Var1 == "High"),]
  pie.low[[i]] <- pie.dat[which(pie.dat$Var1 == "Low"),]
}

# 设置颜色
black  <- "#1E1E1B"
blue   <- "#3C4E98"
yellow <- "#E4DB36"
orange <- "#E19143"
green  <- "#57A12B"
cherry <- "#8D3A86"

# 创建颜色
status.col <- c("grey80",black)
grade.col <- alpha(green, c(0.4, 0.7,0.8, 1))
stage.col <- alpha(blue, c(0.4, 0.6, 0.8,1))
T.col <- alpha(cherry, c(0.4, 0.6, 0.8, 1))
M.col <- alpha(orange,c(0.4,0.7,1))

dev.off()
# 硬核base plot一块一块画，当然也可以把其中的pie chart提取出来后期AI或者PPT拼接也是比较方便的
pdf("Fig 5f pieTable.pdf",width = 7, height = 5)

showLayout <- F # 默认不在最终pdf的首页显示layout结构，不过建议初次绘制的时候改为TRUE看一下，方便理解

# 设置画面布局，相同数字代表同一区块，数字越多代表该区块所占面积越大（一共25个区域）
layout(matrix(c( 1, 1, 1,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5, 5,  6, 6, 6,
                 7, 7, 7,  8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12,
                 7, 7, 7,  8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12,
                 13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                 13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                 19,19,19, 20,20,20, 21,21,21, 22,22,22, 23,23,23, 24,24,24,
                 25,25,25, 25,25,25, 25,25,25, 25,25,25, 25,25,25, 25,25,25),
              byrow = T,nrow = 7))

if(showLayout) {
  layout.show(n = 25) # 直观展示画布分布
}

#-------------------------#
# 画布区域1-6：绘制图抬头 #
#-------------------------#

par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) # 基础参数，各边界距离为0
plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "KIRC",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Status",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Grade",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Stage",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "T",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "M",cex = 2, col = "white") # 显示图标题

#--------------------------------------#
# 画布区域7-12：绘制High组抬头和扇形图 #
#--------------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "High\n(n = 259)",cex = 2, col = "white") # 显示图标题

# High group
pie(pie.high$Status$Pct, 
    col = status.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$Grade$Pct, 
    col = grade.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$Stage$Pct, 
    col = stage.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$T$Pct, 
    col = T.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$M$Pct, 
    col = M.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#--------------------------------------#
# 画布区域13-18：绘制Low组抬头和扇形图 #
#--------------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Low\n(n = 259)",cex = 2, col = "white") # 显示图标题

# Low group
pie(pie.low$Status$Pct, 
    col = status.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$Grade$Pct, 
    col = grade.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$Stage$Pct, 
    col = stage.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$T$Pct, 
    col = T.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$M$Pct, 
    col = M.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#--------------------------------#
# 画布区域19-24：绘制空抬头和p值 #
#--------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Status"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Grade"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n",# 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Stage"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["T"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["M"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black") # 底部封上黑线
abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#----------------------#
# 画布区域25：绘制图例 #
#----------------------#

plot(0,0,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
legend("topleft",
       legend = c("Alive","Dead",
                  "G1","G2","G3","G4",
                  "i","ii","iii","iv",
                  "T1","T2","T3","T4",
                  "M0","M1","MX"),
       fill = c(status.col,
                grade.col,
                stage.col,
                T.col,
                M.col),
       border = NA, # 图例颜色没有边框
       bty = "n", # 图例没有边框
       cex = 1.2,
       #box.lwd = 3,
       x.intersp = 0.05,
       y.intersp = 1,
       text.width = 0.065, # 图例的间隔
       horiz = T) # 图例水平放置

# 关闭图像句柄
invisible(dev.off())




