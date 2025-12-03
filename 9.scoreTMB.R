
#教程公众号：https://mp.weixin.qq.com/s/BCFZjqIMP6q5vBAQvfQujw


rm(list = ls()) 
options(stringsAsFactors = F)
library(tidyverse)
library(maftools)
#合并maf文件
#dir函数从指定的目录（在本例中为"maf_data/"）中检索所有匹配特定模式的文件路径，并将这些路径作为完整路径返回
maffilepath=dir(path = "maf_data/", pattern="masked.maf.gz$", full.names = T, recursive=T)#，dir(path = "maf_data/") 函数调用会获取指定路径下的所有文件和目录的名称,"masked.maf.gz$"，意味着它将匹配所有以masked.maf.gz结尾的文件
head(maffilepath)

#mafdata 将是一个列表，其中每个元素都是一个 MafFrame 对象（成功读取的 MAF 文件）或 NA（读取失败的文件）
mafdata <- lapply(maffilepath, function(x) {
  tryCatch({
    read.maf(x, isTCGA = TRUE)
  }, error = function(e) {
    message(paste("Error reading file:", x, "\n", e$message))
    return(NA)
  })
})
#去掉未成功读取的MAF文件
mafdata_clean1 <- mafdata[!sapply(mafdata, is.na)]

#合并maf
merger_maf=merge_mafs(mafdata_clean1)
maf_tools=merger_maf
save(maf_tools,file = c("BRCA_maftools.Rdata"))#保存rda文件，方便后面随时取用

#
a <- merger_maf@data %>% 
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>% 
  as.data.frame() %>% 
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

gene <- as.character(unique(a$Hugo_Symbol))
sample <- as.character(unique(a$Tumor_Sample_Barcode))

mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))
for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))

gene_top <- gene_count$gene[1:10] # 修改数字，代表TOP多少
#保存文件
write.table(gene_count, file = "gene-mutcount.txt",sep = "\t",row.names = F,quote = F)
save(mat,mat_0_1,file = "TMB-BRCA.rda")
write.csv(mat,"all_mut_type.csv")
write.csv(mat_0_1,"all_mut_01.csv")

#可视化，绘制瀑布图oncoplot
oncoplot(maf = merger_maf,
         top = 30, #前三十个
         fontSize = 0.6, #设置字体大小
         showTumorSampleBarcodes = F) #不显示病人信息

#计算TMB
maf = tmb(maf = merger_maf,
          captureSize = 50,
          logScale = TRUE)   #绘制
maf$sample <- substr(maf$Tumor_Sample_Barcode,1,16)
rownames(maf) <- maf$sample
write.csv(maf,"tmb_res.csv")



##计算TMB
tmb <- tmb(maf = maf_tools,
           captureSize = 50,           
           logScale = T)
head(tmb)
tmb$Tumor_Sample_Barcode <- substr(tmb$Tumor_Sample_Barcode, 1, 12)
write.csv(tmb,"tmb.csv")


#计算mutant-allele tumor heterogeneity
barcode <- unique(maf_tools@data$Tumor_Sample_Barcode)
head(barcode)
MATH <- data.frame()
for (i in barcode){
  out.math = inferHeterogeneity(maf = maf_tools, tsb = i)
  Tumor_Sample_Barcode=unique(out.math$clusterData$Tumor_Sample_Barcode)
  m = unique(out.math$clusterData$MATH)
  out = data.frame(Tumor_Sample_Barcode, m)
  MATH = rbind(MATH, out)
}
head(MATH)      
MATH$Tumor_Sample_Barcode <- substr(MATH$Tumor_Sample_Barcode, 1, 12)
write.csv(MATH,"MATH.csv")


#TMB绘图
#引用包
library(ggpubr)
library(reshape2)
#读取输入文件
tmb=read.csv("tmb.csv",row.names = 1)  #读取TMB数据文件
tmb <- tmb[, c("Tumor_Sample_Barcode", "total_perMB")]
tmb<- tmb[!duplicated(tmb$Tumor_Sample_Barcode), ]
rownames(tmb) <- tmb$Tumor_Sample_Barcode
tmb <- tmb[, -1, drop = FALSE]
colnames(tmb)[which(colnames(tmb) == "total_perMB")] <- "TMB"
#输出整理后的表达数据
tmb=rbind(ID=colnames(tmb), tmb)
write.table(tmb,file="tmb.txt",sep="\t",quote=F,col.names=F)

score=read.csv("TCGA_risk.csv",header = T,row.names = 1) #读取m6A打分的分组文件
rownames(score) <- sub("-01$", "", rownames(score))

#合并数据
tmb=as.matrix(tmb)
tmb[tmb>quantile(tmb,0.975)]=quantile(tmb,0.975)
sameSample=intersect(row.names(tmb), row.names(score))
tmb=tmb[sameSample,,drop=F]
score=score[sameSample,,drop=F]
data=cbind(score, tmb)
data=data[,c("RS", "group", "TMB")]

#设置比较组
data$group=factor(data$group, levels=c("Low", "High"))
group=levels(factor(data$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#设置颜色
bioCol=c("#00468B","#EC0000","#FF9900","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
#绘制箱线图
boxplot=ggboxplot(data, x="group", y="TMB", fill="group",
		          xlab="",
		          ylab="Tumor Burden Mutation",
		          legend.title="group",
		          palette = bioCol )+ 
	    stat_compare_means(comparisons = my_comparisons)
pdf(file="boxplot.pdf",width=5,height=4.5)
print(boxplot)
dev.off()

#相关性图形
length=length(levels(factor(data$group)))
bioCol=c("#00468B","#EC0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
p1=ggplot(data, aes(RS, TMB)) + 
		  xlab("RS")+ylab("Tumor Burden Mutation")+
		  geom_point(aes(colour=group))+
		  scale_color_manual(values=bioCol[1:length])+ 
		  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =RS, y =TMB))
#相关性图形
pdf(file="cor.pdf", width=6, height=4.5)
print(p1)
dev.off()


rm(list = ls()) 
options(stringsAsFactors = F)
#引用包
library(survival)
library(survminer)
TMBFile="tmb.txt"                  #肿瘤突变负荷文件

#读取输入文件
score=read.csv("TCGA_risk.csv",header = T,row.names = 1) #读取m6A打分的分组文件
rownames(score) <- sub("-01$", "", rownames(score))

#读取TMB数据文件
TMB=read.table(TMBFile, header=T, sep="\t", check.names=F, row.names=1)        

#合并数据
sameSample=intersect(row.names(TMB), row.names(score))
TMB=TMB[sameSample,,drop=F]
score=score[sameSample,,drop=F]
data=cbind(score, TMB)

#获取最优cutoff
res.cut=surv_cutpoint(data, time = "OS.time", event = "OS", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
TMBType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$group=="Low", "L-Risk", "H-Risk")
mergeType=paste0(TMBType, "+", scoreType)

#生存曲线函数
bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(OS.time, OS) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS.time, OS) ~ group, data = surData)
  #print(surv_median(fit))
  
  #绘制生存曲线
  width=8
  height=6
  if(length(levels(factor(surData[,"group"])))>2){
    width=8
    height=6
  }
  bioCol=c("#EC0000","#00468B","#FF9900","#4D85BD","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 5,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=T,
                     cumevents=F,
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, onefile = FALSE, width=width, height=height)
  print(surPlot)
  dev.off()
}

#绘制TMB的生存曲线
data$group=TMBType
bioSurvival(surData=data, outFile="TMB.survival.pdf")

#绘制TMB联合打分的生存曲线
data$group=mergeType
bioSurvival(surData=data, outFile="TMB-score.survival.pdf")




#MATH分析
library(ggpubr)
library(reshape2)

#读取输入文件
MATH=read.csv("MATH.csv",row.names = 1)  #读取MATH数据文件
MATH <- MATH[, c("Tumor_Sample_Barcode", "m")]
MATH<- MATH[!duplicated(MATH$Tumor_Sample_Barcode), ]
rownames(MATH) <- MATH$Tumor_Sample_Barcode
MATH <- MATH[, -1, drop = FALSE]
colnames(MATH)[which(colnames(MATH) == "m")] <- "MATH"
#输出整理后的表达数据
math=rbind(ID=colnames(MATH), MATH)
write.table(math,file="MATH.txt",sep="\t",quote=F,col.names=F)

score=read.csv("TCGA_risk.csv",header = T,row.names = 1) #读取m6A打分的分组文件
rownames(score) <- sub("-01$", "", rownames(score))

#合并数据
MATH=as.matrix(MATH)
MATH[MATH>quantile(MATH,0.975)]=quantile(MATH,0.975)
sameSample=intersect(row.names(MATH), row.names(score))
MATH=MATH[sameSample,,drop=F]
score=score[sameSample,,drop=F]
data=cbind(score, MATH)
data=data[,c("RS", "group", "MATH")]

#设置比较组
data$group=factor(data$group, levels=c("Low", "High"))
group=levels(factor(data$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#设置颜色
bioCol=c("#00468B","#EC0000","#FF9900","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]

#绘制箱线图
boxplot=ggboxplot(data, x="group", y="MATH", fill="group",
                  xlab="",
                  ylab="Mutant-Allele Tumor Heterogeneity ",
                  legend.title="group",
                  palette = bioCol )+ 
  stat_compare_means(comparisons = my_comparisons)
pdf(file="MATH.pdf",width=5,height=4.5)
print(boxplot)
dev.off()



rm(list = ls()) 
options(stringsAsFactors = F)
#引用包
library(survival)
library(survminer)
MATHFile="MATH.txt"                  #肿瘤突变负荷文件

#读取输入文件
score=read.table("risk.train.txt",header = T,row.names = 1) #读取m6A打分的分组文件
rownames(score) <- sub("-01A$", "", rownames(score))

#读取MATH数据文件
MATH=read.table(MATHFile, header=T, sep="\t", check.names=F, row.names=1)        

#合并数据
sameSample=intersect(row.names(MATH), row.names(score))
MATH=MATH[sameSample,,drop=F]
score=score[sameSample,,drop=F]
data=cbind(score, MATH)

#获取最优cutoff
res.cut=surv_cutpoint(data, time = "OS.time", event = "OS", variables =c("MATH"))
cutoff=as.numeric(res.cut$cutpoint[1])
MATHType=ifelse(data[,"MATH"]<=cutoff, "L-MATH", "H-MATH")
scoreType=ifelse(data$Risk=="low", "L-Risk", "H-Risk")
mergeType=paste0(MATHType, "+", scoreType)

#生存曲线函数
bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(OS.time, OS) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS.time, OS) ~ group, data = surData)
  #print(surv_median(fit))
  
  #绘制生存曲线
  width=8
  height=6
  if(length(levels(factor(surData[,"group"])))>2){
    width=8
    height=6
  }
  bioCol=c("#EC0000","#00468B","#FF9900","#4D85BD","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 5,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=T,
                     cumevents=F,
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, onefile = FALSE, width=width, height=height)
  print(surPlot)
  dev.off()
}

#绘制MATH的生存曲线
data$group=MATHType
bioSurvival(surData=data, outFile="MATH.survival.pdf")

#绘制MATH联合打分的生存曲线
data$group=mergeType
bioSurvival(surData=data, outFile="MATH-score.survival.pdf")




#绘制高瀑布图
load("BRCA_maftools.Rdata")
tcga.maf=maf_tools
tcga.maf=tcga.maf@data
tcga.maf$Tumor_Sample_Barcode=substr(tcga.maf$Tumor_Sample_Barcode,1,12)

score=read.csv("TCGA_risk.csv",header = T,row.names = 1) #读取m6A打分的分组文件
rownames(score) <- sub("-01$", "", rownames(score))
score$Samples<-rownames(score)

#high
high.sam=score[which(score$group=='High'),"Samples"]
low.sam=score[which(score$group=='Low'),"Samples"]

high.maf=read.maf(maf = tcga.maf[tcga.maf$Tumor_Sample_Barcode %in% high.sam,])
low.maf=read.maf(maf = tcga.maf[tcga.maf$Tumor_Sample_Barcode %in% low.sam,])

library(RColorBrewer)
vc_cols<-brewer.pal(8,"Paired")
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  #'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  #'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
pdf('fig7elow.pdf', height = 8, width = 8)
oncoplot(maf = low.maf, top = 20,
         colors = vc_cols, 
         draw_titv = TRUE)
dev.off()

# 精确颜色映射（基于图片解析）
vc_cols <- c(
  "Missense_Mutation" = "#f46d43",   
  "Splice_Site" = "#66c2a5",         
  "Frame_Shift_Del" = "#d53e4f",    
  "Frame_Shift_Ins" = "#e5f594",     
  "Nonsense_Mutation" = "#fdae61",   
  "Multi_Hit" = "#fee08b"  ,          
  'In_Frame_Del'= "#3288bd"
)

pdf('fig7dhigh.pdf',height = 8,width = 8)
oncoplot(maf = high.maf, top = 20,
         colors = vc_cols,draw_titv = T)
dev.off()



#共性和互斥分析
pdf('fig7f1high.pdf',height = 7,width = 7)
somaticInteractions(high.maf,top = 20,pvalue = c(0.1,0.05))
dev.off()
pdf('fig7f2low.pdf',height = 7,width = 7)
somaticInteractions(low.maf,top = 20,pvalue = c(0.1,0.05),)
dev.off()


# BiocManager::install("TCGAutils")
#CNV分析
load("BRCA_maftools.Rdata")
tcga.maf=maf_tools
tcga.maf=tcga.maf@data
tcga.maf$Tumor_Sample_Barcode=substr(tcga.maf$Tumor_Sample_Barcode,1,12)

score=read.csv("TCGA_risk.csv",header = T,row.names = 1) #读取m6A打分的分组文件
rownames(score) <- sub("-01$", "", rownames(score))
score$Samples<-rownames(score)

#high
high.sam=score[which(score$group=='High'),"Samples"]
low.sam=score[which(score$group=='Low'),"Samples"]

high.maf=read.maf(maf = tcga.maf[tcga.maf$Tumor_Sample_Barcode %in% high.sam,])
low.maf=read.maf(maf = tcga.maf[tcga.maf$Tumor_Sample_Barcode %in% low.sam,])

top10.gene=as.data.frame(table(high.maf@data$Hugo_Symbol))
top10.gene=top10.gene[order(top10.gene$Freq,decreasing = T),]
top10.gene1=as.character(top10.gene$Var1)[1:10]
#
top10.gene=as.data.frame(table(low.maf@data$Hugo_Symbol))
top10.gene=top10.gene[order(top10.gene$Freq,decreasing = T),]
top10.gene2=as.character(top10.gene$Var1)[1:10]
#
top10.genes=unique(c(top10.gene1,top10.gene2))
top10.genes
get_CNV_Preprocess=function(df_cnv){
  df_cnv$`Gene Symbol`=gsub("\\..*","",df_cnv$`Gene Symbol`)
  rownames(df_cnv)=df_cnv$`Gene Symbol`
  
  library(TCGAutils)
  aliquot_id_to_submitter_id=UUIDtoBarcode(colnames(df_cnv)[-c(1:3)]
                                           ,from_type = 'aliquot_ids')
  colnames(df_cnv)[-c(1:3)]=aliquot_id_to_submitter_id[,2]
  colnames(df_cnv)=substr(colnames(df_cnv),1,15)
  if (length(which(duplicated(colnames(df_cnv))))>0){
    df_cnv=df_cnv[,-which(duplicated(colnames(df_cnv)))]
  }
  
  df_cnv=dplyr::distinct(df_cnv,`Gene Symbol`,.keep_all=TRUE)
  
  ensg=df_cnv$`Gene Symbol`
  
  idmap=clusterProfiler::bitr(ensg,fromType="ENSEMBL",toType="SYMBOL",OrgDb="org.Hs.eg.db")
  idmap=dplyr::distinct(idmap,ENSEMBL,.keep_all=TRUE)## 一个基因匹配到多个geneid,随机取一个
  
  cnv.inds=which(!is.na(idmap$SYMBOL)) ## 去掉没有注释到gene symbol的行
  idmap=idmap[cnv.inds,]
  df_cnv=df_cnv[idmap$ENSEMBL,]
  df_cnv$`Gene Symbol`=idmap$SYMBOL
  
  df_cnv=dplyr::distinct(df_cnv,`Gene Symbol`,.keep_all=TRUE)
  df_cnv=df_cnv[,-c(2:3)]
  return(df_cnv)
}


cnv.all=read.delim('BRCA_merge_copy_number.txt',sep = '\t'
                   ,stringsAsFactors = F,header = T,check.names = F)
cnv.all=get_CNV_Preprocess(cnv.all)
cnv.all[1:4,1:4]
rownames(cnv.all)=cnv.all$`Gene Symbol`
cnv.all=cnv.all[,-1]
cnv.all[1:4,1:4]
cnv.all=na.omit(cnv.all[top10.genes,])

GAIN=rowSums(cnv.all> 0)       #拷贝数增加的样品数目
LOSS=rowSums(cnv.all< 0)       #拷贝数缺失的样品数目
GAIN=GAIN/ncol(cnv.all)*100      #拷贝数增加的百分率
LOSS=LOSS/ncol(cnv.all)*100      #拷贝数缺失的百分率
cnv.all.pre=cbind.data.frame(GAIN, LOSS)
cnv.all.pre=cnv.all.pre[order(cnv.all.pre$GAIN,cnv.all.pre$LOSS,decreasing = T),]


#通路基因cnv的统计
cnv.all.pre.max = apply(cnv.all.pre, 1, max)
pdf(file="fig7g.pdf", width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(cnv.all.pre.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(cnv.all.pre.max)))
points(bar,cnv.all.pre[,"GAIN"], pch=20, col="#EC0000", cex=3)
points(bar,cnv.all.pre[,"LOSS"], pch=20, col="#00468B", cex=3)
legend("top", legend=c('GAIN','LOSS'), col=c("#EC0000","#00468B"), pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(cnv.all.pre), adj=1)
dev.off()
