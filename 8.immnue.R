
dir.create('results')
load("TCGA_exp.rda")
library(ggplot2)
library(ggpubr)
tcga.dat=TCGA
#tcga.dat行为基因名，列为样品名
#行列名统一
colnames(tcga.dat) <- sub(".01A$", "", colnames(tcga.dat))
colnames(tcga.dat) <- gsub("\\.", "-", colnames(tcga.dat))

tcga.risk.cli=read.csv("TCGA_risk.csv",header = T,row.names = 1) #读取m6A打分的分组文件
rownames(tcga.risk.cli) <- sub("-01$", "", rownames(tcga.risk.cli))
tcga.risk.cli$Samples<-rownames(tcga.risk.cli)

# library(IOBR)
# #进行ESTIMATE分析
# tcga.estimate <- deconvo_tme(eset = tcga.dat, method = "estimate")

# save(tcga.estimate,file = 'tcga.estimate.RData')

load('tcga.estimate.RData')
tcga.estimate=as.data.frame(tcga.estimate)
head(tcga.estimate)

tcga.estimate.risk=merge(tcga.risk.cli[,c("Samples","group")],
                         tcga.estimate,
                         by.x='Samples',
                         by.y='ID')
head(tcga.estimate.risk)

library(ggplot2)
data<-tcga.estimate.risk
#设置比较组
type=levels(factor(data[,"group"]))
data$group <- factor(data$group, levels = c("Low", "High"))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#设置颜色
bioCol=c("#3474b5","#ef766d","#FF9900","#FF0000","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(data$group))]

#对肿瘤微环境打分进行循环,绘制箱线图
for(i in colnames(data)[3:(ncol(data)-1)]){
  boxplot=ggboxplot(data, x="group", y=i, fill="group",
                    xlab="",
                    ylab=i,
                    legend.title="group",
                    palette=bioCol
  )+ 
    stat_compare_means(comparisons=my_comparisons)
  
  pdf(file=paste0(i, ".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}



#免疫相关通路
immnu_patwhays<-list.files('immnue_pathway/',pattern = '.gmt$')
immnu_patwhay.gene= as.data.frame(data.table::rbindlist(lapply(immnu_patwhays, function(x){
  print(x)
  a=clusterProfiler::read.gmt(paste0('immnue_pathway/',x))
  return(a)
})))
head(immnu_patwhay.gene)
immnu_patwhay.gene$term=gsub('KEGG_','',immnu_patwhay.gene$term)
immnu_patwhay.gene=split(x = immnu_patwhay.gene$gene,f = immnu_patwhay.gene$term)
immnu_patwhay.gene

#定义一个ssGSEA_score函数
ssGSEA_score<-function(dat, kcdf = c("Gaussian", "Poisson")[1], GeneSet =  NULL) {
  # dat：行为基因，列为样本
  # kcdf：log2(FPKM、TPM)数据使用 Gaussian；Counts 数据使用 Poisson
  library(GSVA)

  if (is.null(GeneSet)) {
    gene_set <- read.table(("PMID_28052254.txt.txt"), header = T, stringsAsFactors = F, sep = '\t')
    gene_set <- split(as.matrix(gene_set)[,1], gene_set[,2])
  } else {
    gene_set <- GeneSet
  }
  
ssgsea_res <- gsva(as.matrix(dat),
                     gene_set,
                     method='ssgsea',
                     kcdf=kcdf,
                     abs.ranking=TRUE)
  ssgsea_res <- as.data.frame(t(ssgsea_res))
  return(ssgsea_res)
}

#运行函数进行分析
imm.riskscore <- ssGSEA_score(dat = tcga.dat,
                              kcdf='Gaussian',
                              GeneSet = immnu_patwhay.gene)
head(imm.riskscore)

#定义diff_pathway函数计算p值
diff_pathway<-function(dat,group){
  dat=data.frame(cluster=group,t(dat))
  gr=as.character(unique(group))
  dat1=t(dat[dat$cluster==gr[1],-1])
  dat2=t(dat[dat$cluster==gr[2],-1])
  pathway=unique(c(rownames(dat1),rownames(dat2)))
  p_vale=data.frame()
  for (i in pathway){
    dd1=t.test(as.numeric(dat1[i,]),as.numeric(dat2[i,]))$p.value
    p_vale=rbind(p_vale,data.frame(pathway=i,p.value=dd1))
  }
  return(p_vale)
}

pathwy_p<-diff_pathway(dat=t(imm.riskscore[tcga.risk.cli$Samples,]),
                       group=tcga.risk.cli$group)
head(pathwy_p)
pathwy_p=pathwy_p[pathwy_p$p.value<0.05,]
pathwy_p$lab=paste0(pathwy_p$pathway,
                    ifelse(pathwy_p$p.value<0.001,'***',
                           ifelse(pathwy_p$p.value<0.01,'**',
                                  ifelse(pathwy_p$p.value<0.05,'*',''))))
pathwy_p$lab
imm.riskscore[1:4,1:4]
imm.riskscore1=merge(pathwy_p[,c("pathway","lab")],t(imm.riskscore),by.y=0,by.x='pathway')
head(imm.riskscore1)
rownames(imm.riskscore1)=imm.riskscore1$lab
imm.riskscore1=imm.riskscore1[,-c(1,2)]
imm.riskscore1[1:4,1:5]
library(pheatmap)
anno_col=data.frame(risk=tcga.risk.cli$group,row.names = tcga.risk.cli$Samples)
anno_col <- anno_col[order(anno_col$risk, decreasing = FALSE), , drop = F]  # Low在左
annotation_colors <- list(
  risk = c("Low" = "#3474b5", "High" = "#ef766d")
)
# 1. 高对比度颜色方案
my_colors <- colorRampPalette(c("#00ccff", "white", "#ff3300"))(100)  # 深蓝-白-鲜红
# 2. 调整颜色断点增强对比
breaks <- seq(-6, 6, length.out = 100)  # 根据实际数据范围调整
pdf('results/fig9b.pdf',height = 7,width = 15)
pheatmap(mat = imm.riskscore1[,rownames(anno_col)],scale = 'row',
         annotation_col = anno_col,
         annotation_colors = annotation_colors,  # 自定义注释颜色
         cluster_cols = F, 
         cluster_rows = F,
         show_rownames = T, 
         show_colnames = F,
         color = my_colors,
         breaks = breaks,
         border_color = NA)      # 去除单元格边框 
dev.off()


#CIBERSORT
# tcga.cibersort <- deconvo_tme(eset = tcga.dat,
#                               method = "cibersort",
#                               arrays = FALSE,
#                               perm = 100) #
# save(tcga.cibersort,file = 'results/tcga.cibersort.RData')

load('results/tcga.cibersort.RData')
tcga.cibersort=as.data.frame(tcga.cibersort)
tcga.cibersort[1:4,1:4]
tcga.cibersort=tcga.cibersort[,1:23]
vioplot_plot=function(rt,normal,tumor,normal.name,tumor.name){
  library(vioplot)
  rt=rt[c(normal,tumor),]
  normal=length(normal)
  tumor=length(tumor)
  par(las=1,mar=c(10,6,3,3))
  x=c(1:ncol(rt))
  y=c(1:ncol(rt))
  plot(x,y,
       xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
       main="",xlab="", ylab="Fraction",
       pch=21,
       col="white",
       xaxt="n")
  for(i in 1:ncol(rt)){
    if(sd(rt[1:normal,i])==0){
      rt[1,i]=0.001
    }
    if(sd(rt[(normal+1):(normal+tumor),i])==0){
      rt[(normal+1),i]=0.001
    }
    normalData=rt[1:normal,i]
    tumorData=rt[(normal+1):(normal+tumor),i]
    vioplot(normalData,at=3*(i-1),lty=1,add = T,col = '#3474b5')
    vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = '#ef766d')
    wilcoxTest=wilcox.test(normalData,tumorData)
    p=round(wilcoxTest$p.value,3)
    mx=max(c(normalData,tumorData))
    lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
    #text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.05, paste0("p=",p),''), cex = 0.8)
    text(x=3*(i-1)+0.5, y=mx+0.02, labels= paste0("p=",p), cex = 0.8)
  }
  text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
  legend(x = "topleft", box.col = "black",
         bg ="white", box.lwd = 0.5 , cex = 0.5,
         legend=c(normal.name, tumor.name), 
         fill = c("#3474b5","#ef766d"))
}
high.sample=tcga.risk.cli[which(tcga.risk.cli$group=='High'),"Samples"]
low.sample=tcga.risk.cli[which(tcga.risk.cli$group=='Low'),"Samples"]

rownames(tcga.cibersort)=tcga.cibersort$ID
colnames(tcga.cibersort)=gsub('_CIBERSORT','',colnames(tcga.cibersort))
pdf('results/fig9c.pdf',height = 8,width = 12)
vioplot_plot(rt = tcga.cibersort[,-1],
             normal = low.sample,
             tumor =high.sample,
             normal.name = 'Risk-Low',
             tumor.name = 'Risk-High')
dev.off()


#28种免疫细胞ssGSEA
#定义一个ssGSEA_score函数
ssGSEA_score<-function(dat, kcdf = c("Gaussian", "Poisson")[1], GeneSet =  NULL) {
  # dat：行为基因，列为样本
  # kcdf：log2(FPKM、TPM)数据使用 Gaussian；Counts 数据使用 Poisson
  library(GSVA)
  
  if (is.null(GeneSet)) {
    gene_set <- read.table(("PMID_28052254.txt.txt"), header = T, stringsAsFactors = F, sep = '\t')
    gene_set <- split(as.matrix(gene_set)[,1], gene_set[,2])
  } else {
    gene_set <- GeneSet
  }
  
  ssgsea_res <- gsva(as.matrix(dat),
                     gene_set,
                     method='ssgsea',
                     kcdf=kcdf,
                     abs.ranking=TRUE)
  ssgsea_res <- as.data.frame(t(ssgsea_res))
  return(ssgsea_res)
}
load("TCGA_exp.rda")
tcga.dat=TCGA
#tcga.dat行为基因名，列为样品名
#行列名统一
colnames(tcga.dat) <- sub(".01A$", "", colnames(tcga.dat))
colnames(tcga.dat) <- gsub("\\.", "-", colnames(tcga.dat))
#读取打分的分组文件
tcga.risk.cli=read.csv("TCGA_risk.csv",header = T,row.names = 1) 
rownames(tcga.risk.cli) <- sub("-01$", "", rownames(tcga.risk.cli))
tcga.risk.cli$Samples<-rownames(tcga.risk.cli)

immnu.28.gene=read.delim('PMID_28052254.txt.txt',sep='\t',header = T)

head(immnu.28.gene)
length(table(immnu.28.gene$Cell.type))
immnu.28.gene=split(immnu.28.gene$Metagene,immnu.28.gene$Cell.type)
tcga.28.score<-ssGSEA_score(dat = tcga.dat,
                            GeneSet = immnu.28.gene,
                            kcdf = c("Gaussian", "Poisson")[1])
head(tcga.28.score)
tcga.28.score.risk=cbind.data.frame(tcga.28.score[tcga.risk.cli$Samples,],
                                    risk=tcga.risk.cli$group)
head(tcga.28.score.risk)
write.table(tcga.28.score.risk,'results/tcga.28.score.risk.txt',quote = F,row.names = F,sep='\t')

#箱线图
library(reshape2)
library(ggpubr)
tcga.28.score.risk=reshape2::melt(tcga.28.score.risk)
head(tcga.28.score.risk)
tcga.28.score.risk$risk=factor(tcga.28.score.risk$risk, levels=c("Low", "High"))

fig9d<-ggboxplot(tcga.28.score.risk, x='variable', y='value', 
                  fill = "risk", color = "black",
                  palette = c("#3474b5","#ef766d"), 
                  ylab="Immune infiltration",xlab='',
                  add = "boxplot")+ 
  stat_compare_means(aes(group=risk),method = 'wilcox.test',
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))
fig9d
ggsave('results/fig9d.pdf',fig9d,height = 7,width = 12)
dev.off()


#22种免疫细胞评分与关键基因的相关性分析
#读取基因表达文件,并对数据进行处理
library(limma)
load("TCGA_exp.rda")
rt=TCGA
rt=as.matrix(rt)
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
exp <- as.data.frame(t(data))
exp <- exp[!grepl(".11A$", rownames(exp)), ]
rownames(exp) <- sub(".01A$", "", rownames(exp))
exp <- as.data.frame(t(exp)) 

#读取评分文件
risk=read.table("risk.train.txt", header=T, sep="\t", check.names=F,row.names=1)
head(risk)
riskGenes <- c("TAGLN2", "PCMT1", "PTMA", "CXCL13", "TUBA3D", "DCTPP1") 
tcga.gene.exp=t(exp[riskGenes,])
tcga.gene.exp[1:4,1:4]
rownames(tcga.gene.exp) <- gsub("\\.", "-", rownames(tcga.gene.exp))

#读取cibersort文件
load('results/tcga.cibersort.RData')

tcga.cibersort=as.data.frame(tcga.cibersort)
tcga.cibersort[1:4,1:4]
tcga.cibersort=tcga.cibersort[,1:23]
row.names(tcga.cibersort)<-tcga.cibersort[,1]
tcga.cibersort$ID<-NULL
colnames(tcga.cibersort) <- gsub("_CIBERSORT$", "", colnames(tcga.cibersort))

#查看相同样品名时一定要注意.和-
sameSample=intersect(row.names(tcga.cibersort), row.names(tcga.gene.exp))
tcga.cibersort=tcga.cibersort[sameSample,,drop=F]
tcga.gene.exp=tcga.gene.exp[sameSample,,drop=F]


outTab=data.frame()
for(immune in colnames(tcga.cibersort)){
  for(gene in colnames(tcga.gene.exp)){
    x = as.numeric(tcga.cibersort[, immune])
    y = as.numeric(tcga.gene.exp[, gene])
    corT = cor.test(x, y, method="pearson")  # Ensure "pearson" is used
    cor = corT$estimate
    pvalue = corT$p.value
    text = ifelse(pvalue < 0.001, "***", 
                  ifelse(pvalue < 0.01, "**", 
                         ifelse(pvalue < 0.05, "*", "")))
    outTab = rbind(outTab, cbind(Gene = gene, Immune = immune, cor, text, pvalue))
  }
}

outTab$cor=as.numeric(outTab$cor)
pdf(file="results/fig9e.pdf", width=10, height=10)
ggplot(outTab, aes(Gene, Immune)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  #scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  scale_fill_gradient2(low = "#3474b5", mid = "white", high = "#E64B35") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +    #ȥ?����?
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   #x??????
        axis.text.y = element_text(size = 10, face = "bold")) +       #y??????
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #????ͼ??
  scale_x_discrete(position = "bottom")      #X????????ʾλ??
dev.off()



#22种免疫细胞与风险评分的相关性分析
#读取打分的分组文件
risk=read.csv("TCGA_risk.csv",header = T,row.names = 1) 
rownames(risk) <- sub("-01$", "", rownames(risk))
risk$Samples<-rownames(risk)

load('results/tcga.cibersort.RData')

tcga.cibersort=as.data.frame(tcga.cibersort)
tcga.cibersort[1:4,1:4]
tcga.cibersort=tcga.cibersort[,1:23]
colnames(tcga.cibersort) <- gsub("_CIBERSORT$", "", colnames(tcga.cibersort))

risk=risk[tcga.cibersort$ID,]
#riskScore与免疫细胞相关性计算
imm.risk.ccor=Hmisc::rcorr(as.matrix(cbind.data.frame(tcga.cibersort[,-1], 
                                                      Riskscore=risk$RS)),
                           type = 'spearman')
imm.risk.ccor.r=reshape2::melt(imm.risk.ccor$r)
imm.risk.ccor.p=reshape2::melt(imm.risk.ccor$P)
imm.risk.ccor.r=imm.risk.ccor.r[which(imm.risk.ccor.r$Var1 =='Riskscore' & imm.risk.ccor.r$Var2 !='Riskscore'),]
imm.risk.ccor.p=imm.risk.ccor.p[which(imm.risk.ccor.p$Var1 =='Riskscore' & imm.risk.ccor.p$Var2 !='Riskscore'),]
#汇总
imm.risk.ccor.res=merge(imm.risk.ccor.r,
                        imm.risk.ccor.p,
                        by=c('Var1','Var2'))
head(imm.risk.ccor.res)
colnames(imm.risk.ccor.res)=c('Riskscore','Cell','cor','pvalue')
#定义圆圈颜色的函数
p.col = c('gold','pink','orange','#32cd32','darkgreen')

p.col = c('#77AECD','#066190','#024163','#D98380','#8E0F31')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],
                 ifelse(x>0.6,p.col[2],
                        ifelse(x>0.4,p.col[3],
                               ifelse(x>0.2,p.col[4], p.col[5])
                        )))
  return(color)
}
#定义设置圆圈大小的函数
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],
               ifelse(x<0.2,p.cex[2],
                      ifelse(x<0.3,p.cex[3],
                             ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}
#根据pvalue定义圆圈的颜色
points.color = fcolor(x=imm.risk.ccor.res$pvalue,p.col=p.col)
imm.risk.ccor.res$points.color = points.color
points.cex = fcex(x=imm.risk.ccor.res$cor)
imm.risk.ccor.res$points.cex = points.cex
imm.risk.ccor.res=imm.risk.ccor.res[order(imm.risk.ccor.res$cor),]
#进行可视化
xlim = ceiling(max(abs(imm.risk.ccor.res$cor))*10)/10         
pdf(file="results/fig9f2.pdf", width=9, height=7)      
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(imm.risk.ccor.res)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(imm.risk.ccor.res),col="white",lty=1,lwd=2)
#绘制图形的线段
segments(x0=imm.risk.ccor.res$cor,y0=1:nrow(imm.risk.ccor.res),x1=0,y1=1:nrow(imm.risk.ccor.res),lwd=4)
#绘制图形的圆圈
points(x=imm.risk.ccor.res$cor,y = 1:nrow(imm.risk.ccor.res),col = imm.risk.ccor.res$points.color,pch=16,cex=imm.risk.ccor.res$points.cex)
#展示免疫细胞的名称
text(par('usr')[1],1:nrow(imm.risk.ccor.res),imm.risk.ccor.res$Cell,adj=1,xpd=T,cex=1.5)
#展示pvalue
pvalue.text=ifelse(imm.risk.ccor.res$pvalue<0.001,'<0.001',sprintf("%.03f",imm.risk.ccor.res$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(imm.risk.ccor.res),pvalue.text,adj=0,xpd=T,col=ifelse(abs(imm.risk.ccor.res$cor)>redcutoff_cor & imm.risk.ccor.res$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)
#绘制圆圈大小的图例
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")
#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()



#22种免疫细胞批量做KM曲线，计算预后的关系
dir.create('results/KM')

load('results/tcga.cibersort.RData')

tcga.cibersort=as.data.frame(tcga.cibersort)
tcga.cibersort[1:4,1:4]
tcga.cibersort=tcga.cibersort[,1:23]
colnames(tcga.cibersort) <- gsub("_CIBERSORT$", "", colnames(tcga.cibersort))
row.names(tcga.cibersort)<-tcga.cibersort[,1]
tcga.cibersort$ID<-NULL

risk=read.csv("TCGA_risk.csv",header = T,row.names = 1) 
rownames(risk) <- sub("-01$", "", rownames(risk))
risk$Samples<-rownames(risk)

sameSample=intersect(row.names(tcga.cibersort), row.names(risk))
tcga.cibersort=tcga.cibersort[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
tcga.cibersort.os=cbind.data.frame(tcga.cibersort,
                                   risk[,c("OS","OS.time")])

bioSurvival=function(OS,OS.time,riskScore,title){
  dat=data.frame(OS=OS,OS.time=OS.time,riskScore=riskScore)
  library(survival)
  library(survminer)
  res.cut <- surv_cutpoint(dat,
                           time = "OS.time", 
                           event = "OS", 
                           variables = c("riskScore"))
  cut_va=as.numeric(res.cut$cutpoint[1])
  dat$risk=ifelse(dat$riskScore>=cut_va,'high','low')
  diff=survdiff(Surv(OS.time, OS) ~risk,data = dat)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS.time, OS) ~ risk, data = dat)
  surPlot=ggsurvplot(fit, 
                     data=dat,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title=title,
                     legend.labs=c("high", "low"),
                     xlab="Time(years)",
                     break.time.by = 2.5,
                     palette=c("#E64B35", "#3474b5"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.height=.25)
  print(cut_va)
  return(surPlot)
}
zero_ratio <- colSums(tcga.cibersort.os == 0) / nrow(tcga.cibersort.os)
valid_cols <- names(zero_ratio[zero_ratio < 0.9])  # Keep columns with <90% zeros

for (i in valid_cols) {
  riskScore <- as.numeric(tcga.cibersort.os[[i]])
  riskScore <- riskScore[!is.na(riskScore)]  # Remove NA values
  
  if (length(riskScore) < 10 || length(unique(riskScore)) < 5 || var(riskScore, na.rm = TRUE) < 1e-5) next
  
  tryCatch({
    pdf(paste0('results/KM/', i, '.pdf'), height = 7, width = 7, onefile = FALSE)
    print(bioSurvival(OS = tcga.cibersort.os$OS,
                      OS.time = tcga.cibersort.os$OS.time,
                      riskScore = riskScore, title = i))
    dev.off()
  }, error = function(e) {
    message("Skipping ", i, " due to error: ", e$message)
  })
}


#venn图
#计算p值
diff_pathway<-function(dat,group){
  dat=data.frame(cluster=group,t(dat))
  gr=as.character(unique(group))
  dat1=t(dat[dat$cluster==gr[1],-1])
  dat2=t(dat[dat$cluster==gr[2],-1])
  pathway=unique(c(rownames(dat1),rownames(dat2)))
  p_vale=data.frame()
  for (i in pathway){
    dd1=t.test(as.numeric(dat1[i,]),as.numeric(dat2[i,]))$p.value
    p_vale=rbind(p_vale,data.frame(pathway=i,p.value=dd1))
  }
  return(p_vale)
}

pathwy_p1<-diff_pathway(dat=t(tcga.cibersort),
                        group=risk$group)
diff=pathwy_p1[which(pathwy_p1$p.value<0.05),"pathway"]
diff
cor.pa=imm.risk.ccor.res[which(imm.risk.ccor.res$pvalue<0.05),"Cell"]
cor.pa
venn_plot<-function(data_list,fill=NULL,ven=T){
  crbind2DataFrame<-function(dat,full=F){
    print(class(dat))
    if(class(dat)=='table'){
      if(!is.na(ncol(dat))){
        dat=apply(dat,2,function(x){
          return(x)
        })
      }
    }
    if(class(dat)!='data.frame'){
      dat1=as.data.frame(as.matrix(dat))
    }else{
      dat1=dat
    }
    dat1.class=apply(dat1, 2, class)
    #which(dat1.class!='numeric')
    #print(head(dat1))
    for(i in which(dat1.class!='numeric')){
      dat1[,i]=as.character(dat1[,i])
      if(full){
        dat1[,i]=as.numeric(dat1[,i])
      }else{
        dt=dat1[which(gsub(' ','',dat1[,i])!=''&!is.na(dat1[,i])),i]
        dt=dt[which(dt!='Inf'&dt!='NaN'&dt!='NA')]
        #dt[which(is.na(as.numeric(dt)))]
        if(sum(is.na(as.numeric(dt)))<length(dt)*0.1){
          #print(dat1[,i])
          dat1[,i]=as.numeric(dat1[,i])
        }
      }
    }
    return(dat1)  
  }
  if(length(data_list)<=4&ven){
    library(VennDiagram)
    if(is.null(fill)){
      fill=c('red', 'blue','green','yellow')[1:length(data_list)]
    }
    g=venn.diagram(data_list, filename = NULL,margin = 0.2, 
                   fill = fill, alpha = 0.50, col = 'black', cex = 1, fontfamily = 'serif'
                   #,cat.col = c('black', 'black', 'black', 'black')
                   ,cat.cex = 1, cat.fontfamily = 'serif')
    grid.draw(g)
  }else if(length(data_list)<=6&ven){
    library(venn)
    if(is.null(fill)){
      fill=mg_colors[1:length(data_list)]
    }
    g=venn(data_list, zcolor = fill,ggplot=F,box=F) 
  }else{
    anm=unique(unlist(data_list))
    tbs=cbind()
    for(i in 1:length(data_list)){
      a1=rep(0,length(anm))
      a1[anm%in%data_list[[i]]]=1
      tbs=cbind(tbs,a1)
    }
    colnames(tbs)=names(data_list)
    row.names(tbs)=anm
    tbs=crbind2DataFrame(tbs)
    g=UpSetR::upset(tbs, nsets = ncol(tbs), nintersects = 30, mb.ratio = c(0.5, 0.5),
                    order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE)
    )
    return(g)
  }
}

#自行选择
surva=colnames(tcga.cibersort)[c(2,4,5,7,9,11)]


pdf('results/fig9g.pdf',height = 6,width = 6)
venn_plot(data_list = list(surv=surva,
                           cor=cor.pa,
                           diff=diff))
dev.off()


