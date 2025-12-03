


options(stringsAsFactors = F) 
{ library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(stringr)
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(SingleR)
  library(CCA)
  library(clustree)
  library(cowplot)
  library(monocle)
  library(tidyverse)
  library(SCpubr)
  library(UCell)
  library(GSVA)
  library(GSEABase)
  library(harmony)
  library(plyr)
  library(randomcoloR)
  library(CellChat)
  library(ggpubr)
  library(DoubletFinder) #双细胞检测使用
  library(stringr)  
  library(tidydr)
  library(Nebulosa)
  mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
               '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
               '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
  )
}

sce=readRDS("GSE161529_raw.rds")
#计算线粒体基因比例mito
mito_genes=rownames(sce)[grep("^MT-", rownames(sce),ignore.case = T)] 
print(mito_genes) #可能是13个线粒体基因
sce=PercentageFeatureSet(sce, features = mito_genes, col.name = "percent_mito")

#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
p1=VlnPlot(sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1 
w=length(unique(sce$orig.ident))/3+5;w
ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height = 5)
p2=FeatureScatter(sce, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p2
ggsave(filename="Scatterplot.pdf",plot=p2)

#进行数据过滤
sce.filt <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA< 7000
                   &percent_mito < 20)
#查看过滤后的情况
dim(sce)
dim(sce.filt)
table(sce.filt$tissue)
sce<-sce.filt
###############################################################################                  
#step3: harmony整合多个单细胞样品
###############################################################################
setwd('../')
dir.create("./2-harmony")
setwd("./2-harmony")
#再次进行标准化,使用LogNormalize方法
sce <- NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 10000) 
## 鉴定高变基因
sce <- FindVariableFeatures(sce, nfeatures = 3000) #高变基因数目根据参考文献进行修改
## 归一化
sce <- ScaleData(sce)


###harmony矫正   
library(harmony)
seuratObj <- RunHarmony(sce,group.by.vars = "orig.ident")
names(seuratObj@reductions)
##########矫正后结果可视化
pdf(file = "03.harmony.pdf",width =7.5,height = 5.5)
DimPlot(seuratObj, reduction = "harmony",pt.size = 1)+
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank())
dev.off()
seuratObj <- Seurat::RunTSNE(seuratObj,dims = 1:30,reduction ='harmony')
pdf(file = "03.tsne_harmony.pdf",width =7.5,height = 5.5)
DimPlot(seuratObj, reduction = "tsne",pt.size = 1)+
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank())
dev.off()
seuratObj <- RunUMAP(seuratObj,  dims = 1:30, reduction = "harmony")
pdf(file = "03.umap_harmony.pdf",width =7.5,height = 5.5)
DimPlot(seuratObj, reduction = "umap",pt.size = 1)+
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank())
dev.off()
saveRDS(seuratObj, "sce_harmony.rds")


##确定使用PC个数
sce=seuratObj
# each PC essentially representing a ‘metafeature’
ElbowPlot(sce,ndims = 50,reduction = "pca")  #看Standard Deviation到第几位开始幅度变化不明显来确定

dev.off()

##########0.3.细胞聚类##########################
#sce<-readRDS("sce_harmony.rds")
#确定PC
PC=30   ### PC聚类的值根据前面04.jackstrawplot.pdf确定P值，04.ElbowPlot.pdf看Standard Deviation到第几位开始幅度变化不明显来确定，
##对细胞聚类
#设置不同的分辨率resolution，观察分群效果，dim为PCA选择的主成分数
sce <- FindNeighbors(sce, dims = 1:PC, reduction = "harmony") 

sce.all=sce  #先复制一个数据集进行测试分辨率
#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all=FindClusters(sce.all, resolution = res, algorithm = 1)}
apply(sce.all@meta.data[,grep("RNA_snn_res",colnames(sce.all@meta.data))],2,table)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"), DimPlot(sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low_tsne.pdf",width = 14,height = 8)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.5") + 
                   ggtitle("louvain_0.5"), DimPlot(sce.all, reduction = "tsne", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all, reduction = "tsne", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high_tsne.pdf",width = 14,height = 8)

p2_tree<-clustree(sce.all@meta.data, prefix = "RNA_snn_res.")
pdf(file = "04.clustertree.pdf",width =12,height =10)
p2_tree
dev.off()
saveRDS(sce.all, "sce.all_int.rds")

#接下来分析，按照分辨率为1进行 
sel.clust = "RNA_snn_res.0.5"
sce.all <- SetIdent(sce.all, value = sel.clust)
# 结果聚成几类，用Idents查看
length(levels(Idents(sce.all)))
# 查看每个类别多少个细胞
head(sce.all@meta.data)
table(sce.all@meta.data$seurat_clusters)
table(sce.all@active.ident) 
# 可视化UMAP/tSNE
p1<-DimPlot(sce.all, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank())
ggsave(filename="05-cluster.UMAP.pdf",plot=p1,width = 8,height = 7)
p1
dev.off()
pdf(file = "05-cluster.TSEN.pdf",width =7,height = 5.5)
p2<-DimPlot(sce.all, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank())
p2
dev.off()



#################################  细胞注释  #####################

setwd('../')
dir.create("3-cell")  #创建文件夹，已经有文件夹的直接跳过进行setwd
setwd("3-cell")  
getwd()
mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)
#自动注释
library(SingleR)
library(celldex)
### 注释数据库加载 
load("BlueprintEncode_bpe.se_human.RData")
load("HumanPrimaryCellAtlas_hpca.se_human.RData")
# hpca.se <- celldex::HumanPrimaryCellAtlasData() 
# bpe.se<-celldex::BlueprintEncodeData() 
sce<-sce.all
#使用BP和HPCA两个数据库综合注释，使用list函数读入多个数据库。
anno <- SingleR(sce@assays$RNA@data,
                ref = list(BP = bpe.se, HPCA = hpca.se),
                labels = list(bpe.se$label.main, hpca.se$label.main),
                clusters = sce@meta.data$seurat_clusters)
sce@meta.data$singleR_label <- unlist(lapply(sce@meta.data$seurat_clusters, function(x){anno$labels[x]}))
DimPlot(sce,reduction = 'umap',cols = mycolors,raster=FALSE,pt.size = 0.4,
        group.by = 'singleR_label',label = T,label.box = T)+
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank())
ggsave('singR_markers_umap.pdf',width = 18,height = 14)
DimPlot(sce,reduction = 'tsne',cols = mycolors,pt.size = 0.4,
        group.by = 'singleR_label',label = T,label.box = T)+
  theme_dr(xlength = 0.22, ylength = 0.22, arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))+
  theme(panel.grid = element_blank())
ggsave('singR_markers_tsne.pdf',width = 18,height = 14)



