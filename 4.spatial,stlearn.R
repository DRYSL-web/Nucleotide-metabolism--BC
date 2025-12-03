## 加载包
library(RColorBrewer) 
library(viridis)
library(wesanderson)

library(Seurat)
library(ggplot2)

col_vector=c('#1F77B4','#F87F13','#359C62','#D32929','#69308E','#8C564C','#F33CA9',
  '#F5E801','#08F7F0','#AFC7E6',"#FFFF99")

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector
col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 12, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:12]
dev.off()

# doi.org/10.1158/2159-8290.CD-21-0316
#读取数据-------------
stRNA <- Load10X_Spatial(data.dir ="./GSM6177603//", 
                         filename = "GSM6177603_NYU_BRCA2_Vis_processed_filtered_feature_bc_matrix.h5",
                         slice ="spatial")
library(Seurat)
library(patchwork)
library(dplyr)


pbc_1 <- SpatialFeaturePlot(stRNA, images = "spatial", 
                            features = NULL, alpha = 0) + 
  theme(legend.position = "none") + 
  labs(x = "BRCA4", title = "BRCA4") +   # 设置标题
  theme(plot.title.position = "plot",   # 将标题放置在图形正上方
        plot.title = element_text(hjust = 0.5, size = 20)) +  # 控制标题的大小和位置
  labs(x = "BRCA4")
# 如果你希望标题更大，可以调整标题的大小
pbc_1
ggsave(pbc_1,filename='FigureBRCA4HE.pdf',height = 8,width = 8)
dev.off()

plot1 <- VlnPlot(stRNA, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(stRNA, features = "nCount_Spatial") + theme(legend.position = "right")
plot3 <- wrap_plots(plot1, plot2)
plot3
ggsave(plot3,filename='fig1.pdf',height = 6,width = 8)
dev.off()


stRNA[["percent.mt"]] <- PercentageFeatureSet(stRNA, pattern = "^MT[-]")
plot1 <- VlnPlot(stRNA, features = "percent.mt", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(stRNA, features = "percent.mt") + theme(legend.position = "right")
plot3 <-wrap_plots(plot1, plot2)
plot3
ggsave(plot3,filename='fig2.pdf',height = 6,width = 8)
dev.off()

mt_genes <- row.names(stRNA)[grepl("^MT-", row.names(stRNA))]
rps_genes <- row.names(stRNA)[grepl("^RPS", row.names(stRNA))]
mrp_genes <- row.names(stRNA)[grepl("^MRP", row.names(stRNA))]
rpl_genes <- row.names(stRNA)[grepl("^RPL", row.names(stRNA))]
rb_genes <- c(rps_genes, mrp_genes, rpl_genes)

stRNA <- stRNA[!rownames(stRNA) %in% c(rb_genes, mt_genes), ]

#Genes expressed in less that 10 spots are filtered
stRNA <- stRNA[rowSums(GetAssayData(stRNA, assay = "Spatial") > 0) > 10, ]

## 降维聚类-----
stRNA <- RunPCA(stRNA, assay = "SCT", verbose = T)
ElbowPlot(stRNA, ndims=30, reduction="pca") 
dev.off()

stRNA <- FindNeighbors(stRNA, reduction = "pca", dims = 1:20)
stRNA <- FindClusters(stRNA, verbose = T)
stRNA <- RunUMAP(stRNA, reduction = "pca", dims = 1:20)

## 可视化------
p1 <- DimPlot(stRNA, reduction = "umap", label = TRUE,group.by = 'seurat_clusters',pt.size = 1,cols = col_vector[1:20])
p1
ggsave(p1,filename='fig3.pdf',height = 6,width = 6)
dev.off()

p2 <- SpatialDimPlot(stRNA,label = T)
p2
ggsave(p2,filename='fig4.pdf',height = 6,width = 6)
dev.off()

## 保存
saveRDS(stRNA,file ='stRNA_BRCA4.RDS')


### stlearn------------------
# 见Python------------------------


## 观察基因表达
stRNA=readRDS('./stRNA_BRCA4.RDS')
gene=read.table("NMRGs.txt", header=T, sep="\t", check.names=F)
top_10_genes <- gene$ID[10:20]

p1<-DotPlot(stRNA,features = top_10_genes)+ RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="")
p1
ggsave(p1,filename='fig5.pdf',height = 6,width = 6)
dev.off()


### 代谢评估
stRNA@assays$RNA=stRNA@assays$SCT

library(scMetabolism)
library(ggplot2)
library(rsvd)

stRNA<-sc.metabolism.Seurat(obj = stRNA, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "REACTOME")


input.pathway <- rownames(stRNA@assays[["METABOLISM"]][["score"]])[1:30]
DotPlot.metabolism(obj =stRNA,
                   pathway = input.pathway, phenotype = "seurat_clusters", norm = "y")


a=grep('N-',rownames(stRNA@assays[["METABOLISM"]][["score"]]))

stRNA$Nucleotides_Metabolism=as.numeric(stRNA@assays$METABOLISM$score['Metabolism of nucleotides',])

p1<-SpatialFeaturePlot(stRNA,features = 'Nucleotides_Metabolism')
p1
ggsave(p1,filename='FigureNM.pdf',height = 6,width = 6)
dev.off()
saveRDS(stRNA,file ='stRNA_BRCA4.RDS')
