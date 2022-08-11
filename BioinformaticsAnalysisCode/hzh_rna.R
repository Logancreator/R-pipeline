library(DESeq2)
library(gplots)
library(RColorBrewer)
setwd("/public/home/changjianye/project/hzh/rawData/rnaData/rnaRes/bam")
countData <- read.table("all_Featurecounts",row.names=1,header=T)[,c(6:32)]
countData <- countData[,c(-7,-14,-21,-22)]
countData <- countData[rowMeans(countData)>1,]
condition <- factor(c(rep("GS2_1",3),rep("GS2_3",3),rep("GS2_gw8_1",3),rep("GS2_gw8_3",3),rep("gw8_1",3),rep("gw8_3",3),rep("W0_1",3),rep("W0_3",2)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
head(dds)
dds <- estimateSizeFactors(dds)
raw <- SummarizedExperiment(counts(dds, normalized=FALSE),
                               colData=colData(dds))
nor <- SummarizedExperiment(counts(dds, normalized=TRUE),
                               colData=colData(dds))
vsd <- vst(dds)
rld <- rlog(dds)
pdf("PCA.pdf")
plotPCA( DESeqTransform(raw), intgroup=c("condition") )
plotPCA( DESeqTransform(nor), intgroup=c("condition") )
plotPCA(vsd, intgroup=c("condition"))
plotPCA(rld, intgroup=c("condition"))
dev.off()



ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countData,
colData = colData, design= ~ condition)
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
# 生成颜色
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# 计算相关性pearson correlation
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))

# 层级聚类
hc <- hcluster(t(rlogMat), method="pearson")

# 热图绘制
## 在命令行下运行时，需要去除下面开头的#号，把输出的图保存到文件中
## 输出到文件时，dev.off()命令是关闭输出，从而完成图片的写入。如果不做这一步，则图片则不能写入，从而不能打开
## 如果在Rstudio或其它可视化界面时，可以直接把图输出到屏幕
pdf("qc_pca.pdf")
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
col=hmcol, margins=c(11,11), main="pearson correlation of each
sample")

#样品PCA分析
pca_data <- plotPCA(rld, intgroup=c("condition"), returnData=T, ntop=5000)
pca_data 
dev.off()

