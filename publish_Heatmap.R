library(pheatmap)
library(ComplexHeatmap)
library(corrplot)
library(preprocessCore)
data=read.table("SD28_LD3_LD7_LD28_DESEQ2_merged_count.txt")[,c(4:17)]
colnames(data)=c("SD28_1",
"SD28_2",
"SD28_3_1",
"SD28_3_2",
"LD3_1",
"LD3_2_1",
"LD3_2_2",
"LD3_3",
"LD7_1",
"LD7_2",
"LD7_3",
"LD28_1",
"LD28_2",
"LD28_3"
)
pheatmap(cor(data),cluster_rows=F,cluster_cols=F)
data1=normalize.quantiles(as.matrix(data))
colnames(data1)=c("SD28_1",
"SD28_2",
"SD28_3_1",
"SD28_3_2",
"LD3_1",
"LD3_2_1",
"LD3_2_2",
"LD3_3",
"LD7_1",
"LD7_2",
"LD7_3",
"LD28_1",
"LD28_2",
"LD28_3"
)
bk <- c(seq(0,1,by=0.01))
pheatmap(cor(data1),cluster_rows=F,cluster_cols=F,
color = c(colorRampPalette(colorRampPalette(colors = c("white","red"))(length(bk)),
legend_breaks=seq(0.8,1,0.01),
         breaks=bk
	)



##Corralation计算
	library(DESeq2)
	library(ggrepel)
	library(RColorBrewer)
	library(amap)
	library(gplots)
	library(ggplot2)
	data <- read.csv("SD28_LD3_LD7_LD28_DESEQ2_merged_count.txt", header=F, com='', quote='',
	     check.names=F, sep="\t")[,c(4:17)]
	rownames(data)=paste("peak_",c(1:dim(data)[1]),sep="")
	# 撇掉在多于两个样本中count<1的值，如果样本数多，这个数值可以适当增加
	# 排除极低表达基因的干扰
	data <- data[rowSums(data)>2,]
	head(data)
	colnames(data)=c("SD28_1",
	"SD28_2",
	"SD28_3_1",
	"SD28_3_2",
	"LD3_1",
	"LD3_2_1",
	"LD3_2_2",
	"LD3_3",
	"LD7_1",
	"LD7_2",
	"LD7_3",
	"LD28_1",
	"LD28_2",
	"LD28_3"
	)
	# 读入分组信息
	sample <- read.table("group_info.txt", header=T,sep="\t")

	# sample <- sample[match(colnames(data), rownames(sample)),, drop=F]
	# sample_rowname <- rownames(sample)

	# # 下面的可以忽略，如果没遇到错误不需要执行
	# # 目的是做因子转换
	# sample <- data.frame(lapply(sample, function(x) factor(x, levels=unique(x))))
	# rownames(sample) <- sample_rowname
	# sample
	ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data,
	        colData = sample,  design= ~ conditions)

	dds <- DESeq(ddsFullCountTable)

	## estimating size factors

	## estimating dispersions

	## gene-wise dispersion estimates

	## mean-dispersion relationship

	## final dispersion estimates

	## fitting model and testing
	# normalized=T, 返回标准化的数据
	normalized_counts <- counts(dds, normalized=TRUE)
	head(normalized_counts)

	#根据基因在不同的样本中表达变化的差异程度mad值对数据排序，差异越大的基因排位越前。
	normalized_counts_mad <- apply(normalized_counts, 1, mad)
	normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

	# 标准化后的数据输出
	#write.table(normalized_counts, file="ehbio_trans.Count_matrix.xls.DESeq2.normalized.xls",
	#quote=F, sep="\t", row.names=T, col.names=T)

	# 只在Linux下能运行，目的是填补表格左上角内容
	#system(paste("sed -i '1 s/^/ID\t/'", "ehbio_trans.Count_matrix.xls.DESeq2.normalized.xls"))

	# log转换后的结果
	rld <- rlog(dds, blind=FALSE)
	rlogMat <- assay(rld)
	rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]


	#write.table(rlogMat, file="ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.xls",
	#quote=F, sep="\t", row.names=T, col.names=T)

	# 只在Linux下能运行，目的是填补表格左上角内容
	#system(paste("sed -i '1 s/^/ID\t/'", "ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.xls"))


	#/样品层级聚类分析，判断样品的相似性和组间组内差异
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
	pdf("quail_ATAC.DESeq2.normalized.rlog.pearson.pdf", pointsize=10)
	heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
	col=hmcol, dendrogram="none",margins=c(11,11), main="The pearson correlation",,rowV="as-is",rowInd=c("SD28_1",
	"SD28_2",
	"SD28_3_1",
	"SD28_3_2",
	"LD3_1",
	"LD3_2_1",
	"LD3_2_2",
	"LD3_3",
	"LD7_1",
	"LD7_2",
	"LD7_3",
	"LD28_1",
	"LD28_2",
	"LD28_3"
	))

	heatmap.2(pearson_cor, symm=T, trace="none",Rowv=FALSE,
	col=hmcol, dendrogram="row",margins=c(11,11), main="The pearson correlation")
	dev.off()
##depth
	depth=c(33899777,
	37317363,
	26709212,
	48605397,
	 46400040,
	 39183220,
	38174251,
	 43225089,
	46522112,
	35639439,
	29704280,
	51899068,
	 32220173,
	36360567)
##DESEQ2
	library(DESeq2)
	countsTable <- read.table( "SD28_LD3_LD7_LD28_DESEQ2_merged_count.txt", stringsAsFactors=TRUE )
	rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
	colnames(countsTable)<-c("chr","start","end","SD28_1","SD28_2","SD28_3","SD28_4","LD3_1","LD3_2","LD3_3","LD3_4","LD7_1","LD7_2","LD7_3","LD28_1","LD28_2","LD28_3")
	colData <- data.frame(condition=factor(c(rep("SD28",4),rep("LD3",4),rep("LD7",3),rep("LD28",3))))
	#countsTable_select=countsTable[,c(8,9,10,11,12,13,14)]
	dds<-DESeqDataSetFromMatrix(countsTable[,c(4:17)],colData, formula(~condition)) 
	dds <- estimateSizeFactors(dds)
	dds<- estimateDispersions(dds)
	dds<- DESeq(dds)
##差异分析
	library(DESeq2)
	library(ggplot2)
	my_peak=c()
	countsTable <- read.table( "SD28_LD3_LD7_LD28_DESEQ2_merged_count.txt", stringsAsFactors=TRUE )
	rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
	colnames(countsTable)<-c("chr","start","end","SD28_1",
															"SD28_2",
															"SD28_3_1",
															"SD28_3_2",
															"LD3_1",
															"LD3_2_1",
															"LD3_2_2",
															"LD3_3",
															"LD7_1",
															"LD7_2",
															"LD7_3",
															"LD28_1",
															"LD28_2",
															"LD28_3")
	colData <- data.frame(condition=factor(c("SD28","SD28","SD28","SD28",
		"LD3","LD3","LD3","LD3","LD7","LD7","LD7","LD28","LD28","LD28")))
	dds<-DESeqDataSetFromMatrix(countsTable[,c(4:17)],colData, formula(~condition)) 
	dds <- estimateSizeFactors(dds)
	dds<- estimateDispersions(dds)
	dds<- DESeq(dds)
	normalized_counts <- counts(dds, normalized=TRUE)
	rld <- rlog(dds, blind = FALSE)
    #返回样本名和分组
    pcaData <- plotPCA(rld, intgroup=c("condition"),returnData = T)
	res=results(dds,contrast=c("condition","LD3","SD28"))
	summary(res)
	res <- res[order(res$padj),]
	all_info=merge(countsTable,as.data.frame(res),by="row.names",all.x=TRUE)[,c(2,3,4,19,20,21,22,23,24)]
	head(all_info)
	write.table(all_info,"Information_SD28_LD3_DESEQ2.txt",quote=F,sep="\t",row.names=F,col.names=F)
	diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange >  0.8| log2FoldChange < -0.8))
	write.table(cbind(countsTable[rownames(as.data.frame(diff_gene_deseq2)),c(1,2,3)],
		as.data.frame(diff_gene_deseq2)),
		"SD28_LD3_DESEQ2.txt",quote=F,sep="\t",row.names=F,col.names=F)
	my_peak=c(my_peak,rownames(diff_gene_deseq2))

	res=results(dds,contrast=c("condition","LD7","SD28"))
	summary(res)
	res <- res[order(res$padj),]
	all_info=merge(countsTable,as.data.frame(res),by="row.names",all.x=TRUE)[,c(2,3,4,19,20,21,22,23,24)]
	head(all_info)
	write.table(all_info,"Information_SD28_LD7_DESEQ2.txt",quote=F,sep="\t",row.names=F,col.names=F)	
	diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange >  0.8| log2FoldChange < -0.8))
	write.table(cbind(countsTable[rownames(as.data.frame(diff_gene_deseq2)),c(1,2,3)],
		as.data.frame(diff_gene_deseq2 )),
		"SD28_LD7_DESEQ2.txt",quote=F,sep="\t",row.names=F,col.names=F)
	my_peak=c(my_peak,rownames(diff_gene_deseq2))

	res=results(dds,contrast=c("condition","LD28","SD28"))
	summary(res)
	res <- res[order(res$padj),]
	all_info=merge(countsTable,as.data.frame(res),by="row.names",all.x=TRUE)[,c(2,3,4,19,20,21,22,23,24)]
	head(all_info)
	write.table(all_info,"Information_SD28_LD28_DESEQ2.txt",quote=F,sep="\t",row.names=F,col.names=F)
	diff_gene_deseq2 <-subset(res,padj < 0.05& (log2FoldChange >  0.8| log2FoldChange < -0.8))
	write.table(cbind(countsTable[rownames(as.data.frame(diff_gene_deseq2)),c(1,2,3)],
		as.data.frame(diff_gene_deseq2 )),
		"SD28_LD28_DESEQ2.txt",quote=F,sep="\t",row.names=F,col.names=F)
	my_peak=c(my_peak,rownames(diff_gene_deseq2))

	res=results(dds,contrast=c("condition","LD7","LD3"))
	summary(res)
	res <- res[order(res$padj),]
	diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange >  0.8| log2FoldChange < -0.8))
	write.table(cbind(countsTable[rownames(as.data.frame(diff_gene_deseq2)),c(1,2,3)],
		as.data.frame(diff_gene_deseq2)),
		"LD3_LD7_DESEQ2.txt",quote=F,sep="\t",row.names=F,col.names=F)
	my_peak=c(my_peak,rownames(diff_gene_deseq2))

	res=results(dds,contrast=c("condition","LD28","LD3"))
	summary(res)
	res <- res[order(res$padj),]
	diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange >  0.8| log2FoldChange < -0.8))
	write.table(cbind(countsTable[rownames(as.data.frame(diff_gene_deseq2)),c(1,2,3)],
		as.data.frame(diff_gene_deseq2)),
		"LD3_LD28_DESEQ2.txt",quote=F,sep="\t",row.names=F,col.names=F)
	my_peak=c(my_peak,rownames(diff_gene_deseq2))


	res=results(dds,contrast=c("condition","LD28","LD7"))
	summary(res)
	res <- res[order(res$padj),]
	diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange >  0.8| log2FoldChange < -0.8))
	my_peak=c(my_peak,rownames(diff_gene_deseq2))
	write.table(cbind(countsTable[rownames(as.data.frame(diff_gene_deseq2)),c(1,2,3)],
		as.data.frame(diff_gene_deseq2 )),
		"LD7_LD28_DESEQ2.txt",quote=F,sep="\t",row.names=F,col.names=F)


##热图绘制
	library(ComplexHeatmap)
	pdf("Heatmap_kmeans.pdf")
	# set.seed(111)
	# 	p=Heatmap(matrix=t(scale(t(normalized_counts[unique(my_peak),]))),
	# 	             show_row_names=F,
	# 	             show_row_dend = F,
	# 	            row_labels = F,
	# 	             row_km=4,
	# 	             cluster_columns=F,
	# 	             column_title_gp = gpar(fontsize=15, fontface="bold"),
	# 	             show_parent_dend_line = FALSE,
	# 	             name = "Z-score")
	# 	p
	# set.seed(222)
	# 	p=Heatmap(matrix=t(scale(t(normalized_counts[unique(my_peak),]))),
	# 	             show_row_names=F,
	# 	             show_row_dend = F,
	# 	             row_labels = F,
	# 	             row_km=4,
	# 	             cluster_columns=F,
	# 	             column_title_gp = gpar(fontsize=15, fontface="bold"),
	# 	             show_parent_dend_line = FALSE,
	# 	             name = "Z-score")
	# 	p
	# set.seed(0)
	# 	p=Heatmap(matrix=t(scale(t(normalized_counts[unique(my_peak),]))),
	# 	             show_row_names=F,
	# 	             show_row_dend = F,
	# 	            row_labels = F,
	# 	             row_km=4,
	# 	             cluster_columns=F,
	# 	             column_title_gp = gpar(fontsize=15, fontface="bold"),
	# 	             show_parent_dend_line = FALSE,
	# 	             name = "Z-score")
	# 	p
	# set.seed(1)
	# 	p=Heatmap(matrix=t(scale(t(normalized_counts[unique(my_peak),]))),
	# 	             show_row_names=F,
	# 	             show_row_dend = F,
	# 	            row_labels = F,
	# 	             row_km=4,
	# 	             cluster_columns=F,
	# 	             column_title_gp = gpar(fontsize=15, fontface="bold"),
	# 	             show_parent_dend_line = FALSE,
	# 	             name = "Z-score")
	# 	p
	# set.seed(2)
	# 	p=Heatmap(matrix=t(scale(t(normalized_counts[unique(my_peak),]))),
	# 	             show_row_names=F,
	# 	             show_row_dend = F,
	# 	            row_labels = F,
	# 	             row_km=4,
	# 	             cluster_columns=F,
	# 	             column_title_gp = gpar(fontsize=15, fontface="bold"),
	# 	             show_parent_dend_line = FALSE,
	# 	             name = "Z-score")
	# 	p
	# set.seed(3)
	# 	p=Heatmap(matrix=t(scale(t(normalized_counts[unique(my_peak),]))),
	# 	             show_row_names=F,
	# 	             show_row_dend = F,
	# 	            row_labels = F,
	# 	             row_km=4,
	# 	             cluster_columns=F,
	# 	             column_title_gp = gpar(fontsize=15, fontface="bold"),
	# 	             show_parent_dend_line = FALSE,
	# 	             name = "Z-score")
			
	# 	p
	# set.seed(123)
	# 	p=Heatmap(matrix=t(scale(t(normalized_counts[unique(my_peak),]))),
	# 	             show_row_names=F,
	# 	             show_row_dend = F,
	# 	            row_labels = F,
	# 	             row_km=4,
	# 	             cluster_columns=F,
	# 	             column_title_gp = gpar(fontsize=15, fontface="bold"),
	# 	             show_parent_dend_line = FALSE,
	# 	             name = "Z-score")
			
	# 	p
	# set.seed(321)
	# 	p=Heatmap(matrix=t(scale(t(normalized_counts[unique(my_peak),]))),
	# 	             show_row_names=F,
	# 	             show_row_dend = F,
	# 	            row_labels = F,
	# 	             row_km=4,
	# 	             column_split = factor(c("SD28", "SD28","SD28","SD28","LD3","LD3","LD3","LD3","LD7","LD7","LD7","LD28","LD28","LD28"),levels=c("SD28","LD3","LD7","LD28")),
	# 	             cluster_columns=F,
	# 	             column_title_gp = gpar(fontsize=15, fontface="bold"),
	# 	             show_parent_dend_line = FALSE,
	# 	             name = "Z-score")
			
	# 	p

	# dev.off()
##提取信息
	set.seed(321)
	list<-draw(Heatmap(matrix=t(scale(t(normalized_counts[unique(my_peak),]))),
	             show_row_names=F,
	             show_row_dend = F,
	            row_labels = F,
	             row_km=4,
	             cluster_columns=F,
	             column_title_gp = gpar(fontsize=15, fontface="bold"),
	             show_parent_dend_line = FALSE,
	             name = "Z-score"))
	row_order(list)
	diff_peak=countsTable[unique(my_peak),]
	rownames(diff_peak)=c(1:dim(diff_peak)[1])
	cluster1 = diff_peak[row_order(list)[[1]],column_order(list)] #提取第一簇的相关信息
	cluster2 = diff_peak[row_order(list)[[2]],column_order(list)] #同上
	cluster3 = diff_peak[row_order(list)[[3]],column_order(list)] #同上
	cluster4 = diff_peak[row_order(list)[[4]],column_order(list)] #同上

	rownames(cluster1)<-row_order(list)[[1]]
	rownames(cluster2)<-row_order(list)[[2]]
	rownames(cluster3)<-row_order(list)[[3]]
	rownames(cluster4)<-row_order(list)[[4]]

	cluster1=cluster1[order(cluster1$chr,cluster1$start),]
	cluster2=cluster2[order(cluster2$chr,cluster2$start),]
	cluster3=cluster3[order(cluster3$chr,cluster3$start),]
	cluster4=cluster4[order(cluster4$chr,cluster4$start),]                                                                                                                                                                                                                                                                 

	dim(cluster1)
	dim(cluster2)
	dim(cluster3)
	dim(cluster4)
	write.table(cluster1[,c(1,2,3)],"kmeans/Kmeans_cluster1_DESEQ2.txt"
		,quote=F,sep="\t",row.names=F,col.names=F)
	write.table(cluster2[,c(1,2,3)],"kmeans/Kmeans_cluster2_DESEQ2.txt"
		,quote=F,sep="\t",row.names=F,col.names=F)
	write.table(cluster3[,c(1,2,3)],"kmeans/Kmeans_cluster3_DESEQ2.txt"
		,quote=F,sep="\t",row.names=F,col.names=F)
	write.table(cluster4[,c(1,2,3)],"kmeans/Kmeans_cluster4_DESEQ2.txt"
		,quote=F,sep="\t",row.names=F,col.names=F)


##火山图绘制 
	res=results(dds,contrast=c("condition","LD3","SD28"))
	summary(res)
	res <- res[order(res$padj),]
	diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange >  0.8| log2FoldChange < -0.8))
	my_peak=rownames(diff_gene_deseq2)
	res[which(res$log2FoldChange >= 0.8 & res$padj < 0.05),'sig'] <- 'UP'
	res[which(res$log2FoldChange <= -0.8 & res$padj< 0.05),'sig'] <- 'DOWN'
	res[which(res$padj > 0.05),'sig'] <- 'NO DIFF'
	this_tile <- paste0("LD3 vs SD28",
	    '\ncutoff for abs(logFC) and FDR is 0.8 and 0.05',
	    '\nThe number of up gene is ',nrow(data[data$sig =='UP',]) ,
	    '\nThe number of down gene is ',nrow(data[data$sig =='DOWN',]))
	volcano <- ggplot(as.data.frame(res), aes(log2FoldChange, -log(padj, 10))) +
	  ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
	  geom_point(aes(color = sig), alpha = 0.8, size = 2) +
	  scale_color_manual(values = c('#000080', '#C0C0C0', '#8B0000')) +
	  theme(panel.grid = element_blank(), 
	      panel.background = element_rect(color = 'black', fill = 'transparent'), 
	      ) +
	  theme(legend.title = element_blank(), 
	        legend.key = element_rect(fill = 'transparent'), 
	        legend.background = element_rect(fill = 'transparent')) +
	  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) +
	  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
	  labs(x = 'log2 Fold Change', y = '-log10 padj', color = NA)
	volcano
##Homer
for i in {1..4};do
findMotifsGenome.pl \
	/md01/changjy/data/all/quail_diff/kmeans/Kmeans_cluster$i\_DESEQ2.txt \
	/md01/changjy/data/quail/ref/ncbi_dataset/GCF_001577835.2_Coturnix_japonica_2.1_genomic.fna \
	/md01/changjy/data/all/quail_diff/TF/cluster$i/   -p 40
done


library(ggrepel)
ggplot(data = pcaData,aes(x=PC1,y=PC2,,colour=condition))+
geom_point(size=3)+
geom_text_repel(aes(label = name), size = 3)


setwd("/md01/changjy/data/all/quail_diff/TF/cluster1/knownResults")
for (i in  c(1:10)){
	input=paste("known",i,".motif",sep="")
	output=paste("/md01/changjy/rgtdata/motifs/quail/","C1_",input,sep="")
	data=read.table(input,skip=1)
	change=t(data)*1000
	print(change)
	write.table(change,output,sep="\t",quote=F,col.names=F,row.names=F)}

setwd("/md01/changjy/data/all/quail_diff/TF/cluster2/knownResults")
for (i in  c(1:10)){
	input=paste("known",i,".motif",sep="")
	output=paste("/md01/changjy/rgtdata/motifs/quail/","C2_",input,sep="")
	data=read.table(input,skip=1)
	change=t(data)*1000
	print(change)
	write.table(change,output,sep="\t",quote=F,col.names=F,row.names=F)}


setwd("/md01/changjy/data/all/quail_diff/TF/cluster3/knownResults")
for (i in  c(1:10)){
	input=paste("known",i,".motif",sep="")
	output=paste("/md01/changjy/rgtdata/motifs/quail/","C3_",input,sep="")
	data=read.table(input,skip=1)
	change=t(data)*1000
	print(change)
	write.table(change,output,sep="\t",quote=F,col.names=F,row.names=F)}


setwd("/md01/changjy/data/all/quail_diff/TF/cluster4/knownResults")
for (i in  c(1:10)){
	input=paste("known",i,".motif",sep="")
	output=paste("/md01/changjy/rgtdata/motifs/quail/","C4_",input,sep="")
	data=read.table(input,skip=1)
	change=t(data)*1000
	print(change)
	write.table(change,output,sep="\t",quote=F,col.names=F,row.names=F)}

 	
	data=read.table("/md01/changjy/data/quail/ref/ncbi_dataset/GCF_001577835.2_Coturnix_japonica_2.1_genomic.gtf",sep="\t")
 	data=data[data$V3=="gene",]
	b <- data %>% separate(V9, c("gene_id","other"), "[;]")
	c <-b %>% separate(gene_id, c("gene_id","name"), "[ ]")
	d=c[,c(1,4,5,6,7,10)]
	write.table(d,"gencode.v22.annotation.gtf.gtf",sep="\t",quote=F,row.names=F,col.names=F)
	write.table(d[,c(1,2,3,)],"gencode.v22.annotation.gtf.gtf",sep="\t",quote=F,row.names=F,col.names=F)



ls /md01/changjy/data/all/quail_diff/bam/*.bam|while read bam; do rgt-hint footprinting --organism cj22 --atac-seq $bam /md01/changjy/data/all/quail_diff/SD28_LD3_LD7_LD28_peak_merged.txt  --output-location /md01/changjy/data/all/quail_diff/footprint/ --output-prefix=${bam%.bam*} --paired-end; done
rgt-motifanalysis matching --organism=cj22 --input-files SD28.bed LD3.bed LD7.bed LD28.bed --motif-dbs /md01/changjy/rgtdata/motifs/quail/
cd match/
rgt-hint  differential --organism=cj22 --bc --nc 20 --mpbs-files=SD28_mpbs.bed,LD3_mpbs.bed,LD7_mpbs.bed,LD28_mpbs.bed  --reads-files=/md01/changjy/data/all/quail_diff/bam/SD28.bam,/md01/changjy/data/all/quail_diff/bam/LD3.bam,/md01/changjy/data/all/quail_diff/bam/LD7.bam,/md01/changjy/data/all/quail_diff/bam/LD28.bam --conditions=SD28,LD3,LD7,LD28  --output-location=./


library(DESeq2)
library(ggplot2)
library(ggrepel)
countsTable <- read.table( "/md01/changjy/data/quail/RNAseq/Count_matrix.txt", 
	stringsAsFactors=TRUE,row.names=1 )[,c(4,5,6,7,8)]

colData <- data.frame(condition=factor(c("SD28","SD28","SD28","LD7","LD7")))
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds)
# res <- lfcShrink(dds,
#     contrast = c('condition','LD7','SD28'), res=res, type = 'normal')
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))
write.table(as.data.frame(diff_gene_deseq2),"/md01/changjy/data/quail/RNAseq/SD28_LD7_UP_DEP.txt", row.name=T,sep="\t",quote=F)
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange < -1 ))
write.table(as.data.frame(diff_gene_deseq2),"/md01/changjy/data/quail/RNAseq/SD28_LD7_DOWN_DEP.txt", row.name=T,sep="\t",quote=F)


res=as.data.frame(res)
res=na.omit(res)
res[which(res$log2FoldChange >= 1 & res$padj < 0.05),'sig'] <- 'UP'
res[which(res$log2FoldChange <= -1 & res$padj< 0.05),'sig'] <- 'DOWN'
res[which(!(res$sig%in%c("UP","DOWN"))),'sig'] <- 'NO DIFF'


res$label <-''
res[res$padj < 0.01 & res$log2FoldChange >= 2,]
res[rownames(res)%in%c("DIO3","GPR20","TSHB","CGA","AGRP","SLC16A2","GHRH","DIO2","POMC","PER2"),]$label=c("DIO3","GPR20","TSHB","CGA","AGRP","SLC16A2","GHRH","DIO2","POMC","PER2")
this_tile <- paste0('cutoff for abs(logFC) and FDR is 1 and 0.05',
    '\nThe number of up gene is ',nrow(res[res$sig =='UP',]) ,
    '\nThe number of down gene is ',nrow(res[res$sig =='DOWN',]))
volcano <- ggplot(res, aes(log2FoldChange, -log(padj, 10))) +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
  geom_point(aes(color = sig), alpha = 0.8, size = 2) +
  scale_color_manual(values = c('#000080', '#C0C0C0', '#8B0000')) +
  theme(panel.grid = element_blank(), 
      panel.background = element_rect(color = 'black', fill = 'transparent'), 
      ) +
  #geom_text(aes(label = label), size = 3,vjust=-0.5, alpha=0.8) +  
  theme(legend.title = element_blank(), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 padj', color = NA) +
  xlim(-8, 8)

 p=volcano+ geom_text_repel(data = res, aes(x = res$log2FoldChange, 
                                      y = -log10(res$padj), 
                                      label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black")
 p
 ggsave("/md01/changjy/data/quail/RNAseq/SD28_LD7_volcano.pdf",dpi=900)



library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
#     V1       V2
# SD28 vs LD3  32
# SD28 vs LD7 210
# SD28 vs LD28 509
# LD3 vs LD7  22
# LD3 vs LD28  58
# LD7 vs LD28   2
col=brewer.pal(6,"Dark2")
#col[4]="#99d8c9"
data=read.csv("/md01/changjy/data/all/quail_diff/DEP_number.csv",na.strings=" ",header=F)
data$V1=c("SD28 vs LD3","SD28 vs LD7","SD28 vs LD28","LD3 vs LD7","LD3 vs LD28","LD7 vs LD28")
colnames(data)=c("conditions","values")

 p<-ggplot(data,aes(x = factor(conditions,levels = c("SD28 vs LD3","SD28 vs LD7","SD28 vs LD28","LD3 vs LD7","LD3 vs LD28","LD7 vs LD28"),ordered=T),
                               y =values, 
                               fill=factor(conditions,levels = c("SD28 vs LD3","SD28 vs LD7","SD28 vs LD28","LD3 vs LD7","LD3 vs LD28","LD7 vs LD28"),ordered=T)))+
  geom_bar(stat = "identity" ,width=0.6,position=position_dodge(width = 0.8))+
  geom_text(aes(label=values),position=position_dodge(width = 0.6),size = 6,vjust = -0.25)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.key = element_blank(),legend.title = element_blank(),
        axis.text.x = element_text(size=12,angle=60, vjust=0.5,hjust=0.5,
        	)
  )+  xlab("")+ylab("Number of DARs")+ scale_fill_npg()
p
ggsave("/md01/changjy/data/all/quail_diff/DARs_number.pdf")

conda activate base
nohup Rscript ArchR_1E_ME.sh >>scATAC_1E_ME.log 2>&1 &
nohup Rscript ArchR_2_ME.sh >>scATAC_2_ME.log 2>&1 &
nohup Rscript ArchR_3E_ME.sh >>scATAC_3E_ME.log 2>&1 &
nohup Rscript ArchR_4ME.sh >>scATAC_4ME.log 2>&1 &



nohup Rscript ArchR_1P_plus.sh >>scATAC_1P_plus.log 2>&1 &
nohup Rscript ArchR_2_P.sh >>scATAC_2_P.log 2>&1 &
nohup Rscript ArchR_3E_P.sh >>scATAC_3E_P.log 2>&1 &
nohup Rscript ArchR_4P.sh >>scATAC_4P.log 2>&1 &



library(ArchR)
projHemeNew<-loadArchRProject("/md01/changjy/data/Goose/snATAC_1P_plus")
AddPrefix=FALSE
projHemeNew <- addDoubletScores(
    input = projHemeNew,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1)
saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)
projHemeNew <- addIterativeLSI(ArchRProj = projHemeNew,
                            useMatrix = "TileMatrix",
                            name = "IterativeLSI",
                            iterations = 2,
                            clusterParams = list(resolution = c(0.2),
                                                 sampleCells =13117,
                                                 n.start = 10),
                            varFeatures = 15000,
                            dimsToUse = 1:25,
                            force = TRUE,
                            threads=40
)


chr1=read.table("/md01/changjy/data/Goose/1E-ME/outs/peak_annotation_chrName.txt",header=T)
chr2=read.table("/md01/changjy/data/Goose/2-ME/outs/peak_annotation_chrName.txt",header=T)
chr3=read.table("/md01/changjy/data/Goose/3E-ME/outs/peak_annotation_chrName.txt",header=T)
chr4=read.table("/md01/changjy/data/Goose/4ME/outs/peak_annotation_chrName.txt",header=T)
aaa=intersect(chr1$chrom,chr2$chrom)
bbb=intersect(chr3$chrom,chr4$chrom)
ccc=intersect(aaa,bbb)
setdiff(chromsize$V1,ccc)



ha = HeatmapAnnotation(df = data.frame("Group" = condition, 
	show_annotation_name = T,   
	col = list(Group = c("Squamous cell carcinoma, NOS" =  "red","Adenocarcinoma, NOS" = "blue")),
 	show_legend = T,   
 	annotation_name_side = "left")
 



 plot.atac <- assay(esca.rse)[result$FDR < fdr.cut.off & abs(result$ESCC_minus_ESAD) > diff.cut.off,]
 col <- colorRamp2(seq(min(plot.atac), max(plot.atac),                       by = (max(plot.atac) - min(plot.atac))/99), pal_atac)
 rows.annot <- rowAnnotation(foo = anno_mark(at = c(1,18), labels = rownames(plot.atac)[c(1,18)]))
 ht_list <-   Heatmap(plot.atac,          name = "ATAC-seq log2(counts)",           col = col,          column_names_gp = gpar(fontsize = 8),          show_column_names = F,          heatmap_legend_param = list(legend_direction = "horizontal",                                      labels_gp = gpar(fontsize = 12),                                       title_gp = gpar(fontsize = 12)),          show_row_names = FALSE,          cluster_columns = TRUE,          use_raster = TRUE,          raster_device = c("png"),          raster_quality = 2,          cluster_rows = T,          right_annotation = rows.annot,          row_title = paste0(sum(result$FDR < fdr.cut.off &                                    abs(result$ESCC_minus_ESAD) > diff.cut.off),                             " ATAC-seq peaks"),          #column_order = cols.order,          row_names_gp = gpar(fontsize = 4),          top_annotation = ha,          #width = unit(15, "cm"),          #column_title = paste0("RNA-seq z-score (n = ", ncol(plot.exp),")"),           column_title_gp = gpar(fontsize = 12),           row_title_gp = gpar(fontsize = 12)) 
 draw(ht_list,newpage = TRUE,      column_title = paste0("ATAC-seq ESCC vs ESAD (FDR < ", fdr.cut.off,                           ",  Diff mean log2 Count > ",diff.cut.off,")"),     column_title_gp = gpar(fontsize = 12, fontface = "bold"),     heatmap_legend_side = "bottom",     annotation_legend_side = "right")





 library(BSgenome.goose.ENSEMBLE.ac22)
 eval(parse(text = "BSgenome.goose.ENSEMBLE.ac22"))
 