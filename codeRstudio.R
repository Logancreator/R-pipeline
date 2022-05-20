
    library("ggplot2")
    library("ggh4x")
    library("cowplot")
    data<-read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\cluster_motif\\cluster_motif_top15.csv",
                   head=T)

    colnames(data)=c("TF","p","target","backgroud","enrichment","class")
    data$TF=factor(data$TF,levels=c("RFX(HTH)","RFX2(HTH)","RFX1(HTH)","X-BOX(HTH)","RFX5(HTH)","BAPX1(HOMEOBOX)","NKX2-2(HOMEOBOX)",
                                    "NKX2-1(HOMEOBOX)","NKX2-5(HOMEOBOX)","NKX3-1(HOMEOBOX)","SOX10(HMG)","SOX3(HMG)","SOX6(HMG)","SOX21(HMG)","MEF2C(MADS)",
                                    "MEF2A(MADS)","MEF2B(MADS)","MEF2D(MADS)","EKLF(ZF)","KLF3(ZF)","SPIB(ETS)",
                                    "PGR(NR)","GRE(NR)v1","GRE(NR)v2","ARE(NR)","THRA(NR)v2","THRB(NR)v3","JUN-AP1(BZIP)","HLF(BZIP)","NFIL3(BZIP)","PR(NR)","AR-HALFSITE(NR)","SPL9(SBP)","SPL11(SBP)",
                                    "LXRE(NR)","RORGT(NR)","THRA(NR)v1","THRB(NR)v1","THRB(NR)v2","ELF-1(E2F)","SOX15(HMG)","TGA1(BZIP)","EAR2(NR)","COUP-TFII(NR)"))
    p<-ggplot(data,aes(x=class,y=TF))+
      geom_point(aes(size=`enrichment`,color=-log(p)),pch=16,alpha = 0.6)+
      theme_bw()+
      scale_size_continuous(range=c(2,8))+
      theme(panel.grid = element_blank(), axis.text.x=element_text(angle=30,hjust = 0.5,vjust=0))+
      labs(x=NULL,y=NULL)+
      scale_colour_gradient(low = 'blue', high = 'red')+
      scale_y_discrete(position = "left")+ #y轴文字放左侧
      scale_x_discrete(name=NULL,position = "top",expand = c(0, 1))+
      theme(text =element_text(size=12, color="black", family = "sans"),
            axis.ticks = element_blank(), axis.line = element_blank(),
            axis.text.x=element_text(size=10, angle = 45, hjust=0, color="black", family="sans"),
            axis.text.y=element_text(size=10, family="sans", color="black"))
      
        # theme(legend.position = "top",
        #     legend.text=element_text(size=12, family="sans"),
        #     legend.title=element_text(size=12, family= "sans"),
        #     legend.key.width=unit(0.5,"inches"),
        #     legend.background = element_rect(fill="white", color="white"),
        #     panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        #     legend.key = element_rect(fill="white"))+rremove("ylab") 
    p
    # guide_fill <- get_legend(p+ guides(color = "none") + theme(legend.position = "bottom"))
    # plot_grid(p +
    #             guides(fill = "none") + 
    #             theme(legend.position = "bottom"), 
    #           guide_fill, nrow = 2, rel_heights = c(10, 1))
    ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\test.pdf",dpi=800,height = 10,width = 6)

    
    
    

    
library("Biostrings")    
seqdata <- read.fasta('/md01/changjy/data/quail/ref/ncbi_dataset/GCF_001577835.2_Coturnix_japonica_2.1_rna_from_genomic.fna')
seqdata
df=read.csv("/md01/changjy/data/all/quail_diff/kmeans/Kmeans_cluster1_DESEQ2.txt",sep="\t",header=F)
df$PositionID=rownames(df)
data=read.csv("/md01/changjy/data/all/quail_diff/TF/cluster1/position/motif1.position",header=T,sep="\t")
data=merge(df,data,by="PositionID")
which(data$V1=="NC_029545.1")


for (i in c(1:dim(data)[1])){
  hind <- DNAString(data[i,6])
  
}
  


library(DESeq2)
library(ggplot2)
library(ggrepel)
countsTable <- read.table( "C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD28\\Count_matrix.txt", 
                           stringsAsFactors=TRUE,row.names=1 )[,c(4,5,6,7,8,9,10,11)]
head(countsTable)
colData <- data.frame(condition=factor(c("SD28","SD28","SD28","SD28","LD28","LD28","LD28","LD28")))
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds,contrast = c("condition","LD28","SD28"))
summary(res)
res <- res[order(res$padj),]
#diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))
#write.table(as.data.frame(diff_gene_deseq2),"C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD28\\SD28_LD28_UP_DEP.txt", row.name=T,sep="\t",quote=F)
#diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange < -1 ))
#write.table(as.data.frame(diff_gene_deseq2),"C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD28\\SD28_LD28_DOWN_DEP.txt", row.name=T,sep="\t",quote=F)
write.table(res,"C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD28\\All_results.csv", row.name=T,quote=F,sep=",")



rm(list=ls())
res=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD28\\All_results.csv",row.names = 1)
res=as.data.frame(res)
res=na.omit(res)
res[which(res$log2FoldChange >= 1 & res$padj < 0.05),'sig'] <- 'UP'
res[which(res$log2FoldChange <= -1 & res$padj< 0.05),'sig'] <- 'DOWN'
res[which(!(res$sig%in%c("UP","DOWN"))),'sig'] <- 'NO DIFF'


res$label <-''
res[res$padj < 0.01 & res$log2FoldChange >= 2,]
res[rownames(res)%in%c("DIO3","KLF9","TSHB","DIO2"),]$label=c("DIO3","KLF9","TSHB","DIO2")
this_tile <- paste0('cutoff for abs(logFC) and FDR is 1 and 0.05',
                    '\nThe number of up gene is ',nrow(res[res$sig =='UP',]) ,
                    '\nThe number of down gene is ',nrow(res[res$sig =='DOWN',]))
volcano <- ggplot(res, aes(log2FoldChange, -log(padj, 10))) +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
  geom_point(aes(color = sig), alpha = 0.6, size = 2) +
  #scale_color_manual(values = c('#000080', '#C0C0C0', '#8B0000')) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
  ) +
  geom_text(aes(label = label), size = 3,vjust=2) +  
  theme(legend.title = element_blank(), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05,10), color = 'gray', size = 0.25)+                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ) +
  labs(x = 'log2 Fold Change', y = '-log10 padj')+xlim(-8,8)
volcano

ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD28\\SD28_LD28_volcano.pdf",dpi=900)





library(ggplot2)
library(dplyr)
atac_data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Quadrant\\Information_SD28_LD28_DESEQ2_anno.txt",sep="\t",header=F)
gene=as.character(unique(atac_data$V18))
atac=data.frame(target="",atac_log2FoldChange="",atac_p="",atac_padj="")
atac[,"target"]=as.character(atac[,"target"])
atac[,"atac_log2FoldChange"]=as.character(atac[,"atac_log2FoldChange"])
atac[,"atac_p"]=as.character(atac[,"atac_p"])
atac[,"atac_padj"]=as.character(atac[,"atac_padj"])
for (i in c(1:length(gene))){
  tmp=atac_data[atac_data$V18==gene[i],c(7,10,11)]
  sum_data=sum(tmp[,1])
  atac[i,"target"]=gene[i]
  atac[i,"atac_log2FoldChange"]=sum_data
  atac[i,"atac_p"]=min(tmp[,2])
  atac[i,"atac_padj"]=min(tmp[,3])
}
rownames(atac)=atac[,1]
atac=atac[,-1]
RNA_data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Quadrant\\All_results_SD28_LD28.csv",row.names=1)
RNA=RNA_data[gene,]
colnames(RNA)[2]="RNA_log2FoldChange"
res=cbind(atac,RNA)
res$label <- case_when(abs(as.numeric(res$RNA_log2FoldChange)) >= 1 & abs(as.numeric(res$atac_log2FoldChange)) >= 1 ~ "part1379",
                       abs(as.numeric(res$atac_log2FoldChange)) < 1 & abs(as.numeric(res$RNA_log2FoldChange)) > 1 ~ "part28",
                       abs(as.numeric(res$atac_log2FoldChange)) > 1 & abs(as.numeric(res$RNA_log2FoldChange)) < 1 ~ "part46",
                       abs(as.numeric(res$atac_log2FoldChange)) < 1 & abs(as.numeric(res$RNA_log2FoldChange)) < 1 ~ "part5")          

res$RNA_log2FoldChange=as.numeric(res$RNA_log2FoldChange)
res$atac_log2FoldChange=as.numeric(res$atac_log2FoldChange)
res$atac_p=as.numeric(res$atac_p)
res$atac_padj=as.numeric(res$atac_padj)
str(res)
#开始尝试绘图；
p0 <-ggplot(res,aes(atac_log2FoldChange,RNA_log2FoldChange,color=label))
#添加散点；
p1 <- p0+geom_point(size=2)+guides(color="none")
p1

##自定义半透明颜色
mycolor <- c("#FF9999","#99CC00","#c77cff","gray80")
p2 <- p1 + scale_colour_manual(name="",values=alpha(mycolor,0.7))
p2  

#添加辅助线；
p3 <- p2+geom_hline(yintercept = c(-1,1),
                    size = 0.8,
                    color = "grey40",
                    lty = "dashed")+
  geom_vline(xintercept = c(-1,1),
             size = 0.8,
             color = "grey40",
             lty = "dashed")
p3

#调整横轴和纵轴绘图区域的范围；
#设置y轴范围（上下两端的空白区域设为1），修改刻度标签；
#expansion函数设置坐标轴范围两端空白区域的大小；mult为“倍数”模式，add为“加性”模式；
p4<-p3+
  scale_y_continuous(expand=expansion(add = c(0.5, 0.5)),
                     limits = c(-8, 8),
                     breaks = c(-6,-3,0,3,6),
                     label = c("-6","-3","0","3","6"))+
  scale_x_continuous(expand=expansion(add = c(0.5, 0.5)),
                     limits = c(-8, 8),
                     breaks = c(-6,-3,0,3,6),
                     label = c("-6","-3","0","3","6"))
p4


#自定义图表主题，对图表做精细调整；
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
#隐藏纵轴，并对字体样式、坐标轴的粗细、颜色、刻度长度进行限定；
mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="gray30",size = 15),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
#应用自定义主题；
p5 <- p4+mytheme
p5+ggtitle("SD28_LD7_Quadrant")


# cor = cor(res[,c("RNA_log2FoldChange","atac_log2FoldChange")],use="complete.obs")
# 
# cor=round(cor[2,1],digits = 2)
# #准备作为图形的标题;
# lab = paste("correlation=",cor,sep="")
# lab
# #[1] "correlation=0.35"
# #在图上添加文字标签；
# p5+geom_text(x = -2, y = 5, label = lab)
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Quadrant\\SD28_LD28_quadrant.pdf")


boxplot(c(1.892413159,
          0.443164693,
          2.001325987))


dat = read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\GO\\GO_plot.txt",fileEncoding="UCS-2LE",sep="\t")
library(ggplot2)#没有自己安装 install.package("ggplot2")
library(dplyr)
colnames(dat)
#按照gene number从低到高排序[升序]
dat<- arrange(dat,dat[,5])
#Pathway列最好转化成因子型，否则作图时ggplot2会将所有Pathway按字母顺序重排序
#将Pathway列转化为因子型
dat_cluster1=dat[dat$cluster=="cluster1",][c(1:10),]
dat_cluster2=dat[dat$cluster=="cluster2",][c(1:10),]
dat_cluster3=dat[dat$cluster=="cluster3",][c(1:10),]
dat_cluster4=dat[dat$cluster=="cluster4",][c(1:10),]
dat_cluster1$Description <- factor(dat_cluster1$Description,levels = unique(dat_cluster1$Description))
dat_cluster2$Description <- factor(dat_cluster2$Description,levels = unique(dat_cluster2$Description))
dat_cluster3$Description <- factor(dat_cluster3$Description,levels = unique(dat_cluster3$Description))
dat_cluster4$Description <- factor(dat_cluster4$Description,levels = unique(dat_cluster4$Description))
tmp=rbind(dat_cluster4,dat_cluster3)
tmp=rbind(tmp,dat_cluster2)
tmp=rbind(tmp,dat_cluster1)
tmp$LogP=-tmp$LogP




p <- ggplot(tmp,aes(y=Description,x=cluster,color=LogP,size=GeneInGOAndHitList)) + 
  geom_point() +
  theme_bw() +
  scale_color_gradient(low='blue',high='Red')+ 
  theme(plot.title = element_text(hjust = 0.8),
        strip.text.y = element_text(size = 8),
        legend.position="right",
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        axis.text.x = element_text(size=8,vjust = 0.1,angle =  30),
        axis.text.y = element_text(size=11),
        axis.title.x =NULL ,
        axis.title.y =NULL )+
    scale_y_discrete(position = "right")
p
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\GO_cluster.pdf",height = 9,width=9)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
BiocManager::install(c("TFBSTools"))
#install.packages("Rcpp")

library(TFBSTools)
library(Biostrings)
library(ggplot2)
rm(list=ls())
data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\library_info.csv")
data$class="ATAC"
ggplot(data = data) + geom_boxplot(aes(x =class, y = mapping.ratio., 
                                              fill = class))

library(ggplot2)
library(pheatmap)
data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\wsm.csv",header=T,row.names = 1)
#colnames(data)=c("19.6°C","17.8°C","14.2°C")
data
data1=t(data)
windowsFonts(Times_New_Roman=windowsFont("Times New Roman"))

element_text(family='Times_New_Roman', size = 13, face='bold')
pdf("C:\\Users\\Jeff\\OneDrive\\桌面\\wsm.pdf")
annotation_col = data.frame(
  treatment = factor(c("CCP","CCP","CCP", "Nylon","Nylon","Nylon")), 
  time = factor(c("3d","7d","14d"),levels =c("3d","7d","14d") )
)
rownames(annotation_col) = colnames(data1)
head(annotation_col)



p=pheatmap(data1,scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         annotation_col = annotation_col,
         border_color = "black",
         angle_col = 45,
         main="Genes Expression Heatmap",
         cellwidth = 20, cellheight = 16, 
         fontsize = 8,#fontfamily= "Times New Roman",
         color = colorRampPalette(c("lightblue", "white", "red"))(50))

p
dev.off()

BiocManager::install("sangerseqR")
library(sangerseqR)
library(stringr)
seq<-readsangerseq("C:\\Users\\Jeff\\OneDrive\\桌面\\2-5AOX1_D08.ab1")
seq<-makeBaseCalls(seq,ratio=0.2)#将高度比（某峰比上该处最高峰）的值大于0.2，认为是真正的信号
peakAmp<-peakAmpMatrix(seq)
peakAmp<-as.data.frame(peakAmp)
colnames(peakAmp)<-c("A","C","G","T")

peakAmp$ratio<-apply(peakAmp,1,function(x){a=sort(x,decreasing = T);a[2]/a[1]})
peakAmp$sig<-ifelse(peakAmp$ratio>0.2,T,F)
peakAmp$primaryseq<-t(str_split(primarySeq(seq,string = T),pattern = "",simplify = T))
peakAmp$secondary<-t(str_split(secondarySeq(seq,string = T),pattern = "",simplify = T))

pdf("C:\\Users\\Jeff\\OneDrive\\桌面\\sanger_2-5AOX1_D08_statics.pdf")
chromatogram(seq)
dev.off()
write.csv(peakAmp,"C:\\Users\\Jeff\\OneDrive\\桌面\\sanger_2-5AOX1_D08_statics.csv",quote=F)







library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gg.gap)
library(ggplot2)
library(patchwork)
change=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\SD_steroid(all).csv")
#change$hormone=rownames(change)
#rownames(change)=c(1:dim(change)[1])

tmp1=data.frame(group=rep("SD28d",17),value=t(change[1,c(2:18)]))
tmp1$hormone=rownames(tmp1)
rownames(tmp1)=c(1:dim(tmp1)[1])
tmp2=data.frame(group=rep("SD28d",17),value=t(change[2,c(2:18)]))
tmp2$hormone=rownames(tmp2)
rownames(tmp2)=c(1:dim(tmp2)[1])
tmp3=data.frame(group=rep("SD28d",17),value=t(change[3,c(2:18)]))
tmp3$hormone=rownames(tmp3)
rownames(tmp3)=c(1:dim(tmp3)[1])
tmp4=data.frame(group=rep("LD1d",17),value=t(change[4,c(2:18)]))
tmp4$hormone=rownames(tmp4)
rownames(tmp4)=c(1:dim(tmp4)[1])
tmp5=data.frame(group=rep("LD1d",17),value=t(change[5,c(2:18)]))
tmp5$hormone=rownames(tmp5)
rownames(tmp5)=c(1:dim(tmp5)[1])
tmp6=data.frame(group=rep("LD1d",17),value=t(change[6,c(2:18)]))
tmp6$hormone=rownames(tmp6)
rownames(tmp6)=c(1:dim(tmp6)[1])
tmp7=data.frame(group=rep("LD2d",17),value=t(change[7,c(2:18)]))
tmp7$hormone=rownames(tmp7)
rownames(tmp7)=c(1:dim(tmp7)[1])
tmp8=data.frame(group=rep("LD2d",17),value=t(change[8,c(2:18)]))
tmp8$hormone=rownames(tmp8)
rownames(tmp8)=c(1:dim(tmp8)[1])
tmp9=data.frame(group=rep("LD2d",17),value=t(change[9,c(2:18)]))
tmp9$hormone=rownames(tmp9)
rownames(tmp9)=c(1:dim(tmp9)[1])
tmp10=data.frame(group=rep("LD3d",17),value=t(change[10,c(2:18)]))
tmp10$hormone=rownames(tmp10)
rownames(tmp10)=c(1:dim(tmp10)[1])
tmp11=data.frame(group=rep("LD3d",17),value=t(change[11,c(2:18)]))
tmp11$hormone=rownames(tmp11)
rownames(tmp11)=c(1:dim(tmp11)[1])
tmp12=data.frame(group=rep("LD3d",17),value=t(change[12,c(2:18)]))
tmp12$hormone=rownames(tmp12)
rownames(tmp12)=c(1:dim(tmp12)[1])
tmp13=data.frame(group=rep("LD7d",17),value=t(change[13,c(2:18)]))
tmp13$hormone=rownames(tmp13)
rownames(tmp13)=c(1:dim(tmp13)[1])
tmp14=data.frame(group=rep("LD7d",17),value=t(change[14,c(2:18)]))
tmp14$hormone=rownames(tmp14)
rownames(tmp14)=c(1:dim(tmp14)[1])
tmp15=data.frame(group=rep("LD7d",17),value=t(change[15,c(2:18)]))
tmp15$hormone=rownames(tmp15)
rownames(tmp15)=c(1:dim(tmp15)[1])
tmp16=data.frame(group=rep("LD14d",17),value=t(change[16,c(2:18)]))
tmp16$hormone=rownames(tmp16)
rownames(tmp16)=c(1:dim(tmp16)[1])
tmp17=data.frame(group=rep("LD14d",17),value=t(change[17,c(2:18)]))
tmp17$hormone=rownames(tmp17)
rownames(tmp17)=c(1:dim(tmp17)[1])
tmp18=data.frame(group=rep("LD14d",17),value=t(change[18,c(2:18)]))
tmp18$hormone=rownames(tmp18)
rownames(tmp18)=c(1:dim(tmp18)[1])
tmp19=data.frame(group=rep("LD28d",17),value=t(change[19,c(2:18)]))
tmp19$hormone=rownames(tmp19)
rownames(tmp19)=c(1:dim(tmp19)[1])
tmp20=data.frame(group=rep("LD28d",17),value=t(change[20,c(2:18)]))
tmp20$hormone=rownames(tmp20)
rownames(tmp20)=c(1:dim(tmp20)[1])
tmp21=data.frame(group=rep("LD28d",17),value=t(change[21,c(2:18)]))
tmp21$hormone=rownames(tmp21)
rownames(tmp21)=c(1:dim(tmp21)[1])
colnames(tmp1)=c("group","value","hormone")
colnames(tmp2)=c("group","value","hormone")
colnames(tmp3)=c("group","value","hormone")
colnames(tmp4)=c("group","value","hormone")
colnames(tmp5)=c("group","value","hormone")
colnames(tmp6)=c("group","value","hormone")
colnames(tmp7)=c("group","value","hormone")
colnames(tmp8)=c("group","value","hormone")
colnames(tmp9)=c("group","value","hormone")
colnames(tmp10)=c("group","value","hormone")
colnames(tmp11)=c("group","value","hormone")
colnames(tmp12)=c("group","value","hormone")
colnames(tmp13)=c("group","value","hormone")
colnames(tmp14)=c("group","value","hormone")
colnames(tmp15)=c("group","value","hormone")
colnames(tmp16)=c("group","value","hormone")
colnames(tmp17)=c("group","value","hormone")
colnames(tmp18)=c("group","value","hormone")
colnames(tmp19)=c("group","value","hormone")
colnames(tmp20)=c("group","value","hormone")
colnames(tmp21)=c("group","value","hormone")

tmp=rbind(tmp1,tmp2)
tmp=rbind(tmp,tmp3)
tmp=rbind(tmp,tmp4)
tmp=rbind(tmp,tmp5)
tmp=rbind(tmp,tmp6)
tmp=rbind(tmp,tmp7)
tmp=rbind(tmp,tmp8)
tmp=rbind(tmp,tmp9)
tmp=rbind(tmp,tmp10)
tmp=rbind(tmp,tmp11)
tmp=rbind(tmp,tmp12)
tmp=rbind(tmp,tmp13)
tmp=rbind(tmp,tmp14)
tmp=rbind(tmp,tmp15)
tmp=rbind(tmp,tmp16)
tmp=rbind(tmp,tmp17)
tmp=rbind(tmp,tmp18)
tmp=rbind(tmp,tmp19)
tmp=rbind(tmp,tmp20)
tmp=rbind(tmp,tmp21)
#tmp$hormone=row.names(tmp)
install.packages("gg.gap")
library(ggplot2)
library(reshape2)
library(RColorBrewer)

## 载入数据
df <- tmp[tmp$group%in%c("SD28d","LD3d","LD7d","LD28d"),]
#df <- melt(df, id="Species", variable.name="Attribute", value.name = "Size")
#mycol= brewer.pal(n = 12, name = "Set3")

## 统计 3种鸢尾花形态数据数据均值、标准差、标准误
mean <- aggregate(df$value, by=list(df$group, df$hormone), FUN=mean)
sd <- aggregate(df$value, by=list(df$group, df$hormone), FUN=sd)
len <- aggregate(df$value, by=list(df$group, df$hormone), FUN=length)
df_res <- data.frame(mean, sd=sd$x, len=len$x)
colnames(df_res) = c("Species", "Attribute", "Mean", "Sd", "Count")
str(df_res)
df_res$Se <- df_res$Sd/sqrt(df_res$Count) ### 计算标准差
df_res$Attribute


df_res$Species=factor(df_res$Species,levels = c("SD28d","LD3d","LD7d","LD28d"))
res=ggplot(df_res, aes(x=Attribute, y=Mean, fill=Species)) +
  geom_bar(stat="identity", position=position_dodge(),color="black", width=.6) +
  scale_fill_manual(values = mycol) +
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean +Sd),position=position_dodge(.6), width=.2) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),panel.background =NULL)
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\hormone_level.pdf",width=12,height = 8)





library(ggsignif)
library(ggpubr)
library(RColorBrewer)
mycol= brewer.pal(n = 5, name = "YlGnBu")[c(2:5)]
mycol
compare_means(value ~ Attribute, data = df_res, group.by = "Species")
data1=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\SD28_fragmentSize.txt",header=T)
data2=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\LD3_fragmentSize.txt",header=T)
data3=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\LD7_fragmentSize.txt",header=T)
data4=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\LD28_fragmentSize.txt",header=T)
data=rbind(data1,data2)
data=rbind(data,data3)
data=rbind(data,data4)
data$Sample=factor(data$Sample,levels=c("SD28","LD3","LD7","LD28"))
lwd_pt <- .pt*72.27/96
p1 <- ggplot(data,aes(x =Size,y = Occurrences ,color=Sample))+
  geom_line(size=1)+
  xlab("Fragment length(bp)")+
  xlim(0,800)+
  scale_color_manual(values = mycol)+
  ylab(expression(Normalized ~ read ~ count ~ 10^2))+
  ggtitle("Sample Fragment sizes")+
  theme_set(theme_bw(base_size = 20, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))
p1
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\fragment_size.pdf",width = 10,height = 8)








rm(list=ls())
dat1 = read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\GO_AllLists_Up_top10.csv")
dat1$class="Up"
#按照gene number从低到高排序[升序]
dat1=arrange(dat1,dat1[,12])
dat2 = read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\GO_AllLists_Down_top10.csv")
dat2$class="Down"
#按照gene number从低到高排序[升序]
dat2=arrange(dat2,dat2[,12])
dat=rbind(dat1[c(1:10),],dat2[c(1:10),])
library(ggplot2)#没有自己安装 install.package("ggplot2")
library(dplyr)
colnames(dat)
#dat<- arrange(dat,dat[,11])
#Pathway列最好转化成因子型，否则作图时ggplot2会将所有Pathway按字母顺序重排序
#将Pathway列转化为因子型
dat$Description=factor(dat$Description,levels = unique(dat$Description))
dat$LogP=-dat$LogP
dat$group="group"
p <- ggplot(dat,aes(y=Description,x=group,color=LogP,size=X.GeneInGOAndHitList)) + 
  geom_point() +
  theme_bw() +
  scale_color_gradient(low='blue',high='Red')+ 
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(size = 8),
        legend.position="right",
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        axis.text.x = element_text(size=8,vjust = 0.1,angle =  30),
        axis.text.y = element_text(size=8),
        axis.title.x =NULL ,
        axis.title.y =NULL )+
  scale_y_discrete(position = "right")
p

library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(stringr)
library(tidyr)
#color
mycol= brewer.pal(n = 7, name = "YlGnBu")
##data1
data1=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\bedops_SD28_peak_anno.txt",sep="\t",header=F)
data1=as.data.frame(table(data1$V6))
colnames(data1)=c("position","count")
data1$class="SD28"
##data2
data2=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\bedops_LD3_peak_anno.txt",sep="\t",header=F)
data2=as.data.frame(table(data2$V6))
colnames(data2)=c("position","count")
data2$class="LD3"
##data3
data3=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\bedops_LD7_peak_anno.txt",sep="\t",header=F)
data3=as.data.frame(table(data3$V6))
colnames(data3)=c("position","count")
data3$class="LD7"
##data4
data4=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\bedops_LD28_peak_anno.txt",sep="\t",header=F)
data4=as.data.frame(table(data4$V6))
colnames(data4)=c("position","count")
data4$class="LD28"
##data merge
data=rbind(data1,data2)
data=rbind(data,data3)
data=rbind(data,data4)
data$class=factor(data$class,levels = c("SD28","LD3","LD7","LD28"))
data$position=factor(data$position,levels = unique(data$position))
##plot
ggplot(data=data, aes(x=class, y=count,fill = position)) + 
  geom_bar(stat = "identity")+
  scale_fill_npg()+
  #scale_fill_manual(values=mycol)+
  ylab("Number andgenomic distribution of peaks identified by ATAC-seq")+
    xlab("Sample")+theme(panel.background = element_blank(),panel.border = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.key = element_blank(),legend.title = element_blank(),
                         axis.text.x = element_text(size=12,angle=60, vjust=0.5,hjust=0.5,
                         ))
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\peak_distribution_static.pdf")



rm(list=ls())
library(Sushi)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

## 参数设置


chrom <- "chr12"
chromstart <- 75477413 #全局起始位点
chromend <- 75653942 #全局结束位点
#例图zoom in三次，在这里设置放大区域的起止位点
region <- c(75633289,75644072) #zoom in的范围
region2 <- c(75635224,75636574) #进一步zoom in的范围
region3 <- c(75635776,75635931) #再进一步zoom in的范围


## 输入文件

#Sushi接受输入bedgragh（.bdg）、bedpe、bed、interaction matrix文件。

#此处需要bedgragh和bed文件。如果你的数据类型不符，可以用以下方法做转换：

# - bam文件，可以用bedtools转成bedgragh、bedpe或bed文件
# - gff文件，可以用bedtools转成bed文件
# - bw文件，可以到<http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/>下载“bigWigToBedGraph”，转换为.bdg文件；或者用FigureYa3genomeView来画。
# 
# ChIP-seq、ATAC-seq、DNase-seq、RNA-seq、GRO-seq的信号用bedgragh文件。
# 
# 基因结构、peak或motif所在的位置用bed文件。
# 
# 此处展示了多种测序数据类型，实际操作中根据你需要的类型截取相应的代码来画。
# 
# ```{r}
a <- read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\igvZoom\\gene_exon.bed",header=TRUE,sep="\t")
head(a)
b <- read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\igvZoom\\ATACseq.bdg",header=FALSE,sep="\t")
head(b)
c <- read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\igvZoom\\H3K27ac.bdg",header=FALSE,sep="\t")
d <- read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\igvZoom\\RNA_Pol2.bdg",header=FALSE,sep="\t")
e <- read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\igvZoom\\GROseq+.bdg",header=FALSE,sep="\t")
f <- read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\igvZoom\\GROseq-.bdg",header=FALSE,sep="\t")
g <- read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\igvZoom\\P300.bdg",header=FALSE,sep="\t")
h <- read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\igvZoom\\ATACseq_peaks.bed",header=FALSE,sep="\t")
head(h)
win.graph(width=10, height=10,pointsize=10) 
# ```
layout(matrix(1:9,9,1), heights=c(1,2,1,1,1,1,2,1,1))
par(mgp=c(3,.3,0))
#1 画基因结构
pdf("C:\\Users\\Jeff\\OneDrive\\桌面\\test.pdf")
par(mar=c(0,20,1.5,20)) #bottom, left, top, right
plotGenes(a,chrom,chromstart,chromend,
          labeltext=FALSE,
          bentline=FALSE,bheight=.1,col="black",
          plotgenetype="box",
          maxrows=1)

#画基因组标尺
#sushi没提供画比例尺的功能，如果想要例文那种比例尺，可以每次zoom in都画出这种从头到尾的标尺，然后根据标尺用AI画出比例尺
labelgenome(chrom,chromstart,chromend,
            side=3, #标尺画在上面，或1把标尺画在下面
            n=3, #刻度线的数量
            scale="Mb", #'bp','Kb','Mb'
            line=0, #刻度标签向下移动的距离
            chromcex=0.7,chromadjust=-0.02,scalecex=0.7)
mtext("Genes",side=2,line=7.5,cex=1,font=1,las=2,adj=0,col="black")
#2 画RNA_Pol2的ChIP-seq
par(mar=c(4,20,0,20))
plotBedgraph(d,chrom,chromstart,chromend,
             color="black") #用黑色填充profile
mtext("RNA Pol2",side=2,line=7.5,cex=1,font=1,las=2,adj=0,col="black")
#画被zoom in区域的虚线方框
zoomsregion(region=region,extend=c(0,0.5),wideextend=0.7,offsets=c(0,0),lty=3)
# 第一次zoom in
#3 画GRO-seq的正链，适用于链特异性的RNA-seq数据
par(mar=c(0.2,20,0,20))
plotBedgraph(e,chrom,region[1],region[2],
             color="grey") #用灰色填充profile
mtext("GROseq+",side=2,line=7.5,cex=1,font=1,las=2,adj=0,col="black")
#画zoom in区域的虚线方框
zoombox(passthrough=TRUE,lty=3)
#画被zoom in区域的虚线方框
zoomsregion(region2,extend=c(0,1),wideextend=0,lty=3)
#4 画GRO-seq的负链，适用于链特异性的RNA-seq数据
par(mar=c(1,20,0,20))
plotBedgraph(f,chrom,region[1],region[2],flip=TRUE,color="grey")
mtext("GROseq-",side=2,line=7.5,cex=1,font=1,las=2,adj=0,col="black")
zoombox(passthrough=TRUE,lty=3)
zoomsregion(region2,extend=c(1,1),wideextend=0,lty=3)
#5 画H3K27ac的ChIP-seq
par(mar=c(1,20,0,20))
plotBedgraph(c,chrom,region[1],region[2],color="grey")
mtext("H3K27ac",side=2,line=7.5,cex=1,font=1,las=2,adj=0,col="black")
zoombox(passthrough=TRUE,lty=3)
zoomsregion(region2,extend=c(1,1),wideextend=0,lty=3)
#6 画P300的ChIP-seq
par(mar=c(1,20,0,20))
plotBedgraph(g,chrom,region[1],region[2],color="grey")
mtext("P300",side=2,line=7.5,cex=1,font=1,las=2,adj=0,col="black")
zoombox(passthrough=TRUE,lty=3)
zoomsregion(region2,extend=c(1,1),wideextend=0,lty=3)
#7 画ATACseq
par(mar=c(4,20,0,20))
plotBedgraph(b,chrom,region[1],region[2],color="grey")
mtext("ATACseq",side=2,line=7.5,cex=1,font=1,las=2,adj=0,col="black")
zoombox(lty=3)
zoomsregion(region2,extend=c(1,0.5),wideextend=0.7,offsets=c(0,0),lty=3)
# 第二次zoom in
#8 画ATACseq
par(mar=c(1,20,0,20))
plotBedgraph(b,chrom,region2[1],region2[2],range=c(0,0.5),color="black")
mtext("ATACseq",side=2,line=7.5,cex=1,font=1,las=2,adj=0,col="black")
zoombox(passthrough=TRUE,lty=3,topextend = 1)
zoomsregion(region3,extend=c(0,1),wideextend=0,lty=1,zoomborder=rgb(0,0,0,20,maxColorValue=255),col=rgb(0,0,0,20,maxColorValue=255))
#9 画ATACseq的峰，同样适用于motif位置的展示
par(mar=c(3,20,0,20))
plotBed(beddata=h,chrom=chrom,region2[1],region2[2],row=1,color="black")
mtext("",side=2,line=7.5,cex=1,font=1,las=2,adj=0,col="black")
zoombox(lty=3)
# 画第三次zoom in的框
zoomsregion(region3,extend=c(1,2),wideextend=4,offsets=c(0,0),lty=3,zoomborder=rgb(0,0,0,20,maxColorValue=255),col=rgb(0,0,0,20,maxColorValue=255))
#后面用AI加上motif
dev.off()


library(seqLogo)
library(ggseqlogo)
data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\sanger_1-5AOX1_C08_statics.csv")
data1=data[c(310:349),c(2,3,4,5)]
data2=data[c(350:389),c(2,3,4,5)]
data3=data[c(390:429),c(2,3,4,5)]
#data=data[c(310:369),c(2,3,4,5)]
datNormed1 <- t(apply(data1, 1, function(x) x / (x[1]+x[2]+x[3]+x[4])))
datNormed2 <- t(apply(data2, 1, function(x) x / (x[1]+x[2]+x[3]+x[4])))
datNormed3 <- t(apply(data3, 1, function(x) x / (x[1]+x[2]+x[3]+x[4])))
pdf("C:\\Users\\Jeff\\OneDrive\\桌面\\xmk_seqLogo.pdf",height = 4,width = 15)
pwm1 <- makePWM(t(datNormed1))
pwm2 <- makePWM(t(datNormed2))
pwm3 <- makePWM(t(datNormed3))
seqLogo(pwm1,ic.scale = F)
seqLogo(pwm2,ic.scale = F)
?seqLogo(pwm3,ic.scale = F)
dev.off()

p1=seqLogo(pwm)
plot_grid(p,p1,labels = "AUTO")
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\111.pdf")
ggplot()+geom_logo(pwm)+theme_logo()










library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(dplyr)
mycol= brewer.pal(n = 7, name = "YlGnBu")
data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\sort_merge_diff_peakAnno.txt",sep="\t",header=F)
data=as.data.frame(table(data$V6))
colnames(data)=c("position","count")
data$class="Sample"
data<- arrange(data,data[,2])
data$position=factor(data$position,levels = unique(data$position))
ggplot(data=data, aes(x=class, y=count,fill = position)) + 
  geom_bar(stat = "identity")+
  scale_fill_npg()+
  coord_polar(theta = "y")+blank_theme
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\diff_peak_distribution_static.pdf")

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

ggplot(data=data, mapping=aes(x=class,fill=position))+
  geom_bar(stat="count",width=0.5,position='stack',size=5)+
  coord_polar("y", start=0)+
  scale_fill_npg()

  #geom_text(stat="count",aes(label = scales::percent(..count../100)), size=4, position=position_stack(vjust = 0.5))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("motifmatchr")


rm(list=ls())
library("JASPAR2020")
library(TFBSTools)
library(motifmatchr)
library(Biostrings)
#get JASPAR2020 data
species <- "Homo sapiens"
collection <- "CORE"
opts <- list()
opts["species"] <- species
opts["collection"] <- collection
out  <- TFBSTools::getMatrixSet(JASPAR2020, opts)
if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
  names(out) <- paste(names(out), TFBSTools::name(out), 
                      sep = "_")
##transfer to pwm matrix
pwm <- toPWM(out)
## read my Fa
ht_data=readDNAStringSet("C:\\Users\\Jeff\\OneDrive\\桌面\\Hub_gene_promoter.fa")

## set base path 
my_path="C:\\Users\\Jeff\\OneDrive\\桌面\\"

## Motif scan and save my result to csv file
for (i in c(1:length(ht_data@ranges))){
  seqname=ht_data@ranges@NAMES[i]
  dna <- DNAString(ht_data[[i]])
  siteset <- searchSeq(pwm, dna, seqname=seqname, min.score="95%", strand="*")
  output=paste(my_path,seqname,i,"_TFBS.csv",sep="")
  write.csv(writeGFF3(siteset),output,quote=F)
}

library(gg.gap)
library(ggplot2)
library(patchwork)
#数据
data <-
  data.frame(x = c("Alpha", "Bravo", "Charlie", "Delta"),
             y = c(200, 20, 10, 15))
#画图
p1 = ggplot(data, aes(x = x, y = y, fill = x)) +
  geom_bar(stat = 'identity', position = position_dodge(),show.legend = FALSE) +
  theme_bw() +
  labs(x = NULL, y = NULL)

p1
p2 =gg.gap(plot = res,
           segments = c(5, 40),
           tick_width = 10,
           rel_heights = c(0.25, 0, 0.1),
           ylim = c(0, 50)
)
p2
p1+p2





library(circlize)
circos.initializeWithIdeogram(plotType = c("labels", "axis"))
circos.track(ylim = c(0, 1))
circos.genomicIdeogram() # put ideogram as the third track
circos.genomicIdeogram(track.height = 0.2)
circos.clear()

circos.initializeWithIdeogram()
bed <- generateRandomBed(nr = 100, nc = 4)
col_fun <- colorRamp2(
  c(-1, 0, 1), 
  c("#ef8a62", "#f7f7f7", "#67a9cf")
)
circos.genomicHeatmap(
  bed, col = col_fun, side = "inside", 
  border = "white"
)
circos.clear()
BiocManager::install("clusterprofiler")
library(clusterProfiler)
BiocManager::install("DO.db")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(SeuratData)
library(cowplot)


#InstallData("stxBrain")
install.packages("esquisse")
library(esquisse)

#加载esquisse插件：
esquisse :: esquisser()
diamonds=diamonds   
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DiffBind")

rm(list=ls())
library(DESeq2)
data1=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\Count_matrix.txt",header=T,row.names = 1)
data2=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\Count_matrix (1).txt",header=T,row.names = 1)
data=data.frame(data1[,c(6:10)],data2[,c(6:13)])
colnames(data)=c("SD28_1","SD28_2","SD28_3","LD7_1","LD7_2","SD28_4","SD28_5","SD28_6","SD28_7","LD28_1","LD28_2","LD28_3","LD28_4")
data["DIO2",]
colData <- data.frame(condition=factor(c("SD28","SD28","SD28",
                                      "LD7","LD7","SD28","SD28","SD28","SD28","LD28","LD28","LD28","LD28")))
colData
dds<-DESeqDataSetFromMatrix(data,colData, formula(~condition))
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
rld <- rlog(dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup=c("condition"),returnData = T)
pcaData
plotPCA(rld)




rm(list=ls())
library(ggplot2)
library(ggsci)
library(ggthemes)
SD28=read.csv(
  "C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\QC\\SD28_tss.tsv",
            header = F,fileEncoding = "UTF8",sep="\t")
SD28$V3="SD28"
LD3=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\QC\\LD3_tss.tsv",
              header = F,fileEncoding = "UTF8",sep="\t")
LD3$V3="LD3"
LD7=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\QC\\LD7_tss.tsv",
             header = F,fileEncoding = "UTF8",sep="\t")
LD7$V3="LD7"
LD28=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\QC\\LD28_tss.tsv",
              header = F,fileEncoding = "UTF8",sep="\t")
LD28$V3="LD28"
data=rbind(SD28,LD3)
data=rbind(data,LD7)
data=rbind(data,LD28)
colnames(data)[3]=c("Sample")
data$Sample=factor(data$Sample,levels = c("SD28","LD3","LD7","LD28"))
p=ggplot(data = data, mapping = aes(x = V1, y = V2, colour = Sample)) + 
  geom_line(size=1.5,alpha = 1)+
  scale_color_tableau()+
  theme_bw()+
  xlab("Distance to Transcription Start Site (bp) ")+
  ylab("Normalized Insertion")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank())+
  scale_y_continuous(breaks=c(1:6))
p

ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\TSS_plot.pdf"
       ,height = 3,width=4)

rm(list=ls())
library(dplyr)
library(ggplot2)
library(pheatmap)
res1=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\All_results.csv")
res2=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD28\\All_results.csv")
res1=res1[,c("X","log2FoldChange")]
res2=res2[,c("X","log2FoldChange")]
row.names(res1)=res1[,1]
row.names(res2)=res2[,1]
res=merge(res1,res2,by="row.names")
res=na.omit(res)
data=data.frame(res$Row.names,res$log2FoldChange.x,res$log2FoldChange.y)
rownames(data)=data[,1]
data=data[,-1]
colnames(data)=c("SD28_LD7","SD28_LD28")
data$label=""
data$label <- case_when(as.numeric(data$SD28_LD7) >= 1 & as.numeric(data$SD28_LD28) >= 1 ~ "Up Up",
                       as.numeric(data$SD28_LD7) < -1 & as.numeric(data$SD28_LD28) > 1 ~ "Down Up",
                       as.numeric(data$SD28_LD7) > 1 & as.numeric(data$SD28_LD28) < -1 ~ "Up Down",
                       as.numeric(data$SD28_LD7) < -1 & as.numeric(data$SD28_LD28) < -1 ~ "Down Down")
ggplot() + 
  geom_bar(data = na.omit(data), aes(x = label,fill=label), stat = "count")+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#去除背景
        panel.border = element_blank())+xlab(NULL)#去除边框
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\Figure\\DEG_static.pdf"
       ,height = 4,width = 4)




library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
mycol= rev(brewer.pal(n = 5, name = "Set1")[c(2:5)])
mycol[1]="red"
mycol
compare_means(value ~ Attribute, data = df_res, group.by = "Species")
data1=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\QC\\SD28_fragmentSize.txt",header=T)
data2=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\QC\\LD3_fragmentSize.txt",header=T)
data3=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\QC\\LD7_fragmentSize.txt",header=T)
data4=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\QC\\LD28_fragmentSize.txt",header=T)
data=rbind(data1,data2)
data=rbind(data,data3)
data=rbind(data,data4)
data$Sample=factor(data$Sample,levels=c("SD28","LD3","LD7","LD28"))
lwd_pt <- .pt*72.27/96
p1 <- ggplot(data,aes(x =Size,y = Occurrences ,color=Sample))+
  geom_line(size=0.6)+
  xlab("Fragment length(bp)")+
  xlim(0,1000)+
  scale_color_manual(values = mycol)+
  ylab(expression(Normalized ~ read ~ count ~ 10^2))+
  ggtitle("Sample Fragment sizes")+
  theme_set(theme_bw(base_size = 15, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))
p1
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\fragment_size.pdf",width = 10,height = 8)


args <- list( Species=9606,collection = "CORE", all_versions = TRUE)
motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, args)
motifs


rm(list=ls())
library(DESeq2)
library(ggplot2)
library(ggrepel)
countsTable <- read.table( "C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\Count_matrix.txt",
                           stringsAsFactors=TRUE,row.names=1 )[,c(4,5,6,7,8)]

colData <- data.frame(condition=factor(c("SD28","SD28","SD28","LD7","LD7")))
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)

res=results(dds,contrast = c("condition","LD7","SD28"))
summary(res)
res <- res[order(res$padj),]
write.table(as.data.frame(res),"C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\All_results.csv",quote=F,sep=",")

res=as.data.frame(res)
res=na.omit(res)
res[which(res$log2FoldChange >= 1 & res$padj < 0.05),'sig'] <- 'UP'
res[which(res$log2FoldChange <= -1 & res$padj< 0.05),'sig'] <- 'DOWN'
res[which(!(res$sig%in%c("UP","DOWN"))),'sig'] <- 'NO DIFF'


res$label <-''
res[rownames(res)%in%c("DIO3","GPR20","TSHB","CGA","AGRP","SLC16A2","GHRH","DIO2","POMC","PER2"),]$label=c("DIO3","GPR20","TSHB","CGA","AGRP","SLC16A2","GHRH","DIO2","POMC","PER2")
this_tile <- paste0('cutoff for abs(logFC) and FDR is 1 and 0.05',
                    '\nThe number of up gene is ',nrow(res[res$sig =='UP',]) ,
                    '\nThe number of down gene is ',nrow(res[res$sig =='DOWN',]))
volcano <- ggplot(res, aes(log2FoldChange, -log(padj, 10))) +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
  geom_point(aes(color = sig), size = 1.2) +
  theme(panel.background=element_rect(fill='transparent', color='gray'),
        legend.key=element_rect(fill='transparent', color='transparent'),
        axis.text=element_text(color='black')) +
  theme(legend.title = element_blank(),text=element_text(size=16,  family="serif")) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 padj', color = NA) +
  xlim(-8, 8)

volcano
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\Volcano.pdf",height = 6,width = 6)

rm(list=ls())
library(DESeq2)
library(ggplot2)
library(ggrepel)
countsTable <- read.table( "C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\Count_matrix.txt",
                           stringsAsFactors=TRUE,row.names=1 )[,c(4,5,6,7,8)]

colData <- data.frame(condition=factor(c("SD28","SD28","SD28","LD7","LD7")))
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds,contrast = c("condition","LD7","SD28"))
diff_gene=subset(res,abs(log2FoldChange)>=1 & padj <=0.05)
normalized_counts <- counts(dds, normalized=TRUE)
data=data[which(rowSums(data)>0),]
pheatmap::pheatmap(data[rownames(diff_gene),],
                   scale="row",
                   cluster_cols = F,
                   show_rownames = F,
                   cellwidth = 10, 
                   cellheight = 3)

rm(list=ls())
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(stringr)
library(dplyr)
mycol= brewer.pal(n = 7, name = "YlGnBu")
data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\bedops\\peak_merged.anno",sep="\t",header=F)
colnames(data)=c("chr","start","end","width","strand","annotation","strands","gene_start","gene_end","gene_width","gene_strand","gene","strand")
split_b<-str_split(data$annotation," ")
b<-sapply(split_b,"[",1)
data$annotation=b
data[data$annotation=="3'","annotation"]="3' UTR"
data[data$annotation=="5'","annotation"]="5' UTR"
data=as.data.frame(table(data$annotation))
colnames(data)=c("position","count")
data$class="Sample"
data$ratio=round(data$count/sum(data$count),3)
data<- arrange(data,data[,2])
data$position=factor(data$position,levels = unique(data$position))
blank_theme <- theme_minimal()+
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
ggplot(data=data, aes(x=class, y=ratio,fill = position)) + 
  geom_bar(stat = "identity")+
  scale_fill_npg(labels = paste(data$position," (",scales::percent(data$ratio),") ",sep=""))+
  coord_polar(theta = "y")+blank_theme+
  ylab("Number andgenomic distribution of peaks identified by ATAC-seq")+
  xlab("Sample")
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\peak_distribution_static.pdf")
scales::percent(data$ratio)



data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\All_results.csv",row.names = 1)
data=data[abs(data$log2FoldChange)>=1&data$padj<=0.05,]
write.csv(data,"C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\early_target_gene.csv",quote=F,sep="\t")



library(ggplot2)
library(ggthemes)
library(ggpubr)
data<-read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\cluster_motif\\motif_rpkm.csv",
               header=T)
data$group=factor(data$group,levels = c("SD28_LD7","SD28_LD28"))
data=na.omit(data)
ggplot(data,aes(x=log2FoldChange,y=TF,color=group))+
  xlab("Gene expression log2FoldChange")+
  ylab("TF encode gene")+ 
  geom_point(alpha = 0.8,size=5)+
  scale_color_manual(values=c("#2878B5","#C82423"))+
  theme(axis.title.x =element_text(size=14), 
        legend.background = element_blank(),
        panel.background =  element_blank(), panel.grid.major = element_blank(),
        legend.key = element_rect(fill="white"),
        axis.title.y=element_text(size=14))+
        coord_cartesian(xlim = c(-4,2))
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\motif_rpkm.pdf",height = 6,width = 6)

data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\GO_AllLists_Up_top10.csv")
data_BP=data[data$Category=="GO Biological Processes",]
data_MF=data[data$Category=="GO Molecular Functions",]
data_CC=data[data$Category=="GO Cellular Components",]
colnames(data)
pathway=rbind(data_BP[c(1:5),],data_MF)
pathway=rbind(pathway,data_CC)
ggplot(pathway,aes(-(LogP),Description))+
  geom_point(aes(size=X.GeneInGOAndHitList,color=LogP))+
  scale_color_gradient(low="green",high = "red")+ #
  labs(color=expression(-log[10](Qvalue)),size="Gene",  ##expression函数定义函数样式 []添加下标，^添加上标
       x="Pvalue",      ##自定义标轴
       y="Pathway name",
       title="Pathway enrichment")+
  scale_size(c(1:10))##自定义坐标轴



library(ggplot2)
library(dplyr)
#按照前面的成品预览图，先进行一些背景设定
mytheme <- theme(axis.title=element_text(face="bold", size=10,colour = 'gray25'), #轴标题
                 axis.text=element_text(face="bold", size=10,colour = 'gray25'), #轴标签
                 axis.line = element_line(size=0.5, colour = 'black'), #轴线
                 axis.line.y = element_blank(),  #关闭Y轴线
                 axis.ticks.y = element_blank(), #关闭Y轴刻度线
                 panel.background = element_rect(fill="white"), #背景色
                 panel.grid.major.y=element_blank(), #关闭Y轴主网格线
                 panel.grid.minor.y=element_blank(), #关闭Y轴次网格线
                 panel.grid.minor.x=element_blank()) #关闭X轴次网格线
data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\SD28_LD7_KEGG.txt",sep="\t")
data=data[data$P.Value<=0.05,]
data<- arrange(data,-data[,4])
p<-ggplot(data,aes(y=Term,x=Input.number,fill=P.Value)) + #以Set为横坐标，FE为纵坐标画柱状图，并用P填充颜色
  geom_bar(stat='identity',color='black',width = 0.65) +
  scale_fill_gradient(low='red', high='darkgoldenrod1')  #设定柱子颜色变化范围，随着p值从低到高，柱子的颜色从红色向黄色渐变

p + labs(x='Gene Number',y='KEGG Pathway',fill='p value') +
  mytheme
ggsave(filename = 
         "C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\SD28_LD7_KEGG.pdf",
       height = 3,
       width = 10)





rm(list=ls())
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(stringr)
library(dplyr)
mycol= brewer.pal(n = 7, name = "YlGnBu")
data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\bedops\\peak_merged.anno",sep="\t",header=F)
colnames(data)=c("chr","start","end","width","strand","annotation","strands","gene_start","gene_end","gene_width","gene_strand","gene","strand")
split_b<-str_split(data$annotation," ")
b<-sapply(split_b,"[",1)
data$annotation=b
data[data$annotation=="3'","annotation"]="3' UTR"
data[data$annotation=="5'","annotation"]="5' UTR"

# data=as.data.frame(table(data$annotation))
# colnames(data)=c("position","count")
# data$class="Sample"
# data$ratio=round(data$count/sum(data$count),3)
# data<- arrange(data,data[,2])
# data$position=factor(data$position,levels = unique(data$position))
rm(list=ls())
peak_SD28=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\bedops\\bedops_SD28_peak_anno.txt",sep="\t",header=F)
peak_LD3=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\bedops\\bedops_LD3_peak_anno.txt",sep="\t",header=F)
peak_LD7=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\bedops\\bedops_LD7_peak_anno.txt",sep="\t",header=F)
peak_LD28=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\bedops\\bedops_LD28_peak_anno.txt",sep="\t",header=F)
gtf=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\GCF_001577835.2_Coturnix_japonica_2.1_genomic.gtf",sep="\t",header=F)
genewithpeak_SD28=unique(peak_SD28$V12)
genewithpeak_LD3=unique(peak_LD3$V12)
genewithpeak_LD7=unique(peak_LD7$V12)
genewithpeak_LD28=unique(peak_LD28$V12)
length(genewithpeak)
head(gtf)
gtf=gtf[gtf$V3=="gene",]
data=data.frame(class=c("SD28","LD3","LD7","LD28"),ratio=c(length(genewithpeak_SD28)/dim(gtf)[1],length(genewithpeak_LD3)/dim(gtf)[1],length(genewithpeak_LD7)/dim(gtf)[1],length(genewithpeak_LD28)/dim(gtf)[1]))
data=data.frame(class=c("gene with peak","gene without peak"),ratio=c(length(genewithpeak_SD28)/dim(gtf)[1],1-length(genewithpeak_SD28)/dim(gtf)[1]))

blank_theme <- theme_minimal()+
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
ggplot(data=data, aes(x=class, y=ratio,fill = class)) + 
  geom_bar(stat = "identity")+
  scale_fill_npg(labels = paste(data$class," (",scales::percent(data$ratio),") ",sep=""))+
  coord_polar(theta = "y")+blank_theme+
  ylab("Number andgenomic distribution of peaks identified by ATAC-seq")+
  xlab("Sample")
setwd("C:\\Users\\Jeff\\OneDrive\\桌面")
png(file='runoob-pie.png', height=300, width=300)
cols = c("#22B14C","#FFC90E")
piepercent = paste(data$class ,round(data$ratio*100,2), "% ")
pie(data$ratio,col = cols,labels=piepercent)
legend("topright", names, cex=0.8, fill=cols)
dev.off()

library(dplyr)
library(stringr)
rm(list=ls())
data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\bedops\\bedops_SD28_peak_anno.txt",header=F,sep="\t")
colnames(data)=c("chr","start","end","width","strand","annotation","strands","gene_start","gene_end","gene_width","gene_strand","gene","strand")
split_b<-str_split(data$annotation," ")
b<-sapply(split_b,"[",1)
data$annotation=b
data[data$annotation=="3'","annotation"]="3' UTR"
data[data$annotation=="5'","annotation"]="5' UTR"
data=as.data.frame(table(data$annotation))
colnames(data)=c("position","count")
data$class="Sample"
data$ratio=round(data$count/sum(data$count),3)
data<- arrange(data,data[,2])
data$position=factor(data$position,levels = unique(data$position))
blank_theme <- theme_minimal()+
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
ggplot(data=data, aes(x=class, y=ratio,fill = position)) + 
  geom_bar(stat = "identity")+
  scale_fill_npg(labels = paste(data$position," (",scales::percent(data$ratio),") ",sep=""))+
  coord_polar(theta = "y")+blank_theme+
  ylab("Number andgenomic distribution of peaks identified by ATAC-seq")+
  xlab("Sample")
library("ChIPseeker")
library(ChIPseeker)
library(ggplot2)

files <- getSampleFiles()
f2=getSampleFiles()[[5]]
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.downstreamDistance = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)
x=annotatePeak(f2)
plotAnnoPie(x)



data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\PRJEB12420.tsv",sep="\t")
data=as.data.frame(data[c(1:dim(data)[1]),"fastq_aspera"])
colnames(data)="url"
library(stringr)
#分列
split_b<-str_split(data$url,";")
b<-sapply(split_b,"[",1)
b
c<-sapply(split_b,"[",2)
c
data=c(b,c)
split_b<-str_split(data,":")
b<-sapply(split_b,"[",1)
b
c<-sapply(split_b,"[",2)
c
write.table(c,"C:\\Users\\Jeff\\OneDrive\\桌面\\PRJEB12420_del.tsv",sep="\t",row.names = F,col.names = F,quote=F)






options("repos" = c(CRAN="https://mirror.lzu.edu.cn/CRAN/")) 
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor")
install.packages("RIdeogram")
install.packages("getopt")
BiocManager::install("ggbio")
BiocManager::install("Gviz")
library(RIdeogram)
library(RColorBrewer)
library(ggbio)
library(gggenes)
library(Gviz)
library(ggplot2)
suppressMessages(library(dplyr))
suppressMessages(library(getopt))
setwd("C:\\Users\\Jeff\\OneDrive\\桌面")
install.packages('RIdeogram')
require(RIdeogram)
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")
head(human_karyotype)
head(Random_RNAs_500)
head(gene_density)
ideogram(karyotype = human_karyotype)
convertSVG(svg = "chromosome.svg", device = "png")
human_karyotype <- human_karyotype[,1:3]
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500, label_type = "marker")
id=read.table("gene_family_info.txt")
gene_family=data.frame(Type="gene_family",Shape="triangle",Chr=id$V1,Start=id$V4,End=id$V5,color="6a3d9a")
gene=read.table("gene_info.txt")
gene=data.frame(Chr=gene$V1,Start=gene$V4,End=gene$V5)

## import data
data=gene[,c(1:2)]
# function
getdensity = function(data,binsize,key,maxVal){
  tmp1 = data.frame(filter(data,data[,1] == key),
                    count = 1)
  bin = seq(binsize,maxVal,by = binsize)
  tmp2 = data.frame(key = key,
                    bin = bin,
                    count = 0)
  for (i in 1:length(bin)) {
    x = (filter(tmp1,tmp1[,2]>bin[i]-binsize & tmp1[,2]<bin[i])[,-1]%>%apply(.,2,sum))[2]
    tmp2[i,3] =  x
  }
  tmp2[i+1,1] = key
  tmp2[i+1,2] = maxVal
  tmp2[i+1,3] = (filter(tmp1,tmp1[,2]>bin[i])[,-1]%>%apply(.,2,sum))[2]
  tmp2$test=tmp2$bin-binsize+1
  tmp2=tmp2[,c(1,4,2,3)]
  colnames(tmp2)=c("Chr","Start","End","Value")
  return(tmp2)
}
## excute
out=data.frame(Chr=NULL ,  Start=NULL, End=NULL,   Value=NULL)
for (i in c(1:28)){
  data = data
  binsize = 1000000
  key = chr$Chr[i]
  maxVal = 43627644
  tmp = getdensity(data = data,
                   binsize = binsize,
                   key = key,
                   maxVal = maxVal)

  out=rbind(out,tmp)
}
out
out=out[-which(out$Value==0),]

write.table(out,paste(key,"_",binsize,".Dens.xls",sep = ""),
            row.names = F,
            col.names = F,
            sep = "\t",
            quote = F)
gene_family
ideogram(karyotype = chr, overlaid = out, width = 100 ,label =gene_family, label_type = "marker")
convertSVG(svg = "chromosome.svg", device = "png",width = 20,height = 40)

rm(list=ls())
library(ggplot2)
library(dplyr)
atac_data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Quadrant\\Information_SD28_LD7_DESEQ2_anno.txt",sep="\t",header=F)
split_b<-str_split(atac_data$V12," ")
b<-sapply(split_b,"[",1)
atac_data$V12=b
atac_data=atac_data[atac_data$V12=="Promoter",]
gene=as.character(unique(atac_data$V18))
atac=data.frame(target="",atac_log2FoldChange="",atac_p="",atac_padj="")
atac[,"target"]=as.character(atac[,"target"])
atac[,"atac_log2FoldChange"]=as.character(atac[,"atac_log2FoldChange"])
atac[,"atac_p"]=as.character(atac[,"atac_p"])
atac[,"atac_padj"]=as.character(atac[,"atac_padj"])
for (i in c(1:length(gene))){
  tmp=atac_data[atac_data$V18==gene[i],c(7,10,11)]
  sum_data=sum(tmp[,1])
  atac[i,"target"]=gene[i]
  atac[i,"atac_log2FoldChange"]=sum_data
  atac[i,"atac_p"]=min(tmp[,2])
  atac[i,"atac_padj"]=min(tmp[,3])
}
rownames(atac)=atac[,1]
atac=atac[,-1]
RNA_data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Quadrant\\All_results_SD28_LD7.csv",row.names=1)
RNA=RNA_data[gene,]
colnames(RNA)[2]="RNA_log2FoldChange"
res=cbind(atac,RNA)
res$label <- case_when(abs(as.numeric(res$RNA_log2FoldChange)) >= 1 & abs(as.numeric(res$atac_log2FoldChange)) >= 1 ~ "part1379",
                       abs(as.numeric(res$atac_log2FoldChange)) < 1 & abs(as.numeric(res$RNA_log2FoldChange)) > 1 ~ "part28",
                       abs(as.numeric(res$atac_log2FoldChange)) > 1 & abs(as.numeric(res$RNA_log2FoldChange)) < 1 ~ "part46",
                       abs(as.numeric(res$atac_log2FoldChange)) < 1 & abs(as.numeric(res$RNA_log2FoldChange)) < 1 ~ "part5")          

res$RNA_log2FoldChange=as.numeric(res$RNA_log2FoldChange)
res$atac_log2FoldChange=as.numeric(res$atac_log2FoldChange)
res$atac_p=as.numeric(res$atac_p)
res$atac_padj=as.numeric(res$atac_padj)
str(res)
#开始尝试绘图；
p0 <-ggplot(res,aes(atac_log2FoldChange,RNA_log2FoldChange,color=label))
#添加散点；
p1 <- p0+geom_point(size=2)+guides(color="none")
p1

##自定义半透明颜色
mycolor <- c("#FF9999","#99CC00","#c77cff","gray80")
p2 <- p1 + scale_colour_manual(name="",values=alpha(mycolor,0.7))
p2  

#添加辅助线；
p3 <- p2+geom_hline(yintercept = c(-1,1),
                    size = 0.8,
                    color = "grey40",
                    lty = "dashed")+
  geom_vline(xintercept = c(-1,1),
             size = 0.8,
             color = "grey40",
             lty = "dashed")
p3

#调整横轴和纵轴绘图区域的范围；
#设置y轴范围（上下两端的空白区域设为1），修改刻度标签；
#expansion函数设置坐标轴范围两端空白区域的大小；mult为“倍数”模式，add为“加性”模式；
p4<-p3+
  scale_y_continuous(expand=expansion(add = c(0.5, 0.5)),
                     limits = c(-8, 8),
                     breaks = c(-6,-3,0,3,6),
                     label = c("-6","-3","0","3","6"))+
  scale_x_continuous(expand=expansion(add = c(0.5, 0.5)),
                     limits = c(-4, 4))
                     #breaks = c(-6,-3,0,3,6),
                     #label = c("-4","-3","0","3","6"))
p4


#自定义图表主题，对图表做精细调整；
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
#隐藏纵轴，并对字体样式、坐标轴的粗细、颜色、刻度长度进行限定；
mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="gray30",size = 15),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
#应用自定义主题；
p5 <- p4+mytheme
p5


cor = cor(res[,c("RNA_log2FoldChange","atac_log2FoldChange")],use="complete.obs")

cor=round(cor[2,1],digits = 2)
#准备作为图形的标题;
lab = paste("correlation=",cor,sep="")
lab
#[1] "correlation=0.35"
#在图上添加文字标签；
p5+geom_text(x = -2, y = 5, label = lab)



rm(list=ls())
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(stringr)
library(tidyr)
setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\bedops")
#color
mycol= brewer.pal(n = 7, name = "YlGnBu")
##data1
data1=read.csv("bedops_SD28_peak_anno.txt",sep="\t",header=F)
split_b<-str_split(data1$V6," ")
b<-sapply(split_b,"[",1)
data1$V6=b
data1=as.data.frame(table(data1$V6))
colnames(data1)=c("position","count")
data1$class="SD28"

##data2
data2=read.csv("bedops_LD3_peak_anno.txt",sep="\t",header=F)
split_b<-str_split(data2$V6," ")
b<-sapply(split_b,"[",1)
data2$V6=b
data2=as.data.frame(table(data2$V6))
colnames(data2)=c("position","count")
data2$class="LD3"
##data3
data3=read.csv("bedops_LD7_peak_anno.txt",sep="\t",header=F)
split_b<-str_split(data3$V6," ")
b<-sapply(split_b,"[",1)
data3$V6=b
data3=as.data.frame(table(data3$V6))
colnames(data3)=c("position","count")
data3$class="LD7"
##data4
data4=read.csv("bedops_LD28_peak_anno.txt",sep="\t",header=F)
split_b<-str_split(data4$V6," ")
b<-sapply(split_b,"[",1)
data4$V6=b
data4=as.data.frame(table(data4$V6))
colnames(data4)=c("position","count")
data4$class="LD28"
##data merge
data=rbind(data1,data2)
data=rbind(data,data3)
data=rbind(data,data4)
str(data)
data$position=as.character(data$position)
data$position[which(data$position=="3'")]=c(rep("3' UTR",4))
data$position[data$position=="5'"]=c(rep("5' UTR",4))
data$class=factor(data$class,levels = c("SD28","LD3","LD7","LD28"))
data$position=factor(data$position,levels = unique(data$position))

sum(data1$count)
sum(data2$count)
sum(data3$count)
sum(data4$count)
##plot
ggplot(data=data, aes(x=class, y=count,fill = position)) + 
  geom_bar(stat = "identity")+
  scale_fill_npg()+
  #scale_fill_manual(values=mycol)+
  ylab("Number and genomic distribution of peaks identified by ATAC-seq")+
  xlab("Sample")+theme(panel.background = element_blank(),panel.border = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.key = element_blank(),legend.title = element_blank(),
                       axis.text.x = element_text(size=12,angle=60, vjust=0.5,hjust=0.5,
                       ))
ggsave("merged_peak number and distribution.pdf",height = 4,width = 3.5)



data=as.data.frame(table(data$annotation))
colnames(data)=c("position","count")
data$class="Sample"
data$ratio=round(data$count/sum(data$count),3)
data<- arrange(data,data[,2])
data$position=factor(data$position,levels = unique(data$position))
blank_theme <- theme_minimal()+
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
ggplot(data=data, aes(x=class, y=ratio,fill = position)) + 
  geom_bar(stat = "identity")+geom_text(aes(label=..count..),stat = 'count')+
  scale_fill_npg(labels = paste(data$position," (",scales::percent(data$ratio),") ",sep=""))+
  coord_polar(theta = "y")+blank_theme+
  ylab("Number andgenomic distribution of peaks identified by ATAC-seq")+
  xlab("Sample")
ggsave(
  "C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\peak_distribution_static.pdf",
  height = 5,width = 3.5)

rm(list=ls())
dat1 = read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\early_target_gene.csv",row.names = 1)
dat2 = read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD28\\SD28_LD28_DOWN_DEP.txt",sep="\t")
dat3 = read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD28\\SD28_LD28_UP_DEP.txt",sep="\t")
dat2=rbind(dat2,dat3)
setdiff(rownames(dat1),rownames(dat2))
BiocManager::install("VennDiagram")
write.table(intersect(rownames(dat1),rownames(dat2)),
            "C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\LD7_LD28_common.txt",
            quote=F,sep="\t",row.names = F)


library(stringr)
library(ggplot2)
library(dplyr)
dat = read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\SD28_LD7_GO.txt",header=T,sep="\t")
colnames(dat)
#按照gene number从低到高排序[升序]
dat<- arrange(dat,dat[,8])
#Pathway列最好转化成因子型，否则作图时ggplot2会将所有Pathway按字母顺序重排序
#将Pathway列转化为因子型
dat_cluster1=dat[dat$cluster=="cluster1",]
dat_cluster1<- arrange(dat_cluster1,dat_cluster1[,5])
dat_cluster2=dat[dat$cluster=="cluster2",]
dat_cluster2<- arrange(dat_cluster2,dat_cluster2[,5])
dat_cluster3=dat[dat$cluster=="cluster3",]
dat_cluster3<- arrange(dat_cluster3,dat_cluster3[,5])
dat_cluster4=dat[dat$cluster=="cluster4",]
dat_cluster4<- arrange(dat_cluster4,dat_cluster4[,5])
dat_cluster1$Description <- factor(dat_cluster1$Description,levels = unique(dat_cluster1$Description))
dat_cluster2$Description <- factor(dat_cluster2$Description,levels = unique(dat_cluster2$Description))
dat_cluster3$Description <- factor(dat_cluster3$Description,levels = unique(dat_cluster3$Description))
dat_cluster4$Description <- factor(dat_cluster4$Description,levels = unique(dat_cluster4$Description))
tmp=rbind(dat_cluster4,dat_cluster3)
tmp=rbind(tmp,dat_cluster2)
tmp=rbind(tmp,dat_cluster1)

p <- ggplot(tmp,aes(y=Description,x=cluster,color=-LogP,size=GeneInGOAndHitList)) + 
  geom_point() +
  theme_bw() +
  scale_color_gradient(low='blue',high='Red')+ 
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(size = 12),
        legend.position="right",
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=15,vjust = 0.05,angle =  30),
        axis.text.y = element_text(size=15),
        axis.title.x =NULL ,
        axis.title.y =NULL,
        legend.key.size = unit(20, "pt"))+
        xlab(NULL)+
  scale_size("size_area", range = c(0, 8))+
  scale_y_discrete(position = "right")
p

ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\ATAC\\Figure\\GO_cluster.pdf",height = 6,width=9)


rm(list=ls())
library(stringr)
library(ggplot2)
library(dplyr)
dat = read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\SD28_LD7_GO.txt",header=T,sep="\t",row.names = 1)
colnames(dat)
#按照gene number从低到高排序[升序]
dat<- arrange(dat,dat[,1])
dat$Name=rownames(dat)
dat$Number=factor(dat$Number,levels = c(4,5,7,8,10,14))
dat$Name=factor(dat$Name,levels = dat$Name)
dat$raw_P_value=-log(dat$raw_P_value)
p <- ggplot(dat,aes(x =dat$Name,y=dat$Number)) + 
  coord_flip()+
  geom_col(aes(fill=raw_P_value))+
  scale_fill_gradient(low = "blue", high = "red")+
  xlab(NULL)+
  ylab("Gene Number")+
  theme_bw() +
  scale_color_gradient(low='blue',high='Red')+ 
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(size = 12),
        legend.position="right",
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=15,angle =  30),
        axis.text.y = element_text(size=15),
        axis.title.x =NULL ,
        axis.title.y =NULL,
        legend.key.size = unit(20, "pt"))

p
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\quail_ATAC_info\\RNA\\SD28_LD7\\SD28_LD7_GO.pdf",
       height = 3,
       width=10)


library(ggplot2)
library(ggsci)
data=data.frame(gene=c("DIO2","DIO2","DIO2","DIO2","DIO2","DIO3","DIO3","DIO3","DIO3","DIO3"),
                count=c(19.12239037,	9.764849314,	8.312591303,	40.8844515,	54.62884788,47.0476496,	27.01958982	,40.5521114	,1.463050463,	1.352067711),
            class=c("SD28","SD28","SD28","LD7","LD7","SD28","SD28","SD28","LD7","LD7"))


data$class=factor(data$class,levels = c("SD28","LD7","LD28"))
# 基本箱线图

p = ggplot(data, aes(x=gene, y=count,fill=class)) + 
  geom_boxplot()+theme(axis.line.x=element_line(linetype=1,color="black",size=1),
                       axis.line.y=element_line(linetype=1,color="black",size=1),
                       panel.grid = element_blank(), 
                       axis.text.x=element_text(angle=30,hjust = 0.5,vjust=0),
                       panel.background =  element_rect(fill="white"))+scale_fill_npg()
p
BiocManager::install("DiffBind")
library("DiffBind")
  setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\")
  dbObj <- dba(sampleSheet="sample_info.csv")
  dbObj_count <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
  pdf("Goose_PCA_Pearson.pdf")
  dba.plotPCA(dbObj,  attributes=DBA_CONDITION,vcolors=c("blue"))
  plot(dbObj)
  dev.off()
  dbObj_contrast <- dba.contrast(dbObj_count, categories=DBA_CONDITION,minMembers = 4)
  dbObj_diff <- dba.analyze(dbObj_contrast, method=DBA_ALL_METHODS,)
  comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
  comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
out <- as.data.frame(comp1.edgeR)  
write.table(out, file="Nanog_vs_Pou5f1_edgeR.txt", sep="\t", quote=F, col.names = NA)
out <- as.data.frame(comp1.deseq) 
write.table(out, file="Nanog_vs_Pou5f1_deseq.txt", sep="\t", quote=F, col.names = NA)
# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
dim(edge.bed[which(edge.bed$Fold<=0),])
dim(edge.bed[which(edge.bed$Fold>=0),])
write.table(edge.bed, file="Nanog_vs_Pou5f1_edgeR_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
dim(deseq.bed[which(deseq.bed$Fold<=0),])
dim(deseq.bed[which(deseq.bed$Fold>=0),])
write.table(deseq.bed, file="Nanog_vs_Pou5f1_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

setwd("E:\\")
library(Seurat)
library(SeuratData)
BiocManager::install("SeuratDisk")

remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
InstallData("stxBrain")
InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final
SaveH5Seurat(pbmc3k.final, filename = "pbmc3k.h5Seurat")
Convert("pbmc3k.h5Seurat", dest = "h5ad")


seurat_obj <- LoadH5Seurat("E:\\pbmc3k.h5Seurat")
#替换列名

colnames(seurat_obj@meta.data)[7]<-"nCount_RNA"
colnames(seurat_obj@meta.data)[8]<-"nFeature_RNA"
colnames(seurat_obj@meta.data)[9]<-"percent.mt"