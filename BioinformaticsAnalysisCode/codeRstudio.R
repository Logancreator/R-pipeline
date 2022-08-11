#install.packages("FactoMineR")
#鑾峰緱FPKM/RPKM璁＄畻闇�瑕佺殑鍩哄洜闀垮害锛堣�冭檻exon涔嬮棿鐨刼verlap锛?
library(data.table)
library("IRanges")
require("rtracklayer")

gtf <- readGFF("C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gtf")

anno <- setDT(gtf)

anno <- anno[type=="exon",]

setnames(anno,c("seqid","start","end","gene","exon_number"),c("Chr","ExonStart","ExonEnd","Gene","Exon_number"))  ##杩欓噷瑕佹敞鎰忎竴涓?

Exon_region <- unique(anno[,.(Chr,ExonStart,ExonEnd,Exon_number,Gene)])

Exon_region <- Exon_region[,{x <- IRanges(ExonStart,ExonEnd);y <- reduce(x); list(ExonStart=y@start,ExonEnd=y@start+y@width-1)},by=.(Gene,Chr)]

Exon_region[,Exon_num:=1:.N,by=Gene]

Exon_region <- Exon_region[,.(Chr,ExonStart,ExonEnd,Exon_num,Gene)]

Exon_len <- Exon_region[,.(ExonLen = ExonEnd - ExonStart + 1),by=.(Exon_num,Gene)]

gene_len <- as.data.frame(Exon_len[,.(Length = sum(ExonLen)),by=Gene])

##fwrite(Exon_region,file="/md01/changjy/data/quail/ref/ncbi_dataset/All_quailgene_exon.bed", sep = "\t", col.names = T)
##fwrite(gene_len, file = "/md01/changjy/data/quail/ref/ncbi_dataset/All_quailgene_len.txt", sep = "\t", col.names = T)

#璁＄畻RPKM

data=read.table("C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\RNA-seq\\count.txt",row.names=1,header=T)

data=data[,c(6:dim(data)[2])]

colnames(data)= c("B", "B", "B", "L", "L", "L")

dim(data)

dim(gene_len)

length(intersect(rownames(data),gene_len$Gene))

gene_len=gene_len[gene_len$Gene %in% intersect(rownames(data),gene_len$Gene),]

dim(gene_len)

head(gene_len)

head(data)

fpkm<-array(0,dim=c(length(rownames(data)),length(colnames(data))))

for (i in c(1:length(rownames(data)))){
  for (j in c(1:length(colnames(data)))){
    fpkm[i,j]=data[i,j]/gene_len[gene_len$Gene==rownames(data)[i],2][1]/sum(data[,j])*10^9
  }
}

colnames(fpkm)=colnames(data)

rownames(fpkm)=rownames(data)

fpkm 

write.csv(fpkm,"C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\RNA-seq\\Count_matrix_rpkm.csv",quote=F,row.names=T)


library(DiffBind)
library(ChIPpeakAnno)
setwd("/public/home/changjianye/project/magang/ATAC")
ATAC <- dba(sampleSheet = "sample_info.csv")
ATAC1 <- dba.count(ATAC,bUseSummarizeOverlaps=TRUE)

dba.plotPCA(ATAC1)
ATAC1
counts <- dba.peakset(ATAC1, bRetrieve=TRUE, writeFile="count.csv")
counts

ATAC2 <- dba.contrast(ATAC1, categories = DBA_CONDITION, minMembers = 4)
ATAC2
ATAC3 <- dba.analyze(ATAC2)


	library(DESeq2)
	library(ggplot2)
 												"L2P1_L1_809F10",
                            "L2P2_L1_810F10",
                            "L5P1_L1_811F10",
                            "L5P2_L1_812F10",
                            "B1P1_L2_801D73",
                            "B1P2_L1_802D73",
                            "B5P1_L4_806D73",
                            "B5P2_L2_808D73")
	colData <- data.frame(condition=factor(c("L","L","L","L",
"B","B","B","B")))
	dds<-DESeqDataSetFromMatrix(countsTable[,c(4:11)],colData, formula(~condition))  
  dds<- estimateSizeFactors(dds)
	dds<- estimateDispersions(dds)
	dds<- DESeq(dds)
	normalized_counts <- counts(dds, normalized=TRUE)
	rld <- rlog(dds, blind = FALSE)
    #杩斿洖鏍锋湰鍚嶅拰鍒嗙粍
    pcaData <- plotPCA(rld, intgroup=c("condition"),returnData = T)
	res=results(dds,contrast=c("condition","B","L"))
	summary(res)
	res <- res[order(res$padj),]
	all_info=merge(countsTable,as.data.frame(res),by="row.names",all.x=TRUE)
	head(all_info)
	write.table(all_info,
              "Information_Lay_Breed_DESEQ2.txt",
              quote=F,
              sep="\t",
              row.names=F) setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\deseq2")
	my_peak=c()
	countsTable <- read.table( "Lay_Breed_merged.count", stringsAsFactors=TRUE )
	rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
	colnames(countsTable)<-c("chr","start","end",
		
	sig=all_info[ which(all_info$padj <=0.05),-1]
	write.table(sig,"Lay_Breed_DESEQ2.txt",quote=F,
	            sep="\t",row.names=F)
  sigUp=sig[ which(all_info$log2FoldChange >= 0.8 ),-1]
  write.table(sigUp,"Lay_Breed_DESEQ2_up.txt",quote=F,
              sep="\t",row.names=F)
  sigDown=sig[ which(all_info$log2FoldChange <= -0.8 ),-1]
  write.table(sigDown,"Lay_Breed_DESEQ2_down.txt",quote=F,
              sep="\t",row.names=F)
#BiocManager::install("DiffBind",version = "3.11")

library(DiffBind)
library(ChIPpeakAnno)
library(pheatmap)
setwd("/public/home/changjianye/project/magang/DiffBindv4/")

dbObj <- dba(sampleSheet="/public/home/changjianye/project/magang/ATAC/sample_info.csv")

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

pdf("DiffBindv4.pdf")

#plot
dba.plotPCA(dbObj, attributes=DBA_CONDITION,label=DBA_ID)

#DBA_SCORE_READS

#DBA_SCORE_NORMALIZED

#DBA_SCORE_CONTROL_READS

#DBA_SCORE_READS_MINUS

#DBA_SCORE_READS_FOLD

#DBA_SCORE_RPKM

#DBA_SCORE_RPKM_FOLD

#DBA_SCORE_RPKM_MINUS

plot(dbObj)

# Establishing a contrast 
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION,minMembers = 4)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)

#  summary of results
dba.show(dbObj, bContrasts=T)

#  overlapping peaks identified by the two different tools (DESeq2 and edgeR)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
write.table(as.data.frame(comp1.deseq), file="Lay_vs_Breed_deseq2_all.txt", sep="\t", quote=F, row.names = F , col.names = F)
write.table(as.data.frame(comp1.edgeR), file="Lay_vs_Breed_edgeR_all.txt", sep="\t", quote=F, row.names = F , col.names = F)


# edgeR
out <- as.data.frame(comp1.edgeR)
out[which(out$start <=0),"start"] = 1
out <- out[which(abs(out$Fold)>=1&out$FDR<=0.05),]
write.table(out, file="Lay_vs_Breed_edgeR.txt", sep="\t", quote=F, row.names = F , col.names = F)

# Create bed files for each keeping only significant peaks (p < 0.05)
# edgeR
out <- as.data.frame(comp1.edgeR)
out[which(out$start <=0),"start"] = 1
out <- out[which(abs(out$Fold)>=1&out$FDR<=0.05),
      c("seqnames", "start", "end", "strand", "Fold")]
edgeR.bed <- out
write.table(edgeR.bed, file="Lay_vs_Breed_edgeR_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(edgeR.bed[which(deseq.bed$Fold>=1),], file="Lay_vs_Breed_edgeR_up.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(edgeR.bed[which(deseq.bed$Fold<=-1),], file="Lay_vs_Breed_edgeR_down.bed", sep="\t", quote=F, row.names=F, col.names=F)


# DESeq2
out <- as.data.frame(comp1.deseq)
out[which(out$start <=0),"start"] = 1
out <- out[which(abs(out$Fold)>=1&out$FDR<=0.05),]
write.table(out, file="Lay_vs_Breed_deseq2.txt", sep="\t", quote=F, row.names = F , col.names = F)

# Create bed files for each keeping only significant peaks (p < 0.05)
# DESeq2
out <- as.data.frame(comp1.deseq)
out[which(out$start <=0),"start"] = 1
out <- out[which(abs(out$Fold)>=1&out$FDR<=0.05),
			c("seqnames", "start", "end", "strand", "Fold")]
deseq.bed <- out
write.table(deseq.bed, file="Lay_vs_Breed_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(deseq.bed[which(deseq.bed$Fold>=1),], file="Lay_vs_Breed_deseq2_up.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(deseq.bed[which(deseq.bed$Fold<=-1),], file="Lay_vs_Breed_deseq2_down.bed", sep="\t", quote=F, row.names=F, col.names=F)


df=as.data.frame(dbObj$binding)
deseq.bed=as.data.frame(deseq.bed)
head(rownames(df))
head(rownames(deseq.bed))
data=df[rownames(df)%in%rownames(deseq.bed),]

dim(data)
pheatmap(data[,c(4:11)],scale="row",show_rownames=F)
dev.off()









#BiocManager::install("ggthemes")
rm(list=ls())
library(ggplot2)
library(ggsci)
library(ggthemes)
Lay=read.csv(
  "C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\ATAC-seq\\QC\\Lay_tss.tsv",
            header = F,fileEncoding = "UTF8",sep="\t")
Lay$Sample= c("Lay")
Breed=read.csv(
  "C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\ATAC-seq\\QC\\Breed_tss.tsv",
            header = F,fileEncoding = "UTF8",sep="\t")
Lay$Sample= c("Lay")
Breed$Sample= c("Breed")
data=rbind(Lay,Breed)

data$Sample=factor(data$Sample,levels = c("Lay","Breed"))
p=ggplot(data = data, mapping = aes(x = V1, y = V2, colour = Sample)) + 
  geom_line(size=1,alpha = 1)+
  scale_color_tableau()+
  theme_bw()+
  xlab("Distance to Transcription Start Site (bp) ")+
  ylab("Normalized Insertion")+
  xlim(-1000,1000)+
  ylim(1,5)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank())
  #scale_y_continuous(breaks=c(1:6))
p

ggsave("C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\ATAC-seq\\QC\\TSS_plot.pdf"
       ,height = 3,width=4)



library(dplyr)
library(grid)
library(ggplot2)
library(gridExtra)
library(Rsamtools)
library(GenomicRanges)
options(stringsAsFactors=FALSE)
rm(list = ls())
gtf <- read.table("/public/home/changjianye/project/duck/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gtf", sep="\t")
gtf_sub <- gtf[which(gtf[, 3] == "gene"),]
#genes <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
#gene_id=gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf_sub[, 9])))
x=strsplit(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),split=";")
genes_id=as.character(lapply(x, function(x) x[1]))
genes <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
                 gene_id=genes_id,symbol=genes_id)

genes
gtf_sub <- gtf[which(gtf[, 3] == "exon"),]
#exons <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
#gene_id=gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf_sub[, 9])))

x=strsplit(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),split=";")
genes_id=as.character(lapply(x, function(x) x[1]))
exons <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
                 gene_id=genes_id,symbol=genes_id)
exons


gtf_sub <- gtf[which(gtf[, 3] == "start_codon"),]
#TSS <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 4]), strand=Rle(strand(gtf_sub[, 7])),
#gene_id=gsub("\\..*", "", gsub(".*transcript_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*transcript_name ", "", gtf_sub[, 9])))
x=strsplit(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),split=";")
genes_id=as.character(lapply(x, function(x) x[1]))
TSS<-GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
           gene_id=genes_id,symbol=genes_id)

geneAnnotation <- SimpleList(genes=genes, exons=exons, TSS=TSS)
saveRDS(geneAnnotation, "/public/home/changjianye/project/duck/gencode.ct22.annotation.rds")

geneAnnotation$genes
geneAnnotation$exons
geneAnnotation$TSS


#library(stringr)
#str_sub(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),1,str_locate(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),"\\;")[2])


get_polygon <- function(x, y, w, h, arrow="*", rate=1) {
  if (arrow == "-")
  {
    res <- data.frame(x=x, y=y)
    res <- rbind(res, c(x+min(w, rate*h/2), y+h/2))
    res <- rbind(res, c(x+w, y+h/2))
    res <- rbind(res, c(x+w, y-h/2))
    res <- rbind(res, c(x+min(w, rate*h/2), y-h/2))
    return(res)
  }
  res <- data.frame(x=x, y=y+h/2)
  if (arrow == "+") {
    res <- rbind(res, c(x+max(0, w-rate*h/2), y+h/2))
    res <- rbind(res, c(x+w, y))
    res <- rbind(res, c(x+max(0, w-rate*h/2), y-h/2))
  } else {
    res <- rbind(res, c(x+w, y+h/2))
    res <- rbind(res, c(x+w, y-h/2))
  }
  res <- rbind(res, c(x, y-h/2))
  return(res)
}


geneAnnotation <- readRDS("/public/home/changjianye/project/duck/gencode.ct22.annotation.rds")
file_list <- gsub(".bam", "", list.files(path="/public/home/changjianye/project/magang/publish/", pattern="*.bam$"))
file_list

read_count <- data.frame(Sample=file_list, Count=0)
read_count
for (i in 1:length(file_list)) read_count$Count[i] <- sum(read.delim(paste0(file_list[i], "_chr.tsv"), h=F)[, 2])
read_base <- min(read_count$Count)
read_count$Rate <- floor(log10(read_count$Count/read_base))
read_count$Rate[which(read_count$Rate == 0)] <- 1
read_count
read_base

gene_list <- c("PRL","TSHB","NR5A1","CGA","VIT","CETP","TSHR")
data_info <- data.frame(ID=gene_list, Symbol=gene_list)

tileSize <- 5
expand <- 3000
imgWidth <- 12
for (geneSymbol in gene_list)
{
  print(geneSymbol)
  region <- geneAnnotation$genes[which(mcols(geneAnnotation$genes)$gene_id == geneSymbol)]
  chr <- names(table(as.character(seqnames(region))))
  strand <- names(table(as.character(strand(region))))
  if (length(chr) != 1 | length(strand) != 1) print("multiple terms")
  if(strand == "+") start(region) <- start(region)-expand else end(region) <- end(region)+expand
  pos <- c(trunc(max(min(start(region))-tileSize, 0)/tileSize)*tileSize, trunc((max(end(region))+tileSize)/tileSize)*tileSize)
  subrng <- GRanges(seqnames=chr, ranges=IRanges(pos[1], end=pos[2]))
  regionTiles <- seq(trunc(pos[1]/tileSize), trunc(pos[2]/tileSize))*tileSize
  
  frag_list <- SimpleList()
  for (file in paste0(file_list, ".bam")) frag_list <- c(frag_list, scanBam(BamViews(file, bamRanges=subrng)))
  group_mat <- matrix(0, nrow=length(regionTiles), ncol=length(file_list))
  for (i in 1:length(file_list))
  {
    ts <- match(trunc(frag_list[[i]][[1]]$pos/tileSize)*tileSize, regionTiles, nomatch = 0)
    ids <- which(ts > 0)
    te <- match(trunc((frag_list[[i]][[1]]$pos+frag_list[[i]][[1]]$isize)/tileSize)*tileSize, regionTiles, nomatch = 0)
    ide <- which(te > 0)
    for (nb in c(ts[ids], te[ide])) group_mat[nb, i] <- group_mat[nb, i] + 1
  }

  df <- data.frame(which(group_mat > 0, arr.ind=TRUE))
  df$y <- group_mat[cbind(df[,1], df[,2])]
  dfm1 <- df
  dfm1$row <- dfm1$row - 1
  dfm1$y <- 0
  dfp1 <- df
  dfp1$row <- dfp1$row + 1
  dfp1$y <- 0
  df <- rbind(df, dfm1, dfp1, data.frame(row=rep(c(1, length(regionTiles)), each=ncol(group_mat)), col=rep(1:ncol(group_mat), 2), y=0))
  df <- df[!duplicated(df[,1:2]),]
  df <- df[df$row > 0 & df$row < (length(regionTiles)+1),]
  df$x <- regionTiles[df$row]
  df$group <- file_list[df$col]
  df <- df[order(df$group, df$x),]
  df <- df[,c("x", "y", "group")]
  df$y <- round(df$y/read_count$Rate[match(df$group, read_count$Sample)])
  df$group <- factor(df$group, levels=file_list)
  #ylim <- c(0, quantile(df$y, probs=c(0.999)))
  #df$y[df$y < ylim[1]] <- ylim[1]
  #df$y[df$y > ylim[2]] <- ylim[2]
  #write.csv(df, paste0("frag_test_", groupBy, "_", geneSymbol, ".csv"))

  gene_info <- geneAnnotation$genes[which(mcols(geneAnnotation$genes)$gene_id == geneSymbol)]
  mcols(gene_info)$gene_id <- paste(mcols(gene_info)$gene_id, 1:length(gene_info), sep=".")
  exon_info <- data.frame()
  for (i in 1:length(gene_info))
  {
    exon_sub <- reduce(subsetByOverlaps(geneAnnotation$exons, gene_info[i]))
    exon_sub <- data.frame(start=start(exon_sub), end=end(exon_sub), id=mcols(gene_info)$gene_id[i], strand="*")
    if (strand == "+") exon_sub$strand[which.max(exon_sub$end)] <- "+"
    if (strand == "-") exon_sub$strand[which.min(exon_sub$start)] <- "-"
    exon_info <- rbind(exon_info, exon_sub)
  }
  gene_info <- data.frame(start=start(gene_info), end=end(gene_info), id=mcols(gene_info)$gene_id, strand=strand)

  term <- data_info$Symbol[match(geneSymbol, data_info$ID)]
  tss_min <- 0
  tss_max <- 0
  if (strand == "+") tss_min <- min(gene_info$start)-1000 else tss_min <- min(gene_info$end)-1000
  if (strand == "+") tss_max <- max(gene_info$start)+1000 else tss_max <- max(gene_info$end)+1000
  tss_min <- max(tss_min, pos[1])
  tss_max <- min(tss_max, pos[2])
  pa <- ggplot(df, aes(x=x, y=y))+geom_area(stat="identity", colour="gray20", fill="gray20")+facet_wrap(facets= ~ group, strip.position='left', ncol=1)+
    scale_x_continuous(expand=expansion(), limits=c(pos[1], pos[2]))+guides(fill=FALSE, colour=FALSE)+
    geom_vline(xintercept=tss_min, colour="gray60", linetype="dashed", size=1)+
    
    geom_vline(xintercept=tss_max, colour="gray60", linetype="dashed", size=1)+
    theme(panel.background=element_blank(), panel.spacing=unit(0, "lines"), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
          strip.text.y.left=element_text(size=12, color="black", angle=0, hjust=1), strip.background=element_blank())

  bg <- data.frame(get_polygon(pos[1], 0, pos[2]-pos[1], nrow(gene_info)), group="Gene")
  pc <- ggplot(bg, aes(x=x, y=y))+geom_polygon(color="white", fill="white", alpha=0)+facet_wrap(facets= ~ group, strip.position='left', ncol=1)+
    scale_x_continuous(expand=expansion(), limits=c(pos[1], pos[2]), breaks=c(tss_min, tss_max))+
    labs(x=paste0("\nATAC Signal Coverage (0-", max(df$y), ")\nRelated to ", term, " (", chr, ": ", pos[1], "-", pos[2], ")"))+
    geom_vline(xintercept=tss_min, colour="red", linetype="dashed", size=1)+
    #geom_vline(xintercept=6187792, colour="red", linetype="dashed", size=1)
  #geom_vline(xintercept=6179270, colour="red", linetype="dashed", size=1)+
    geom_vline(xintercept=tss_max, colour="red",linetype="dashed", size=1)+
    theme(panel.background=element_blank(), panel.spacing=unit(0, "lines"), axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_text(size=15, color="black"), axis.line.x=element_line(linetype=1,colour="black"), axis.text.x=element_text(size=12, colour="black"),
          strip.text.y.left=element_text(size=12, color="black", angle=0, hjust=1), strip.background=element_blank())
  for (i in 1:nrow(gene_info)) pc <- pc+geom_polygon(data=get_polygon(gene_info$start[i], i-nrow(gene_info)/2-0.5, gene_info$end[i]-gene_info$start[i], 0.1), 
                                                     aes(x=x, y=y), color="gray20", fill="gray20")
  eid <-match(exon_info$id, gene_info$id)
  xy_rate <- (pos[2] - pos[1])/50
  for (i in 1:nrow(exon_info)) pc <- pc+geom_polygon(data=get_polygon(exon_info$start[i], eid[i]-nrow(gene_info)/2-0.5, exon_info$end[i]-exon_info$start[i], 0.8, 
                                                                      exon_info$strand[i], xy_rate), aes(x=x, y=y), color="gray20", fill="gray20")


  imageHeight <- imgWidth*(length(file_list)+1+nrow(gene_info)/2)/20
  ggsave(plot=patchwork::wrap_plots(A=pa, B=pc, ncol=1, heights=c(length(file_list)*2, nrow(gene_info))), 
         width=imgWidth, height=imageHeight, dpi=600, filename=paste0("bulk_sample_", term, "_", geneSymbol, ".pdf"), limitsize=F)

}



library(ggplot2)
library("ggsci")
res=read.table("C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\ATAC-seq\\diffbind\\Lay_vs_Breed_deseq2_all.txt",header=T)
head(res,1)
res=na.omit(res)
res[which(res$Fold >= 1 & res$FDR <= 0.01),'sig'] <- 'UP'
res[which(res$Fold <= -1 & res$FDR< 0.01),'sig'] <- 'DOWN'
res[which(!(res$sig%in%c("UP","DOWN"))),'sig'] <- 'NO DIFF'
head(res,2)
volcano <- ggplot(res, aes(Fold, -log(FDR, 10))) +
 theme(plot.title = element_text(size=2,hjust = 0.5))+
  geom_point(aes(color = sig), size = 0.6) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'white'), 
  ) +
  theme(legend.title = element_blank(),text=element_text(size=8,  family="serif")) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 padj', color = NA) +
  scale_color_manual(
    # 鎸囧畾鏍囪壊
    values =
      c("steelblue", "grey60", "orangered")
    )

volcano
ggsave("C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\ATAC-seq\\diffbind\\Volcano.pdf")

BiocManager::install("JASPAR2020")
library(motifmatchr)
library(GenomicRanges)
library(JASPAR2020)
library(ChIPseeker)
species <- "Homo sapiens"
collection <- "CORE"
opts <- list()
opts["species"] <- species
opts["collection"] <- collection
out  <- TFBSTools::getMatrixSet(JASPAR2020, opts)
 
if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
  names(out) <- paste(names(out), TFBSTools::name(out), 
                      sep = "_")
motif <- out
motif$MA1540.1_NR5A1
peak <- readPeakFile("C:\\Users\\Jeff\\OneDrive\\妗岄潰\\peak.bed")
motif_ix <- matchMotifs(motif$MA1540.1_NR5A1, peak)
?matchMotifs
a=Biostrings::readDNAStringSet("C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.fa")



library(ggplot2)
library(ggpubr)
library(RColorBrewer)
lwd_pt <- .pt*72.27/96
mycol= brewer.pal(n = 5, name = "YlGnBu")[c(2:5)]
setwd("C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\ATAC-seq\\insertFrag")
Breed=read.table("Breed_fragmentSize_shift.txt",header=T)
Lay=read.table("Lay_fragmentSize_shift.txt",header=T)
head(Breed)
head(Lay)
data=rbind(Breed,Lay)
head(data)
ggplot(data,aes(x =Size,y = Occurrences ,color=Sample))+
  geom_line(size=1)+
  xlab("Fragment length(bp)")+
  xlim(0,800)+
  scale_color_manual(values = mycol)+
  ylab(expression(Normalized ~ read ~ count ~ 10^2))+
  ggtitle("Sample Fragment sizes")+
  theme_set(theme_bw(base_size = 10, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))
ggsave("Fragment_size.pdf")

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
mycol= brewer.pal(n = 5, name = "YlGnBu")[c(2:5)]
setwd("C:\\Users\\Jeff\\OneDrive\\妗岄潰\\Goose_ATAC_info\\ATAC-seq\\bedops")
distance_breed=read.csv("Breed.narrowPeak.anno",header=F,sep="\t")
distance_lay=read.csv("Lay.narrowPeak.anno",header=F,sep="\t")

distance_breed=data.frame(distance=distance_breed[c(1:dim(distance_breed)[1]),c("V13")],sample=rep("Breed",dim(distance_breed)[1]))
head(distance_breed)
dim(distance_breed)

distance_lay=data.frame(distance=distance_lay[c(1:dim(distance_lay)[1]),c("V13")],sample=rep("Lay",dim(distance_lay)[1]))
head(distance_lay)
dim(distance_lay)
dat=rbind(distance_breed,distance_lay)
colnames(dat)=c("distance","sample")


ggplot(dat,aes(x=abs(distance)/1000,color=sample,fill=sample))+
geom_histogram(bins = 50)+
xlab('Distance to TSS')+
ylab('Density')+
xlim(0,500)+
ylim(0,2000)+
theme_classic()
ggsave('distanceToTss.density.pdf',dpi = 1080,width=4,height=2.5)





ggplot(dat,aes(x=abs(distance),color=sample))+
geom_line()+
xlab('Distance to TSS (Lay)')+
ylab('Density')
ggsave('distanceToTss.density.pdf',dpi = 1080)


  rm(list=ls())
	library(DESeq2)
	library(ggplot2)
  library(RColorBrewer)
  library(amap)
  #library(gplot)
  #BiocManager::install("gplot")
  library(pheatmap)
  setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\7.10run\\")
	my_peak=c()
	countsTable <- read.table( "count.txt", stringsAsFactors=TRUE )
	countsTable[,c(4:11)]=round(countsTable[,c(4:11)]/colSums(countsTable[,c(4:11)])*mean(colSums(countsTable[,c(4:11)])))
  
	countsTable[,c(4:11)]
	colSums(countsTable[,c(4:11)])
	mean(colSums(countsTable[,c(4:11)])))
	dim(countsTable)
	head(countsT)
	rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
	colnames(countsTable)<-c("chr","start","end",
"L2P1","L2P2","L5P1","L5P2","B1P1","B1P2","B5P1","B5P2")
	
	countsTable["peak_1174",]
	
	colData <- data.frame(condition=factor(c("L","L","L","L","B","B","B","B")))
	dds<-DESeqDataSetFromMatrix(countsTable[,c(4:11)],colData, formula(~condition)) 
	dds <- estimateSizeFactors(dds)
	dds<- estimateDispersions(dds)
	dds<- DESeq(dds)
	normalized_counts <- counts(dds, normalized=TRUE)
	rld <- rlog(dds, blind = FALSE)

  pcaData <- plotPCA(rld, intgroup=c("condition"),returnData = T)
  
  normalized_counts <- counts(dds, normalized=TRUE)
  head(normalized_counts)
  normalized_counts_mad <- apply(normalized_counts, 1, mad)
  normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

  write.table(normalized_counts, file="LayBreed_trans.Count_matrix.xls.DESeq2.normalized.xls",
  quote=F, sep="\t", row.names=T, col.names=T)
  res <- results(dds,c("condition","B","L")) 
  res
  res[which(res$log2FoldChange >= 1 & res$padj <= 0.05),]
  res[which(res$log2FoldChange <=- 1 & res$padj <= 0.05),]
  boxplot(res$log2FoldChange)
  plotMA(res)
  write.table(res, file="LayBreed_res.Count_matrix.xls.DESeq2.normalized.xls",
              quote=F, sep="\t", row.names=T, col.names=T)
  
  res=as.data.frame(res)
  res=na.omit(res)
  res[which(res$log2FoldChange >= 1 & res$padj <= 0.05),]
  res[which(res$log2FoldChange  <= -1 & res$padj <= 0.05),]
  res[which(res$log2FoldChange >= 1 & res$padj <= 0.05),'sig'] <- 'UP'
  res[which(res$log2FoldChange <= -1 & res$padj<= 0.05),'sig'] <- 'DOWN'
  res[which(!(res$sig%in%c("UP","DOWN"))),'sig'] <- 'NO DIFF'
  dim(res[res$sig=="UP",])
  dim(res[res$sig=="DOWN",])
  
  
  volcano <- ggplot(as.data.frame(res), aes(log2FoldChange, -log(padj, 10))) +
    theme(plot.title = element_text(size=8,hjust = 0.5))+
    geom_point(aes(color = sig), size = 2) +
    theme_classic() +
    scale_color_manual(values = c("lightblue","gray","red"))+
    geom_vline(xintercept = c(-0.8, 0.8), color = 'gray', size = 0.25) +
    geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
    labs(x = 'log2 Fold Change', y = '-log10 padj', color = NA)
  volcano
  ggsave("LayBreed_res_Volcano.pdf",height = 4,width = 4)
  
  
  
  
  
  
  
  
  rld <- rlog(dds, blind=TRUE)
  rlogMat <- assay(rld)
  rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]

  write.table(rlogMat, file="LayBreed_trans.Count_matrix.xls.DESeq2.normalized.rlog.xls",
  quote=F, sep="\t", row.names=T, col.names=T)
  
  res <- results(dds)
  
  resSig <- subset(res, padj <= 0.05)

  resSig_up <- subset(resSig, log2FoldChange > 1)
  
  resSig_down <- subset(resSig, log2FoldChange < -1)
  write.table(cbind(countsTable[rownames(resSig),],resSig), file="LayBreed_resSig.Count_matrix.xls.DESeq2.normalized.rlog.xls",
              quote=F, sep="\t", row.names=T, col.names=T)
  write.table(cbind(countsTable[rownames(resSig_up),],resSig_up), file="LayBreed_resSig_up.Count_matrix.xls.DESeq2.normalized.rlog.xls",
              quote=F, sep="\t", row.names=T, col.names=T)
  write.table(cbind(countsTable[rownames(resSig_down),],resSig_down), file="LayBreed_resSig_down.Count_matrix.xls.DESeq2.normalized.rlog.xls",
              quote=F, sep="\t", row.names=T, col.names=T)
  
  
  pdf("LayBreed_res_Heatmap.pdf",height = 4,width = 3)
  Heatmap_Data=normalized_counts[rownames(resSig),]
  annotation_col = data.frame(
    Conditions  = factor(c(rep("Laying", 4),rep("Breeding", 4)) ))
  head(annotation_col)
  
  
  colnames(Heatmap_Data)=factor(colnames(Heatmap_Data))
  pheatmap(Heatmap_Data,
           scale = "row",
           show_rownames = F,
           cluster_cols = F,
           main = "DAR Heatmap",
           cutree_rows = 2,
           gaps_col  = 4,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           ColSideColors = rep(c("pink", "purple"), each = 4))
  
  dev.off()
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)


  pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))

  
  hc <- hcluster(t(rlogMat), method="pearson")


  pdf("LayBreed_trans.Count_matrix.xls.DESeq2.normalized.rlog.pearson.pdf", pointsize=10)
  heatmap.2(pearson_cor, Colv = as.dendrogram(hc), Rowv=as.dendrogram(hc), symm=T, trace="none",
  col=hmcol, main="The pearson correlation of each
  sample")

  dev.off()
  pdf("LayBreed_trans.Count_matrix.xls.DESeq2.normalized.rlog.PCA.pdf", pointsize=10)
  ?plotPCA(rld, intgroup=c("condition"), ntop=100)
  plotPCA(rld, intgroup=c("condition"), ntop=200)
  plotPCA(rld, intgroup=c("condition"), ntop=500)
  plotPCA(rld, intgroup=c("condition"), ntop=1000)
  plotPCA(rld, intgroup=c("condition"), ntop=2000)
  plotPCA(rld, intgroup=c("condition"), ntop=3000)
  plotPCA(rld, intgroup=c("condition"), ntop=5000)
  plotPCA(rld, intgroup=c("condition"), ntop=10000)
  dev.off()
  #Inputisamatrixoflogtransformedvalues 
  rld<-rlog(dds,blind=T) 
  rld_mat<-assay(rld) 
  pca<-prcomp(t(rld_mat))
  
  #CreatedataframewithmetadataandPC3andPC4valuesforinputtoggplot 
  df<-cbind(pca$x)
  df=as.data.frame(df)
  df$sample=c("L","L","L","L","B","B","B","B")
  ggplot(df)+geom_point(aes(x=PC3,y=PC4,color=sample)) 

  

rm(list=ls())
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(stringr)
library(dplyr)
mycol= brewer.pal(n = 7, name = "YlGnBu")
data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\bedops\\Breed.narrowPeak.anno",sep="\t",header=F)
colnames(data)=c("chr","start","end","width","strand","annotation","strands","gene_start","gene_end","gene_width","gene_strand","gene","strand")
split_b<-str_split(data$annotation," ")
b<-sapply(split_b,"[",1)
data$annotation=b
data[data$annotation=="5'","annotation"]="Proximal"
data[data$annotation=="3'","annotation"]="Distal"
data[data$annotation=="Exon","annotation"]="Genic"
data[data$annotation=="Intron","annotation"]="Genic"
data[data$annotation=="Distal","annotation"]="Distal"
data[data$annotation=="Downstream","annotation"]="Genic"
data[data$annotation=="Distal","annotation"]="Distal"
data[data$annotation=="Promoter","annotation"]="Proximal"
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
    plot.title=element_text(size=5, face="bold")
  )
data %>%ggplot( aes(x=class, y=count,fill = position)) + 
  geom_bar(stat = "identity")+
  #scale_fill_npg(labels = paste(data$position," (",scales::percent(data$ratio),") ",sep=""))+
  blank_theme+
  ylab("Number andgenomic distribution of peaks identified by ATAC-seq")+
  xlab("Sample")
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\bedops\\Breed.narrowPeak_distribution.pdf",width=4,height = 4)




data=read.csv("C:/Users/Jeff/OneDrive/桌面/peak_number.csv")
data
library(ggplot2)
library(RColorBrewer)
col=brewer.pal(12,"YlGnBu")
p<-ggplot(data,aes(x = factor(sample,levels = c("SD28_1","SD28_2","SD28_3_1","SD28_3_2","LD3_1","LD3_2_1","LD3_2_2","LD3_3","LD7_1","LD7_2","LD7_3","LD28_1","LD28_2","LD28_3")),y =number, fill=factor(condition,levels = c("SD28","LD3","LD7","LD28"))))+theme_set(theme_bw())+
  geom_bar(stat = "identity" ,width=0.6,position=position_dodge(width = 0.8))+
  theme(panel.background = element_blank(),panel.border = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.key = element_blank(),legend.title = element_blank(),
        axis.text.x = element_text(size=6,angle=60, vjust=0.5,hjust=0.5,
        ))+  xlab("")+ylab("Number of Peak")+ scale_y_continuous(breaks=seq(0, 90000, 10000))

p+ geom_errorbar(aes(y=merged, ymax=merged, ymin=merged),colour="black",width=0.6,size=1)
ggsave("Peak_number.pdf")




library(ggplot2)
lwd_pt <- .pt*72.27/96
peak=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\peak\\Peaknumber.txt",header=T)
peak=peak[-9,]
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 10, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))
ggplot(peak, aes(x=class, y=number/1000,fill=class)) +
  geom_boxplot()+
  #geom_text(x = 1, y = 7, label = "fontsize = 10",  size = 10/.pt) +
  labs(title="Peak Number Statics",x="Condition", y = "Number/1000")

ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\peak\\Peaknumber_Statics.pdf",height = 3,width = 3)





rm(list=ls())
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(amap)
library(gplots)
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
install.packages('gplots')
library(pheatmap)
setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\RNA-seq\\deseq2")
countsTable <- read.table( "count.txt", stringsAsFactors=TRUE,header=T,row.names = 1 )

#rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
colnames(countsTable)<-c("chr","start","end","strand","length",
                         "B4","B5","B6","L2","L3","L5")
colData <- data.frame(condition=factor(c("B","B","B","L","L","L")))
dds<-DESeqDataSetFromMatrix(countsTable[,c(6:11)],colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
rld <- rlog(dds, blind = FALSE)

pcaData <- plotPCA(rld, intgroup=c("condition"),returnData = T)
pdf("LayBreed_res.DESeq2.normalized.PCA.pdf",height = 4,width = 4)
pcaData <- plotPCA(rld, intgroup=c("condition"))
pcaData
dev.off()
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

write.table(normalized_counts, file="LayBreed_trans.Count_matrix.xls.DESeq2.normalized.xls",
            quote=F, sep="\t", row.names=T, col.names=T)
res <- results(dds,c("condition","B","L")) 
res
write.table(res, file="LayBreed_res.Count_matrix.xls.DESeq2.normalized.xls",
            quote=F, sep="\t", row.names=T, col.names=T)

res=as.data.frame(res)
res=na.omit(res)
res[which(res$log2FoldChange >= 1 & res$padj < 0.05),'sig'] <- 'UP'
res[which(res$log2FoldChange <= -1 & res$padj< 0.05),'sig'] <- 'DOWN'
res[which(!(res$sig%in%c("UP","DOWN"))),'sig'] <- 'NO DIFF'
volcano <- ggplot(as.data.frame(res), aes(log2FoldChange, -log(padj, 10))) +
  theme(plot.title = element_text(size=5,hjust = 0.5))+
  geom_point(aes(color = sig), size = 1) +
  theme_classic() +
  scale_color_manual(values = c("lightblue","gray","red"))+
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 padj', color = NA) 

volcano
ggsave("LayBreed_res_Volcano.pdf",height = 4,width = 4)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
rld <- rlog(dds, blind=TRUE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]

pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))


hc <- hcluster(t(rlogMat), method="pearson")


pdf("LayBreed_trans.Count_matrix.xls.DESeq2.normalized.rlog.pearson.pdf", pointsize=10)
heatmap.2(pearson_cor, Colv = as.dendrogram(hc), Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, main="The pearson correlation of each
  sample")

dev.off()

pdf("LayBreed_res_Heatmap.pdf",height = 4,width = 3)
res <- results(dds,c("condition","B","L"))

resSig <- subset(res, padj <= 0.05)

resSig_up <- subset(resSig, log2FoldChange > 1)
resSig_down <- subset(resSig, log2FoldChange < -1)
Heatmap_Data=normalized_counts[rownames(resSig),]
annotation_col = data.frame(
  Conditions  = factor(c(rep("Laying", 3),rep("Breeding", 3)) ))
head(annotation_col)


colnames(Heatmap_Data)=factor(colnames(Heatmap_Data))
pheatmap(Heatmap_Data,
         scale = "row",
         show_rownames = F,
         cluster_cols = F,
         main = "DEG Heatmap",
         cutree_rows = 2,
         gaps_col  = 3,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         ColSideColors = rep(c("pink", "purple"), each = 4))

dev.off()

rm(list=ls())
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(tidyverse)
library(stringr)
library(tidyr)
#install.packages("tidyverse")
setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\peak")
#color
mycol= brewer.pal(n = 7, name = "YlGnBu")
##data1
data1=read.csv("B1P1_L2_801D73_peaks.narrowPeak.anno",sep="\t",header=F)
data1$V13=separate(data=data1,col=V13,into=c("a","b"),sep=" ")$a
data1=as.data.frame(table(data1$V13))
colnames(data1)=c("position","count")
data1$class="B1P1"
##data2
data2=read.csv("B1P2_L1_802D73_peaks.narrowPeak.anno",sep="\t",header=F)
data2$V13=separate(data=data2,col=V13,into=c("a","b"),sep=" ")$a
data2=as.data.frame(table(data2$V13))
colnames(data2)=c("position","count")
data2$class="B1P2"
##data3
data3=read.csv("B5P1_L4_806D73_peaks.narrowPeak.anno",sep="\t",header=F)
data3$V13=separate(data=data3,col=V13,into=c("a","b"),sep=" ")$a
data3=as.data.frame(table(data3$V13))
colnames(data3)=c("position","count")
data3$class="B5P1"
##data4
data4=read.csv("B5P2_L2_808D73_peaks.narrowPeak.anno",sep="\t",header=F)
data4$V13=separate(data=data4,col=V13,into=c("a","b"),sep=" ")$a
data4=as.data.frame(table(data4$V13))
colnames(data4)=c("position","count")
data4$class="B5P2"
##data1
data5=read.csv("L2P1_L1_809F10_peaks.narrowPeak.anno",sep="\t",header=F)
data5$V13=separate(data=data5,col=V13,into=c("a","b"),sep=" ")$a
data5=as.data.frame(table(data5$V13))
colnames(data5)=c("position","count")
data5$class="B1P2"
##data2
data6=read.csv("L2P2_L1_810F10_peaks.narrowPeak.anno",sep="\t",header=F)
data6$V13=separate(data=data6,col=V13,into=c("a","b"),sep=" ")$a
data6=as.data.frame(table(data6$V13))
colnames(data6)=c("position","count")
data6$class="L2P2"
##data3
data7=read.csv("L5P1_L1_811F10_peaks.narrowPeak.anno",sep="\t",header=F)
data7$V13=separate(data=data7,col=V13,into=c("a","b"),sep=" ")$a
data7=as.data.frame(table(data7$V13))
colnames(data7)=c("position","count")
data7$class="L5P1"
##data4
data8=read.csv("L5P2_L1_812F10_peaks.narrowPeak.anno",sep="\t",header=F)
data8$V13=separate(data=data8,col=V13,into=c("a","b"),sep=" ")$a
data8=as.data.frame(table(data8$V13))
colnames(data8)=c("position","count")
data8$class="L5P2"

##data merge
data=rbind(data1,data2)
data=rbind(data,data3)
data=rbind(data,data4)
data=rbind(data,data5)
data=rbind(data,data6)
data=rbind(data,data7)
data=rbind(data,data8)

data$class=factor(data$class,levels = c("B1P1","B1P2","B5P1","B5P2",
                                        "L2P1","L2P2","L5P1","L5P2"))
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
ggsave("peak_distribution_static.pdf")



library(pheatmap)
setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\")
dap=read.csv("DAP-seq\\NR5A1.annotation.csv",header=F)
rna=read.csv("RNA-seq\\deseq2\\count__rpkm.csv",header=T,row.names = 1)
head(dap)
head(rna)

up=read.csv("RNA-seq\\deseq2\\degup.csv",header=T,row.names = 1)
down=read.csv("RNA-seq\\deseq2\\degdown.csv",header=T,row.names = 1)
res=c(rownames(up),rownames(down))

data=na.omit(rna[dap$V29,])

data=data[which(rowSums(data) > 0),]

pheatmap(data,scale="row",cluster_cols = F)

intersect(res,dap$V29)

BiocManager::install("ggh4x")
library("ggplot2")
library("ggh4x")
library("cowplot")
library("cowplot")
data<-read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\DiffBindv4\\reshape\\TF_info.txt",sep="\t",header = T)
data$class=factor(data$class,levels = c("Up","Down"))

data$pvalue=factor(data$pvalue,levels = unique(data$pvalue[order(data$pvalue)]))

data$motif=toupper(data$motif)

data$motif=factor(data$motif,levels = rev(data$motif))
p<-ggplot(data,aes(x=class,y=motif))+
  geom_point(aes(size=enrichment,color=-logPvalue))+
  theme_bw()+
  xlab(NULL)+
  ylab(NULL)+
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

p+theme(legend.position = "top",
    legend.text=element_text(size=6, family="sans"),
    legend.title=element_text(size=6, family= "sans"),
    legend.background = element_rect(fill="white", color="white"),
    panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
    legend.key = element_rect(fill="white"))
p
ggsave("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\DiffBindv4\\reshape\\Motif_bubble_plot.pdf",height=5,width=4)








rm(list=ls())
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(tidyverse)
library(stringr)
library(tidyr)
#install.packages("tidyverse")
setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\peak\\peak_0.05/")
#color
mycol= brewer.pal(n = 7, name = "YlGnBu")
##data1
data1=read.csv("B1P1_L2_801D73_peaks.narrowPeak",sep="\t",header=F)

##data2
data2=read.csv("B1P2_L1_802D73_peaks.narrowPeak",sep="\t",header=F)

##data3
data3=read.csv("B5P1_L4_806D73_peaks.narrowPeak",sep="\t",header=F)

##data4
data4=read.csv("B5P2_L2_808D73_peaks.narrowPeak",sep="\t",header=F)

##data1
data5=read.csv("L2P1_L1_809F10_peaks.narrowPeak",sep="\t",header=F)

##data2
data6=read.csv("L2P2_L1_810F10_peaks.narrowPeak.anno",sep="\t",header=F)

##data3
data7=read.csv("L5P1_L1_811F10_peaks.narrowPeak",sep="\t",header=F)

##data4
data8=read.csv("L5P2_L1_812F10_peaks.narrowPeak",sep="\t",header=F)


##data merge
data=data.frame(number=c(dim(data1)[1],dim(data2)[1],dim(data3)[1],dim(data4)[1],dim(data5)[1],dim(data6)[1],dim(data7)[1],dim(data8)[1]),
                file=c("B1P1","B1P2","B5P1","B5P2","L2P1","L2P2","L5P1","L5P2"),
                class=c("B","B","B","B","L","L","L","L"))

data$class=factor(data$class,levels = c("B1P1","B1P2","B5P1","B5P2",
                                        "L2P1","L2P2","L5P1","L5P2"))
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


rm(list=ls())
library(openxlsx)
library(readxl)
library(ggplot2)
library(pheatmap)
setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\RNA-seq\\deseq2")
data = read.csv("count__rpkm.csv",row.names = 1)
data=data[order(data$B,decreasing = T),]
data = data[c(1:10),c(4,5,6,1,2,3)]
L=data.frame(row.names = rownames(data),sample="L",number=data[,"L"])
L1=data.frame(row.names = rownames(data),sample="L",number=data[,"L"])
L2=data.frame(row.names = rownames(data),sample="L",number=data[,"L"])
B=data.frame(row.names = rownames(data),sample="B",number=data[,"L"])
B1=data.frame(row.names = rownames(data),sample="B",number=data[,"L"])
B2=data.frame(row.names = rownames(data),sample="B",number=data[,"L"])

res=rbind(L,L1)
res=rbind(res,L2)
res=rbind(res,B)
res=rbind(res,B1)
res=rbind(res,B2)
rownames(res) = factor(rownames(res),levels = rownames(res))
ggplot(data=res)+
  geom_boxplot(aes(x = sample, y = number))

diff1=read.csv("degdown.csv",row.names = 1)
diff2=read.csv("degup.csv",row.names = 1)
names=c(rownames(diff1),rownames(diff2))
names
data=data[names,]
colnames(data)=c("B1","B2","B3","L1","L2","L3")
data=data[,c(4,5,6,1,2,3)]
pdf("RNAseq_Heatmap.pdf",width = 4,height = 7)
pheatmap(data[names,],
         scale = "row", border = F,
         color=colorRampPalette(rev(c("red","white","blue")))(102),
         cutree_rows = 2,
         cutree_cols = 2,
         ColSideColors = rep(c("pink", "purple"), each = 4),
         show_rownames = F)
dev.off()

library(ggpubr)
data=read.csv("PRL.csv")
data=data[which(data$GENE=="PRL"),]
data$CLASS=factor(data$CLASS,levels = c("Lay","Broody"))
data$FPKM=c(11879.90394,36686.0746,32988.86584,90702.22405,45488.79565	,85167.57483	)
data
pdf("PRL_boxplot_color.pdf",height = 2.5,width = 2.5)
ggplot(data=data)+
  geom_boxplot(aes(x = GENE, y = FPKM,fill= CLASS))+
  theme_classic()+
  scale_y_continuous(position = "right")
data
dev.off()
ggsave("PRL_boxplot_color.pdf",height = 2,width = 2)



library(openxlsx)
library(ggplot2)
data=read.csv("LayBreed_res.Count_matrix.xls.DESeq2.normalized.csv",row.names = 1)
head(data)
data=data[-which(rownames(data)%in%c("COX1","ND1","ND4","ND5","COX2","ATP6","COX3","CYTB")),]
data=data[c(1:10),]
gene=rownames(data)
gene

data$X = factor(data$X,levels = rev(data$X))
data$baseMean=round(data$baseMean)
#data$baseMean = factor(data$baseMean,levels = data$baseMean)
ggplot(data = data, mapping = aes(x = rownames(data), y = baseMean)) + 
  geom_bar(stat = 'identity') + 
  theme_classic()+scale_y_continuous(position = "right")+
  theme(axis.text.x = element_text(angle = 70,vjust = 0.85,hjust = 0.85))+
  scale_color_gradient(low="black",high = "red")
ggsave("top10_rpkm.pdf",height = 3,width = 4)

##PCA
library(ggrepel)
library(ggplot2)
data=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\DiffBindv4\\dba_binding.txt")
head(data)
p <- ggplot(data,aes(x=PC1,y=PC2,color=class))+ geom_point(size=5)+
  theme_bw() +xlim(-100,100)+ylim(-100,100)+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"))+
  xlab("Principal Component #1 [59%]") +
  ylab("Principal Component #2 [21%]") +
  geom_text_repel(aes(PC1, PC2, label = rownames(data)))
p

library(ggplot2)
library(factoextra)
library(FactoMineR)
head(data)
data=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\RNA-seq\\deseq2\\LayBreed_trans.Count_matrix.xls.DESeq2.normalized.xls",header=T)
dim(data)
df=t(data)

dim(df)
iris.pca<- PCA(df,graph = FALSE)
iris.pca$var
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             pointsize =4,
             mean.point=F,
             pointshape = 19,
             col.ind = c("Broody","Broody","Broody","Lay","Lay","Lay"), # color by groups
             palette = c("#00AFBB", "#E7B800"),# 
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             title="")+
  theme_bw() +
  theme(text=element_text(size=10,face="plain",color="black"),
        axis.title=element_text(size=10,face="plain",color="black"),
        axis.text = element_text(size=10,face="plain",color="black"),
        legend.title = element_text(size=10,face="plain",color="black"),
        legend.text = element_text(size=10,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1)
  )

col=c(rep("red",3),rep("blue",3))
pca.plot = function(data,col){
  
  df.pca <- PCA(t(data), graph = FALSE)
  fviz_pca_ind(df.pca,
               mean.point=F,
               geom.ind = "point",
               col.ind = c("B","B","B","L","L","L"),
               legend.title = "Groups"
  )
}

pdf("Deseq2-PCA-2.pdf",height = 4,width = 4)
dim(df)
pca.plot(data,col)
dev.off()
data=read.table("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\DiffBindv4\\Lay")





rm(list=ls())
library(tidyverse)
library(dplyr)

setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\DiffBindv4\\reshape\\")
data=read.csv("Lay_vs_Breed_deseq2_sig_reshape.anno",header=F,sep="\t")
score1=separate(data = data, col = V8, into = c("score1", "score2"), sep = "\\(")["score1"]
head(score1)
data$V8=score1
head(data)
table(data$V8)/dim(data)[1]



library(pheatmap)
setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\DiffBindv4\\")

pdf("DiffBind_deseq_Heatmap.pdf",height = 4,width = 3)
Heatmap_Data=read.table("Lay_vs_Breed_deseq2_pheatmap_inR.txt",row.names = 1,header=T)


annotation_col = data.frame(
  Conditions  = factor(c(rep("Lay", 4),rep("Breeding", 3)) ))
head(annotation_col)
rownames(annotation_col)=colnames(Heatmap_Data)
class(Heatmap_Data)
colnames(Heatmap_Data)=factor(colnames(Heatmap_Data))

pheatmap(Heatmap_Data,
scale = "row",
show_rownames = F,
cluster_cols = F,
main = "DAR Heatmap",
cutree_rows = 2,
gaps_col  = 4 ,
color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
ColSideColors = c("pink","pink","pink","pink", "purple", "purple", "purple"))

dev.off()




library(openxlsx)
library(ggplot2)
setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\RNA-seq\\deseq2")
# data=read.csv("count__rpkm.csv",row.names = 1)
# head(data)
# data=data[which(rownames(data)%in%
#                   c("PRL","CHGA","LOC106046928","EEF1A1","POMC","FTH1","DKK3","RPL23","SCG2","HSP90AA1",
#                     "LOC106044772","RPL23")),]
# head(data)
# dim(data)
# gene=rownames(data)
# gene
# head(data)
# data$name=rownames(data)
# data$name=factor(data$name,levels = rev(data$name))
# rownames(data) = factor(rownames(data),levels = rownames(data))
# #data$baseMean=round(data$baseMean)
# #data$baseMean = factor(data$baseMean,levels = data$baseMean)
# write.csv(data,"rpkm_top10.csv")
position_dodge(width = NULL, preserve = c("total", "single"))

position_dodge2(
  width = NULL,
  preserve = c("total", "single"),
  padding = 0.1,
  reverse = FALSE
)
data=read.csv("rpkm_top10.csv")
data
data$gene=factor(data$gene,levels = c("PRL","CHGA","LOC106046928","EEF1A1","POMC","DKK3","FTH1","RPL23","SCG2","LOC106044772","HSP90AA1"))
ggplot(data = data, mapping = aes(x = gene,y = rpkm,fill=class)) + 
  geom_bar(position=position_dodge(0),stat="identity")+ 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 70,vjust = 0.85,hjust = 0.85))+
  ylab(expression(FPKM~ value~ 10^3))
#dev.off()
ggsave("top10_rpkm.pdf",height = 3,width = 4)



library(ggplot2)
library(dplyr)
atac_data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\ATAC-seq\\DiffBindv4\\diffbind_deseq2_quantant.anno",sep="\t",header=F)
head(atac_data)
gene=as.character(unique(atac_data$V20))
gene
atac=data.frame(target="",atac_log2FoldChange="",atac_p="",atac_padj="")
atac[,"target"]=as.character(atac[,"target"])
atac[,"atac_log2FoldChange"]=as.character(atac[,"atac_log2FoldChange"])
atac[,"atac_p"]=as.character(atac[,"atac_p"])
atac[,"atac_padj"]=as.character(atac[,"atac_padj"])
atac_data[atac_data$V20==gene[1],c(11,12,13)]
for (i in c(1:length(gene))){
  tmp=atac_data[atac_data$V20==gene[i],c(11,12,13)]
  sum_data=sum(tmp[,1])
  atac[i,"target"]=gene[i]
  atac[i,"atac_log2FoldChange"]=sum_data
  atac[i,"atac_p"]=min(tmp[,2])
  atac[i,"atac_padj"]=min(tmp[,3])
}
rownames(atac)=atac[,1]
atac=atac[,-1]
RNA_data=read.csv("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\RNA-seq\\deseq2\\LayBreed_res.Count_matrix.xls.DESeq2.normalized.csv",row.names=1)
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
p5


cor = cor(res[,c("RNA_log2FoldChange","atac_log2FoldChange")],use="complete.obs")

cor=round(cor[2,1],digits = 2)
cor
#准备作为图形的标题;
#lab = paste("correlation=",cor,sep="")
#lab
#[1] "correlation=0.35"
#在图上添加文字标签；
#p5+geom_text(x = -2, y = 5, label = lab)
ggsave(width = 4,height = 4,"C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\quadrant.pdf")


##tobias
rm(list=ls())
library(ggplot2)
library(ggrepel)
setwd("C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\")
data=read.table("bindetect_results.txt",header=T)
colnames(data)
data$label=""
data[which(data$Lay_Broody_change >= 0.1 & -log10(data$Lay_Broody_pvalue) >= 60),"label"]="high score in Lay condition"
data[which(data$Lay_Broody_change <= -0.1 & -log10(data$Lay_Broody_pvalue) >= 60),"label"]="high score in Broody condition"

data$rep <- ifelse(-log10(data$Lay_Broody_pvalue) >= 60 & abs(data$Lay_Broody_change) >= 0.1,
                        as.character(data$name), "")


p <- ggplot(
  # 数据、映射、颜色
  data, aes(x = Lay_Broody_change, y = -log10(Lay_Broody_pvalue),colour = label)) +
  geom_point(size=1) +
  scale_color_manual(values=c( "#d2dae2",
                               "#546de5",
                               "#ff4757"))+
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_classic()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
pdf("aaaaaaa.pdf")
p + geom_label_repel(data = data, aes(x = data$Lay_Broody_change, 
                                         y = -log10(data$Lay_Broody_pvalue), 
                                         label = rep),
                     size = 3, box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE)

# 注释
p+geom_text_repel(
  data = data,
  aes(label = rep),
  size = 5,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"))
dev.off()
