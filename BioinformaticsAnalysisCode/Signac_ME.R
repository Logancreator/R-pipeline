library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Dmelanogaster.Ensemble.dm6)
library(future)
library(ggplot2)
library(clustree)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2) # for 100 Gb RAM
setwd('/public/home/changjianye/project/Signac/another/ME')
# #########################################################
#prepare
fai=read.table("/public/home/changjianye/project/Signac/ref/Sichuan_goose_rename_addN/fasta/genome.fa.fai")
start=fai[,1]
length=fai[,2]
genome <- Seqinfo(seqnames=start,seqlengths=length)

raw <- rtracklayer::import("/public/home/changjianye/project/Signac/ref/rename_addN.gff")
gtf = raw[raw$type%in%c("gene","mRNA","exon","CDS"),c("source","type","gene","gene_biotype")]
gtf$transcript_source = 'ensemble'
gtf$gene_biotype = 'protein_coding'
gtf$gene_id = gtf$gene
gtf[is.na(gtf$gene_id)]$gene_id=paste('Mt_exon_',c(1:24),sep="")
gtf$gene_name = gtf$gene_id
fa <- Rsamtools::FaFile("/public/home/changjianye/project/Signac/ref/Sichuan_goose_rename_addN/fasta/genome.fa")
# #########################################################
# read in peak sets
signac_1E <- read.table(
  file = "/public/home/changjianye/project/Signac/another/1E-ME/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
signac_2E <- read.table(
  file = "/public/home/changjianye/project/Signac/another/2E-ME/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
signac_3E <- read.table(
  file = "/public/home/changjianye/project/Signac/another/3E-ME/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)
signac_4E <- read.table(
  file = "/public/home/changjianye/project/Signac/another/4E-ME/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.signac_1E <- makeGRangesFromDataFrame(signac_1E)
gr.signac_2E <- makeGRangesFromDataFrame(signac_2E)
gr.signac_3E <- makeGRangesFromDataFrame(signac_3E)
gr.signac_4E <- makeGRangesFromDataFrame(signac_4E)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.signac_1E, gr.signac_2E, gr.signac_3E, gr.signac_4E))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#############################################################
##create object
# load metadata
signac_1E_cell <- read.table(
  file = "/public/home/changjianye/project/Signac/another/1E-ME/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

signac_2E_cell <- read.table(
  file = "/public/home/changjianye/project/Signac/another/2E-ME/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

signac_3E_cell <- read.table(
  file = "/public/home/changjianye/project/Signac/another/3E-ME/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

signac_4E_cell <- read.table(
  file = "/public/home/changjianye/project/Signac/another/4E-ME/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
#############################################################
# perform an initial filtering of low count cells
signac_1E_cell <- signac_1E_cell[signac_1E_cell$passed_filters > 2000, ]
signac_2E_cell <- signac_2E_cell[signac_2E_cell$passed_filters > 2000, ]
signac_3E_cell <- signac_3E_cell[signac_3E_cell$passed_filters > 2000, ]
signac_4E_cell <- signac_4E_cell[signac_4E_cell$passed_filters > 2000, ]
#############################################################
# create fragment objects
signac_1E_frags <- CreateFragmentObject(
  path = "/public/home/changjianye/project/Signac/another/1E-ME/outs/fragments.tsv.gz",
  cells = rownames(signac_1E_cell)
)
signac_2E_frags <- CreateFragmentObject(
  path = "/public/home/changjianye/project/Signac/another/2E-ME/outs/fragments.tsv.gz",
  cells = rownames(signac_2E_cell)
)
signac_3E_frags <- CreateFragmentObject(
  path = "/public/home/changjianye/project/Signac/another/3E-ME/outs/fragments.tsv.gz",
  cells = rownames(signac_3E_cell)
)
signac_4E_frags <- CreateFragmentObject(
  path = "/public/home/changjianye/project/Signac/another/4E-ME/outs/fragments.tsv.gz",
  cells = rownames(signac_4E_cell)
)
#############################################################
# Quantify the peaks in each data set
signac_1E_counts <- FeatureMatrix(
  fragments = signac_1E_frags ,
  features = combined.peaks,
  cells = rownames(signac_1E_cell)
)

signac_2E_counts <- FeatureMatrix(
  fragments = signac_2E_frags,
  features = combined.peaks,
  cells = rownames(signac_2E_cell)
)

signac_3E_counts <- FeatureMatrix(
  fragments = signac_3E_frags,
  features = combined.peaks,
  cells = rownames(signac_3E_cell)
)

signac_4E_counts <- FeatureMatrix(
  fragments = signac_4E_frags,
  features = combined.peaks,
  cells = rownames(signac_4E_cell)
)
# #############################################################
# create object
signac_1E_assay <- CreateChromatinAssay(counts=signac_1E_counts, fragments = signac_1E_frags , min.cells = 10,min.features = 2000 , sep = c(":","-") , genome = genome , annotation = gtf)
signac_1E <- CreateSeuratObject(counts=signac_1E_assay, assay = "ATAC", meta.data=signac_1E_cell)
DefaultAssay(signac_1E) <- "ATAC"
Annotation(signac_1E) <- gtf

signac_2E_assay <- CreateChromatinAssay(counts=signac_2E_counts, fragments = signac_2E_frags , min.cells = 10,min.features = 2000 , sep = c(":","-") , genome = genome , annotation = gtf)
signac_2E <- CreateSeuratObject(counts=signac_2E_assay, assay = "ATAC", meta.data=signac_2E_cell)
DefaultAssay(signac_2E) <- "ATAC"
Annotation(signac_2E) <- gtf

signac_3E_assay <- CreateChromatinAssay(counts=signac_3E_counts, fragments = signac_3E_frags , min.cells = 10,min.features = 2000 , sep = c(":","-") , genome = genome , annotation = gtf)
signac_3E <- CreateSeuratObject(counts=signac_3E_assay, assay = "ATAC", meta.data=signac_3E_cell)
DefaultAssay(signac_3E) <- "ATAC"
Annotation(signac_3E) <- gtf

signac_4E_assay <- CreateChromatinAssay(counts=signac_4E_counts, fragments = signac_4E_frags , min.cells = 10,min.features = 2000 , sep = c(":","-") , genome = genome , annotation = gtf)
signac_4E <- CreateSeuratObject(counts=signac_4E_assay, assay = "ATAC", meta.data=signac_4E_cell)
DefaultAssay(signac_4E) <- "ATAC"
Annotation(signac_4E) <- gtf
############################################################
# merge object
# add information to identify dataset of origin
signac_1E$dataset <- '1E-ME'
signac_2E$dataset <- '2E-ME'
signac_3E$dataset <- '3E-ME'
signac_4E$dataset <- '4E-ME'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = signac_1E,
  y = list(signac_2E, signac_3E, signac_4E),
  add.cell.ids = c("1E", "2E", "3E", "4E")
)
#combined[["ATAC"]]

pdf("/public/home/changjianye/project/Signac/another/ME/combined.pdf")

DefaultAssay(combined) <- 'ATAC'
## UMAP
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 30)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 1:10, reduction = 'lsi')
combined  <- FindNeighbors(object = combined , reduction = 'lsi', dims = 1:10)
combined <- FindClusters(object = combined ,resolution = 0.1)
DimPlot(object = combined, label = TRUE) + NoLegend()

###############################################################
#FindClusters
combined <- FindClusters(object = combined ,resolution = c(seq(.4,1.6,.2)))
clustree(combined@meta.data,prefix='ATAC_snn_res.')


#peakPlot
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
CoveragePlot(
  object = combined,
  group.by = 'dataset',
  region = "scaffold0001-500-5500"
)
dev.off()

##fliter
# compute nucleosome signal score per cell
combined <- NucleosomeSignal(object = combined)

# compute TSS enrichment score per cell 
# fast = TRUE
combined <- TSSEnrichment(object = combined, fast = TRUE)

# add blacklist ratio and fraction of reads in peaks
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

# combined <- subset(
#   x = combined,
#   subset = peak_region_fragments > 1000 &
#     peak_region_fragments < 30000 &
#     pct_reads_in_peaks > 10 &
#     blacklist_ratio < 0.05 &
#     nucleosome_signal < 3 &
#     TSS.enrichment > 1
# )
combined

saveRDS(object = combined,file = "/public/home/changjianye/project/Signac/another/ME/combined_subset_cluster.rds")

# ##########################################################################
# gene_1E=unique(read.csv('/public/home/changjianye/project/Signac/another/1E-ME/outs/peak_annotation.tsv',sep="\t",row.names=NULL,header=T)$end)
# gene_2E=unique(read.csv('/public/home/changjianye/project/Signac/another/2E-ME/outs/peak_annotation.tsv',sep="\t",row.names=NULL,header=T)$end)
# gene_3E=unique(read.csv('/public/home/changjianye/project/Signac/another/3E-ME/outs/peak_annotation.tsv',sep="\t",row.names=NULL,header=T)$end)
# gene_4E=unique(read.csv('/public/home/changjianye/project/Signac/another/4E-ME/outs/peak_annotation.tsv',sep="\t",row.names=NULL,header=T)$end)

# gene=unique(Reduce(intersect,list(gene_1E,gene_2E,gene_3E,gene_4E)))
# gene=intersect(gene,gtf$gene)
# #RNA assay
# DefaultAssay(combined) <- 'ATAC'
# gene.activities <- GeneActivity(combined,features=gene)
# # add the gene activity matrix to the Seurat object as a new assay and normalize it
# combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
# combined <- NormalizeData(
#   object = combined,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(combined$nCount_RNA)
# )

# pdf("genePlot.pdf")
# DefaultAssay(combined) <- 'RNA'

# FeaturePlot(
#   object = combined,
#   features = c("TSHB"),
#   pt.size = 0.01,
#   max.cutoff = 'q95',
#   min.cutoff = "q10",
#   ncol = 1)

# dev.off()
# saveRDS(object = combined,file = "combined_subset_cluster_RNA.rds")
# ################################################################
# ##add Motif
# # Get a list of motif position frequency matrices from the JASPAR database
# # 使用getMatrixSet函数从JASPAR数据库中提取Motif的PFM矩阵信息
# pfm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

# # Scan the DNA sequence of each peak for the presence of each motif
# # 使用CreateMotifMatrix函数构建Motif矩阵对象
# motif.matrix <- CreateMotifMatrix(
#   features = granges(combined),
#   pwm = pfm,
#   genome = fa,
#   use.counts = FALSE)

# # # Create a new Mofif object to store the results
# # # 使用CreateMotifObject函数构建Motif对象
# motif <- CreateMotifObject(
#   data = motif.matrix,
#   pwm = pfm)

# combined <- SetAssayData(
#   object = combined,
#   slot = 'motifs',
#   new.data = motif
# )

# # Add the Motif object to the assay
# # 使用AddMotifObject函数将Motif类添加到Seurat对象中
# combined <- AddMotifs(combined, genome = fa, pfm = pfm,assay='ATAC')

# pdf("motif_CTCF_MA0139.pdf")
# combined <- Footprint(combined, genome = fa, in.peaks = TRUE, motif.name = "MA0139.1")
# p2 <- PlotFootprint(combined, features = "MA0139.1")
# p2
# dev.off()

##############################################################
#celltype identity
# '''

# #Astrocyte
# c("AGT","AQP4")

# #Glutamatergic neurons
# c("SLC17A6")

# #GABA neurons
# c("GAD1","TH","SLC18A2","SLC17A6","DPP10")

# #Neuron
# c("CACNA1E","CACNA2D1")

# #Oligodendrocyte
# c("ST18","ERMN","ABCA2","SOX10")

# #Oligodendrocyte Precursor Cell
# "MEGF11"

# #Microglial
# c("C1QA","RUNX1","CX3CR1")

# markerGenes=c("AGT","AQP4","SLC17A6","GAD1","TH","SLC18A2","SLC17A6","DPP10","CACNA1E","CACNA2D1","ST18","ERMN","ABCA2","SOX10","MEGF11","C1QA","RUNX1","CX3CR1")
# '''

# new.cluster.ids=c('GABA neurons','Oligodendrocyte','Astrocyte','Oligodendrocyte Precursor Cell','Unknown','Microglial','Glutamatergic neurons','Unknown','Unknown','Unknown')
# combined@meta.data$celltype <- combined@meta.data$seurat_clusters
# levels(combined@meta.data$celltype) <- new.cluster.ids
# #Dotplot
# DotPlot(combined, features = unique(markerGenes),group.by = "celltype")+RotatedAxis()+
#   scale_x_discrete("")+scale_y_discrete("")
# ggsave("celltype_markerGenes_dot.pdf",width = 9.5,height = 6)

###############################################################
##something to run
# markerGenes=c("AGT","AQP4","SLC17A6","GAD1","TH","SLC18A2","SLC17A6","DPP10","CACNA1E","CACNA2D1","ST18","ERMN","ABCA2","SOX10","MEGF11","C1QA","RUNX1","CX3CR1")
# combined <- readRDS("combined_subset_cluster_RNA.rds")
# new.cluster.ids=c('GABA neurons','Oligodendrocyte','Astrocyte','Oligodendrocyte Precursor Cell','Unknown','Microglial','Glutamatergic neurons','Unknown','Unknown','Unknown')
# combined@meta.data$celltype <- combined@meta.data$seurat_clusters
# levels(combined@meta.data$celltype) <- new.cluster.ids
# #Dotplot
# DotPlot(combined, features = unique(markerGenes),group.by = "celltype")+RotatedAxis()+
#   scale_x_discrete("")+scale_y_discrete("")
# ggsave("celltype_markerGenes_dot.pdf",width = 9.5,height = 6)


# pdf("TSHB.pdf")
# DefaultAssay(combined) <- 'RNA'
# FeaturePlot(
#   object = combined,
#   features = c("TSHB"),
#   pt.size = 0.01,
#   max.cutoff = 'q95',
#   min.cutoff = "q10",
#   ncol = 1)
# dev.off()
# saveRDS(object = combined,file = "combined_subset_cluster_RNA_rename.rds")
###############################################################
##call peak
# DefaultAssay(combined) <- 'ATAC'
# peaks <- CallPeaks(
#   object = combined,
#   group.by = "celltype",
#   macs2.path = "macs2"
# )
# saveRDS(object = peaks,file = "peaks.rds")
# #############################################################
# # set path
# setwd("/public/home/changjianye/project/Signac/")
# '''
# chr=unique(read.csv("/public/home/changjianye/project/Signac/2E-ME-P/2E-ME/outs/peaks.bed",comment.char="#",header=F,sep="\t")$V1)[1]
# write.table(chr,"chrUse.txt",quote=F,row.names=F,col.names=F)
# fai=read.table("/public/home/changjianye/project/Signac/ref/Sichuan_goose/fasta/genome.fa.fai")
# start=intersect(fai[,1],chr)
# length=fai[which(fai$V1%in%start),2]
# genome <- Seqinfo(seqnames=start,seqlengths=length)
# '''
# fai=read.table("/public/home/changjianye/project/Signac/ref/Sichuan_goose_rename_addN/fasta/genome.fa.fai")
# start=fai[,1]
# length=fai[,2]
# genome <- Seqinfo(seqnames=start,seqlengths=length)
# '''
# gtf=read.csv("/public/home/changjianye/project/Signac/ref/Sichuan_goose/genes/genes.gtf.gz",sep="\t",comment.char="#",header=F)
# gtf=gtf[which(gtf$V1%in%chr),]
# write.table(gtf,"/public/home/changjianye/project/Signac/ref/Sichuan_goose/genes/new_genes.gtf",sep="\t",quote=F,row.names=F,col.names=F)
# '''
# raw <- rtracklayer::import("/public/home/changjianye/project/Signac/ref/rename_addN.gff")
# gtf = raw[raw$type%in%c("gene","pseudogene","ncRNA","snRNA","C_gene_segment","D_loop","snoRNA","rRNA","V_gene_segment","lnc_RNA"),c("source","type","gene","gene_biotype")]
# gtf$gene_id = gtf$gene
# gtf$gene_name = gtf$gene
# fa <- Rsamtools::FaFile("/public/home/changjianye/project/Signac/ref/Sichuan_goose_rename_addN/fasta/genome.fa")

# counts <- Read10X_h5(filename = "/public/home/changjianye/project/Signac/another/2E-ME/outs/filtered_peak_bc_matrix.h5")
# single <- "/public/home/changjianye/project/Signac/another/2E-ME/outs/singlecell.csv"
# fragments <- '/public/home/changjianye/project/Signac/another/2E-ME/outs/fragments.tsv.gz'
# metadata <- read.csv(
#   file = single,
#   header = TRUE,
#   row.names = 1
# )




# chrom_assay <- CreateChromatinAssay(
#   counts = counts,
#   sep = c(":","-"),
#   genome = genome,
#   annotation = gtf,
#   fragments = fragments,
#   min.cells = 10,
#   min.features = 2000
# )

# pbmc <- CreateSeuratObject(
#   counts = chrom_assay,
#   assay = "peaks",
#   meta.data = metadata
# )

# DefaultAssay(pbmc) <- "peaks"
# Annotation(pbmc) <- gtf
# # compute nucleosome signal score per cell
# pbmc <- NucleosomeSignal(object = pbmc)

# # compute TSS enrichment score per cell
# pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# # add blacklist ratio and fraction of reads in peaks
# pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
# pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# # pdf("QC.pdf")
# # pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
# # TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

# # pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
# # FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

# # VlnPlot(
# #   object = pbmc,
# #   features = c('pct_reads_in_peaks', 'peak_region_fragments',
# #                'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
# #   pt.size = 0.1,
# #   ncol = 5
# # )

# # dev.off()


# pbmc <- subset(
#   x = pbmc,
#   subset = peak_region_fragments > 3000 &
#     peak_region_fragments < 20000 &
#     pct_reads_in_peaks > 15 &
#     blacklist_ratio < 0.05 &
#     nucleosome_signal < 4 &
#     TSS.enrichment > 2
# )
# pbmc

# pdf("DepthCor.pdf")
# pbmc <- RunTFIDF(pbmc)
# pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
# pbmc <- RunSVD(pbmc)
# DepthCor(pbmc)
# dev.off()

# pdf("DimPlot.pdf")
# pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
# pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
# pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
# DimPlot(object = pbmc, label = TRUE) + NoLegend()
# dev.off()

# gene.activities <- GeneActivity(pbmc)
# # add the gene activity matrix to the Seurat object as a new assay and normalize it
# pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
# pbmc <- NormalizeData(
#   object = pbmc,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(pbmc$nCount_RNA)
# )

# pdf("genePlot.pdf")
# DefaultAssay(pbmc) <- 'RNA'

# FeaturePlot(
#   object = pbmc,
#   features = c('DIO2',"CGA","CTCF"),
#   pt.size = 0.01,
#   max.cutoff = 'q95',
#   min.cutoff = "q10",
#   ncol = 2
# )

# dev.off()


# ##rename
# # pbmc <- RenameIdents(
# #   object = pbmc,
# #   '0' = 'Epithelial',
# #   '1' = 'Epithelial',
# #   '2' = 'Basophil',
# #   '3' = 'Myeloid_1',
# #   '4' = 'Myeloid_2',
# #   '5' = 'Tcell'
# # )


# # change back to working with peaks instead of gene activities
# DefaultAssay(pbmc) <- 'peaks'

# da_peaks <- FindMarkers(
#   object = pbmc,
#   ident.1 = "1",
#   ident.2 = "2",
#   min.pct = 0.05,
#   test.use = 'LR',
#   latent.vars = 'peak_region_fragments'
# )

# head(da_peaks)

# ##add Motif
# # Get a list of motif position frequency matrices from the JASPAR database
# # 使用getMatrixSet函数从JASPAR数据库中提取Motif的PFM矩阵信息
# pfm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

# # Scan the DNA sequence of each peak for the presence of each motif
# # 使用CreateMotifMatrix函数构建Motif矩阵对象
# motif.matrix <- CreateMotifMatrix(
#   features = granges(pbmc),
#   pwm = pfm,
#   genome = fa,
#   use.counts = FALSE)

# # # Create a new Mofif object to store the results
# # # 使用CreateMotifObject函数构建Motif对象
# motif <- CreateMotifObject(
#   data = motif.matrix,
#   pwm = pfm)


# # motif <- AddMotifs(
# #     object = granges(x = pbmc),
# #     genome = fa,
# #     pfm = pfm,
# #     verbose = TRUE
# #   )
# pbmc <- SetAssayData(
#   object = pbmc,
#   slot = 'motifs',
#   new.data = motif
# )

# # Add the Motif object to the assay
# # 使用AddMotifObject函数将Motif类添加到Seurat对象中
# pbmc <- AddMotifs(pbmc, genome = fa, pfm = pfm,assay='peaks')

# pdf("MA0139.pdf")
# pbmc <- Footprint(pbmc, genome = fa, in.peaks = TRUE, motif.name = "MA0139.1")
# p2 <- PlotFootprint(pbmc, features = "MA0139.1")
# p2
# dev.off()


# '''
# fai=read.table("/public/home/changjianye/project/Signac/ref/Sichuan_goose_rename/fasta/genome.fa.fai")
# start=fai[,1]
# length=fai[,2]
# genome <- Seqinfo(seqnames=start,seqlengths=length)

# gtf <- rtracklayer::import("/public/home/changjianye/project/Signac/ref/Sichuan_goose_rename/genes/genes.gtf.gz")
# gtf = gtf[which(gtf$type=="transcript"),c("type","gene_id")]
# gtf = gtf[which(gtf$type=="transcript"),c("type","gene_id")]
# gtf$gene_biotype = "protein_coding"
# gtf$gene_name = gtf$gene_id

# fa <- Rsamtools::FaFile("/public/home/changjianye/project/Signac/ref/Sichuan_goose_rename/fasta/genome.fa")

# dic=read.table("/public/home/changjianye/project/Signac/ref/Sichuan_goose/rename/chr_rename.txt")
# row=rownames(counts)
# chr=strsplit(row[1],split = ":") [[1]][1]
# num=strsplit(row[1],split = ":") [[1]][2]
# chrUse=NULL
# numUse=NULL
# for (i in c(1:length(row))){
#   chr=strsplit(row[i],split = ":") [[1]][1]
#   chrUse=append(chrUse,chr)
#   num=strsplit(row[i],split = ":") [[1]][2]
#   numUse=append(numUse,num)
# }



# id=match(chrUse,dic$V1)
# name=NULL
# for (i in c(1:length(id))){

#   if(id[i]< 10 & id[i] >= 0){
#   name <- append(name,paste("000", as.character(id[i]),sep=""))
#   }else if(id[i] < 100 & id[i] >= 10){
#   name <- append(name,paste("00", as.character(id[i]),sep=""))
#   }else if(id[i] < 1000 & id[i] >= 100){
#   name <- append(name,paste("0", as.character(id[i]),sep=""))
#   }else{
#   name <- append(name,as.character(id[i]))
#   }

# }
# chrName=paste("scaffold",name,sep="")
# rowNew=paste(chrName,":",numUse,sep="")
# rownames(counts)=rowNew
# '''
