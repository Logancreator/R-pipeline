library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(patchwork)
library(AnnotationHub)
library(ELISAtools)

#creat reference 
# hub <- AnnotationHub()
# query(hub, "Sus scrofa")
# pig <- hub[["AH98171"]]

annotations <- readRDS("/home/chenruipu/intestine/scATAC/pig_ref.rds")
## create seurat object
counts <- Read10X_h5(filename = "/home/chenruipu/snp/scATAC/cellranger/scATAC-3/outs/filtered_peak_bc_matrix.h5")


# keep peaks having at least 100 fragments
# rs <- rowSums(counts)
# peaks.keep <- rownames(counts[rs > 100, ])
# # counts <- counts[peaks.keep, ]


# chromatinassay <- CreateChromatinAssay(counts = counts, genome = "susScr11")

frag.path <- "/home/chenruipu/snp/scATAC/cellranger/scATAC-3/outs/fragments.tsv.gz"

chrom_assay <- CreateChromatinAssay(
  sep = c(":", "-"),
  counts = counts,
  genome = "susScr11",
  fragments = frag.path,
  min.cells = 0,
  min.features = 200,
  validate.fragments = F
)
metadata <- read.csv(
  file = "/home/chenruipu/snp/scATAC/cellranger/scATAC-3/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  project = "ATAC",
  min.cells = 1,
  meta.data = metadata,
)
# 使用NucleosomeSignal函数计算Nucleosome banding pattern
atac <- NucleosomeSignal(object = atac)
# 计算peaks中reads的比例
atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
# 计算比对到“blacklist”区域的reads比率
# atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments

# 质控指标的可视化
VlnPlot(atac, c("pct_reads_in_peaks", "peak_region_fragments", "nucleosome_signal"), pt.size = 0.1)
ggsave('QC_vol.pdf',width = 15,height = 5)
# VlnPlot(atac, c("nucleosome_signal"), pt.size = 0.1) & scale_y_log10()
# ggsave("nucleosome_signal.pdf", width = 5, height = 5)

# 根据核小体条带信号强度进行分组
atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 10, "NS > 10", "NS < 10")
# 使用FragmentHistogram函数可视化核小体条带分布模式
FragmentHistogram(object = atac, group.by = "nucleosome_group")


#annotations <- GetGRangesFromEnsDb(ensdb = pig)
# change to UCSC style since the data was mapped to hg19
#seqlevelsStyle(annotations) <- "UCSC"
#genome(annotations) <- "susScr11"
# add the gene information to the object

Annotation(atac) <- annotations

# to save time use the first 2000 TSSs
# 使用TSSEnrichment函数计算TSS富集得分
# compute nucleosome signal score per cell
atac <- NucleosomeSignal(object = atac)
# compute TSS enrichment score per cell
atac <- TSSEnrichment(object = atac, fast = FALSE)
# add blacklist ratio and fraction of reads in peaks
atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
#atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments
atac$high.tss <- ifelse(atac$TSS.enrichment > 10, "High", "Low")
TSSPlot(atac, group.by = "high.tss") + NoLegend()
ggsave('TSS_enricment.pdf',width = 10,height = 5)

atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 10, "NS > 10", "NS < 10")
# 使用FragmentHistogram函数可视化核小体条带分布模式
FragmentHistogram(object = atac, group.by = "nucleosome_group")
ggsave("FragmentHistogram.pdf", width = 10, height = 5)
# 根据不同QC指标进行数据过滤
atac <- subset(
  x = atac,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    nucleosome_signal < 10 &
    TSS.enrichment > 2
)

# 使用RunTFIDF函数进行数据归一化
atac <- RunTFIDF(atac)
# 使用FindTopFeatures函数选择可变的features
atac <- FindTopFeatures(atac, min.cutoff = "q0")
# 使用RunSVD函数进行线性降维
atac <- RunSVD(
  object = atac,
  assay = "peaks",
  reduction.key = "LSI_",
  reduction.name = "lsi"
)
DepthCor(atac)
ggsave('lsi.pdf')
# 使用RunUMAP函数进行UMAP非线性降维
atac <- RunUMAP(object = atac, reduction = "lsi", dims = 2:30)
# 对细胞执行基于图的聚类
atac <- FindNeighbors(object = atac, reduction = "lsi", dims = 2:30)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3)
# 使用DimPlot函数进行数据可视化
DimPlot(object = atac, label = TRUE) + NoLegend()
ggsave("umap.pdf", width = 10, height = 10)

#save merged object
saveRDS(atac, "scATAC-3.Rds")
#load merged object
atac<-readRDS('')
