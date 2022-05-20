##prefix

##Make the BSgenome model 
Package: BSgenome.Anser_cygnoides.UCSC
Title: Full genome sequences for Anser_cygnoides (UCSC version)
Description: Full genome sequences for Anser_cygnoides as provided by UCSC and stored in Biostrings objects.
Version: 1.4.2
Suggests: TxDb.Anser_cygnoides.UCSC.knownGene
organism: Anser_cygnoides
common_name: Goose
provider: UCSC
provider_version: goose_v1.0
release_date: Feb. 2009
release_name: Genome Reference Consortium Anser_cygnoides
organism_biocview: Homo_sapiens
BSgenomeObjname: Anser_cygnoides
    ## ---------------------------------------------------------------------
    ## Upstream sequences
    ## ---------------------------------------------------------------------
    ## Starting with BioC 3.0, the upstream1000, upstream2000, and
    ## upstream5000 sequences for hg19 are not included in the BSgenome data
    ## package anymore. However they can easily be extracted from the full
    ## genome sequences with something like:
    .
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    gn <- sort(genes(txdb))
    up1000 <- flank(gn, width=1000)
    up1000seqs <- getSeq(genome, up1000)
    .
    ## IMPORTANT: Make sure you use a TxDb package (or TxDb object),
    ## that contains a gene model based on the exact same reference genome
    ## as the BSgenome object you pass to getSeq(). Note that you can make
    ## your own custom TxDb object from various annotation resources.
    ## See the makeTxDbFromUCSC(), makeTxDbFromBiomart(),
    ## and makeTxDbFromGFF() functions in the GenomicFeatures
    ## package.
seqs_srcdir: /fh/fast/morgan_m/BioC/BSgenomeForge/srcdata/BSgenome.Hsapiens.UCSC.hg19/seqs



##Make the cellranger-atac reference from Goose's Assembled genome 
{
    organism: "Anser_cygnoides"
    genome: ["Anser_cygnoides_chr"]
    input_fasta: ["/md01/changjy/data/Goose_1.0/ref/Anser_cygnoides/fasta/genome.fa""/md01/changjy/data/Goose_1.0/ref/Anser_cygnoides/fasta/genome.fa"]
    input_gtf: ["/md01/changjy/data/Goose_1.0/ref/Anser_cygnoides/genes/genes.gtf.gz"]
    non_nuclear_contigs: ["Scaffold_310"]
}

##cellranger-atac mapping 
cellranger-atac count \
--id=1P-plus \
--reference=/md01/changjy/data/Goose/ref/Anser_cygnoides_chr \
--fastqs=/md01/changjy/data/Goose/LBFC20210423/file_1 \
--sample=1P-plus \
--localcores=20 \
--localmem=128

##sort tsv file
gunzip fragments.tsv.gz
##del some annotation
sed -i '1,161d' fragments.tsv
# del some peak
python /md01/changjy/software/del.py fragments.tsv new.fragment.tsv
##merge file
cat annotation new.fragment.tsv > new.fragments.tsv
##compress the tsv file
bgzip new.fragments.tsv 
#make the index for tsv.gz file
tabix -p bed new.fragments.tsv.gz


# ####Creat ArchRProject
# rm(list = ls()) 
# Sys.setenv(R_MAX_NUM_DLLS=999)
# library('ArchR')
# library('pheatmap')
# library('Rsamtools')
# library('scran')
# library(GenomicFeatures)
# library(BSgenome.Acygnoides.Assembly.go12)

# ##Import the previous documents
# #projHemeNew<-loadArchRProject("/md01/changjy/data/Goose/snATAC_ME_new")
# projHemeNew<-loadArchRProject("/md01/changjy/data/Goose/snATAC_PT_new_filter")


# gtf <- read.table("/md01/changjy/data/Goose/assemble_ref/ref_reduce_chrM/gff_9.txt", sep="\t")
# gtf_sub <- gtf[which(gtf[, 3] == "gene"),]
# genes <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])), 
#  gene_id=gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf_sub[, 9])))
# gtf_sub <- gtf[which(gtf[, 3] == "exon"),]
# exons <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])), 
#  gene_id=gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf_sub[, 9])))
# gtf_sub <- gtf[which(gtf[, 3] == "transcript"),]
# TSS <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 4]), strand=Rle(strand(gtf_sub[, 7])), 
#  gene_id=gsub("\\..*", "", gsub(".*transcript_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*transcript_name ", "", gtf_sub[, 9])))
# geneAnnotation <- SimpleList(genes=genes, exons=exons, TSS=TSS)
# # Scaffold_310    .   transcript  22399   24700   .   +   .   transcript_id "G14617.1.1"; gene_id "G14617.1.1";
# # Scaffold_310    .   exon    22399   22459   .   +   .   transcript_id "G14617.1.1"; gene_id "G14617.1.1";
# # Scaffold_310    .   CDS 22343   22464   .   +   0   transcript_id "G14617.1.1"; gene_id "G14617.1.1";

# ###Creating the geneAnnotation
# genes<-read.csv("/md01/changjy/data/Goose/ref/Anser_cygnoides_chr/change_rename_sorted.genes.txt",sep="\t",header=F)
# TSS<-read.csv("/md01/changjy/data/Goose/ref/Anser_cygnoides_chr/change_rename_sorted.TSS.txt",sep="\t",header=F)
# exons<-read.csv("/md01/changjy/data/Goose/ref/Anser_cygnoides_chr/change_rename_sorted.exons.txt",sep="\t",header=F)

# genes <- GRanges(seqnames=Rle(genes[,1]), 
#             ranges=IRanges(genes[, 2], end=genes[, 3]),
#             strand=Rle(strand(genes[, 5])), 
#             gene_id=genes[,6], 
#             symbol=genes[,6])

# TSS <- GRanges(seqnames=Rle(TSS[,1]), 
#             ranges=IRanges(TSS[, 2], end=TSS[, 3]),
#             strand=Rle(strand(TSS[, 5])), 
#             gene_id=TSS[,6], 
#             symbol=TSS[,6])

# exons <- GRanges(seqnames=Rle(exons[,1]), 
#             ranges=IRanges(exons[, 2], end=exons[, 3]),
#             strand=Rle(strand(exons[, 5])), 
#             gene_id=exons[,6], 
#             symbol=exons[,6])


# geneAnnotation <- SimpleList(genes=genes, exons=exons, TSS=TSS)

# ####Creating the genomeAnnotation
# library(BSgenome.Acygnoides.Assembly.go12)
# chromsize<-read.table("/md01/changjy/data/Goose/Goose_bulk/ref/ensemble/Anser_cygnoides_ensemble/fasta/genome.fa.fai")
# #chromsize<-chromsize$V2
# gr <- GRanges(chromsize$V1,IRanges(start=rep(1,dim(chromsize)[1]),end=chromsize$V2))
# genomeAnnotation<-createGenomeAnnotation(BSgenome.Acygnoides.Assembly.go12,gr)
# genomeAnnotation$chromSizes<-gr
# genomeAnnotation$chromSizes@seqinfo@seqlengths<-chromsize$V2


rm(list = ls()) 
library('ArchR')
library('pheatmap')
library('Rsamtools')
library('scran')
library(GenomicFeatures)
library(BSgenome.Acygnoides.Assembly.go12)

chr1=read.table("/md01/changjy/data/Goose/1E-ME/outs/peak_annotation_chrName.txt",header=T)
chr2=read.table("/md01/changjy/data/Goose/2-ME/outs/peak_annotation_chrName.txt",header=T)
chr3=read.table("/md01/changjy/data/Goose/3E-ME/outs/peak_annotation_chrName.txt",header=T)
chr4=read.table("/md01/changjy/data/Goose/4ME/outs/peak_annotation_chrName.txt",header=T)
aaa=intersect(chr1$chrom,chr2$chrom)
bbb=intersect(chr3$chrom,chr4$chrom)
ccc=intersect(aaa,bbb)
gtf <- read.table("/md01/changjy/data/Goose/Goose_bulk/ref/ensemble/new_Anser_cygnoides.GooseV1.0.101.gtf", sep="\t")
gtf=gtf[gtf[,1]%in%ccc,]

gtf_sub <- gtf[which(gtf[, 3] == "gene"),]
genes <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])), 
 gene_id=gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf_sub[, 9])))
gtf_sub <- gtf[which(gtf[, 3] == "exon"),]
exons <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])), 
 gene_id=gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf_sub[, 9])))
gtf_sub <- gtf[which(gtf[, 3] == "transcript"),]
TSS <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 4]), strand=Rle(strand(gtf_sub[, 7])), 
 gene_id=gsub("\\..*", "", gsub(".*transcript_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*transcript_name ", "", gtf_sub[, 9])))
geneAnnotation <- SimpleList(genes=genes, exons=exons, TSS=TSS)

library(BSgenome.Acygnoides.Assembly.go12)
chromsize<-read.table("/md01/changjy/data/Goose/Goose_bulk/ref/ensemble/Anser_cygnoides_ensemble/fasta/genome.fa.fai")
chromsize<-chromsize[chromsize[,1]%in%ccc,]
chromsize$V2=chromsize$V2+1
gr <- GRanges(as.character(chromsize$V1),IRanges(start=rep(0,dim(chromsize)[1]),end=chromsize$V2))
gr$seqinfo=gr$seqinfo[gr$seqinfo$seqnames%in%ccc,]
genomeAnnotation<-createGenomeAnnotation(BSgenome.Acygnoides.Assembly.go12,gr)
genomeAnnotation$chromSizes<-gr 
genomeAnnotation$chromSizes@seqinfo@seqlengths<-as.factor(chromsize$V2)



#Import the file
inputfiles<-c("/md01/changjy/data/Goose/1E-ME/outs/fragments.tsv.gz",

 "/md01/changjy/data/Goose/2-ME/outs/fragments.tsv.gz",

 "/md01/changjy/data/Goose/3E-ME/outs/fragments.tsv.gz",

 "/md01/changjy/data/Goose/4ME/outs/fragments.tsv.gz")
samplenames<-c("1E-ME","2-ME","3E-ME","4ME")


##Creating Arrow files
AddPrefix=FALSE
ArrowFiles <- createArrowFiles(
    inputFiles=inputfiles,
    sampleNames=samplenames,
    geneAnnotation=geneAnnotation,
    genomeAnnotation=genomeAnnotation,
    minTSS = 3, 
    minFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    force=TRUE,
    subThreading=FALSE,
    threads=1,
    cleanTmp = FALSE,
    nChunk = 1,
    )


##Creating ArchRProject
projHemeNew <- ArchRProject(ArrowFiles = ArrowFiles,
                         outputDirectory = "/md01/changjy/data/Goose/snATAC",
                         copyArrows = TRUE,
                         geneAnnotation=geneAnnotation,
                         genomeAnnotation=genomeAnnotation)












                         
##plot the FragmentSize and TSS enrichment
p1 <- plotFragmentSizes(ArchRProj = projHemeNew)
p2 <- plotTSSEnrichment(ArchRProj = projHemeNew)
plotPDF(p1,p2, name = "unprocess_QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)
##save ArchRProject
saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)###Creating the geneAnnotation

p1 <- plotFragmentSizes(ArchRProj = projHemeNew)
p2 <- plotTSSEnrichment(ArchRProj = projHemeNew)
plotPDF(p1,p2, name = "unprocess_QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)
##save ArchRProject
saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)###Creating the geneAnnotation


##Add the DoubletScores
AddPrefix=FALSE
projHemeNew <- addDoubletScores(
    input = projHemeNew,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1,
    threads=1)
saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)


##plot log10(nFrags) and TSSEnrichment
df <- getCellColData(projHemeNew, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = 3.74, lty = "dashed")
plotPDF(p, name = "unprocess_TSS-vs-Frags.pdf", ArchRProj = projHemeNew, addDOC = FALSE)

##plot the Tss enrichment rideges 
p3 <- plotGroups(
    ArchRProj = projHemeNew, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
plotPDF(p3, name = "unprocess_QC-Sample-Ridges.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)

#plot violin
p4 <- plotGroups(
    ArchRProj = projHemeNew, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p4, name = "unprocess_QC-Sample-Boxplot.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)



##filter low TSS cel
TSSEnrichmentScore=4
idxPass <- which(projHemeNew$TSSEnrichment >= TSSEnrichmentScore )
cellsPass <- projHemeNew$cellNames[idxPass]
projHemeNew1=projHemeNew[cellsPass, ]

##filter low fragment cel
nFrags=1000
idxPass <- which(projHemeNew$nFrags>=nFrags )
cellsPass <- projHemeNew$cellNames[idxPass]
projHemeNew=projHemeNew[cellsPass, ]



saveArchRProject(ArchRProj = projHemeNew,load = TRUE,threads=40,overwrite=TRUE)




##Add the DoubletScores
AddPrefix=FALSE
projHemeNew <- addDoubletScores(
    input = projHemeNew,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1,force=TRUE)
saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)


##plot log10(nFrags) and TSSEnrichment
df <- getCellColData(projHemeNew, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = 3.74, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHemeNew, addDOC = FALSE)

##plot the Tss enrichment rideges 
p3 <- plotGroups(
    ArchRProj = projHemeNew, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
plotPDF(p3, name = "QC-Sample-Ridges.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)

#plot violin
p4 <- plotGroups(
    ArchRProj = projHemeNew, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p4, name = "QC-Sample-Boxplot.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)

###filter the binucleated cell
projHemeNew <- filterDoublets(projHemeNew)
saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)
##plot them rerun


###Use the LSI algorithm to reducedims
##sampleCells
##varFeatures
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

#saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)

####add Clusters 
projHemeNew <- addClusters(input = projHemeNew,
                        reducedDims = "IterativeLSI",
                        method = "Seurat",
                        name = "Cluster",
                        resolution = 0.4 ,force = TRUE)


########## Add embedding using UMAP
projHemeNew <- addUMAP(ArchRProj = projHemeNew,
                    reducedDims = "IterativeLSI",
                    name = "UMAP",
                    nNeighbors = 25, minDist = 0.6 ,force = TRUE,
                    metric = "cosine",
                    threads=45)
#coloring by 'Sample', we can also  color by 'Cluster'
p1 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData",name = "Cluster",embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData", 
                    name = "Sample", embedding = "UMAP")

plotPDF(p1,name = "define_Embedding_colorby_cellColData_UMAP_(cellColData and Sample).pdf", 
        ArchRProj = projHemeNew, 
        addDOC = FALSE, 
        width = 4, 
        height = 4)
saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)

########## Add embedding using TSNE
projHemeNew <- addTSNE(
    ArchRProj = projHemeNew, 
    reducedDims = "IterativeLSI", 
    name = "TSNE", 
    perplexity = 30,
    force=TRUE
)
p1 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData", name = "Cluster", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)



###scran
projHemeNew <- addClusters(
    input = projHemeNew,
    reducedDims = "IterativeLSI",
    method = "scran",
    name = "ScranClusters",
    k = 15
)
projHemeNew <- addUMAP(ArchRProj = projHemeNew,
                    reducedDims = "IterativeLSI",
                    name = "ScranClusters",
                    nNeighbors = 25, minDist = 0.6 ,force = TRUE,
                    metric = "cosine",
                    threads=45)


##Harmony detect and process batch effects
projHemeNew <- addHarmony(
    ArchRProj = projHemeNew,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force=TRUE
)

##UMAP use the Harmony matrix
projHemeNew <- addUMAP(
    ArchRProj = projHemeNew, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
p1 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData",name = "Cluster",embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData", 
                    name = "Sample", embedding = "UMAPHarmony")

plotPDF(p1,name = "define_Embedding_colorby_cellColData_UMAPHarmony_(cellColData and Sample).pdf", 
        ArchRProj = projHemeNew, addDOC = FALSE, width = 4, height = 4)
saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)



### Creating a confusion matrix which represent cluster structure
cM_Cluster<- confusionMatrix(i = projHemeNew$Cluster,j = projHemeNew$Sample)
cM_Cluster <- cM_Cluster / Matrix::rowSums(cM_Cluster)
p1 <- pheatmap::pheatmap(
    mat = as.matrix(cM_Cluster), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
plotPDF(p1, name = "Celltype_Sample.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 4, height = 4)


##change the gene name
# genename=read.csv("/md01/changjy/data/Goose_1.0/GeneName.txt",sep="\t",header=F)
# genename$V2=toupper(genename$V2)
# celltype=read.csv("/md01/changjy/data/Goose_1.0/celltype_gene.csv",sep=",",header=F)

##Get the Marker Features
markersGS <- getMarkerFeatures(
    ArchRProj = projHemeNew, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Cluster",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
write.table(as.data.frame(markerList$C4),"/md01/changjy/data/Goose/snATAC_PT_new_filter/markerList$C4.txt",sep="\t",quote=F)




markerGenes=c(
"PRL", "CARTPT","TSHB","POMC-A","TGAS113E22.1","TNFAIP8","HSPB1","RLN3","HPGDS","S100A10","S100A6","WNT5A","CHGA",
  "CHGB","DIO2","CALB2","NDRG1","KRT7","NPBWR2","IGFBP7","CFD","PODXL","TFPI2","RAMP2","JCHAIN","TXNDC5","CRIP1",
  "HBAD","GPX1")
markerGenes=c( "LHX3","NKX2-2","NR4A1","NR4A2", "HES1","ZFP36L1", "NR3C1","NFIB", "ZNF521","NR2F2","SOX2","S100B","POMC-A","PAX7","POU1F1","PRL","NR5A1","SOX11","ZBTB20","TSHB","GATA2","ISL1","ESR1","CEBPD","AR","ASCL1","NEUROD1","E2F8")



markerGenes=c("ALDH1A2","SOX2") #Stem
markerGenes=c("DIO2","TSHB","TRHR")#Thyrotrope
markerGenes=c("GH","DLK1") ##Somatotrope
markerGenes=c("PRL","OLFM1") ##Lactotrope
markerGenes=c("POMC-A","NR3C1","AR","HMGA2")##Corticotrope
markerGenes=c("S100A6", "S100A10", "WNT5A", "S100A1","IL6")

markerGenes=data[which(data[,7]=="Gona"),1]## Gona
markerGenes=data[which(data[,7]=="Lac"),1]##Lac
markerGenes=data[which(data[,7]=="Thy"),1]##Thy
markerGenes=data[which(data[,7]=="Cort"),1]##Cort
markerGenes=data[which(data[,7]=="Endo"),1]##Endo
markerGenes=data[which(data[,7]=="WBC"),1]##WBC
markerGenes=data[which(data[,7]=="RBC"),1]##RB
p=plotEmbedding(
    ArchRProj = projHemeNew, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes[markerGenes%in%gene], 
    embedding = "UMAPHarmony",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(projHemeNew))


# p2 <- lapply(p, function(x){
#     x + guides(color = FALSE, fill = FALSE) + 
#     theme_ArchR(baseSize = 6.5) +
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#     theme(
#         axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(), 
#         axis.text.y=element_blank(), 
#         axis.ticks.y=element_blank()
#    )
#})
#do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(p, name = "Plot-WBC_marker_gene.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)



p=plotEmbedding(
    ArchRProj = projHemeNew, 
    colorBy = "GeneScoreMatrix", 
    name = "NKX2-2", 
    embedding = "UMAPHarmony",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(projHemeNew))



















cell_cls=projHemeNew$Cluster
gene_list <- c("SYT1","SNAP25","SLC32A1","SLC17A6")
#"ERMN","CD9","ELF1","ITM2A","RX1","CCDC153","OLIG2","PDGFRA","PECAM1","FLT1","C1QC","AGT","GFAP")
featureDF <- ArchR:::.getFeatureDF(getArrowFiles(projHemeNew), "GeneScoreMatrix")
gsm <- t(ArchR:::.getMatrixValues(projHemeNew, gene_list, "GeneScoreMatrix",T))
p_list <- lapply(gene_list, function(x){
	data <- data.frame(Cluster=factor(cell_cls,
  levels=rev(c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15"))), Value=gsm[, x])
	p <- ggplot(data, aes(x=Cluster, y=Value, fill=Cluster))+geom_boxplot(width=15)+coord_flip()+
	scale_y_continuous(breaks=10)+
	labs(title=NULL, x=x, y=NULL)+
	theme(legend.position="none", plot.margin=unit(c(-1, 0, -1, 0), "cm"), 
	axis.text.y=element_blank(), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank())+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
	return(p)})
#p_list[[1]] <- p_list[[1]]+theme(axis.text.y=element_text(size=10, colour="black"))
ggsave(plot=patchwork::wrap_plots(plotlist=p_list, nrow=1),
	width=10, height=6, dpi=200, "marker_cls_violin.pdf")








P=plotEmbedding(
    ArchRProj = projHemeNew, 
    colorBy = "GeneScoreMatrix", 
    name =c("FOLR1") , 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(projHemeNew))



plotPDF(P, name = "Plot-FLT1_marker_gene.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)

markerGenes=c("S100B", "S100A6", "S100A1","S100A11","IGFBP7","CFD","HBAD","CRIP1","IFI30")

markerGenes=c("Pou1f1","Prl","Gh","Ghrhr","Tshb","Top2a","Mki67","Pomc","Pax7","Crhr1" ,"Sox2","Vim","Fstl1","Col1a1","Ctss")
markerGenes=c( "RLN3","HPGDS", "CPLX1",  "HEXB", "NR5A1")
markerGenes=c("SOX2","POMC-A","PRL","TSHB","ALDH1A2","NR3C1","DLK1","OLFM1","DIO2")
p=plotEmbedding(
    ArchRProj = projHemeNew, 
    colorBy = "GeneScoreMatrix", 
    name =  toupper(markerGenes), 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(projHemeNew))



plotPDF(p, name = "Plot-TEST_marker_gene.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)


##plot the Heatmap of markerGenes
##T cell
   labelMarkers = markerGenes
heatmapGS <- plotMarkerHeatmap(seMarker = markerGenes,
                           cutOff = "FDR <= 0.01& Log2FC >= 1",
                            labelMarkers=labelMarkers,  ###label the gene name
                           transpose = TRUE)
ComplexHeatmap::draw(heatmapGS,heatmap_legend_side = 'bot', annotation_legend_side = 'bot')
plotPDF(heatmapGS, name = 'GeneScores-marker-heatmap2',width = 8, height = 6, ArchRProj = projHemeNew, addDOC = F)


#Impute weights to our ArchRProject
projHemeNew <- addImputeWeights(projHemeNew)  
p <-plotEmbedding(ArchRProj = projHemeNew,
                  colorBy = 'GeneScoreMatrix',
                  name = markerGenes,
                  embedding = 'UMAP',
                  imputeWeights = getImputeWeights(projHemeNew))
p2 <- lapply(p, function(x){x + 
    guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0,0,0,0), 'cm')) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
})
do.call(cowplot::plot_grid, c(list(ncol =3), p2))    #to plot all genes using cowplot into a single plot.
plotPDF(plotList = p, name = 'plot-UMAP-Marker-genes2-W-Imputation.pdf', ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height =5)


##Find marker peaks
# making pseudo-bulk replicate
projHemeNew <- addGroupCoverages(ArchRProj = projHemeNew, groupBy = 'Cluster') 

# call peak with macs2
pathToMacs2 <- '~/software/miniconda2/bin/macs2'

#Creating the PeakCalls
#Creating the PeakCalls
#Creating the PeakCalls
##Add the Reproducible PeakSet
projHemeNew <- addReproduciblePeakSet(ArchRProj = projHemeNew, 
    groupBy = 'Cluster',
    pathToMacs2 = pathToMacs2,
   #genomeSize=1137194038,
    peakMethod="Tiles",
    method="p")

#Add a new matrix
projHemeNew <- addPeakMatrix(projHemeNew) 
getAvailableMatrices(projHemeNew)



#Get marker peaks
markersPeaks <- getMarkerFeatures(ArchRProj = projHemeNew,
                                  useMatrix ='PeakMatrix',
                                  groupBy = 'Cluster',
                                  bias = c('TSSEnrichment','log10(nFrags)'),
                                  testMethod = 'wilcoxon') 
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")


##plot the Marker peak heatmaps
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks,
                              cutOff = 'FDR <= 0.05 & Log2FC >= 1',
                              transpose = T)
draw(heatmapPeaks, heatmap_legend_side = 'bot', annotation_legend_side = 'bot')
plotPDF(heatmapPeaks, name = 'Peak-Marker-heatmap', width =8, height =6, ArchRProj = projHemeNew, addDOC = F )


###plot the DEP by MA or Volcano
pma <- markerPlot(seMarker = markersPeaks, 
    name = "B",                             ###B cell,we can change the celltype to plot
    cutOff = "FDR <= 0.1 & Log2FC >= 1", 
    plotAs = "MA")

pv <- markerPlot(seMarker = markersPeaks, 
    name = "B",                             ###B cell,we can change the celltype to plot
    cutOff = "FDR <= 0.1 & Log2FC >= 1", 
    plotAs = "Volcano")
plotPDF(pma, pv, name = "B-Markers-peak-MA-Volcano", width = 5, height = 5, ArchRProj = projHemeNew, addDOC = FALSE)

##plot the Browser Track
p <- plotBrowserTrack(
    ArchRProj = projHemeNew, 
    groupBy = "Cluster", 
    geneSymbol = "VCAN",          ## Interested gene list
    upstream = 50000,
    downstream = 50000,
    ylim=c(0, 2)
)
#grid::grid.newpage()
#grid::grid.draw(p$GATA1)

plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projHemeNew, addDOC = FALSE)

###Find marker Features between B cells and T cells
markerTest <- getMarkerFeatures(
  ArchRProj = projHemeNew, 
  useMatrix = "PeakMatrix",
  groupBy = "subtype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Unknown",
  bgdGroups = "T"
)
##plot above by MA and Vocalno
pma <- markerPlot(seMarker = markerTest, name = "B", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "B", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = "Tcells-vs-Bcells-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHemeNew, addDOC = FALSE)
===============================================================================================================
##Add Motif Annotation to projHemeNew
projHemeNew <- addMotifAnnotations(ArchRProj = projHemeNew, motifSet = "cisbp", name = "Motif")
===============================================================================================================
##Find Tcell different between the Tcell and Bcell DEP motif enrichment
motifsUp <- peakAnnoEnrichment(
    seMarker = markerList,
    ArchRProj = projHemeNew,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
##For plotting,create a dataframe including motif,padj,rank
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
#ggplot绘制motif enrichment and use ggrepel to label
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 2,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

##Find Bcell different between the Tcell and Bcell DEP motif enrichment
motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHemeNew,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"   ###Attention: Log2FC<=-0.5
  )
##For plotting,create a dataframe including motif,padj,rank
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
#ggplot绘制motif enrichment and use ggrepel to label
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 2,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
plotPDF(ggUp, ggDo, name = "Unknown-vs-T-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projHemeNew, addDOC = FALSE)
===============================================================================================================
##markersPeaks motif enrichment
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHemeNew,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

##Heatmap of enrichmentMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHemeNew2, addDOC = FALSE)
===============================================================================================================
###Encode
# projHemeNew2 <- addArchRAnnotations(ArchRProj = projHemeNew2, collection = "EncodeTFBS")
# enrichEncode <- peakAnnoEnrichment(
#     seMarker = markersPeaks,
#     ArchRProj = projHemeNew,
#     peakAnnotation = "EncodeTFBS",
#     cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
#   )
# heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
# ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHemeNew2, addDOC = FALSE)
===============================================================================================================
# projHemeNew2 <- addArchRAnnotations(ArchRProj = projHemeNew2, collection = "ATAC")
# enrichATAC <- peakAnnoEnrichment(
#     seMarker = markersPeaks,
#     ArchRProj = projHemeNew2,
#     peakAnnotation = "ATAC",
#     cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
#   )
# heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 7, transpose = TRUE)
# ComplexHeatmap::draw(heatmapATAC, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHemeNew2, addDOC = FALSE)
===============================================================================================================
# projHemeNew <- addArchRAnnotations(ArchRProj = projHemeNew2, collection = "Codex")
# enrichCodex <- peakAnnoEnrichment(
#     seMarker = markersPeaks,
#     ArchRProj = projHemeNew2,
#     peakAnnotation = "Codex",
#     cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
#   )
=============================================================================================================================
##(13)Motif deviation
projHemeNew <- addMotifAnnotations(ArchRProj = projHemeNew, motifSet = "cisbp",name = "Motif")
if('Motif' %ni% names(projHemeNew@peakAnnotation)){
  projHemeNew <- addMotifAnnotations(ArchRProj = projHemeNew, motifSet = 'cisbp', name = 'Motif')
} 
projHemeNew <- addBgdPeaks(projHemeNew) # to add a set of backgroud peaks which are used in computing deviations.
projHemeNew <- addDeviationsMatrix(ArchRProj = projHemeNew,
                                peakAnnotation = 'Motif',
                                #method="ArchR",
                                force = T)
plotVarDev <- getVarDeviations(projHemeNew, name = 'MotifMatrix', plot = T)
plotVarDev # plot these variable deviations
plotPDF(plotVarDev, name = 'Variable-Motif-Deviation-Scores', width = 5, height =5, ArchRProj = projHemeNew, addDOC = F)

## extract a subset of motifs for downstream analysis.
motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- getFeatures(projHemeNew, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep('z:', markerMotifs, value = T)
markerMotifs <- markerMotifs[markerMotifs %ni% 'z:SREBF1_22']
markerMotifs

## plot the distribution of chromVAR deviation scores for each cluster.
p <- plotGroups(ArchRProj = projHemeNew,
                groupBy = 'subtype',
                colorBy = 'MotifMatrix',embedding = 'UMAPHarmony',
                name = markerMotifs,plotAs="violin",
                imputeWeights = getImputeWeights(projHemeNew))
p <- plotGroups(ArchRProj = projHemeNew,
                groupBy = 'cluster',
                colorBy = 'MotifMatrix',embedding = 'UMAPHarmony',
                name = markerMotifs)

## overlay the z-scores on our UMAP embedding as we've done previously for gene scores
p <- plotEmbedding(ArchRProj = projHemeNew, 
                   colorBy = 'MotifMatrix',
                   name = sort(markerMotifs),
                   embedding = 'UMAPHarmony',
                   imputeWeights = getImputeWeights(projHemeNew))
p <- plotEmbedding(ArchRProj = projHemeNew, 
                   colorBy = 'MotifMatrix',
                   name = sort(markerMotifs),
                   embedding = 'UMAPHarmony')
## plot all of these motif UMAPs using cowplot
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(p, name = 'Plot-Groups-Deviations-w-Imputation-UMAP', width = 5, height = 5, ArchRProj = projHemeNew, addDOC = F)
============================================================================================================
###(14)Motif footprinting

motifPositions <- getPositions(projHemeNew) # this creates a GRangesList object where each TF motif is represented by a separate GRanges object.

motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs

projHemeNew <- addGroupCoverages(ArchRProj = projHemeNew, groupBy = "subtype",force=TRUE)

seFoot <- getFootprints(
  ArchRProj = projHemeNew2, 
  positions = motifPositions[markerMotifs], 
  groupBy = "subtype")  #then we can plot them.

##Subtracting by the Tn5 Bias: a second strategy for normalization divides the footprinting signal by the Tn5 bias signal
p<-plotFootprints(seFoot = seFoot,
               ArchRProj = projHemeNew,
               normMethod = "Subtract",
               plotName = 'Plot-Footprints-Subtract-Bias',
               addDOC = F,
               smoothWindow = 5)
##Dividing by the Tn5 Bias: a second strategy for normalization divides the footprinting signal by the Tn5 bias signal

p<-plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHemeNew, 
  normMethod = "Divide",
  plotName = "Plot-Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

seTSS <- getFootprints(
  ArchRProj = projHemeNew2, 
  positions = GRangesList(TSS = getTSS(projHemeNew)), 
  groupBy = "subtype",
  flank = 2000
)

plotFootprints(
  seFoot = seTSS,
  ArchRProj = projHemeNew, 
  normMethod = "None",
  plotName = "TSS-No-Normalization",
  addDOC = FALSE,
  flank = 2000,
  flankNorm = 100
)
#plotPDF(p, name = "Plot-Footprints-Divide-Bias-PorHcelltype.pdf", ArchRProj = projHemeNew2, addDOC = FALSE, width = 6, height = 8)

================================================================================================================================
####Plotting browser tracks of cis-element coaccessibility
projHemeNew <- addCoAccessibility(ArchRProj = projHemeNew, 
                               reducedDims = 'IterativeLSI')
cA <- getCoAccessibility(ArchRProj = projHemeNew,
                         corCutOff = 0.5,
                         resolution = 1,
                         returnLoops = F)
cA

metadata(cA)[[1]]
cA <- getCoAccessibility(ArchRProj = projHemeNew,
                         corCutOff = 0.5,
                         resolution = 1,
                         returnLoops = T)
cA[[1]]

cA <- getCoAccessibility(ArchRProj = projHemeNew,
                         corCutOff = 0.5,
                         resolution = 1000, #when decrease the resolution of our loops to resolution = 1000, this can help with over-plotting of co-accessibility interactions
                         returnLoops = T)

cA[[1]]

cA <- getCoAccessibility(ArchRProj = projHemeNew,
                         corCutOff = 0.5,
                         resolution = 10000,  returnLoops = T) #when decrease the resolution of our loops to resolution = 10000, we identify even fewer co-accessibility interactions

cA[[1]]

##Plotting browser tracks of Co-accessibility
# once we have added co-accessibility information to our ArchRProject, we can use this as a loop track when plotting browser tracks.

markerGenes <- c('CD19','MS4A1','IL4R','CD27','CD38','XBP1','TCL1A','BCL11A','BAFF','LEF1','TBX21','REL','ERG','EBF1','CSMD1','CTNNA2','OPCML','GAS7','CXCR4','CD37','MEF2C','CCR6','BHLHE41','TNFSF9','POU2AF1',
                 'TNFRSF13C','PLCG2','PAX5','KCNN4','EBF1','TGFBR3','ID2','GATA3','RORA','IL7R','PRKCA')
markerGenes <- c('TNFSF13B')
p <- plotBrowserTrack(ArchRProj = projHemeNew,
                      groupBy = "Clusters_1",
                      geneSymbol = markerGenes,
                      upstream = 50000,
                      downstream = 50000,
                      loops = getCoAccessibility(projHemeNew))

p <- plotBrowserTrack(ArchRProj = projHemeNew,
                      groupBy = "sampletype2",
                      geneSymbol = markerGenes,
                      upstream = 50000,
                      downstream = 50000,
                      loops = getCoAccessibility(projHemeNew))

#ArchRBrowser(projHemeNew)
grid::grid.newpage()
plotPDF(plotList = p,
        name = 'Plot-B-sampletype2-Tracks-Marker-Genes-with-CoAccessibility-BAFF.pdf',
        ArchRProj = projHemeNew,
        addDOC = F, width = 5, height = 5)

                returnLoops = T)


##chromVAR bias和z-score bias
seGroupMotif <- getGroupSE(ArchRProj = projHemeNew2, useMatrix = "MotifMatrix", groupBy = "subtype")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
    ArchRProj = projHemeNew,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)


corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]#corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]


corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])



p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )


===========================================================================================================================
####Trajectory construction
trajectory <- c("T", "NK",'Unknown','plasma')
projHemeNew <- addTrajectory(
  ArchRProj = projHemeNew, 
  name = "Trajectory", 
  groupBy = "subtype",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)
head(projHemeNew$Trajectory[!is.na(projHemeNew$Trajectory)])
p1 <- plotTrajectory(projHemeNew, trajectory = "Trajectory", groupBy = "subtype",colorBy = "cellColData", name = "Trajectory",alpha=0.5)
p1[[1]]
plotPDF(p1, name = "Plot-COVHD-B1-Traj-UMAP.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)


trajMM  <- getTrajectory(ArchRProj = projHemeNew, name = "Trajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))

trajGSM <- getTrajectory(ArchRProj = projHemeNew, name = "Trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))

plotPDF(p1, p2, name = "Plot-Trajectory-Traj-Heatmaps.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 6, height = 8)



corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]

trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ] #subset our corresponding trajectory SummarizedExperiment objects to only contain the elements that passed significance above.
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]

#To best order these features, we can create a new trajectory where the values of these two trajectories are multiplied. 
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

#gene score trajectory
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
#the motif trajectory
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
plotPDF(ht1,ht2, name = "Plot-COVHD-Traj-integration-GSM-MM-2.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 6, height = 8)













motifs <-c('REL','RELA','RELB','BCL11A','IRF4','EBF1','POU2F2','BATF','TBX21','LEF1','NFKB1','SPIB','TCF4','NFKB2',
           'JUN','JUND','FOSB','RUNX1','TP73','POU2F2','MEF2C','EOMES','RUNX2')
markerMotifs <- getFeatures(projHemeNew, select = paste(motifs, collapse = '|'), useMatrix = 'Motif-celltypeMatrix')
markerMotifs <- grep('z:', markerMotifs, value = T)
markerMotifs <- markerMotifs[markerMotifs %ni% 'z:SREBF1_22']
markerMotifs

p <- plotGroups(ArchRProj = projHemeNew2, 
  groupBy = "subtype", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projHemeNew2)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = projHemeNew, addDOC = FALSE)


p <- plotEmbedding(
    ArchRProj = projHemeNew, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHemeNew))

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(p, name = "Plot-UMAP-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = projHemeNew, addDOC = FALSE)

markerRNA <- getFeatures(projHemeNew, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
markerRNA

p <- plotEmbedding(
    ArchRProj = projHemeNew, 
    colorBy = "GeneScoreMatrix", 
    name = sort(markerRNA), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHemeNew)
)


p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(p, name = "Plot-UMAP-mRNA-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = projHemeNew2, addDOC = FALSE)





=================================================================================================================================
============================================================================================================






ego <- enrichGO(gene  = a$1,keyType = "SYMBOL",universe = names(geneList),OrgDb= org.Gg.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE)




p <- plotBrowserTrack(
    ArchRProj = projHemeNew, 
    groupBy = "Cluster", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)


pdf("GO_cluster.pdf")
erich.go.BP <- enrichGO(gene=as.factor(markerList$Plasma$name),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05)
ego <- enrichKEGG(
          gene = markerList$Unknown$name,
          keyType = "kegg",
          organism  = 'hsa',
          pvalueCutoff  = 0.05,
          pAdjustMethod  = "BH",
          qvalueCutoff  = 0.05
)

p1=dotplot(erich.go.BP)
p1
erich.go.BP <- enrichGO(gene=as.factor(markerList$C2$name),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05)
p2=dotplot(erich.go.BP)
p2
erich.go.BP <- enrichGO(gene=as.factor(markerList$C3$name),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05)
p3=dotplot(erich.go.BP)
p3
erich.go.BP <- enrichGO(gene=as.factor(markerList$C4$name),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05)
p4=dotplot(erich.go.BP)
p4
erich.go.BP <- enrichGO(gene=as.factor(markerList$C5$name),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05)
p5=dotplot(erich.go.BP)
p5
erich.go.BP <- enrichGO(gene=as.factor(markerList$C6$name),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05)
p6=dotplot(erich.go.BP)
p6
erich.go.BP <- enrichGO(gene=as.factor(markerList$C7$name),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05)
p7=dotplot(erich.go.BP)
p7

plotPDF(p1,p2,p3,p4,p5,p6,p7, name = "GO_cluster.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)






#LMPP
#CLP
#proB
#preB
#Naive B
#Memory B
#plasma
#Naive CD4+ T
#Naive CD8+ T
#Memory CD4' T
#Memory CD8+T
#NK
#CMP
#GMP
#MDP
#pDC
#Early.ery
#Late.ery
#Early.baso
#CD14mono
#CD16mono




markerGenes=c(
    "FCN2",
    "MME",
    "CD19","IL7R",
    "CD79B","MS4A1",
    "IL4R","PAX5","MS4A1",
    "PAX5","MS4A1",
    "SDC1",
    "ITGA6","CCR7",
    "CD8A","LEF1",
    "CD4","CD52",
    "CD8A","EOMES",
    "GNLY","NKG7","FCGR3B","FASLG",
    "TAL1","GATA1",
    "SPII","MPO",
    "FLT3","MPO",
    "DERL3","FLT3",
    "HBB","GATA1",
    "GATA1",
    "PF4","ITGA2B",
    "NCR1","CEBPB",
    "FCGR3A","SIGLEC10")
DNAJB8
c("ENPP3","CD63","KIT", "ITGAM", "IL3RY", "FUT4" ) ##Basophil
c("PF4", "PLA2G12A", "PPBP","CD133", "KDR", "CD31", "CD34", "CD45") ##progenitor
c("CD14", "CD163", "CD68", "CSF1R", "FCGR3A")   ##Macrophage
c("CD34","HPRT1","CD105", "CD29", "NANOG")#stem cell
c("CD3G","CD8A","CD3E","CCR7", "CD27", "IL7R","ENTPD1","CD4", "CTLA4", "FOXP3", "IL2RA") ##T
pdf("Myeloid.pdf")

plotEmbedding(
    ArchRProj = projHemeNew, 
    colorBy = "GeneScoreMatrix", 
    name =  "CD163"  ,
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(projHemeNew))
)
dev.off()

pdf("UMAP-pbmc-cell-type.pdf")
df = data.frame(projHemeNew@embeddings$UMAP$df, projHemeNew$subtype)
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = projHemeNew.subtype,size=10)) +
  geom_point(size=0.2,alpha=0.5) + 
  theme_bw() +
  #scale_color_manual(values = c("#ffcbcb","#0072B5FF","#FCF5B3","#E18727FF","#21DB00","#66E1E6")) + 
  theme(legend.title=element_text(size=10),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Celltype")
  dev.off()


eg = bitr(markerList$Unknown$name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  kk <- enrichKEGG(gene = eg$ENTREZID, 
                 organism ='hsa',
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 minGSSize = 1)
                 #readable = TRUE ,
                 #use_internal_data =FALSE)
  barplot(kk)






featureDF <- ArchR:::.getFeatureDF(getArrowFiles(projHemeNew), "GeneScoreMatrix")
gsm <- t(ArchR:::.getMatrixValues(projHemeNew, featureDF$name[match(markerGenes, gsub("\\..*", "", featureDF$name))], "GeneScoreMatrix", TRUE))
write.csv(gsm, "archr_gs_mat.csv")





name=read.table("/md01/changjy/data/Goose_1.0/ref/Anser_cygnoides_chr/New_GeneName.txt")
genes=as.data.frame(genes)
for (i in c(1:dim(genes)[1])){
res=name[which(name[,1]==genes[i,6]),2]
genes[i,6]=res
genes[i,7]=res
}
 write.table(genes,"/md01/changjy/data/Goose_1.0/ref/Anser_cygnoides_chr/rename_sorted.genes.txt",sep="\t",quote=F,row.names=F)

exons=as.data.frame(exons)
for (i in c(1:dim(exons)[1])){
res=name[which(name[,1]==exons[i,6]),2]
exons[i,6]=res
exons[i,7]=res
} write.table(exons,"/md01/changjy/data/Goose_1.0/ref/Anser_cygnoides_chr/rename_sorted.exons.txt",sep="\t",quote=F,row.names=F)

TSS=as.data.frame(TSS)
for (i in c(1:dim(TSS)[1])){
res=name[which(name[,1]==TSS[i,6]),2]

TSS[i,6]=res
TSS[i,7]=res
}
 write.table(TSS,"/md01/changjy/data/Goose_1.0/ref/Anser_cygnoides_chr/rename_sorted.TSS.txt",sep="\t",quote=F,row.names=F)

