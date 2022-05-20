rm(list = ls()) 
library('ArchR')
library('pheatmap')
library('Rsamtools')
library('scran')
library("harmony")
library(GenomicFeatures)
library(BSgenome.goose.ENSEMBLE.ac23)
# chr1=read.table("/md01/changjy/data/Goose/1E-ME/outs/peak_annotation_chrName.txt",header=T)
# chr2=read.table("/md01/changjy/data/Goose/2-ME/outs/peak_annotation_chrName.txt",header=T)
# chr3=read.table("/md01/changjy/data/Goose/3E-ME/outs/peak_annotation_chrName.txt",header=T)
# chr4=read.table("/md01/changjy/data/Goose/4ME/outs/peak_annotation_chrName.txt",header=T)
# aaa=intersect(chr1$chrom,chr2$chrom)
# bbb=intersect(chr3$chrom,chr4$chrom)
# ccc=intersect(aaa,bbb)
# gtf <- read.table("/md01/changjy/data/Goose/Goose_bulk/ref/ensemble/new_Anser_cygnoides.GooseV1.0.101.gtf", sep="\t")
# gtf=gtf[gtf[,1]%in%ccc,]

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

# library(BSgenome.Acygnoides.Assembly.go12)
# chromsize<-read.table("/md01/changjy/data/Goose/Goose_bulk/ref/ensemble/Anser_cygnoides_ensemble/fasta/genome.fa.fai")
# chromsize<-chromsize[chromsize[,1]%in%ccc,]
# chromsize$V2=chromsize$V2+1
# gr <- GRanges(as.character(chromsize$V1),IRanges(start=rep(0,dim(chromsize)[1]),end=chromsize$V2))
# gr$seqinfo=gr$seqinfo[gr$seqinfo$seqnames%in%ccc,]
# genomeAnnotation<-createGenomeAnnotation(BSgenome.Acygnoides.Assembly.go12,gr)
# genomeAnnotation$chromSizes<-gr 
# genomeAnnotation$chromSizes@seqinfo@seqlengths<-as.factor(chromsize$V2)



# #Import the file
# inputfiles<-c("/md01/changjy/data/Goose/1E-ME/outs/fragments.tsv.gz",

#  "/md01/changjy/data/Goose/2-ME/outs/fragments.tsv.gz",

#  "/md01/changjy/data/Goose/3E-ME/outs/fragments.tsv.gz",

#  "/md01/changjy/data/Goose/4ME/outs/fragments.tsv.gz")
# samplenames<-c("1E-ME","2-ME","3E-ME","4ME")


# ##Creating Arrow files
# AddPrefix=FALSE
# ArrowFiles <- createArrowFiles(
#     inputFiles=inputfiles,
#     sampleNames=samplenames,
#     geneAnnotation=geneAnnotation,
#     genomeAnnotation=genomeAnnotation,
#     minTSS = 1, 
#     minFrags = 100, 
#     addTileMat = TRUE,
#     addGeneScoreMat = TRUE,
#     force=TRUE,
#     subThreading=FALSE,
#     threads=1,
#     cleanTmp = FALSE,
#     nChunk = 1,
#     )
# projHemeNew <- ArchRProject(ArrowFiles = ArrowFiles,
#                          outputDirectory = "/md01/changjy/data/Goose/snATAC_noprocess",
#                          copyArrows = TRUE,
#                          geneAnnotation=geneAnnotation,
#                          genomeAnnotation=genomeAnnotation)

# saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)###Creating the geneAnnotation
projHemeNew<-loadArchRProject("/md01/changjy/data/Goose/snATAC")
##Add the DoubletScores
AddPrefix=FALSE
# projHemeNew <- addDoubletScores(
#     input = projHemeNew,
#     k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#     knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#     LSIMethod = 1,
#     threads=1)
# saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)


##plot log10(nFrags) and TSSEnrichment
# df <- getCellColData(projHemeNew, select = c("log10(nFrags)", "TSSEnrichment"))
# p <- ggPoint(
#     x = df[,1], 
#     y = df[,2], 
#     colorDensity = TRUE,
#     continuousSet = "sambaNight",
#     xlabel = "Log10 Unique Fragments",
#     ylabel = "TSS Enrichment",
#     xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
#     ylim = c(0, quantile(df[,2], probs = 0.99))
# ) + geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = 3.74, lty = "dashed")
# plotPDF(p, name = "unprocess_TSS-vs-Frags.pdf", ArchRProj = projHemeNew, addDOC = FALSE)

# ##plot the Tss enrichment rideges 
# p3 <- plotGroups(
#     ArchRProj = projHemeNew, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "TSSEnrichment",
#     plotAs = "ridges"
#    )
# plotPDF(p3, name = "QC-Sample-Ridges.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)

# #plot violin
# p4 <- plotGroups(
#     ArchRProj = projHemeNew, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "TSSEnrichment",
#     plotAs = "violin",
#     alpha = 0.4,
#     addBoxPlot = TRUE
#    )
# plotPDF(p4, name = "QC-Sample-Boxplot.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)
# projHemeNew<-loadArchRProject("/md01/changjy/data/Goose/snATAC/")
# AddPrefix=FALSE
# # projHemeNew <- addDoubletScores(
# #     input = projHemeNew,
# #     k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
# #     knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
# #     LSIMethod = 1,force=TRUE)

# # projHemeNew <- filterDoublets(projHemeNew)
# # saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)

# projHemeNew <- addIterativeLSI(ArchRProj = projHemeNew,
#                             useMatrix = "TileMatrix",
#                             name = "IterativeLSI",
#                             iterations = 2,
#                             clusterParams = list(resolution = c(0.2),
#                                                  sampleCells =13117,
#                                                  n.start = 10),
#                             varFeatures = 15000,
#                             dimsToUse = 1:25,
#                             force = TRUE,
#                             threads=40)

# saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)

# ####add Clusters 
# projHemeNew <- addClusters(input = projHemeNew,
#                         reducedDims = "IterativeLSI",
#                         method = "Seurat",
#                         name = "Cluster",
#                         resolution = 0.4 ,force = TRUE)


# ########## Add embedding using UMAP
# projHemeNew <- addUMAP(ArchRProj = projHemeNew,
#                     reducedDims = "IterativeLSI",
#                     name = "UMAP",
#                     nNeighbors = 25, minDist = 0.6 ,force = TRUE,
#                     metric = "cosine",
#                     threads=45)
# #coloring by 'Sample', we can also  color by 'Cluster'
# p1 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData",name = "Cluster",embedding = "UMAP")
# p2 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData", 
#                     name = "Sample", embedding = "UMAP")

# plotPDF(p1,name = "define_Embedding_colorby_cellColData_UMAP_(cellColData and Sample).pdf", 
#         ArchRProj = projHemeNew, 
#         addDOC = FALSE, 
#         width = 4, 
#         height = 4)
# saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)

# ########## Add embedding using TSNE
# projHemeNew <- addTSNE(
#     ArchRProj = projHemeNew, 
#     reducedDims = "IterativeLSI", 
#     name = "TSNE", 
#     perplexity = 30,
#     force=TRUE
# )
# p1 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
# p2 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData", name = "Cluster", embedding = "TSNE")
# ggAlignPlots(p1, p2, type = "h")
# plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)



# ###scran
# projHemeNew <- addClusters(
#     input = projHemeNew,
#     reducedDims = "IterativeLSI",
#     method = "scran",
#     name = "ScranClusters",
#     k = 15
# )
# projHemeNew <- addUMAP(ArchRProj = projHemeNew,
#                     reducedDims = "IterativeLSI",
#                     name = "ScranClusters",
#                     nNeighbors = 25, minDist = 0.6 ,force = TRUE,
#                     metric = "cosine",
#                     threads=45)


# ##Harmony detect and process batch effects
# projHemeNew <- addHarmony(
#     ArchRProj = projHemeNew,
#     reducedDims = "IterativeLSI",
#     name = "Harmony",
#     groupBy = "Sample",
#     force=TRUE
# )

# ##UMAP use the Harmony matrix
# projHemeNew <- addUMAP(
#     ArchRProj = projHemeNew, 
#     reducedDims = "Harmony", 
#     name = "UMAPHarmony", 
#     nNeighbors = 30, 
#     minDist = 0.5, 
#     metric = "cosine"
# )
# p1 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData",name = "Cluster",embedding = "UMAPHarmony")
# p2 <- plotEmbedding(ArchRProj = projHemeNew, colorBy = "cellColData", 
#                     name = "Sample", embedding = "UMAPHarmony")

# plotPDF(p1,name = "define_Embedding_colorby_cellColData_UMAPHarmony_(cellColData and Sample).pdf", 
#         ArchRProj = projHemeNew, addDOC = FALSE, width = 4, height = 4)




# ### Creating a confusion matrix which represent cluster structure
# cM_Cluster<- confusionMatrix(i = projHemeNew$Cluster,j = projHemeNew$Sample)
# cM_Cluster <- cM_Cluster / Matrix::rowSums(cM_Cluster)
# p1 <- pheatmap::pheatmap(
#     mat = as.matrix(cM_Cluster), 
#     color = paletteContinuous("whiteBlue"), 
#     border_color = "black"
# )
# plotPDF(p1, name = "Celltype_Sample.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 4, height = 4)



# ##Get the Marker Features
# markersGS <- getMarkerFeatures(
#     ArchRProj = projHemeNew, 
#     useMatrix = "GeneScoreMatrix", 
#     groupBy = "Cluster",
#     bias = c("TSSEnrichment", "log10(nFrags)"),
#     testMethod = "wilcoxon"
# )
# markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)
# projHemeNew <- addImputeWeights(projHemeNew)
# markerGenes=c("PRL","DIO2")
# p=plotEmbedding(
#     ArchRProj = projHemeNew, 
#     colorBy = "GeneScoreMatrix", 
#     name = markerGenes, 
#     embedding = "UMAP",
#     quantCut = c(0.01, 0.95),
#     imputeWeights = getImputeWeights(projHemeNew))
# plotPDF(p, name = "Plot-11_marker_gene.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)
# saveArchRProject(ArchRProj = projHemeNew, load = TRUE,threads=40,overwrite=TRUE)
# # markerGenes=c("Pou1f1","Prl","Gh","Ghrhr","Tshb","Top2a","Mki67","Pomc","Pax7","Crhr1" ,"Sox2","Vim","Fstl1","Col1a1","Ctss")
# # p=plotEmbedding(
# #     ArchRProj = projHemeNew, 
# #     colorBy = "GeneScoreMatrix", 
# #     name =  toupper(markerGenes), 
# #     embedding = "UMAP",
# #     quantCut = c(0.01, 0.95),
# #     imputeWeights = getImputeWeights(projHemeNew))
# # plotPDF(p, name = "Plot-11_marker_gene.pdf", ArchRProj = projHemeNew, addDOC = FALSE, width = 5, height = 5)
addArchRThreads(1)
projHemeNew <- addGroupCoverages(ArchRProj = projHemeNew, groupBy = 'Cluster',force=TRUE) 

# call peak with macs2
pathToMacs2 <- findMacs2()

#Creating the PeakCalls
#Creating the PeakCalls
#Creating the PeakCalls
##Add the Reproducible PeakSet
projHemeNew <- addReproduciblePeakSet(
    ArchRProj = projHemeNew, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2,
    force=TRUE
)
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
