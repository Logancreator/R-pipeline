rm(list = ls()) 
library('ArchR')
library('pheatmap')
library('Rsamtools')
library('scran')
library("harmony")
library(GenomicFeatures)
library(BSgenome.Zmays.Ensemble.zmv4)
setwd("/public/home/changjianye/project/scatac/P/")

addArchRThreads(8)
projHeme2 <- loadArchRProject("/public/home/changjianye/project/scatac/P/Save-ProjHeme_P-filterDoublets/")
# ##双胞计算
# projHeme2 <- addDoubletScores(
#     input = projHeme2,
#     k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#     knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#     LSIMethod = 1
# )
# ##去除双胞
# projHeme2 <- filterDoublets(projHeme2)

# ##save
# saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme_P-v1.0/",threads=40,overwrite=TRUE,load = TRUE)

##LSI
projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = length(projHeme2$cellNames), 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    projectCellsPre=TRUE,
    saveIterations=TRUE,
    threads=5,
    seed=1,
    force=TRUE)

projHeme2 <- addIterativeLSI( 
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix",
    name = "IterativeLSI2", 
    iterations = 2, 
    clusterParams = list(
        resolution = c(1),
        sampleCells = length(projHeme2$cellNames), 
        n.start = 10), 
    force=TRUE,
    varFeatures = 20000, 
    dimsToUse = 1:15)





##Batch Effect Correction wtih Harmony
projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)

projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI2",
    name = "Harmony2",
    groupBy = "Sample",
    force = TRUE
)






##Seurat降维
projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = TRUE
)

projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI2",
    method = "Seurat",
    name = "Clusters2",
    resolution = 0.8,
    force = TRUE
)




# ##细胞组成
# library(pheatmap)
# cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
# cM <- cM / Matrix::rowSums(cM)
# p <- pheatmap::pheatmap(
#     mat = as.matrix(cM), 
#     color = paletteContinuous("whiteBlue"), 
#     border_color = "black"
# )

#Clustering using Seurat’s FindClusters() function
projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)

projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI2", 
    name = "UMAP2", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)



# Clustering using scran
projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI",
    method = "scran",
    name = "ScranClusters",
    k = 15,
    force = TRUE
)


projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI2",
    method = "scran",
    name = "ScranClusters2",
    k = 15,
    force = TRUE
)



p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP2")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP2")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters2.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)



p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "ScranClusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-ScranClusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP2")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "ScranClusters2", embedding = "UMAP2")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-ScranClusters2.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)



#t-Stocastic Neighbor Embedding (t-SNE)
projHeme2 <- addTSNE(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI", 
    name = "TSNE", 
    perplexity = 30,
    force = TRUE
)

projHeme2 <- addTSNE(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI2", 
    name = "TSNE2", 
    perplexity = 30,
    force = TRUE
)



p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample2", embedding = "TSNE2")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters2", embedding = "TSNE2")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters2.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)


p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "ScranClusters", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-tSNE-Sample-ScranClusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)


p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE2")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "ScranClusters2", embedding = "TSNE2")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-tSNE-Sample-ScranClusters2.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)



#Dimensionality Reduction After Harmony
projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)


p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters2", embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters2.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)



projHeme2 <- addTSNE(
    ArchRProj = projHeme2, 
    reducedDims = "Harmony", 
    name = "TSNEHarmony", 
    perplexity = 30,
    force = TRUE
)

p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Plot-TSNE2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters2", embedding = "TSNEHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Plot-TSNE2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme_P-filterDoublets/",threads=40,overwrite=TRUE,load = TRUE)


# ##addImputeWeights
# projHeme2 <- addImputeWeights(projHeme2)

# ##Making Pseudo-bulk Replicates
# projHeme2 <- addGroupCoverages(ArchRProj = projHeme2, groupBy = "Clusters")
# pathToMacs2 <- findMacs2()
# projHeme2 <- addReproduciblePeakSet(
#     ArchRProj = projHeme2, 
#     groupBy = "Clusters", 
#     pathToMacs2 = pathToMacs2
# )
# getPeakSet(projHeme2)

# ##save
# saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme_P-filterDoublets/",threads=40,overwrite=TRUE,load = TRUE)

# #Add Peak Matrix
# projHeme2 <- addPeakMatrix(projHeme2)
# getAvailableMatrices(projHeme2)

# #Identifying Marker Peaks with ArchR
# markersPeaks <- getMarkerFeatures(
#     ArchRProj = projHeme2, 
#     useMatrix = "PeakMatrix", 
#     groupBy = "Clusters",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markersPeaks
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
# markerList


# ##Marker Peak MA and Volcano Plots
# pma <- markerPlot(seMarker = markersPeaks, name = "C1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
# pma

# pv <- markerPlot(seMarker = markersPeaks, name = "C1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
# pv
# plotPDF(pma, pv, name = "C1-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

# #Pairwise Testing Between Groups
# markerTest <- getMarkerFeatures(
#   ArchRProj = projHeme2, 
#   useMatrix = "PeakMatrix",
#   groupBy = "Clusters",
#   testMethod = "wilcoxon",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   useGroups = "C1",
#   bgdGroups = "C2"
# )
# pma <- markerPlot(seMarker = markerTest, name = "C1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
# pma
# pv <- markerPlot(seMarker = markerTest, name = "C1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
# pv
# plotPDF(pma, pv, name = "C1-vs-C2-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

# saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme_P-filterDoublets/",threads=40,overwrite=TRUE,load = TRUE)
# dev.off()
