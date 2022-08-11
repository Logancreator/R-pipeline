library(ArchR)
library(ggplot2)
set.seed(1)

addArchRGenome("mm10")

###add new cell type
all_new =loadArchRProject("/md01/nieyg/scATAC-ArchR/Save-allcelltype-rescue2")
all_new =all_new[which(all_new$subtypename =="Endothelial" |all_new$subtypename =="Fibroblast" | 
  all_new$subtypename =="Epicardial" |all_new$subtypename =="Smoothmuscle" |
  all_new$subtypename =="B_cell" |all_new$subtypename =="Glial" |
  all_new$subtypename =="granulocyte" |all_new$subtypename =="Macrophages" |
  all_new$subtypename =="Pericyte" |all_new$subtypename =="T_cell"),]

#plot cell type in umap
df = data.frame(projHeme2@embeddings$UMAP$df, projHeme2$celltype)
celltype = projHeme2$celltype

pdf("UMAP-all-NCM-cell-type-split-IC-size0.4-alpha0.5.pdf")
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = celltype)) +
  geom_point() + 
  theme_bw() + 
  # scale_color_manual(values = c("#ffcbcb","#0072B5FF","#FCF5B3","#E18727FF","#21DB00","#66E1E6",
  #   "#20854EFF","#1f441e","#7876B1FF","#0F3057")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Celltype")
dev.off()

#plot genotype in umap
meta<-matrix(data=NA,nrow=57899,ncol=2);
colnames(meta)=c("Sample","Sample2");
meta[,1] = all_new$Sample
cl <- meta[,1] 
for (i in 1:57899){
#  if(cl[i] == "AR3" |cl[i] ==  "AR7" |cl[i] == "AR14")
#    meta[i,2]="AR"
#  if(cl[i] == "NAR3" |cl[i] ==  "NAR7" |cl[i] == "NAR14")
#      meta[i,2]="NAR"
  if(cl[i] == "P3" |cl[i] ==  "P7" |cl[i] == "P14")
      meta[i,2]="P"
  else 
    meta[i,2]="Other"
}
meta = data.frame(meta)
table(meta$Sample2)
all_new$treatment = meta$Sample2
table(all_new$treatment)

df = data.frame(all_new@embeddings$UMAP$df, all_new$treatment)
Sample = all_new$treatment

pdf("./plot/10celltype/P-UMAP-emphasize3.pdf")
ggplot(df, aes(x = IterativeLSI.dim8.UMAP_Dimension_1 , y = IterativeLSI.dim8.UMAP_Dimension_2,color = Sample)) +
  geom_point(size=0.05,alpha=0.5) + 
  theme_bw() + 
  scale_color_manual(values = c("#D5D2CD" ,"#65B953")) + 
#  scale_color_manual(values = c("#E64B35B2","#00468BB2","#42B540B2" )) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Genotype")
dev.off()
#P #65B953 AR #D23831 #NAR #172168 #other #D5D2CD

#plot timepoint in umap
meta<-matrix(data=NA,nrow=57899,ncol=2);
colnames(meta)=c("Sample","Sample2");
meta[,1] = all_new$Sample
cl <- meta[,1]
for (i in 1:57899){
#  if(cl[i] == "AR3" |cl[i] ==  "NAR3" |cl[i] == "P3")
#    meta[i,2]="Day3"
#  if(cl[i] == "AR7" |cl[i] ==  "NAR7" |cl[i] == "P7")
#      meta[i,2]="Day7"
  if(cl[i] == "AR14" |cl[i] ==  "NAR14" |cl[i] == "P14")
      meta[i,2]="Day14"
  else 
    meta[i,2]="Other"
}
meta = data.frame(meta)
table(meta$Sample2)
all_new$timepoint = meta$Sample2
table(all_new$timepoint)

df = data.frame(all_new@embeddings$UMAP$df, all_new$timepoint)
Timepoint = all_new$timepoint

pdf("./plot/10celltype/Day14-UMAP-emphasize.pdf")
ggplot(df, aes(x = IterativeLSI.dim8.UMAP_Dimension_1 , y = IterativeLSI.dim8.UMAP_Dimension_2,color = Timepoint)) +
  geom_point(size=0.05,alpha=0.5) + 
  theme_bw() + 
  scale_color_manual(values = c("#79BFB5","#D5D2CD")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Timepoint")
dev.off()
#3 #E7922A 7 #B98ABB 14 #79BFB5

#plot track of specific gene of cell types 
markergene <- c("Cd74", #B cell
	"Adgrf5", # EC
	"Bnc1",#Epi
	"Pdgfra", #FB
	"Plp1",# Glial
	"S100A9", #Granulocy
	"Adgre1",# MP
	"Fgfr1", #Pericyte
	"Abcc9", #SMC
	"Cd3e" #T cell
	)
p <- plotBrowserTrack(
    ArchRProj = all_new, 
    groupBy = "subtypename", 
    geneSymbol = markergene, 
    upstream = 2000,
    downstream = 2000,
    pal =c("#48453d","#0072B5FF","#fbaba9","#E18727FF","#b92419","#66E1E6",
    "#20854EFF","#8edcaf","#7876B1FF","#0F3057")
    )

#heatmap of specific gene in single cell
###extract matrix
GSM<-getMatrixFromProject(ArchRProj = all_new, useMatrix = "GeneScoreMatrix")
head(GSM@assays@data$GeneScoreMatrix)[,1:10]
genename<-rowData(GSM)$name ####�?24333个基�?

EC <- all_new[which(all_new$subtypename =="Endothelial"),]
FB <- all_new[which(all_new$subtypename =="Fibroblast"),]
Epi <- all_new[which(all_new$subtypename =="Epicardial"),]
SMC <- all_new[which(all_new$subtypename =="Smoothmuscle"),]
B <- all_new[which(all_new$subtypename =="B_cell"),]
Gl <- all_new[which(all_new$subtypename =="Glial"),]
Gran <- all_new[which(all_new$subtypename =="granulocyte"),]
MP <- all_new[which(all_new$subtypename =="Macrophages"),]
Per <- all_new[which(all_new$subtypename =="Pericyte"),]
T <- all_new[which(all_new$subtypename =="T_cell"),]

ECcount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(EC))]
FBcount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(FB))]
Epicount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(Epi))]
SMCcount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(SMC))]
Bcount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(B))]
Glcount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(Gl))]
Grancount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(Gran))]
MPcount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(MP))]
Percount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(Per))]
Tcount<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(T))]

ECcount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(EC))]
FBcount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(FB))]
Epicount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(Epi))]
SMCcount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(SMC))]
Bcount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(B))]
Glcount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(Gl))]
Grancount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(Gran))]
MPcount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(MP))]
Percount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(Per))]
Tcount<-GSM@assays@data$GeneIntegrationMatrix[,which(colnames(GSM) %in% rownames(T))]

ECcount<-as.matrix(ECcount)
FBcount<-as.matrix(FBcount)
Epicount<-as.matrix(Epicount)
SMCcount<-as.matrix(SMCcount)
Bcount<-as.matrix(Bcount)
Glcount<-as.matrix(Glcount)
Grancount<-as.matrix(Grancount)
MPcount<-as.matrix(MPcount)
Percount<-as.matrix(Percount)
Tcount<-as.matrix(Tcount)

ECcell =colnames(ECcount)
FBcell =colnames(FBcount)
Epicell =colnames(Epicount)
SMCcell =colnames(SMCcount)
Bcell =colnames(Bcount)
Glcell =colnames(Glcount)
Grancell =colnames(Grancount)
MPcell =colnames(MPcount)
Percell =colnames(Percount)
Tcell =colnames(Tcount)

ECcount1<-sample(ECcell,100)
FBcount1<-sample(FBcell,100)
Epicount1<-sample(Epicell,100)
SMCcount1<-sample(SMCcell,100)
Bcount1<-sample(Bcell,100)
Glcount1 =sample(Glcount,100)
Grancount1<-sample(Grancell,100)
MPcount1<-sample(MPcell,100)
Percount1<-sample(Percell,100)
Tcount1<-sample(Tcell,100)


head(ECcount2[,1:5])
ECcount2<-ECcount[,which(colnames(ECcount) %in% ECcount1)]
FBcount2<-FBcount[,which(colnames(FBcount) %in% FBcount1)]
Epicount2<-Epicount[,which(colnames(Epicount) %in% Epicount1)]
SMCcount2<-SMCcount[,which(colnames(SMCcount) %in% SMCcount1)]
Bcount2<-Bcount[,which(colnames(Bcount) %in% Bcount1)]
Glcount2<-Glcount[,which(colnames(Glcount) %in% Glcount1)]
Grancount2<-Grancount[,which(colnames(Grancount) %in% Grancount1)]
MPcount2<-MPcount[,which(colnames(MPcount) %in% MPcount1)]
Percount2<-Percount[,which(colnames(Percount) %in% Percount1)]
Tcount2<-Tcount[,which(colnames(Tcount) %in% Tcount1)]

ECcount2 <-  data.frame(cbind(genename,ECcount2))
FBcount2 <-  data.frame(cbind(genename,FBcount2))
Epicount2 <-  data.frame(cbind(genename,Epicount2))
SMCcount2 <-  data.frame(cbind(genename,SMCcount2))
Bcount2 <-  data.frame(cbind(genename,Bcount2))
Glcount2 <-  data.frame(cbind(genename,Glcount2))
Grancount2 <-  data.frame(cbind(genename,Grancount2))
MPcount2 <-  data.frame(cbind(genename,MPcount2))
Percount2 <-  data.frame(cbind(genename,Percount2))
Tcount2 <-  data.frame(cbind(genename,Tcount2))

rownames(ECcount2) = ECcount2$genename
rownames(FBcount2) = FBcount2$genename
rownames(Epicount2) = Epicount2$genename
rownames(SMCcount2) = SMCcount2$genename
rownames(Bcount2) = Bcount2$genename
rownames(Glcount2) = Glcount2$genename
rownames(Grancount2) = Grancount2$genename
rownames(MPcount2) = MPcount2$genename
rownames(Percount2) = Percount2$genename
rownames(Tcount2) = Tcount2$genename

ECcount2 <- ECcount2[,-1]
FBcount2 <- FBcount2[,-1]
Epicount2 <- Epicount2[,-1]
SMCcount2 <- SMCcount2[,-1]
Bcount2 <- Bcount2[,-1]
Glcount2 <- Glcount2[,-1]
Grancount2 <- Grancount2[,-1]
MPcount2 <- MPcount2[,-1]
Percount2 <- Percount2[,-1]
Tcount2 <- Tcount2[,-1]

SUM<-cbind(ECcount2,FBcount2)
SUM<-cbind(SUM,SMCcount2) 
SUM<-cbind(SUM,Epicount2)
SUM<-cbind(SUM,Percount2) 
SUM<-cbind(SUM,Glcount2) 
SUM<-cbind(SUM,MPcount2) 
SUM<-cbind(SUM,Tcount2) 
SUM<-cbind(SUM,Bcount2) 
SUM<-cbind(SUM,Grancount2) 

#marker gene
Endo <- c("Cdh5","Kdr","Adgrf5","Cldn5")
Fibr<- c("Col3a1","Col6a3","Col6a1","Ddr2","Col5a1", "Pdgfra")
Smc <- c("Rgs5","Pdgfrb","Accsl","Abcc9","Rgs4") 
epi <- c("Ezr","Efna5","Bnc1","Mbp", "Prr15","Slc9a3r1","Lrrn4")
per <- c("Gfpt2","Kcnj8","Slc1a7","Fgfr1")
glcell <- c("Chl1","Kcna6","Slc35f1","Plp1")
mpcell <- c("Cd86","Mpeg1","Cd68","Ccl3","Cd44")
tcell <- c("Cd3e","Cd3d", "Cd8a", "Cd8b1","nkg7","Igfbp4","Lat","Cd53","Itk")
bcell <- c("Cd74","Cd22","Cd79b","Cxcr4")
grancell <- c("Samsn1","Clec4d","S100A9", "S100A8")

markergene <- c(Endo,Fibr,Smc,epi,per,glcell,mpcell,tcell,bcell,grancell)
SUM_table <- SUM

name <- rownames(heatmapGS) 

SUM_sel1 =SUM_table[which(rownames(SUM_table)%in% Endo),]
SUM_sel2 =SUM_table[which(rownames(SUM_table)%in% Fibr),]
SUM_sel3 =SUM_table[which(rownames(SUM_table)%in% Smc),]
SUM_sel4 =SUM_table[which(rownames(SUM_table)%in% epi),]
SUM_sel5 =SUM_table[which(rownames(SUM_table)%in% per),]
SUM_sel6 =SUM_table[which(rownames(SUM_table)%in% glcell),]
SUM_sel7 =SUM_table[which(rownames(SUM_table)%in% mpcell),]
SUM_sel8 =SUM_table[which(rownames(SUM_table)%in% tcell),]
SUM_sel9 =SUM_table[which(rownames(SUM_table)%in% bcell),]
SUM_sel10 =SUM_table[which(rownames(SUM_table)%in% grancell),]

SUM_sel <- rbind(SUM_sel1,SUM_sel2)
SUM_sel <- rbind(SUM_sel,SUM_sel3)
SUM_sel <- rbind(SUM_sel,SUM_sel4)
SUM_sel <- rbind(SUM_sel,SUM_sel5)
SUM_sel <- rbind(SUM_sel,SUM_sel6)
SUM_sel <- rbind(SUM_sel,SUM_sel7)
SUM_sel <- rbind(SUM_sel,SUM_sel8)
SUM_sel <- rbind(SUM_sel,SUM_sel9)
SUM_sel <- rbind(SUM_sel,SUM_sel10)

SUM_sel  <- SUM_table[which(rownames(SUM_table) %in% markergene),]
loc = match(markergene,rownames(SUM_table))
SUM_sel2 <- SUM_table[loc,]

write.table(SUM_sel,"all-11celltype-sample-100cells_genescoreMatrix.txt",sep = "\t",quote = F)
write.table(SUM_sel,"all-11celltype-sample-100cells_GeneIntegrationMatrix.txt",sep = "\t",quote = F)

SUM_sel=read.table("all-11celltype-sample-100cells_genescoreMatrix.txt", header = T)
SUM_sel=read.table("all-11celltype-sample-100cells_genescoreMatrixtxt",header = T)

head(SUM_sel2[,1:5])
SUM_sel2 =na.omit(SUM_sel2)

write.table(SUM_sel2,"all-11celltype-sample-100cells_genescoreMatrix-removeNA.txt",sep = "\t",quote = F)
write.table(SUM_sel2,"all-11celltype-sample-100cells_GeneIntegrationMatrix-removeNA.txt",sep = "\t",quote = F)

SUM_sel3=read.table("all-11celltype-sample-100cells_genescoreMatrix-removeNA.txt")
SUM_sel3=read.table("all-11celltype-sample-100cells_GeneIntegrationMatrix-removeNA.txt",header = T)

head(SUM_sel3[,1:5])
SUM_sel3 <- as.matrix(SUM_sel3)

SUM_sel4 =log10(SUM_sel3 +1)
SUM_sel4=t(scale(t(SUM_sel4)))
SUM_sel4 =na.omit(SUM_sel4)
library(preprocessCore)
SUM_sel5= normalize.quantiles(SUM_sel4)
rownames(SUM_sel5) =rownames(SUM_sel4)
colnames(SUM_sel5) =colnames(SUM_sel4)
head(SUM_sel5[,1:5])

bk=unique(c(seq(-2.5,2.5, length=100)))
pdf("all-11celltype-Genescore-random-100-select.pdf")
ped =pheatmap(SUM_sel5, cluster_rows=FALSE, 
              clustering_method = "ward.D2",
              show_rownames=TRUE, treeheight_row = 0,
#              breaks=bk,
              color = paletteContinuous(set = "blueYellow", n = 256, reverse = FALSE),
#              color = colorRampPalette((rev(brewer.pal(n = 11, name="RdBu"))))(100), 
              cluster_cols=FALSE,cutree_rows =0,main = "",
              show_colnames = FALSE,fontsize_row = 7.5)
ped
dev.off()

#Endohelial cell 
ArchR.EC <- ArchR[which(ArchR@cellColData$Clusters == "C14" |ArchR@cellColData$Clusters == "C15" |ArchR@cellColData$Clusters == "C16" 
                         |ArchR@cellColData$Clusters == "C17" |ArchR@cellColData$Clusters == "C18" |ArchR@cellColData$Clusters == "C19" |
                           ArchR@cellColData$Clusters == "C20" |ArchR@cellColData$Clusters == "C21" ),]

ArchR.EC <- addIterativeLSI( ArchRProj = ArchR.EC,useMatrix = "TileMatrix",name = "IterativeLSI", iterations = 2, 
  clusterParams = list(resolution = c(1),sampleCells = 10000, n.start = 10), 
  varFeatures = 30000, dimsToUse = 1:15)

####add cluster####
ArchR.EC <- addClusters(input = ArchR.EC, reducedDims = "IterativeLSI",name = "Clusters_3", force = TRUE, resolution = 1, dimsToUse = 1:15)
table(ArchR.EC$Clusters_3)
cM <- confusionMatrix(paste0(ArchR.EC$Clusters_3), paste0(ArchR.EC$Sample))
cM
ArchR.EC <- addUMAP(ArchRProj = ArchR.EC, reducedDims = "IterativeLSI", nNeighbors = 10, minDist = 0.3, metric = "cosine", force = TRUE)

#####define genotype####
meta<-matrix(data=NA,nrow=28055,ncol=2);
colnames(meta)=c("Sample","Genotype");
meta[,1] = ArchR.EC$Sample
for ( i in 1:nrow(meta))
{
  if (meta[i,1]  == "AR3" |meta[i,1]  == "AR7" |meta[i,1]  == "AR14" )
    meta[i,2] = "AR"
  else if (meta[i,1]  == "NAR3"|meta[i,1]  == "NAR7" |meta[i,1]  == "NAR14")
    meta[i,2]  = "NAR"
  else
    meta[i,2]  = "P"
}
head(meta)
meta = data.frame(meta)
ArchR.EC$Genotype  = meta$Genotype
table(ArchR.EC$Genotype)

####Call peak#####
table(ArchR.EC@cellColData$Clusters_3)
pathToMacs2 <- "/public/home/Chenzh275//miniconda3/bin/macs2"
ArchR.EC <- addReproduciblePeakSet(ArchRProj = ArchR.EC, groupBy = "Clusters_3", pathToMacs2 = pathToMacs2, cutOff = 0.01)
getPeakSet(ArchR.EC)

projHemeTmp <- addReproduciblePeakSet(ArchRProj = ArchR.EC, groupBy = "Clusters_3", peakMethod = "Tiles", method = "p", cutOff = 0.01)
getPeakSet(projHemeTmp)

####add peak matrix####
ArchR.EC <- addPeakMatrix(ArchR.EC)
getAvailableMatrices(ArchR.EC)

#define subtype
table(ArchR.EC$Clusters_3)
meta<-matrix(data=NA,nrow=28055,ncol=2);
colnames(meta)=c("cluster","subtype");
meta[,1] = ArchR.EC$Clusters_3
cl <- meta[,1]
for ( i in 1:length(cl))
{
  if (cl[i] == "C4" |cl[i] == "C5" |cl[i] == "C6" |cl[i] == "C7" |cl[i] == "C10"  )
    meta[i,2] = "EC1"
  else if (cl[i] == "C2"|cl[i] == "C3" )
    meta[i,2]  = "EC2"
  else if (cl[i] == "C8" )
    meta[i,2]  = "EC3"
  else if (cl[i] == "C1" | cl[i] == "C12" |cl[i] == "C13")
    meta[i,2]  = "EC4"
  else
    meta[i,2]  = "EC5"
}
ArchR.EC$Subtype  = meta$subtype
table(ArchR.EC$Subtype)

####marker motif 
ArchR.EC <- addMotifAnnotations(ArchRProj = ArchR.EC, motifSet = "homer", name = "Motif", force = TRUE)
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,  ArchRProj = ArchR.EC,
                                   peakAnnotation = "Motif",  cutOff = "FDR <= 0.5 & Log2FC >= 0.5")
enrichMotifs

saveArchRProject(ArchRProj = ArchR.EC2, outputDirectory = "Save-Proj2-EC5-subtype", load = FALSE)

EC4proj<-loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
table(EC4proj$Subtype)

#####re draw umap 
head(df)
df = data.frame(EC4proj@embeddings$UMAP$df, EC4proj$Subtype)
Subtype = EC4proj$Subtype

##subtype umap
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = Subtype)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#cbb01c","#663b0e","#cb761c","#331e07","#f4d2ae")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Subtype")


markersGS <- getMarkerFeatures(ArchRProj = EC4proj, useMatrix = "GeneScoreMatrix", groupBy = "EC4proj$Subtype",
                               bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersGS
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
markerList$EC2
#write.table(markerList,"EC-subtype-markerlist-FDR001-log2FC05.txt",quote = F,sep = "\t")

heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", transpose = F, 
                               labelMarkers = EC_marker, clusterCols = T,returnMat = T)
library(pheatmap)
pheatmap(heatmapGS,cluster_cols = F,cluster_rows = F,
         color = paletteContinuous(set = "blueYellow", n = 256, reverse = FALSE),
         #cellwidth = 15, cellheight = 12,
         show_rownames=F,show_colnames=T,cutree_rows =5)


##markerpeak
markersPeaks <- getMarkerFeatures(ArchRProj = EC4proj, useMatrix = "PeakMatrix", 
                                  groupBy = "EC4proj$Subtype", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.5 & Log2FC >= 0.5")
write.table(markerList,"markerList-Peak-Marker-EC-Subtype-FDR05log2FC05.txt",quote = F,sep = "\t")
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,  ArchRProj = EC4proj,
                                   peakAnnotation = "Motif", cutOff = "FDR <= 0.5 & Log2FC >= 0.5")
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, transpose = F,returnMat = T)
pdf("EC-Motifs-Enriched-Marker-Heatmap-FDR05-log2FC05-n.pdf")
pheatmap(heatmapEM,cluster_cols = F,cluster_rows = F,
         color = paletteContinuous(set = "comet", n = 256, reverse = FALSE),
         #cellwidth = 15, cellheight = 12,
         show_rownames=T,show_colnames=T)
dev.off()

##GO term 

library(clusterProfiler)
library("org.Mm.eg.db")
library("ChIPseeker")
library(AnnotationDbi)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

Allgene <-read.table("markersList-EC-gene-subtype-cluster_3-FC05.txt",header = T)


EC1 <-Allgene[which(Allgene$group_name =="EC1"),]
EC2 <-Allgene[which(Allgene$group_name =="EC2"),]
EC3 <-Allgene[which(Allgene$group_name =="EC3"),]
EC4 <-Allgene[which(Allgene$group_name =="EC4"),]
EC5 <-Allgene[which(Allgene$group_name =="EC5"),]

EC1_gene  = EC1$name
EC1_geneID <- bitr(EC1_gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
EC1_ego <- enrichGO(EC1_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
EC1_ego <- as.data.frame(EC1_ego)
EC1_ego$Description
write.table(EC1_ego,"EC1_ego.txt",sep="\t",quote = FALSE, row.names=F)
dotplot(EC1_ego, showCategory = 30 )

EC2_gene  = EC2$name
EC2_geneID <- bitr(EC2_gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
EC2_ego <- enrichGO(EC2_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
EC2_ego <- as.data.frame(EC2_ego)
EC2_ego$Description
write.table(EC2_ego,"EC2_ego.txt",sep="\t",quote = FALSE, row.names=F)
dotplot(EC2_ego, showCategory = 30 )

EC3_gene  = EC3$name
EC3_geneID <- bitr(EC3_gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
EC3_ego <- enrichGO(EC3_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
EC3_ego <- as.data.frame(EC3_ego)
EC3_ego$Description
write.table(EC3_ego,"EC3_ego.txt",sep="\t",quote = FALSE, row.names=F)
dotplot(EC3_ego, showCategory = 30 )

EC4_gene  = EC4$name
EC4_geneID <- bitr(EC4_gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
EC4_ego <- enrichGO(EC4_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
EC4_ego <- as.data.frame(EC4_ego)
EC4_ego$Description
write.table(EC4_ego,"EC4_ego.txt",sep="\t",quote = FALSE, row.names=F)
dotplot(EC4_ego, showCategory = 30 )

EC5_gene  = EC5$name
EC5_geneID <- bitr(EC5_gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
EC5_ego <- enrichGO(EC5_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
EC5_ego <- as.data.frame(EC5_ego)
EC5_ego$Description
write.table(EC5_ego,"EC5_ego.txt",sep="\t",quote = FALSE, row.names=F)
dotplot(EC5_ego, showCategory = 30 )

#peak link to gene


##MP
ArchR.IC <- ArchR[which(ArchR@cellColData$Clusters == "C12" |ArchR@cellColData$Clusters == "C13"),]
ArchR.IC <- addIterativeLSI( ArchRProj = ArchR.IC,useMatrix = "TileMatrix",name = "IterativeLSI", iterations = 2, 
                             clusterParams = list(resolution = c(1),sampleCells = 10000, n.start = 10), 
                             varFeatures = 30000, dimsToUse = 1:15)

####add cluster####
ArchR.IC <- addClusters(input = ArchR.IC, reducedDims = "IterativeLSI",name = "Clusters_2", force = TRUE, resolution = 1, dimsToUse = 1:15)
table(ArchR.IC$Clusters_2)
cM <- confusionMatrix(paste0(ArchR.IC$Clusters_2), paste0(ArchR.IC$Sample))
cM

ArchR.IC <- addUMAP(ArchRProj = ArchR.IC, reducedDims = "IterativeLSI", nNeighbors = 10, minDist = 0.5, metric = "cosine", force = TRUE)

#####define genotype####
meta<-matrix(data=NA,nrow=4114,ncol=2);
colnames(meta)=c("Sample","Genotype");
meta[,1] = ArchR.IC$Sample
for ( i in 1:nrow(meta))
{
  if (meta[i,1]  == "AR3" |meta[i,1]  == "AR7" |meta[i,1]  == "AR14" )
    meta[i,2] = "AR"
  else if (meta[i,1]  == "NAR3"|meta[i,1]  == "NAR7" |meta[i,1]  == "NAR14")
    meta[i,2]  = "NAR"
  else
    meta[i,2]  = "P"
}
head(meta)
meta = data.frame(meta)
ArchR.IC$Genotype  = meta$Genotype
table(ArchR.IC$Genotype)

#####define subtype####
table(ArchR.IC$Clusters_2)
meta<-matrix(data=NA,nrow=4114,ncol=2);
colnames(meta)=c("cluster","subtype");
meta[,1] = ArchR.IC$Clusters_2
cl <- meta[,1]
for ( i in 1:length(cl))
{
  if (cl[i] == "C1" |cl[i] == "C4" |cl[i] == "C5" | cl[i] == "C6" )
    meta[i,2] = "IC1"
  else if (cl[i] == "C2" )
    meta[i,2]  = "IC2"
  else if (cl[i] == "C3" )
    meta[i,2]  = "IC3"
  else if (cl[i] == "C7")
    meta[i,2]  = "IC4"
  else if (cl[i] == "C9")
    meta[i,2]  = "IC5"
  else
    meta[i,2]  = "IC6"
}
head(meta)
meta = data.frame(meta)
ArchR.IC$Subtype  = meta$subtype
table(ArchR.IC$Subtype)

####label timepoint 
meta<-matrix(data=NA,nrow=4114,ncol=2);
colnames(meta)=c("Sample","timepoint");
meta[,1] = ArchR.IC$Sample
cl <- meta[,1]
for ( i in 1:nrow(meta))
{
  if (meta[i,1]  == "AR3" |meta[i,1]  == "NAR3" |meta[i,1]  == "P3" )
    meta[i,2] = "Day3"
  else if (meta[i,1]  == "AR7"|meta[i,1]  == "NAR7" |meta[i,1]  == "P7")
    meta[i,2]  = "Day7"
  else
    meta[i,2]  = "Day14"
}
head(meta)
meta = data.frame(meta)
ArchR.IC$timepoint  = meta$timepoint
table(ArchR.IC$timepoint)

#####extract macrophage
MP =ArchR.IC[which(ArchR.IC$Subtype =="IC1" |ArchR.IC$Subtype =="IC2" |ArchR.IC$Subtype =="IC3" ),]
MP3 <- addIterativeLSI(ArchRProj = MP3,useMatrix = "TileMatrix",name = "IterativeLSI_MP3", iterations = 2, 
                       clusterParams = list(resolution = c(0.8),sampleCells = 10000, n.start = 10), 
                       varFeatures = 50000, dimsToUse = 1:25,force = TRUE)
###MP3
####re add cluster####
MP3 <- addClusters(input = MP3, reducedDims = "IterativeLSI_MP3",name = "MP3_Clusters", force = TRUE, 
                   resolution = 0.8, dimsToUse = 1:25)
table(MP3$MP3_Clusters)
cM <- confusionMatrix(paste0(MP3$MP3_Clusters), paste0(MP3$Sample))
cM

MP3 <- addUMAP(ArchRProj = MP3, reducedDims = "IterativeLSI_MP3", nNeighbors = 10, minDist = 0.5, metric = "cosine", force = TRUE)

###re define subset
table(MP3$MP3_Clusters)
meta<-matrix(data=NA,nrow=3367,ncol=2);
colnames(meta)=c("cluster","subtype");
meta[,1] = MP3$MP3_Clusters
cl <- meta[,1]
for ( i in 1:length(cl))
{
  if (cl[i] == "C1" )
    meta[i,2] = "MP1"
  else if (cl[i] == "C2" |cl[i] == "C3")
    meta[i,2]  = "MP2"
  else if (cl[i] == "C4" )
    meta[i,2]  = "MP3"
  else 
    meta[i,2]  = "MP4"
}
head(meta)
meta = data.frame(meta)
MP3$Subtype  = meta$subtype
table(MP3$Subtype)

####Making Pseudo-bulk Replicates######
MP3 <- addGroupCoverages(ArchRProj = MP3, groupBy = "MP3$Subtype")

####Call peak#####
pathToMacs2 <- "/public/home/Chenzh275/miniconda3/bin/macs2"
MP3 <- addReproduciblePeakSet(ArchRProj = MP3, groupBy = "MP3$Subtype", pathToMacs2 = pathToMacs2, cutOff = 0.01)
getPeakSet(MP3)

MP3 <- addPeakMatrix(MP3)
getAvailableMatrices(MP3)

##markerPeak
table(MP3$Subtype)
markersPeaks <- getMarkerFeatures(ArchRProj = MP3, useMatrix = "PeakMatrix", 
                                  groupBy = "MP3$Subtype", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks,  cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

##motif
MP3 <- addMotifAnnotations(ArchRProj = MP3, motifSet = "homer", name = "MP_Motif",force = TRUE)

#Motif Enrichment in Marker Peaks
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,  ArchRProj = MP3,
                                   peakAnnotation = "MP_Motif", cutOff = "FDR <= 0.5 & Log2FC >= 0.5")

#saveArchRProject(ArchRProj = MP3, outputDirectory = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype", load = FALSE)
MP3proj =loadArchRProject("/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")


#####re draw umap 
df = data.frame(MP3@embeddings$UMAP$df, MP3$Genotype)
df = data.frame(MP3@embeddings$UMAP$df, MP3$timepoint)
df = data.frame(MP3@embeddings$UMAP$df, MP3$Subtype)
Timepoint = MP3$timepoint
Subtype = MP3$Subtype
Sample = MP3$Genotype
###timepoint
ggplot(df, aes(x = IterativeLSI_MP3.UMAP_Dimension_1 , y = IterativeLSI_MP3.UMAP_Dimension_2, color = Timepoint)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#FFA500","#0073C2FF", "#20854EFF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Timepoint")

##subtype umap
ggplot(df, aes(x = IterativeLSI_MP3.UMAP_Dimension_1 , y = IterativeLSI_MP3.UMAP_Dimension_2, color = Subtype)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#615fa4","#7688b1","#9f76b1","#37365e")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Subtype")

###sample
ggplot(df, aes(x = IterativeLSI_MP3.UMAP_Dimension_1 , y = IterativeLSI_MP3.UMAP_Dimension_2,color = Sample)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#E64B35B2","#00468BB2","#42B540B2" )) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Genotype")

##marker gene vlnplot
markerGenes <- c("Timd4", "Lyve1","Folr2","Mrc1","Cd163","F13a1","Igf1" ,  ###resident macrophages (TIMD+ CCR2-)
                "Il1b","Ccr2","H2-Aa","H2-Eb1","Cd74", ###resident macrophages (TIMD- CCR2+)
                "Cd14","Mmp12","Cx3cr1",###resident macrophages (TIMD- CCR2-)
                "Ly6c1","Isg" 
                )


p2 <- plotGroups(          
  ArchRProj = MP3proj, 
  groupBy = "MP3proj$Subtype", 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes ,
  plotAs = "violin",
  pal =c("#0C321D","#9DDFD3","#FFC93C","#B0C793"),
  alpha = 0.4,
  addBoxPlot = TRUE
)

p3 <- lapply(p2, function(x){x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 3),p3))
plotPDF(plotList = p2, name = "Plot-Voilin-type-Marker-Genes-MP-W-Imputation-4-Subtype.pdf", ArchRProj = MP3proj, addDOC = FALSE, width = 5, height = 5)


#marker gene UMAP
markerGenes <- c("Timd4", "Lyve1","Folr2","Mrc1","Cd163","F13a1","Igf1", ###resident macrophages (TIMD+ CCR2-)
                "Il1b","Ccr2","H2-Aa","H2-Eb1","Cd74",  ###resident macrophages (TIMD- CCR2+)
                "Cd14","Mmp12","Cx3cr1", ###resident macrophages (TIMD- CCR2-)
                "Ly6c1","Ly6c2" ##recruited monocyte
                )
MP3proj =loadArchRProject("/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")

MP3proj <- addImputeWeights(MP3proj)
getAvailableMatrices(MP3proj)

p <- plotEmbedding(
    ArchRProj = MP3proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(MP3proj)
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
plotPDF(plotList = p, name = "Plot-umap-MarkerGenes-sepecficGene-FDR001logFC05-ATAC", ArchRProj = MP3proj, addDOC = FALSE, width = 5, height = 5)

##protein folding TF link to gene promoter region
MP3proj =loadArchRProject("/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")

FB3proj<-loadArchRProject(path = "/md01/nieyg/scATAC-ArchR/Save-fibroblast2")

EC4proj<-loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")

SMC1proj<-loadArchRProject(path="/md01/nieyg/scATAC-ArchR/Save-SMC-merge12")


#查看有gene link 的peak是否有sp1 klf转录因子的结合位�?
motifPositions <- getPositions(MP3proj,name = "MP_Motif") ##peakAnnotate name of project
motifPositions <- getPositions(FB3proj,name = "FB_Motif") 
motifPositions <- getPositions(EC4proj,name = "Motif") 
motifPositions <- getPositions(SMC1proj,name = "Motif") 
motifPositions

##MP3 "Sp1","KLF14","KLF3","KLF5","Klf9","Klf4" 
##FB3 "Sp1","KLF14","KLF3","Klf9","KLF5","Klf4"
##EC4 "Sp1","KLF14","KLF3","Klf9","KLF5","Klf4"
##SMC "Sp1","KLF14","Klf9","KLF3","KLF5","Klf4"
motifs <- c("Sp1","KLF14","KLF3","Klf9","KLF5","Klf4")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

MP3_peak <- readPeakFile("/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/homer/markerPeak_FDR01log2FC05/newMP3-markerlist_homer.bed")
FB3_peak <- readPeakFile("/md01/Chenzh275/ori/Data/scATAC/ArchR/matrix/FB_peak/FB3-markerList-peak-homer.bed")
EC4_peak <- readPeakFile("/md01/Chenzh275/ori/Data/scATAC/ArchR-Endothelial_Cell/homer/EC4-markerList-peak.bed")
SMC1_peak <- readPeakFile("/md01/Chenzh275/ori/Data/scATAC/ArchR/matrix/SMC_peak/SMC1-markerList-peak-FDR01log2FC05-homer.bed")

setwd("/md01/Chenzh275/ori/Data/scATAC/ProteinFolding-TF/MP3/TF_binding_site")
setwd("/md01/Chenzh275/ori/Data/scATAC/ProteinFolding-TF/FB3/TF_binding_site")
setwd("/md01/Chenzh275/ori/Data/scATAC/ProteinFolding-TF/EC4/TF_binding_site")
setwd("/md01/Chenzh275/ori/Data/scATAC/ProteinFolding-TF/SMC1/TF_binding_site")

Sp1_total <- motifPositions$Sp1.Zf_266 ##Sp1 所有的peak
Sp1_overlaps <- findOverlaps(Sp1_total,MP3_peak)  #FB3_peak EC4_peak SMC1_peak
#write.table(Sp1_total,"SP1-binding-site.txt",sep ="\t",quote =F)

KLF14_total <- motifPositions$KLF14.Zf_145 ##KLF14 所有的peak
KLF14_overlaps <- findOverlaps(KLF14_total,MP3_peak)  
#write.table(KLF14_total,"KLF14-binding-site.txt",sep ="\t",quote =F)

KLF3_total <- motifPositions$KLF3.Zf_146 ##KLF14 所有的peak
KLF3_overlaps <- findOverlaps(KLF3_total,MP3_peak)  
#write.table(KLF3_total,"KLF3-binding-site.txt",sep ="\t",quote =F)

Klf9_total <- motifPositions$Klf9.Zf_149 ##KLF14 所有的peak
Klf9_overlaps <- findOverlaps(Klf9_total,MP3_peak)  
#write.table(Klf9_total,"Klf9-binding-site.txt",sep ="\t",quote =F)

KLF5_total <- motifPositions$KLF5.Zf_148 ##KLF14 所有的peak
KLF5_overlaps <- findOverlaps(KLF5_total,MP3_peak)  
#write.table(KLF5_total,"KLF5-binding-site.txt",sep ="\t",quote =F)

Klf4_total <- motifPositions$Klf4.Zf_147 ##KLF14 所有的peak
Klf4_overlaps <- findOverlaps(Klf4_total,MP3_peak)  
#write.table(Klf4_total,"Klf4-binding-site.txt",sep ="\t",quote =F)

##
peakset <- (MP3proj@peakSet)
peakset <- (FB3proj@peakSet)
peakset <- (EC4proj@peakSet)
peakset <- (SMC1proj@peakSet)

Hspa8 <- peakset[which(peakset$nearestGene == "Hspa8")]
Hspa1a <- peakset[which(peakset$nearestGene == "Hspa1a")]
Hspa1b <- peakset[which(peakset$nearestGene == "Hspa1b")]
Hsp90aa1 <- peakset[which(peakset$nearestGene == "Hsp90aa1")]

Hspa8_promoter <- Hspa8[which(Hspa8$peakType =="Promoter")]
Hspa1a_promoter <- Hspa1a[which(Hspa1a$peakType =="Promoter")]
Hspa1b_promoter <- Hspa1b[which(Hspa1b$peakType =="Promoter")]
Hsp90aa1_promoter <- Hsp90aa1[which(Hsp90aa1$peakType =="Promoter")]
 
Hspa8_promoter <- Hspa8_promoter[which(Hspa8_promoter$distToTSS  <= 1000 )] #FB3 EC4 SMC1
Hspa1a_promoter <- Hspa1a_promoter[which(Hspa1a_promoter$distToTSS  <= 1000)]
Hspa1b_promoter <- Hspa1b_promoter[which(Hspa1b_promoter$distToTSS  <= 1000)]
Hsp90aa1_promoter <- Hsp90aa1_promoter[which(Hsp90aa1_promoter$distToTSS  <= 1000)]

promoter <- union(Hspa8_promoter,Hspa1a_promoter)
promoter <- union(promoter,Hspa1b_promoter)
promoter <- union(promoter,Hsp90aa1_promoter)

Sp1_total <- motifPositions$Sp1.Zf_266 ##Sp1 所有的peak
Sub_Sp1 <- subsetByOverlaps(Sp1_total,Hsp90aa1_promoter)
Sub_Sp1
KLF14_total <- motifPositions$KLF14.Zf_145 ##KLF14 所有的peak
Sub_KLF14 <- subsetByOverlaps(KLF14_total,Hsp90aa1_promoter)
Sub_KLF14
KLF3_total <- motifPositions$KLF3.Zf_146 ##KLF14 所有的peak
Sub_KLF3 <- subsetByOverlaps(KLF3_total,Hsp90aa1_promoter)
Sub_KLF3
Klf9_total <- motifPositions$Klf9.Zf_149 ##KLF14 所有的peak
Sub_Klf9 <- subsetByOverlaps(Klf9_total,Hsp90aa1_promoter)
Sub_Klf9
KLF5_total <- motifPositions$KLF5.Zf_148 ##KLF14 所有的peak
Sub_KLF5 <- subsetByOverlaps(KLF5_total,Hsp90aa1_promoter)
Sub_KLF5
Klf4_total <- motifPositions$Klf4.Zf_147 ##KLF14 所有的peak
Sub_Klf4 <- subsetByOverlaps(Klf4_total,Hsp90aa1_promoter)
Sub_Klf4



