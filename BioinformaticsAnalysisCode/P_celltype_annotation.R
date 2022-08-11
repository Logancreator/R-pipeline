library(ArchR)
setwd("/public/home/changjianye/project/scatac/soul_P")
projHeme2 <- loadArchRProject("Save-ProjHeme_P-filterDoublets_filterFrags/")
##cell filter
idxPass <- which(log10(projHeme2$nFrags) <= 4.25 & projHeme2$TSSEnrichment >=2 )
cellsPass <- projHeme2$cellNames[idxPass]
projHeme2=projHeme2[cellsPass, ]


##TSS vs Frag
df <- getCellColData(projHeme2, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 1, lty = "dashed") + geom_vline(xintercept = 3.5, lty = "dashed")

p
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme2, addDOC = FALSE)





projHeme2 <- addImputeWeights(projHeme2)
# markerGenes=c("GFAP","AGT","SLC1A2","AQP4","ETNPPL", #astrocyte

# "SLC17A6","GRIN1",#glutamatergic neurons

# "GAD1","TH","DDC","SLC18A2","SLC32A1","SLC17A6","DPP10",#GABA

# "CACNA1E","CACNA2D1", #neuron

# "MEGF11" ,"VCAN" ,"LHFPL3", #Oligodendrocyte precursor cell

# "ELF1",#ependymocyte

# "RAX",#tanycyte

# "MBP","CNP","ERMN","ABCA2","SOX10","ST18" #oligodendrocyte

# ,"C1QA","C1QB","C1QC","RUNX1","SLC2A5","PTPRC","CX3CR1",#Microglial

# "GATA2")  #endothelial
p=plotEmbedding(
            ArchRProj = projHeme2,
            colorBy = "GeneScoreMatrix",
            name = markerGenes,
            embedding = "UMAPHarmony",
            imputeWeights = getImputeWeights(projHeme2))
plotPDF(p)



projHeme2$celltype=projHeme2$Clusters
projHeme2$celltype[projHeme2$celltype%in%c("C3")]="Endotheliocyte or RBC" ## yes
projHeme2$celltype[projHeme2$celltype%in%c("C8")]="Somatotrophs" ##
projHeme2$celltype[projHeme2$celltype%in%c("C6","C20")]="Thyrotrophs" ##
projHeme2$celltype[projHeme2$celltype%in%c("C23")]="folliculo-stellate cells"   ##
projHeme2$celltype[projHeme2$celltype%in%c("C16","C17","C18","C19")]="Gonadotrophs" ##yes
projHeme2$celltype[projHeme2$celltype%in%c("C12")]="Corticotrophs"  ##
projHeme2$celltype[projHeme2$celltype%in%c("C7","C9","C10")]="Lactotrophs"  ##
projHeme2$celltype[projHeme2$celltype%in%c("C25")]="Red_blood_cell"  ##
projHeme2$celltype[projHeme2$celltype%in%c("C11")]="progenitor or stem cell"
projHeme2$celltype[projHeme2$celltype%in%c("C1")]="posterior pituicyte pituitary"
projHeme2$celltype[projHeme2$celltype%in%c("C2","C4","C5","C13","C14","C15","C21","C22","C24")]="Unknown"

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony",baseSize=5)
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "celltype", embedding = "UMAPHarmony",baseSize=5)
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-celltype.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)



markerGenes=c("SLC17A6", "GRIN1",  ##HT-EX-1,HT-EX-2
    "GAD1", "GAD2" ,##GABA
"TH",  "DDC", "SLC18A2", "LMX1A", ##Neuron
"CACNA1E", "CACNA2D1","KCNQ3", "RYR3","ANO4" ,
"CACNB2","MEIS2",
"P2RY12","CALCR","RUNX1")
p=plotEmbedding(
            ArchRProj = projHeme2,
            colorBy = "GeneScoreMatrix",
            name = markerGenes,
            embedding = "UMAP",
            imputeWeights = getImputeWeights(projHeme2))
plotPDF(p)





df = data.frame(projHeme2@embeddings$UMAP$df, projHeme2$celltype)
celltype = projHeme2$celltype

library(ggsci)
pdf("/public/home/changjianye/project/scatac/soul_ME/Save-ProjHeme_ME-filterDoublets_filterFrags/Plots/UMAP-all-NCM-cell-type-split-IC-size0.4-alpha0.5.pdf")
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = projHeme2.celltype)) +
  geom_point(size= 0.3) + 
  theme_bw() + 
  scale_color_jco()+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Celltype")
dev.off()



markersGS <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)


#Bubble plot for ArchR data script
geneintegration= getMatrixFromProject(projHeme2, useMatrix="GeneScoreMatrix",threads=15)
geneintedata=assay(geneintegration)
genedata=rowData(geneintegration)

# geneintedata_matrix <- matrix(data = 0, nrow=dim(geneintedata)[1], ncol = dim(geneintedata)[2])
# for (i in seq_len(nrow(geneintedata))){
#     for (j in seq_len(ncol(geneintedata))){
#         geneintedata_matrix[i][j] <- geneintedata[i][j]
#     }}


geneinteinfo=cbind(genedata, geneintedata)
geneinteinfo2=as.data.frame(t(as.data.frame(geneinteinfo[,7:length(colnames(geneinteinfo))])))
colnames(geneinteinfo2)=genedata$name
rownames(geneinteinfo2)=colnames(geneinteinfo)[7:length(colnames(geneinteinfo))]


for(i in rownames(geneinteinfo2)){
  geneinteinfo2[i,"celltype"]=projHeme2$Clusters[which(projHeme2$cellNames==i)]
  geneinteinfo2[i,"Samples"]=projHeme2$Sample[which(projHeme2$cellNames==i)]
}
##gene list 
gene_list=c("GFAP","AGT","SLC1A2","AQP4","ETNPPL", #astrocyte

"SLC17A6","GRIN1",#glutamatergic neurons

"GAD1","TH","DDC","SLC18A2","SLC32A1","SLC17A6","DPP10",#GABA

"CACNA1E","CACNA2D1", #neuron

"MEGF11" ,"VCAN" ,"LHFPL3", #Oligodendrocyte precursor cell

"ELF1",#ependymocyte

"RAX",#tanycyte

"MBP","CNP","ERMN","ABCA2","SOX10","ST18" #oligodendrocyte

,"C1QA","C1QB","C1QC","RUNX1","SLC2A5","PTPRC","CX3CR1",#Microglial

"GATA2")  #endothelial
##compute pct_exp and avg_exp
bubble_plot_info=data.frame()
for(i in gene_list){
  for(k in 1:length(unique(geneinteinfo2$Clusters))){
    a=nrow(bubble_plot_info)
    l=unique(geneinteinfo2$Clusters)[k]
    bubble_plot_info[a+1,"gene_name"]=i
    bubble_plot_info[a+1,"Clusters"]=l
    eval(parse(text=(paste("bubble_plot_info[",a,"+1,'pct_exp']=(length(geneinteinfo2[(geneinteinfo2$",i,">1 & geneinteinfo2$Clusters=='",l,"'),'",i,"'])/nrow(geneinteinfo2[geneinteinfo2$Clusters=='",l,"',]))*10", sep=""))))
    eval(parse(text=(paste("bubble_plot_info[",a,"+1,'avg_exp']=mean(geneinteinfo2[geneinteinfo2$",i,">0 & geneinteinfo2$Clusters=='",l,"','",i,"'])", sep=""))))
  }
}
##plot
pdf("/public/home/changjianye/project/scatac/soul_ME/Save-ProjHeme_ME-filterDoublets_filterFrags/Plots/Bubble-plot.pdf", width=7, height=17)

ggplot(data = bubble_plot_info, mapping = aes_string(x = 'gene_name', y = 'Clusters')) +
  geom_point(mapping = aes_string(size = 'pct_exp', color = "avg_exp")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(size = guide_legend(title = 'Percent Expressed')) +
  labs(
    x = 'gene_name',
    y = 'Clusters'
  )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 0))+ scale_colour_gradient(low="white",high="blue")+
  scale_size_area()+
  coord_flip()
dev.off()



##markerGene
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)")
)
markersGS
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList
write.table(as.data.frame(markerList$C1),"markerList_C1.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C2),"markerList_C2.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C3),"markerList_C3.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C4),"markerList_C4.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C5),"markerList_C5.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C6),"markerList_C6.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C7),"markerList_C7.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C8),"markerList_C8.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C9),"markerList_C9.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C10),"markerList_C10.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C11),"markerList_C11.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C12),"markerList_C12.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C13),"markerList_C13.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C14),"markerList_C14.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C15),"markerList_C15.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C16),"markerList_C16.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C17),"markerList_C17.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C18),"markerList_C18.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C19),"markerList_C19.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C20),"markerList_C20.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C21),"markerList_C21.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C22),"markerList_C22.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C23),"markerList_C23.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C24),"markerList_C24.txt",quote=F,sep="\t",row.names=F)
write.table(as.data.frame(markerList$C25),"markerList_C25.txt",quote=F,sep="\t",row.names=F)
