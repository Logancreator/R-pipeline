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
) + geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = 3.5, lty = "dashed")

p
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme2, addDOC = FALSE)





projHeme2 <- addImputeWeights(projHeme2)
markerGenes=c("GFAP","AGT","SLC1A2","AQP4","ETNPPL", #astrocyte

"SLC17A6","GRIN1",#glutamatergic neurons

"GAD1","TH","DDC","SLC18A2","SLC32A1","SLC17A6","DPP10",#GABA

"CACNA1E","CACNA2D1", #neuron

"MEGF11" ,"VCAN" ,"LHFPL3", #Oligodendrocyte precursor cell

"ELF1",#ependymocyte

"RAX",#tanycyte

"MBP","CNP","ERMN","ABCA2","SOX10","ST18" #oligodendrocyte

,"C1QA","C1QB","C1QC","RUNX1","SLC2A5","PTPRC","CX3CR1",#Microglial

"GATA2")  #endothelial

markerGenes=toupper(c("Slc6a13",
  "Acta2",
  "Tbx18",
  "Igfbp2",
  "Rax",

  "Kcnj8",
  "Cspg4",

  "Mal",
  "Apod",



  "Fyn",

  "Syt1",
  "Snap25",
  "Slc17a6",
  "Gad1",
  "Tmem119",
  "P2ry12",
  "Csf1r",
  "C1qa",
  "Cx3cr1",


  "Adgrf5",

  "Cd44",
  "Gja1",
  "Gfap",
  "Aqp4",
  "Agt",
  "col25a1"))
p=plotEmbedding(
            ArchRProj = projHeme2,
            colorBy = "GeneScoreMatrix",
            name = markerGenes,
            embedding = "UMAP",
            imputeWeights = getImputeWeights(projHeme2))
plotPDF(p)


projHeme2$celltype=projHeme2$Clusters
projHeme2$celltype[projHeme2$celltype%in%c("C1","C2","C3","C4","C5","C6")]="astrocyte" ##
projHeme2$celltype[projHeme2$celltype%in%c("C21","C22")]="neuron" ##
projHeme2$celltype[projHeme2$celltype%in%c("C10","C11")]="oligodendrocyte precursor cell" ##
projHeme2$celltype[projHeme2$celltype%in%c("C8","C9")]="oligodendrocyte"   ##
projHeme2$celltype[projHeme2$celltype%in%c("C15","C16","C19")]="microglial/ependymocyte" ##
projHeme2$celltype[projHeme2$celltype%in%c("C20","C7","C17","C18","C19")]="endothelial/tanycyte"  ##
projHeme2$celltype[projHeme2$celltype%in%c("C12","C13","C14")]="unknown"  ##
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP",baseSize=5)
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "celltype", embedding = "UMAP",baseSize=5)
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-celltype.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 3, height = 3)



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