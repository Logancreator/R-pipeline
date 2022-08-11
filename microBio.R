##å¾®ç”Ÿä¿?
#GO
rm(list=ls())
setwd("C:\\Users\\Jeff\\OneDrive\\×ÀÃæ\\quail_ATAC_info\\RNA\\SD28_LD7\\")
go=read.csv("SD28_LD7downGO.csv")
go_bp=go[which(go$Category=="GO Biological Processes"),]
go_bp$class="BP"
go_mf=go[which(go$Category=="GO Molecular Functions"),]
go_mf$class="MF"
go_cc=go[which(go$Category=="GO Cellular Components"),]
go_cc$class="CC"
res=rbind(go_bp[c(1:15),],go_mf[c(1:3),])
res=rbind(res,go_cc[c(1:2),])
res$pvalue=10^(res$LogP)
head(res)
colnames(res)
colnames(res)[4]="pathway"
colnames(res)[7]="enrichment"
colnames(res)[12]="count"
res=res[,c("pathway","enrichment","pvalue","count","class")]
write.csv(res,"metascape_GO_degdown_bpmfcc.csv",quote=F,row.names = F)

##å¾®ç”Ÿä¿?
#KEGG
setwd("C:\\Users\\Jeff\\OneDrive\\æ¡Œé¢\\Goose_ATAC_info\\RNA-seq\\GO\\")
kegg=read.csv("kegg_all.txt",sep="\t",comment.char = "#",header=T)
kegg=kegg[which(kegg$P.Value <= 0.05),]
colnames(kegg)
colnames(kegg)[4]="Count"
colnames(kegg)[6]="pvalue"
kegg=kegg[,c(1,4,6)]
write.csv(kegg,"kegg.csv",quote=F,row.names = F)
##å¾®ç”Ÿä¿?




##å¾®ç”Ÿä¿?
#Heatmap
setwd("C:\\Users\\Jeff\\OneDrive\\æ¡Œé¢\\Goose_ATAC_info\\RNA-seq\\deseq2")
degup=read.csv("degup.csv",row.names = 1)
degdown=read.csv("degdown.csv",row.names = 1)
trans=read.table("LayBreed_trans.Count_matrix.xls.DESeq2.normalized.xls")
heat=trans[c(rownames(degup),rownames(degdown)),]
write.csv(heat,"degHeatmap.csv")

##å¾®ç”Ÿä¿?
library(pheatmap)
setwd("C:\\Users\\Jeff\\OneDrive\\æ¡Œé¢\\Goose_ATAC_info\\ATAC-seq\\DiffBindv4\\reshape")

pdf("LayBreed_res_Heatmap.pdf",height = 4,width = 3)
Heatmap_Data=read.csv("degHeatmap.csv",row.names = 1,header=F)

annotation_col = data.frame(
  Conditions  = factor(c(rep("Laying", 3),rep("Breeding", 3)) ))
head(annotation_col)

class(Heatmap_Data)
colnames(Heatmap_Data)=factor(colnames(Heatmap_Data))

pheatmap(Heatmap_Data),
         scale = "row",
         show_rownames = F,
         cluster_cols = F,
         main = "DEG Heatmap",
         cutree_rows = 2,
         gaps_col  = 3,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         ColSideColors = rep(c("pink", "purple"), each = 4))

dev.off()









setwd("C:\\Users\\Jeff\\OneDrive\\×ÀÃæ\\Goose_ATAC_info\\ATAC-seq\\DiffBindv4\\reshape\\Lay_vs_Breed_deseq2_up_reshape.motifDir")

library("ggplot2")
library("ggh4x")
library("cowplot")
data<-read.csv("tf.txt",sep="\t",encoding = "UTF-8")

colnames(data)=c("TF","p","target","backgroud","enrichment","class")
data$TF=factor(data$TF,levels=c("RFX(HTH)","RFX2(HTH)","RFX1(HTH)","X-BOX(HTH)","RFX5(HTH)","BAPX1(HOMEOBOX)","NKX2-2(HOMEOBOX)",
                                "NKX2-1(HOMEOBOX)","NKX2-5(HOMEOBOX)","NKX3-1(HOMEOBOX)","SOX10(HMG)","SOX3(HMG)","SOX6(HMG)","SOX21(HMG)","MEF2C(MADS)",
                                "MEF2A(MADS)","MEF2B(MADS)","MEF2D(MADS)","EKLF(ZF)","KLF3(ZF)","SPIB(ETS)",
                                "PGR(NR)","GRE(NR)v1","GRE(NR)v2","ARE(NR)","THRA(NR)v2","THRB(NR)v3","JUN-AP1(BZIP)","HLF(BZIP)","NFIL3(BZIP)","PR(NR)","AR-HALFSITE(NR)","SPL9(SBP)","SPL11(SBP)",
                                "LXRE(NR)","RORGT(NR)","THRA(NR)v1","THRB(NR)v1","THRB(NR)v2","ELF-1(E2F)","SOX15(HMG)","TGA1(BZIP)","EAR2(NR)","COUP-TFII(NR)"))
p<-ggplot(data,aes(x=class,y=TF))+
  geom_point(aes(size=`enrichment`,color=-log(p)),pch=16,alpha = 0.7)+
  theme_bw()+
  scale_size_continuous(range=c(2,8))+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=30,hjust = 0.5,vjust=0))+
  labs(x=NULL,y=NULL)+
  scale_colour_gradient(low = 'blue', high = 'red')+
  scale_y_discrete(position = "left")+ #yÖáÎÄ×Ö·Å×ó²à
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
ggsave("C:\\Users\\Jeff\\OneDrive\\×ÀÃæ\\quail_ATAC_info\\ATAC\\Motif_bubble_plot.pdf",height=10,width=6)
# guide_fill <- get_legend(p+ guides(color = "none") + theme(legend.position = "bottom"))
# plot_grid(p +
#             guides(fill = "none") + 
#             theme(legend.position = "bottom"), 
#           guide_fill, nrow = 2, rel_heights = c(10, 1))


data=data.frame(TF=c(133,11),class=c("up","down"))
ggplot() + 
  geom_col(data = data, aes(x = class,y= TF))+theme_classic()




