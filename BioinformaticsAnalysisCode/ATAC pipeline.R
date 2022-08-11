
###统计基因组gtf文件中的基因数量
cut -f3 quail.gtff|grep  'gene'|wc -l


ls *.R1.fastq.gz|while read file
do
fq1=$file
fq2=${file%.R1.fastq.gz*}.R2.fastq.gz
out=${fq1%_L*}
echo " trim_galore cut adapters started at $(date)"
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
            --paired $fq1 $fq2  \
            --gzip -o /md01/changjy/data/Goose/Goose_bulk/trim_galore/$out
echo "trim_galore cut adapters finished at $(date)"
done



ls *1.fq.gz|sed 's/_1.fq.gz//g'|\
awk '{print "~/miniconda3/bin/bwa mem \
/home/share/hzh/genome/HJX74bwaindex/HJX74.fa \
/home/share/hzh/genome/cx27/"$1"_1.fq.gz / \
home/share/hzh/genome/cx27/"$1"_2.fq.gz>"$1".mem"}'>work.sh
split -l 1 work.sh -d -a 2 run_
ls run_0*|awk '{print "mv "$1" "$1".sh"}'|sh
ls run_0*.sh|less|awk '{print "nohup bash "$1" 1>"$1".log 2>"$1".err &"}'|less|sh


#bowtie2
cat all_number.txt|while read line
do
arr=(${line//,/})
fq1=${arr[0]}
fq2=${arr[1]}
out=${arr[3]}
ref=/md01/changjy/data/quail/ref/ncbi_dataset/index/quail
echo "$fq1"
echo "$fq2"
echo "ref index:$ref"
echo "bowtie2 start at $(date +%F%n%T)" 

bowtie2  -p 5  --very-sensitive -X 2000 -x  $ref -1 $fq1 -2 $fq2 |samtools sort  -O bam  -@ 5 -o - > $out.bam  

echo "bowtie2 mapping finished at $(date +%F%n%T)"
done




#awk ChrM rate
ls |while read line
do
sam_file=$line
echo "sam_file is : $sam_file"
echo "caculate ChrM rate at $(date +%F%n%T)"
grep "NC_023832.1" $sam_file|wc -l
wc -l $sam_file
echo "caculate ChrM rate finished ... at $(date +%F%n%T)"
done
>ChrM.log
#unChrM  ##Goose NC_023832.1
ls *.bam|while read line
do
bam_file=$line
out=${bam_file%.bam*}.q10.unchrM.bam
echo "----------------------------------"
echo "sam_file is : $sam_file"
echo "extract unChrM.sam at $(date +%F%n%T)"
#awk '$3!="NC_023832.1"|| NF<10' $bam_file |samtools view -S -b -f 0x2 -q 10 |samtools sort -o   > $out
awk '$3!="NC_023832.1"' $bam_file |samtools view -S -b -f 0x2 -q 10 |samtools sort -o  > $out
echo "extract unChrM  finished ... at $(date +%F%n%T)"
done

#transfer bam and q10 sort
ls *.unchrM.bam|while read line
do
unchrM_file=$line
out=${unchrM_file%.*}
echo "sam_file is : $sam_file"
echo "output is $out.pe.q10.sort."
echo "extract unChrM.sam at $(date +%F%n%T)"
samtools view -S -b -f 0x2 -q 10 $unchrM_file|samtools sort -o $out.pe.q10.sort.bam
echo "extract unChrM  finished ... at $(date +%F%n%T)"
echo "----------------------------------"
echo "\n"
echo "\n"
done




# 标记PCR重复
ls *pe.q10.sort.bam|while read bam
do
time gatk MarkDuplicates -I $bam -O ${bam%.bam*}.markdup.bam -M ${bam%.bam*}.markdup_metrics.txt && echo "** markdup done **"  #去除这些由PCR扩增所形成的duplicates
done

#Picard
ls *.sort.bam|while read line
do
unchrM_file=$line
out=${unchrM_file%.*}
echo "bam_file is : $line"
echo "output is $out.rmdup.bam"
echo "Start: Picard remove PCR duplicate at $(date +%F%n%T)"
picard MarkDuplicates 
INPUT=$line 
OUTPUT=$out.rmdup.bam 
METRICS_FILE=$out.Picard_Metrics_unfiltered_bam.txt 
VALIDATION_STRINGENCY=LENIENT 
ASSUME_SORTED=true 
REMOVE_DUPLICATES=true
echo "Sucessfully::Picard remove PCR duplicate finished ... at $(date +%F%n%T)"
echo "----------------------------------"
echo "\n"
echo "\n"
done> Picard.log


#Preseq (library complexity) 
ls *q10.sort.bam|while read line
do
out=${unchrM_file%.*}
echo "bam_file is : $line"
echo "output is $out.rmdup.bam"
echo "Start: Picard remove PCR duplicate at $(date +%F%n%T)"
echo "Library complexity estimation by Preseq"
echo "Start: Preseq caculate library complexity at $(date +%F%n%T)"
preseq  lc_extrap-B -P -o /md01/changjy/data/quail/qc/library_complexity/$line.library_complexity $line
echo "Sucessfully:Preseq caculate library complexity finished... at $(date +%F%n%T)"
echo "----------------------------------"

#fragment_distribution
perl /md01/jinxu/ori/pipeline_perl/ATAC-seq-perl-jin/Code/fragment_length_dist.pl $INPUT $OUTPUT.fragment_length_dist
sort $OUTPUT.fragment_length_dist >../fragment/$OUTPUT.fragment_length_sort_dist
pdf("SD-P8-M5-2_L1_806D74.fragment_length.pdf")
hist(data$V1,breaks=150,xlab="Insertsize(bp)",ylab="count",xlim=c(1,500> main="SD-P8-M5-2_L1_806D74 Fragment distribution")
dev.off()

ls *q10.sort.rmdup.bam|while read line
do
samtools index $line > $line.bai
done

#bw file convert
ls *q10.sort.rmdup.bam|while read line
do
out=${line%.bam*}
bamCoverage -p 20 --normalizeUsing RPKM -b $line -o $out.RPKM.bw
done
#计算TSS，及图像

ls *.bw|while read line
do 
out=${line%.bw*}
computeMatrix reference-point -p 15 --referencePoint TSS -b 2000 -a 2000 -R /md01/changjy/data/quail/ref/ncbi_dataset/Quail.TSS -S $line --skipZeros  -out $out.TSS.gz  --outFileSortedRegions $out.genes.bed
done


#转化为bed文件
ls *sort.rmdup.bam|while read line
do
out=${line%.bam*}
bedtools bamtobed -i $line  > ./$out.bed
done
##shift
ls *.rmdup.bed|while read line
do
out=${line%.bed*}	
awk  '{{if($6=="+")print $1"\t"$2+4"\t"$3"\t"$4"\t"$5"\t"$6}if($6=="-"){print $1"\t"$2"\t"$3-5"\t"$4"\t"$5"\t"$6}}' $line >$out.shift.bed
done



plotProfile -m matrix.gz -out profile.pdf
plotHeatmap -m matrix.gz -out heatmap.pdf

#Macs2  callpeak
macs2 callpeak -f BED -t mybed -g 900738757 -n peak.names --outdir ./ --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -q 0.05
idr --samples LD-P12-M1-1_L3_801D73.unchrM.pe.q10.sort.rmdup.bed_peaks.narrowPeak LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bed.peaks_peaks.narrowPeak --rank p.value --output-file LD-P12-M1-1_L3andLD-P12-M2-1_L2_801D74.idr  --plot


ggplot(data,aes(x=c1,c2))+geom_point()+theme_classic()+theme(axis.text.x=element_text(angle=90, hjust=1))+scale_y_continuous(limits=c(0,50))+xlab("sample")+ylab("chrM ratio")+title("sample chrM ratio")



####DiffBind
library("DiffBind")
library("AnnotationHub")
library("biomaRt")
library("clusterProfiler")
library("topGO")
library("Rgraphviz")
library("pathview")
library("ChIPseeker")
library("GenomicFeatures")
library("dplyr")
##载入sample文件
sampleSheet="/md01/changjy/data/quail/next_1.0/diffbind/no_shift/noshift_Quail_sd_7_28.csv"
sample_info <-dba(sampleSheet=sampleSheet)
#sample_info <-dba(sampleSheet="/md01/changjy/data/quail/next_1.0/no_shift/diffbind/noshift_Quail_sd_7_28.csv")
sample_info
#dba_meta <- dba(minOverlap = 1, sampleSheet = sample_info)
#dba_meta
dba_count1 <- dba.count(sample_info,bUseSummarizeOverlaps=TRUE,bParallel=40)
dba_count2 <- dba.count(sample_info,summits=100,bParallel=40)
dba_count1
dba_count2

###构建对比集
dba_contrast1 <- dba.contrast(dba_count1, categories=DBA_CONDITION,minMembers=3)
dba_contrast2 <- dba.contrast(dba_count2, categories=DBA_CONDITION,minMembers=3)
##设种
##差异分析
dba_diff1 <- dba.analyze(dba_contrast1,DBA_ALL_METHODS)
dba_diff2 <- dba.analyze(dba_contrast2,DBA_ALL_METHODS)

#dba.plotHeatmap(dba_count,RowAttributes = DBA_TISSUE,ColAttributes = F)
#dba.plotPCA(dba_count,label=DBA_ID)
#dev.off()
#pdf(paste0("plot/",experiment_name,"_MAplot.pdf"))
#dba.plotMA(dba_diff,cex.main=0.8,contrast=?,method=DBA_EDGER,bXY=TRUE,bSmooth=FALSE,bNormalized=TRUE)
#abline(h = c(-?,?),col = "#ec008c", lty = 5)
#dev.off()
##存储差异分析结果
# dba_report_all_1 <- dba.report(dba_diff,th = 1,contrast=1)
# dba_report_all_2<- dba.report(dba_diff,th = 1,contrast=2)
# dba_report_all_3<- dba.report(dba_diff,th = 1,contrast=3)
# dba_report_all_4<- dba.report(dba_diff,th = 1,contrast=4)
# dba_report_all_5<- dba.report(dba_diff,th = 1,contrast=5)
# dba_report_all_6<- dba.report(dba_diff,th = 1,contrast=6)
#experiment为实验组对照组的名字
#dba_report_all$feature_id <- paste0(shortday_longday3day,"_",names(dba_report_all))
dba_report_diff1_1<-dba.report(dba_diff1,method=DBA_DESEQ2,contrast=1)
dba_report_diff1_2<-dba.report(dba_diff1,method=DBA_DESEQ2,contrast=2)
dba_report_diff1_3<-dba.report(dba_diff1,method=DBA_DESEQ2,contrast=3)
dba_report_diff1_4<-dba.report(dba_diff1,method=DBA_DESEQ2,contrast=4)
dba_report_diff1_5<-dba.report(dba_diff1,method=DBA_DESEQ2,contrast=5)
dba_report_diff1_6<-dba.report(dba_diff1,method=DBA_DESEQ2,contrast=6)

dba_report_diff2_1<-dba.report(dba_diff2,method=DBA_DESEQ2,contrast=1)
dba_report_diff2_2<-dba.report(dba_diff2,method=DBA_DESEQ2,contrast=2)
dba_report_diff2_3<-dba.report(dba_diff2,method=DBA_DESEQ2,contrast=3)
dba_report_diff2_4<-dba.report(dba_diff2,method=DBA_DESEQ2,contrast=4)
dba_report_diff2_5<-dba.report(dba_diff2,method=DBA_DESEQ2,contrast=5)
dba_report_diff2_6<-dba.report(dba_diff2,method=DBA_DESEQ2,contrast=6)


#查看peak上调或下调
# sum(dba_report_diff_1$Fold>0)
# sum(dba_report_diff_1$Fold<0)
# sum(dba_report_diff_2$Fold>0)
# sum(dba_report_diff_2$Fold<0)
# sum(dba_report_diff_3$Fold>0)
# sum(dba_report_diff_3$Fold<0)
# sum(dba_report_diff_4$Fold>0)
# sum(dba_report_diff_4$Fold<0)
# sum(dba_report_diff_5$Fold>0)
# sum(dba_report_diff_5$Fold<0)
# sum(dba_report_diff_6$Fold>0)
# sum(dba_report_diff_6$Fold<0)


#DC143C
    #DA70D6
        #9400D3

txdb<-makeTxDbFromGFF("/md01/changjy/data/quail/ref/ncbi_dataset/quail.gff",organism="Coturnix japonica",dataSource="NCBI_ANNOTATION")
peakAnno1_1 <- annotatePeak(dba_report_diff1_1,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno1_2<- annotatePeak(dba_report_diff1_2,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno1_3<- annotatePeak(dba_report_diff1_3,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno1_4<- annotatePeak(dba_report_diff1_4,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno1_5<- annotatePeak(dba_report_diff1_5,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno1_6<- annotatePeak(dba_report_diff1_6,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")

peakAnno2_1 <- annotatePeak(dba_report_diff2_1,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno2_2<- annotatePeak(dba_report_diff2_2,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno2_3<- annotatePeak(dba_report_diff2_3,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno2_4<- annotatePeak(dba_report_diff2_4,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno2_5<- annotatePeak(dba_report_diff2_5,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
peakAnno2_6<- annotatePeak(dba_report_diff2_6,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")

#保存chipseeker注释结果
#peakAnno
#peakAnno@anno

##dba_count1 <- dba.count(sample_info,bUseSummarizeOverlaps=TRUE,bParallel=40) shift

write.table(as_tibble(peakAnno1_1@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/DiffBind_Deseq2_peakAnno_control_3d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_2@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/DiffBind_Deseq2_peakAnno_control_7d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_3@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/DiffBind_Deseq2_peakAnno_control_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_4@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/DiffBind_Deseq2_peakAnno_3d_7d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_5@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/DiffBind_Deseq2_peakAnno_3d_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_6@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/DiffBind_Deseq2_peakAnno_7d_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)


##dba_count2 <- dba.count(sample_info,summits=100,bParallel=40) shift
write.table(as_tibble(peakAnno2_1@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/summit100_DiffBind_Deseq2_peakAnno_control_3d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_2@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/summit100_DiffBind_Deseq2_peakAnno_control_7d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_3@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/summit100_DiffBind_Deseq2_peakAnno_control_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_4@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/summit100_DiffBind_Deseq2_peakAnno_3d_7d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_5@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/summit100_DiffBind_Deseq2_peakAnno_3d_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_6@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/shift/summit100_DiffBind_Deseq2_peakAnno_7d_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)

###dba_count1 <- dba.count(sample_info,bUseSummarizeOverlaps=TRUE,bParallel=40)  no shift
write.table(as_tibble(peakAnno1_1@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_DiffBind_Deseq2_peakAnno_control_3d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_2@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_DiffBind_Deseq2_peakAnno_control_7d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_3@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_DiffBind_Deseq2_peakAnno_control_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_4@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_DiffBind_Deseq2_peakAnno_3d_7d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_5@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_DiffBind_Deseq2_peakAnno_3d_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno1_6@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_DiffBind_Deseq2_peakAnno_7d_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)


##dba_count2 <- dba.count(sample_info,summits=100,bParallel=40)  no shift
write.table(as_tibble(peakAnno2_1@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_summit100_DiffBind_Deseq2_peakAnno_control_3d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_2@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_summit100_DiffBind_Deseq2_peakAnno_control_7d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_3@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_summit100_DiffBind_Deseq2_peakAnno_control_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_4@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_summit100_DiffBind_Deseq2_peakAnno_3d_7d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_5@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_summit100_DiffBind_Deseq2_peakAnno_3d_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(as_tibble(peakAnno2_6@anno), file = "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_summit100_DiffBind_Deseq2_peakAnno_7d_28d.txt",sep = '\t', quote = FALSE, row.names = FALSE)





library(biomaRt)
library(curl)
datasets<-listDatasets(my_mart)
View(datasets)
my_dataset<-useEnsembl("cjaponica_gene_ensembl", mirror = "useast",biomart = "ensembl")
pheatmap(matrix,scale="row",k=15)

#bed文件根据基因组进行排序
sort -k1,1 -k2,2n DiffBind_Deseq2_peakAnno_all.txt > DiffBind_Deseq2_peakAnno_all_sorted.txt
bedtools merge -i DiffBind_Deseq2_peakAnno_all_sorted.txt > DiffBind_Deseq2_peakAnno_all_sorted_merged.bed
bedtools multicov -bams *.bams -bed  DiffBind_Deseq2_peakAnno_all_sorted_merged.bed> DiffBind_Deseq2_peakAnno_all_sorted_merged_count.txt


#矩阵quantile normalization 然后scale 画图

# matrix<-read.table("/md01/changjy/data/quail/next/DiffBind_Deseq2_peakAnno_all_sorted_merged_PeakCount.txt")
# matrix<-matrix[,c(4:17)]
# set.seed(123)
# Heatmap(matrix=t(scale(t(normalize.quantiles(as.matrix(matrix))),center=F,scale=T)),cluster_columns=F,show_row_names=F)
# set.seed(123)
# pic<-Heatmap(matrix=t(scale(t(normalize.quantiles(as.matrix(obj))),center=F,scale=T)),
# show_row_names=F,
# row_km=2,
# cluster_columns=F,
# column_title = "sample", row_title = "peakset",
# column_title_gp = gpar(fontsize=20, fontface="bold"),
# col = c("#003C67B2","white","#ED0000B2"),  
# show_parent_dend_line = FALSE,
# name = "legend",
# row_title_gp = gpar(fontsize=20, fontface="bold"))
library(ComplexHeatmap)
library(preprocessCore)
matrix<-read.table("/md01/changjy/data/quail/next/DiffBind_Deseq2_peakAnno_all_sorted_merged_PeakCount.txt")
matrix<-matrix[,c(4:17)]
matrix<-t(scale(t(normalize.quantiles(as.matrix(matrix))),center=T,scale=T))
colnames(matrix)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
#提取图画信息
set.seed(1)
list<-draw(Heatmap(matrix=matrix,
             show_row_names=F,
            # show_column_names = T,
             show_row_dend = F,
            row_labels = F,
             row_km=4,
             cluster_columns=F,
             #column_title = "sample", row_title = "peakset",
             column_title_gp = gpar(fontsize=15, fontface="bold"),
             #col = c("#FEAC5E","#C779D0","#4bc0c8"), 
#            cluster_column_slices = F,
            # top_annotation =  ha_column1,
             show_parent_dend_line = FALSE,
             name = "Z-score"))
             #row_title_gp = gpar(fontsize=15, fontface="bold"))

row_order(list)
column_order(list)
cluster1 = matrix[row_order(list)[[2]],column_order(list)] #提取第一簇的相关信息
cluster2 = matrix[row_order(list)[[1]],column_order(list)] #同上
cluster3 = matrix[row_order(list)[[3]],column_order(list)] #同上
cluster4 = matrix[row_order(list)[[4]],column_order(list)] #同上


rownames(cluster1)<-row_order(list)[[2]]
rownames(cluster2)<-row_order(list)[[1]]
rownames(cluster3)<-row_order(list)[[3]]
rownames(cluster4)<-row_order(list)[[4]]

dim(cluster1)
dim(cluster2)
dim(cluster3)
dim(cluster4)






peakanno<-read.csv("peak.txt",)
#下种，进行knn聚类进行分簇
set.seed(123)
list<-draw(Heatmap(matrix,
show_row_names=T,
row_km=4,
cluster_columns=F,
column_title = "sample", row_title = "peakset",
column_title_gp = gpar(fontsize=20, fontface="bold"),
col = c("#003C67B2","white","#ED0000B2"),  
show_parent_dend_line = FALSE,
name = "legend",
row_title_gp = gpar(fontsize=20, fontface="bold")))
row_order(list)
column_order(list)
###complexheatmap extract cludter
d

#画MAplot 图，看看应该筛选哪几个？
CM12[["sign"]]= ifelse( abs(CM12$log2FoldChange) >= 1.5 & log2(CM12$baseMean)>=5.5, "DEP", "normal")
aes = aes(x=log2(baseMean),y=log2FoldChange,fill = sign)
ggplot(CM12, aes) + geom_point(aes(color=sign),size=1) + 
  scale_fill_manual(values = c( "DEP" = "red", "normal" = "grey"))


##GO注释以及KEGG注释

library("AnnotationHub")
library("biomaRt")
library("clusterProfiler")
library("topGO")
library("Rgraphviz")
library("pathview")
# 导入文件
data <- read.table(file="data/quail/peak/test/cluster1_anno.txt",sep='\t')
# 搜寻 OrgDb（比如我想查找物种 Theobroma cacao）
hub <- AnnotationHub::AnnotationHub()
query(hub,"Anser cygnoides")
quail_org<-hub[["AH86240"]]
columns(quail_org)
# 查看Gene ID 类型
columns(quail_org)
# ID转换（将我们的ID转换成与Org.Db一致的ID类型）
data <- as.character(data$geneId)
data_id <- mapIds(x = quail_org,keys = data,keytype = "SYMBOL",column = "ENTREZID")
# 去除NA值
na.omit(data_id)
# 富集（其他富集函数参见R包文档）
erich.go.BP <- enrichGO(gene=data,OrgDb = quail_org,keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = T)
# 绘制条形图
barplot(erich.go.BP)
# 绘制气泡图
dotplot(erich.go.BP)

bitr(keys(quail_org)[1], 'ENTREZID', c("REFSEQ", "GO", "ONTOLOGY"), quail_org)
sample_genes
res = enrichGO(sample_genes, OrgDb=quail_org, pvalueCutoff=1, qvalueCutoff=1)



#HOMER用法
for i in {1..10}; \
do findMotifsGenome.pl /md01/changjy/data/quail/peak/DESeq2/test/cluster/cluster1.txt \
 /md01/changjy/data/quail/ref/ncbi_dataset/index/quail.fa ./position -find known$i.motif \
 -p 40 >motif$i.position; \
done

for i in {1..32}
findMotifsGenome.pl /md01/changjy/data/quail/peak/deseq2 ../../ref/ncbi_dataset/quail.fa  ../../motif/cluster_1_motif/ -find ../../motif/cluster_1_motif/knownResults/known$i.motif -p 40 > cluster_1_known$i.motif.position

str_extract_all(b,"[0-9]+[0-9]")
as.data.frame(do.call(rbind,str_extract_all(as.character(cluster_1_motif$PositionID),"[0-9]+[0-9]")))
cluster_1_1_motif_annotation<-peak_annotation[which(rownames(peak_annotation)%in%rownames(cluster_1_motif)),]

#find之后的motif进行注释 
rm(list=ls())
peak_id<-paste("cluster_1_",c(1:1054),sep="")
for (i in 1:70){
	motif_name<-paste("data/quail/motif/motif_position/cluster1_find/cluster_1_known",i,".motif.position",sep="")
	file_name<-paste("data/quail/motif/motif_annotation/cluster_1_motif_annotation/cluster_1_",i,"_motif.anno",sep="")
	motif<-read.csv(motif_name,sep="\t",header=T)
	anno<-read.csv("data/quail/peak/peakAnno_tb.csv",sep="\t",header=T)
	peak_id<-paste("cluster_1_",c(1:1063),sep="")
	anno<-cbind(anno[,-1],peak_id,sep="\t")
	head(anno)
	library(stringr)
#	PositionID<-str_extract_all(motif$PositionID,"[0-9]+[0-9]")
#	motif[,1]<-as.vector(unlist(PositionID))
	head(motif)
	names(motif)[names(motif) == 'PositionID'] <- "peak_id"
	motif_anno<-merge(motif,anno,by="peak_id")
	write.table(motif_anno,file_name,sep = '\t', quote = FALSE, row.names = FALSE)

}

rm(list=ls())
peak_id<-paste("cluster_2_",c(1:1063),sep="")
for (i in 1:28){
	motif_name<-paste("data/quail/motif/motif_position/cluster2_find/cluster_2_known",i,".motif.position",sep="")
	file_name<-paste("data/quail/motif/motif_annotation/cluster_2_motif_annotation/cluster_2_",i,"_motif.anno",sep="")
	motif<-read.csv(motif_name,sep="\t",header=T)
	anno<-read.csv("data/quail/peak/peakAnno_tb.csv",sep="\t",header=T)
	peak_id<-paste("cluster_2_",c(1:1063),sep="")
	anno<-cbind(anno[,-1],peak_id,sep="\t")
	head(anno)
	library(stringr)
#	PositionID<-str_extract_all(motif$PositionID,"[0-9]+[0-9]")
#	motif[,1]<-as.vector(unlist(PositionID))
	head(motif)
	names(motif)[names(motif) == 'PositionID'] <- "peak_id"
	motif_anno<-merge(motif,anno,by="peak_id")
	write.table(motif_anno,file_name,sep = '\t', quote = FALSE, row.names = FALSE)

}


rm(list=ls())
peak_id<-paste("cluster_3_",c(1:1063),sep="")
for (i in 1:28){
	motif_name<-paste("/md01/changjy/data/quail/motif/motif_position/cluster3_find/cluster_3_known",i,".motif.position",sep="")
	file_name<-paste("/md01/changjy/data/quail/motif/motif_annotation/cluster_3_motif_annotation/cluster_3_",i,"_motif.anno",sep="")
	motif<-read.csv(motif_name,sep="\t",header=T)
	anno<-read.csv("/md01/changjy/data/quail/peak/peakAnno_tb.csv",sep="\t",header=T)
	peak_id<-paste("cluster_3_",c(1:1063),sep="")
	anno<-cbind(anno[,-1],peak_id,sep="\t")
	head(anno)
	library(stringr)
#	PositionID<-str_extract_all(motif$PositionID,"[0-9]+[0-9]")
#	motif[,1]<-as.vector(unlist(PositionID))
	head(motif)
	names(motif)[names(motif) == 'PositionID'] <- "peak_id"
	motif_anno<-merge(motif,anno,by="peak_id")
	write.table(motif_anno,file_name,sep = '\t', quote = FALSE, row.names = FALSE)

}


library(tidyverse)
rm(list=ls())
peak_id<-paste("cluster_4_",c(1:1063),sep="")
for (i in 1:25){
	motif_name<-paste("data/quail/motif/motif_position/cluster4_find/cluster_4_known",i,".motif.position",sep="")
	file_name<-paste("data/quail/motif/motif_annotation/cluster_4_motif_annotation/cluster_4_",i,"_motif.anno",sep="")
	motif<-read.csv(motif_name,sep="\t",header=T)
	anno<-read.csv("data/quail/peak/peakAnno_tb.csv",sep="\t",header=T)
	peak_id<-paste("cluster_4_",c(1:1063),sep="")
	anno<-cbind(anno[,-1],peak_id,sep="\t")
	head(anno)
	library(stringr)
#	PositionID<-str_extract_all(motif$PositionID,"[0-9]+[0-9]")
#	motif[,1]<-as.vector(unlist(PositionID))
	head(motif)
	names(motif)[names(motif) == 'PositionID'] <- "peak_id"
	motif_anno<-merge(motif,anno,by="peak_id")
	write.table(motif_anno,file_name,sep = '\t', quote = FALSE, row.names = FALSE)

}


enrich.go.BP = enrichGO(gene       = cluster1_symbol,
                        OrgDb        = quail_org,
                        keyType      = 'SYMBOL',ont= "BP",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.05,
                        readable     = T)


library(ggplot2)
ggplot(a, aes(type,Name) + 
    geom_point(aes(size=Log.Pvalue))) + 
    theme_minimal() + 
    xlab(NULL) + ylab(NULL) + 
    scale_size_continuous(range=c(1, 10))

#motif的简单可视化

library(ggplot2)
p<-ggplot(a,aes(type,TFs_family))+
	ggtitle("HOMER de nove motif")+
	theme(plot.title = element_text(hjust = 0.5))+
	geom_point(aes(size=abs(Log_Pvalue),color=log(enrichment.score)))+
	scale_color_gradient(low="green",high = "red")+
	theme_minimal()+xlab(NULL)+ylab(NULL)+
	scale_size_continuous(range=c(1, 4),breaks=seq(0, 10, 5))+
	theme(axis.text.x=element_text(size=12,angle=-10),axis.text.y=element_text(size=7,family="myFont")) 
	
	 ##自定义坐标轴

p 
 

ggplot(a, aes(x = rank(Log.Pvalue), y = enrichment.score, shape = type, colour = classification)) +

geom_point()

ggplot(a, aes(y = rank(P.value), x = log(enrichment),color=family))+
ggtitle("known_motif")
scale_y_continuous(limits = c(0,140))+
geom_point(size=5)+
scale_x_continuous(limits = c(2,6))+
scale_color_manual(values=c("green", "#E69F00", "#56B4E9","red","blue","grey","black","yellow"))+
theme(panel.grid.major =element_blank(), 
	panel.grid.minor = element_blank(),
	panel.background = element_blank(),
	axis.line = element_line(colour = "black"))



#bowtie2比对
cat all_number.txt|while read line
do
arr=(${line//,/})
fq1=${arr[0]}
fq2=${arr[1]}
out=${arr[3]}
ref=/md01/changjy/data/quail/ref/ncbi_dataset/index/quail
echo "$fq1"
echo "$fq2"
echo "ref index:$ref"
echo "bowtie2 start at $(date +%F%n%T)" 
bowtie2	-p 10   --very-sensitive   -x $ref -1 $fq1 -2 $fq2  -S ../mapping/$out.sam   
echo "bowtie2 mapping finished at $(date +%F%n%T)"
done


#计算ChrM比对率 
ls *.sam|while read line
do
sam_file=$line
echo "sam_file is : $sam_file"
echo "caculate ChrM rate at $(date +%F%n%T)"
grep "LR605957.1" $sam_file|wc -l
wc -l $sam_file
echo "caculate ChrM rate finished ... at $(date +%F%n%T)"
done


#移除chrM的reads
samtools view -h LD-P10-M10 |grep -v NC_003408.1|samtools view -bS -q 10 -F 1804 -f 0x2 - >LD-P10-M1.filter.bam
ls *.bam|while read line
do
bam_file=$line
echo "----------------------------------"
echo "bam_file is : $sam_file"
echo "extract unChrM.bam at $(date +%F%n%T)"
samtools view $line |awk '$3!="NC_003408.1"|| NF<10' $bam_file |samtools view -S -b   > $bam_file.unchrM.bam
echo "extract unChrM  finished ... at $(date +%F%n%T)"
done

#转为bam并留取quality>10的reads
ls *.unchrM.bam|while read line
do
unchrM_file=$line
out=${unchrM_file%.*}
echo "sam_file is : $sam_file"
echo "output is $out.pe.q10.sort."
echo "extract unChrM.sam at $(date +%F%n%T)"
samtools view -S -b -f 0x2 -q 10 $unchrM_file|samtools sort -o $out.pe.q10.sort.bam
echo "extract unChrM  finished ... at $(date +%F%n%T)"
echo "----------------------------------"
echo "\n"
echo "\n"
done


#Picard去pcr重复
ls *filter.bam|while read line
do
unchrM_file=$line
out=${unchrM_file%.*}
echo "bam_file is : $line"
echo "output is $out.rmdup.bam"
echo "Start: Picard remove PCR duplicate at $(date +%F%n%T)"
picard MarkDuplicates INPUT=$line OUTPUT=$out.rmdup.bam METRICS_FILE=$out.Picard_Metrics_unfiltered_bam.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
echo "Sucessfully::Picard remove PCR duplicate finished ... at $(date +%F%n%T)"
echo "----------------------------------"
echo "\n"
echo "\n"
done> Picard.log


#Preseq 检测文库复杂度
ls *q10.sort.bam|while read line
do
out=${unchrM_file%.*}
echo "bam_file is : $line"
echo "output is $out.rmdup.bam"
echo "Start: Picard remove PCR duplicate at $(date +%F%n%T)"
echo "Library complexity estimation by Preseq"
echo "Start: Preseq caculate library complexity at $(date +%F%n%T)"
preseq  lc_extrap-B -P -o /md01/changjy/data/quail/qc/library_complexity/$line.library_complexity $line
echo "Sucessfully:Preseq caculate library complexity finished... at $(date +%F%n%T)"
echo "----------------------------------"


#绘制插入片段长度直方图
perl /md01/jinxu/ori/pipeline_perl/ATAC-seq-perl-jin/Code/fragment_length_dist.pl $INPUT $OUTPUT.fragment_length_dist
sort $OUTPUT.fragment_length_dist >../fragment/$OUTPUT.fragment_length_sort_dist
pdf("SD-P8-M5-2_L1_806D74.fragment_length.pdf")
hist(data$V1,breaks=150,xlab="Insertsize(bp)",ylab="count",xlim=c(1,500> main="SD-P8-M5-2_L1_806D74 Fragment distribution")
dev.off()


#bam转为bigwig文件
ls *q10.sort.rmdup.bam|while read line
do
out=${line%.bam*}
bamCoverage -p 20 --normalizeUsing RPKM -b $line -o ./bw/$out.RPKM.bw
# ls *-M?|while read id 
# do
# input=$id
# bamCoverage \
# --binSize 10 \
# -p 20 \
# --bam $input \
# --normalizeUsing  RPKM \
# --outFileName ./bw/$id.rpkm.bigwig
# done


#computeMatrix
ls *.bw|while read line
do 
out=${line%.bw*}
computeMatrix reference-point -p 25 --referencePoint TSS -b 3000 -a 3000 -R /md01/changjy/data/quail/ref/ncbi_dataset/Quail.Tss -S $line --skipZeros  -out $out.TSS.gz  --outFileSortedRegions $out.genes.bed
done


##computeMatrix得到的矩阵进行画图
plotProfile -m matrix.gz -out profile.pdf
plotHeatmap --plotFileFormat pdf --heatmapHeight 20 --heatmapWidth 5 --colorMap Greens --missingDataColor white --plotTitle LD --matrixFile SD-P8-M5-2_L1_806D74.unchrM.pe.q10.sort.rmdup.RPKM.TSS.gz  --outFileName heatmap.pdf


###bam文件转为bed文件
ls *.bam|while read line
do
out=${line%.bam*}
 bamtobed -i $line  > ./bed/$out.bed
done


##shift   对于正链上的reads需要向右偏移4bp 负链上的reads, 则向左偏移5bp 即++4bp --5bp
ls *.rmdup.bed|while read line
do
out=${line%.bed*}	
awk  '{{if($6=="+")print $1"\t"$2+4"\t"$3"\t"$4"\t"$5"\t"$6}if($6=="-"){print $1"\t"$2"\t"$3-5"\t"$4"\t"$5"\t"$6}}' $line >$out.shift.bed
done
##shift version2

def tn5_shift_ta(ta, out_dir)：
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))
    shifted_ta = '{}.tn5.tagAlign.gz'.format(prefix)
    cmd = 'zcat -f {} | '
    cmd += 'awk \'BEGIN {{OFS = "\\t"}}'
    cmd += '{{ if ($6 == "+") {{$2 = $2 + 4}} '
    cmd += 'else if ($6 == "-") {{$3 = $3 - 5}} print $0}}\' | '
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(
        ta,
        shifted_ta)
    run_shell_cmd(cmd)    
    return shifted_ta
    ''


###此时的bed已移位  call peak
ls *.bed|while read line
do
out=${line%.bed*}
macs2  callpeak \
-t  $line \
-f BED \
-g  900738757 \
--outdir /md01/changjy/data/quail/next_1.0/MACS_peak \
--shift 37 \
--extsize 73 \
--nomodel  \
-q 0.01 \
-n $out 
done

ls ../../mapping/sort.rmdup.bam/*.bam|while read bam; do out=${bam%.bam*}; macs2 callpeak -t $bam  --nomodel --keep-dup=all  -f BAMPE  -g  900738757 --outdir /md01/changjy/data/quail/next_1.0/bam_peak/ -n $out

##idr求得高重复peak 
idr --samples LD-P12-M1-1_L3_801D73.unchrM.pe.q10.sort.rmdup.bed_peaks.narrowPeak LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bed.peaks_peaks.narrowPeak \
--rank p.value \
--output-file LD-P12-M1-1_L3andLD-P12-M2-1_L2_801D74.idr \
--plot


###plot
ggplot(data,aes(x=c1,c2))+
  geom_point()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  scale_y_continuous(limits=c(0,50))+
  xlab("sample")+
  ylab("chrM ratio")+
  title("sample chrM ratio")




####DiffBind
library("DiffBind")
library("AnnotationHub")
library("biomaRt")
library("clusterProfiler")
library("topGO")
library("Rgraphviz")
library("pathview")
##载入sample文件
sample_info <-dba(sampleSheet="Quail_sd_7_28.csv")
sample_info
#dba_meta <- dba(minOverlap = 1, sampleSheet = sample_info)
#dba_meta
dba_count <- dba.count(sample_info,bParallel=40)
all=as.data.frame(dba_count$binding)
dba_count

idr=as.data.frame(dba_count$merged)
idr[which(idr$CHR=="1"),]$CHR="NC_029516.1"
idr[which(idr$CHR=="2"),]$CHR="NC_029517.1"
idr[which(idr$CHR=="3"),]$CHR="NC_029518.1"
idr[which(idr$CHR=="4"),]$CHR="NC_029519.1"
idr[which(idr$CHR=="5"),]$CHR="NC_029520.1"
idr[which(idr$CHR=="6"),]$CHR="NC_029521.1"
idr[which(idr$CHR=="7"),]$CHR="NC_029522.1"
idr[which(idr$CHR=="8"),]$CHR="NC_029523.1"
idr[which(idr$CHR=="9"),]$CHR="NC_029524.1"
idr[which(idr$CHR=="10"),]$CHR="NC_029525.1"
idr[which(idr$CHR=="11"),]$CHR="NC_029526.1"
idr[which(idr$CHR=="12"),]$CHR="NC_029527.1"
idr[which(idr$CHR=="13"),]$CHR="NC_029528.1"
idr[which(idr$CHR=="14"),]$CHR="NC_029529.1"
idr[which(idr$CHR=="15"),]$CHR="NC_029530.1"
idr[which(idr$CHR=="16"),]$CHR="NC_029531.1"
idr[which(idr$CHR=="17"),]$CHR="NC_029532.1"
idr[which(idr$CHR=="18"),]$CHR="NC_029533.1"
idr[which(idr$CHR=="19"),]$CHR="NC_029534.1"
idr[which(idr$CHR=="20"),]$CHR="NC_029535.1"
idr[which(idr$CHR=="21"),]$CHR="NC_029536.1"
idr[which(idr$CHR=="22"),]$CHR="NC_029537.1"
idr[which(idr$CHR=="23"),]$CHR="NC_029538.1"
idr[which(idr$CHR=="24"),]$CHR="NC_029539.1"
idr[which(idr$CHR=="25"),]$CHR="NC_029540.1"
idr[which(idr$CHR=="26"),]$CHR="NC_029541.1"
idr[which(idr$CHR=="27"),]$CHR="NC_029542.1"
idr[which(idr$CHR=="28"),]$CHR="NC_029543.1"
idr[which(idr$CHR=="29"),]$CHR="NC_029544.1"
idr[which(idr$CHR=="30"),]$CHR="NC_029545.1"
idr[which(idr$CHR=="31"),]$CHR="NC_029547.1"

library(ChIPseeker)
library(GenomicFeatures)
txdb<-makeTxDbFromGFF("/md01/changjy/data/quail/ref/ncbi_dataset/chr_split/genomic.gff",organism="Coturnix japonica",dataSource="NCBI_ANNOTATION")
peak=readPeakFile("/md01/changjy/data/quail/peak/DiffBind_IDRmerged_peak.txt")
#peak=readPeakFile(idr)
idr_peakAnno <- annotatePeak(peak,
                         TxDb = txdb,tssRegion=c(-3000,3000),
                         level = "gene")
idr_peakAnno
idr_peakAnno@anno
library(dplyr)
idr_peakAnno_tb <- as_tibble(idr_peakAnno@anno)
write.table(idr_peakAnno_tb, file = "/md01/changjy/data/quail/peak/idr_peakAnno.txt",sep = '\t', quote = FALSE, row.names = FALSE)
write.table(idr,
                        "/md01/changjy/data/quail/peak/DiffBind_IDRmerged_peak.txt",
                        quote=F,
                        row.names=F,
                        sep="\t")

###构建对比集
dba_contrast <- dba.contrast(dba_count, categories=DBA_CONDITION,minMembers=3)

##设种
setseed(123)

##差异分析
dba_diff <- dba.analyze(dba_contrast,DBA_ALL_METHODS)
#dba.plotHeatmap(dba_count,RowAttributes = DBA_TISSUE,ColAttributes = F)
#dba.plotPCA(dba_count,label=DBA_ID)
#dev.off()
#pdf(paste0("plot/",experiment_name,"_MAplot.pdf"))
#dba.plotMA(dba_diff,cex.main=0.8,contrast=?,method=DBA_EDGER,bXY=TRUE,bSmooth=FALSE,bNormalized=TRUE)
#abline(h = c(-?,?),col = "#ec008c", lty = 5)
#dev.off()


##存储差异分析结果
dba_report_all_1 <- dba.report(dba_diff,contrast=1,method=DBA_DESEQ2)
dba_report_all_2 <- dba.report(dba_diff,contrast=2,method=DBA_DESEQ2)
dba_report_all_3 <- dba.report(dba_diff,contrast=3,method=DBA_DESEQ2)
dba_report_all_4 <- dba.report(dba_diff,contrast=4,method=DBA_DESEQ2)
dba_report_all_5 <- dba.report(dba_diff,contrast=5,method=DBA_DESEQ2)
dba_report_all_6 <- dba.report(dba_diff,contrast=6,method=DBA_DESEQ2)

write.table(dba_report_all_1,"control_3d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_2,"control_7d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_3,"control_28d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_4,"3d_7d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_5,"3d_28d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_6,"7d_28d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
cat *.csv|awk '{print $2"\t"$3"\t"$4}'|sort -k1,1 -k2,2n >all.sorted.DEP
bedtools merge -i all.sorted.DEP >all.sorted.merged.DEP
bedtools multicov \ 
-bed all.sorted.merged.DEP \ 
 -bams \
 /md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M1-1_L3_801D73.unchrM.pe.q10.sort.rmdup.bam \
 /md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bam \
 /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-1_L2_801F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-2_L1_802F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M1-1_L1_807D75.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-1_L1_803F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-2_L1_806F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M5-2_L1_806D74.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M1-2_L1_808F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M3-2_L1_808D73.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M6-2_L1_808D74.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M1-1_L1_809D73.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M3-1_L1_809D74.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M6-1_L2_801D75.unchrM.pe.q10.sort.rmdup.bam > all.sorted.merged.DEP.count


dba_report_all_1=a[rownames(as.data.frame(dba_report_all_1)),]
dba_report_all_2=a[rownames(as.data.frame(dba_report_all_2)),]
dba_report_all_3=a[rownames(as.data.frame(dba_report_all_3)),]
dba_report_all_4=a[rownames(as.data.frame(dba_report_all_4)),]
dba_report_all_5=a[rownames(as.data.frame(dba_report_all_5)),]
dba_report_all_6=a[rownames(as.data.frame(dba_report_all_6)),]

data=rbind(dba_report_all_1,dba_report_all_2)
data=rbind(data,dba_report_all_3)
data=rbind(data,dba_report_all_4)
data=rbind(data,dba_report_all_5)


bedtools multicov \ 
-bed all.sorted.merged.DEP \ 
 -bams \
 /md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M1-1_L3_801D73.unchrM.pe.q10.sort.rmdup.bam \
 /md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bam \
 /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-1_L2_801F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-2_L1_802F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M1-1_L1_807D75.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-1_L1_803F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-2_L1_806F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M5-2_L1_806D74.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M1-2_L1_808F10.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M3-2_L1_808D73.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M6-2_L1_808D74.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M1-1_L1_809D73.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M3-1_L1_809D74.unchrM.pe.q10.sort.rmdup.bam \
  /md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M6-1_L2_801D75.unchrM.pe.q10.sort.rmdup.bam > all.sorted.merged.DEP.count



# data=read.table("edgeR/sorted.merged.DEP.count")
# data=data[4:17]
# head(data)
# dim(data)

# #data=log(data+1)
# data=t(scale(t(data),center=T,scale=T))
# data=normalize.quantiles(as.matrix(data))

# colnames(data)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
# pdf("Heatmap.pdf",height=6.5,width=4)
# data[data>1]=1
# data[data<-1]=-1
# Heatmap(data,cluster_columns=F,row_km=4)
# dev.off()


# data=read.table("/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_DiffBind_Deseq2_peakAnno_all_sorted_merged_count.txt")
# data=data[4:17]
# pdf("Heatmap.pdf",height=6.5,width=4)

# data=t(scale(t(data),center=T,scale=T))
# data=normalize.quantiles(as.matrix(data))
# #data=normalize.quantiles(as.matrix(data))
# colnames(data)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
# bk=unique(c(seq(-2,2, length=100)))
# colors <-brewer.pal(9,"YlGnBu")
# set.seed(123)

# # ped =pheatmap(data, cluster_rows=T, 
# #               show_rownames=F,treeheight_row = 0,
# #               cluster_col=F,
# #               #color = colors, 
# #               #breaks = bk,
# #               cutree_rows=5,
# #               main = "", gaps_row = c(10, 14),
# #               show_colnames = T,fontsize_row = 7.5)
# list<-draw(Heatmap(matrix=data,
#              show_row_names=F,
#             # show_column_names = T,
#              show_row_dend = F,
#             row_labels = F,
#              row_km=4,
#              cluster_columns=F,
#              #column_title = "sample", row_title = "peakset",
#              column_title_gp = gpar(fontsize=15, fontface="bold"),
#              #col = c("#FEAC5E","#C779D0","#4bc0c8"), 
# #            cluster_column_slices = F,
#             # top_annotation =  ha_column1,
#              show_parent_dend_line = FALSE,
#              name = "Z-score"))
#              #row_title_gp = gpar(fontsize=15, fontface="bold"))

# Heatmap(matrix=data,
#              show_row_names=F,
#             # show_column_names = T,
#              show_row_dend = F,
#             row_labels = F,
#              row_km=4,
#              cluster_columns=F,
#              #column_title = "sample", row_title = "peakset",
#              column_title_gp = gpar(fontsize=15, fontface="bold"),
#              #col = c("#FEAC5E","#C779D0","#4bc0c8"), 
# #            cluster_column_slices = F,
#             # top_annotation =  ha_column1,
#              show_parent_dend_line = FALSE,
#              name = "Z-score")
#              #row_title_gp = gpar(fontsize=15, fontface="bold"))

# dev.off()
#experiment为实验组对照组的名字
#dba_report_all$feature_id <- paste0(shortday_longday3day,"_",names(dba_report_all))
#dba_report_diff<-dba.report(dba_diff,method=DBA_EDGER)

#查看peak上调或下调
sum(dba_report_diff$Fold>0)
sum(dba_report_diff$Fold<0)


##利用GFF文件个性化制作TxDb文件
library(GenomicFeatures)
txdb<-makeTxDbFromGFF("/md01/changjy/data/quail/ref/ncbi_dataset/chr_split/genomic.gff",organism="Coturnix japonica",dataSource="NCBI_ANNOTATION")
peakAnno <- annotatePeak(peak,
                         TxDb = txdb,tssRegion=c(-3000,3000)
                         level = "gene")
·
peakAnno
peakAnno@anno
library(dplyr)
peakAnno_tb <- as_tibble(peakAnno@anno)
write.table(peakAnno_tb, file = "peakAnno_tb.txt",sep = '\t', quote = FALSE, row.names = FALSE)

library(biomaRt)
library(curl)
datasets<-listDatasets(my_mart)
View(datasets)
my_dataset<-useEnsembl("cjaponica_gene_ensembl", mirror = "useast",biomart = "ensembl")
pheatmap(matrix,scale="row",k=15)


sort -k1,1 -k2,2n in.bed
Heatmap(matrix=t(scale(t(normalize.quantiles(as.matrix(matrix))),center=F,scale=T)),cluster_columns=F,show_row_names=F)
set.seed(123)
pic<-Heatmap(matrix,
show_row_names=F,
row_km=4,
cluster_columns=F,
column_title = "sample", row_title = "peakset",
column_title_gp = gpar(fontsize=20, fontface="bold"),
col = c("#003C67B2","white","#ED0000B2"),  
show_parent_dend_line = FALSE,
name = "legend",
row_title_gp = gpar(fontsize=20, fontface="bold"))


peakanno<-read.csv("peak.txt",)
set.seed(123)
list<-draw(Heatmap(matrix,
show_row_names=T,
row_km=4,
cluster_columns=F,
column_title = "sample", row_title = "peakset",
column_title_gp = gpar(fontsize=20, fontface="bold"),
col = c("#003C67B2","white","#ED0000B2"),  
show_parent_dend_line = FALSE,
name = "legend",
row_title_gp = gpar(fontsize=20, fontface="bold")))
row_order(list)
column_order(list)
###complexheatmap extract cludter
d

CM12[["sign"]]= ifelse( abs(CM12$log2FoldChange) >= 1.5 & log2(CM12$baseMean)>=5.5, "DEP", "normal")
aes = aes(x=log2(baseMean),y=log2FoldChange,fill = sign)
ggplot(CM12, aes) + geom_point(aes(color=sign),size=1) + 
  scale_fill_manual(values = c( "DEP" = "red", "normal" = "grey"))




library("AnnotationHub")
library("biomaRt")
library("clusterProfiler")
library("topGO")
library("Rgraphviz")
library("pathview")
# 导入文件
data <- read.table(file="data/quail/peak/test/cluster1_anno.txt",sep='\t')
# 搜寻 OrgDb（比如我想查找物种 Theobroma cacao）
hub <- AnnotationHub::AnnotationHub()
query(hub,"Coturnix japonica")
quail_org<-hub[["AH86235"]]
columns(quail_org)
# 查看Gene ID 类型
columns(quail_org)
# ID转换（将我们的ID转换成与Org.Db一致的ID类型）
data <- as.character(data$geneId)
data_id <- mapIds(x = quail_org,keys = data,keytype = "SYMBOL",column = "ENTREZID")
# 去除NA值
na.omit(data_id)
# 富集（其他富集函数参见R包文档）
erich.go.BP <- enrichGO(gene=data,OrgDb = quail_org,keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = T)
# 绘制条形图
barplot(erich.go.BP)
# 绘制气泡图
dotplot(erich.go.BP)

bitr(keys(quail_org)[1], 'ENTREZID', c("REFSEQ", "GO", "ONTOLOGY"), quail_org)
sample_genes
res = enrichGO(sample_genes, OrgDb=quail_org, pvalueCutoff=1, qvalueCutoff=1)

ls -1 ../../motif/cluster_1_motif/knownResults/known*.motif|awk -F "/" '{print $6}'
do
findMotifsGenome.pl ../../peak/deseq2/DiffBind_DESEQ2_3day_DP_new.bed ../../ref/ncbi_dataset/quail.fa  ../../motif/cluster_1_motif/ -find $id -p 40 > cluster_1_$id.position
done

for i in {1..32}
findMotifsGenome.pl /md01/changjy/data/quail/peak/deseq2 ../../ref/ncbi_dataset/quail.fa  ../../motif/cluster_1_motif/ -find ../../motif/cluster_1_motif/knownResults/known$i.motif -p 40 > cluster_1_known$i.motif.position

str_extract_all(b,"[0-9]+[0-9]")
as.data.frame(do.call(rbind,str_extract_all(as.character(cluster_1_motif$PositionID),"[0-9]+[0-9]")))
cluster_1_1_motif_annotation<-peak_annotation[which(rownames(peak_annotation)%in%rownames(cluster_1_motif)),]



motif<-read.csv("data/quail/motif/motif_position/cluster1_find/cluster_1_known1.motif.position",sep="\t")
head(motif)
anno<-read.csv("data/quail/peak/peakAnno_tb.csv",sep="\t")
head(anno)
library(stringr)
PositionID<-str_extract_all(motif$PositionID,"[0-9]+[0-9]")
PositionID<-as.vector(unlist(PositionID))
motif$PositionID<- PositionID
head(motif)
names(motif)[names(motif) == 'PositionID'] <- 'peak_id'
motif_anno<-merge(motif,anno,by="peak_id")
write.table(peakAnno_tb, file = "peakAnno_tb.txt",sep = '\t', quote = FALSE, row.names = FALSE)




#split库
retrieve record with 'object[["AH86235"]]'



ls *-M?|while read id 
do
input=$id
bamCoverage \
--binSize 10 \
-p 20 \
--bam $input \
--normalizeUsing  RPKM \
--outFileName ./bw/$id.rpkm.bigwig
done

bamCoverage -p 25 --binSize 10 -b ../LD-P12-M1-1_L3_801D73.unchrM.pe.q10.sort.rmdup.bam --normalizeUsing RPKM -o  rpkm.bigwig


#footprint FootprintScores分析
ls *.unchrM.pe.q10.sort.rmdup_uncorrected.bw|while read file
do
output=${file%.*}
TOBIAS FootprintScores --signal $file --regions ../peak/peak.txt --output $output.footprint.bw --cores 15
done

#footprint BINDetect分析
ls *unchrM.pe.q10.sort.rmdup_uncorrected.footprint.bw |while read bigwig
  do
  output=${bigwig%unchrM*}
  TOBIAS BINDetect --motifs /md01/changjy/data/quail/motif/cluster_1_motif/knownResults/known1.motif \
    --signals $bigwig \
    --genome /md01/changjy/data/quail/ref/ncbi_dataset/index/quail.fa \
    --peaks /md01/changjy/data/quail/peak/peak.txt  \
    --outdir NKX2-2_$output \
    --cond_names quail \
    --cores 8 
  done

  ls *unchrM.pe.q10.sort.rmdup_uncorrected.footprint.bw |while read bigwig
  do
  output=${bigwig%unchrM*}
  TOBIAS BINDetect --motifs /md01/changjy/data/quail/motif/cluster_1_motif/knownResults/known10.motif \
    --signals $bigwig \
    --genome /md01/changjy/data/quail/ref/ncbi_dataset/index/quail.fa \
    --peaks /md01/changjy/data/quail/peak/peak.txt  \
    --outdir NKX3-1_$output \
    --cond_names quail \
    --cores 8 
  done


    ls *unchrM.pe.q10.sort.rmdup_uncorrected.footprint.bw |while read bigwig
  do
  output=${bigwig%unchrM*}
  TOBIAS BINDetect --motifs /md01/changjy/data/quail/motif/cluster_1_motif/knownResults/known3.motif \
    --signals $bigwig \
    --genome /md01/changjy/data/quail/ref/ncbi_dataset/index/quail.fa \
    --peaks /md01/changjy/data/quail/peak/peak.txt  \
    --outdir NKX2-1_$output \
    --cond_names quail \
    --cores 8 
  done
      ls *unchrM.pe.q10.sort.rmdup_uncorrected.footprint.bw |while read bigwig
  do
  output=${bigwig%unchrM*}
  TOBIAS BINDetect --motifs /md01/changjy/data/quail/motif/cluster_1_motif/knownResults/known11.motif \
    --signals $bigwig \
    --genome /md01/changjy/data/quail/ref/ncbi_dataset/index/quail.fa \
    --peaks /md01/changjy/data/quail/peak/peak.txt  \
    --outdir MEF2C_$output \
    --cond_names quail \
    --cores 8 
  done

      ls *unchrM.pe.q10.sort.rmdup_uncorrected.footprint.bw |while read bigwig
  do
  output=${bigwig%unchrM*}
  TOBIAS BINDetect --motifs /md01/changjy/data/quail/motif/cluster_1_motif/knownResults/known12.motif \
    --signals $bigwig \
    --genome /md01/changjy/data/quail/ref/ncbi_dataset/index/quail.fa \
    --peaks /md01/changjy/data/quail/peak/peak.txt  \
    --outdir MEF2A_$output \
    --cond_names quail \
    --cores 8 
  done

        ls *unchrM.pe.q10.sort.rmdup_uncorrected.footprint.bw |while read bigwig
  do
  output=${bigwig%unchrM*}
  TOBIAS BINDetect --motifs /md01/changjy/data/quail/motif/cluster_1_motif/knownResults/known29.motif \
    --signals $bigwig \
    --genome /md01/changjy/data/quail/ref/ncbi_dataset/index/quail.fa \
    --peaks /md01/changjy/data/quail/peak/peak.txt  \
    --outdir FOXK1_$output \
    --cond_names quail \
    --cores 8 
  done



pdf("/md01/changjy/data/quail/next_1.0/diffbind/no_shift/length_boxplot.pdf",width=2.5,height=3.6)
p=ggplot(data,aes(x=name,y=length,fill="white"))+
stat_boxplot(geom = "errorbar",width=0.15)+
ylim(200,3000)+geom_boxplot()+
xlab("peak")+
title("Different peak length after merged")+
theme_classic()
p
dev.off()

pdf("/md01/changjy/data/quail/next_1.0/diffbind/no_shift/length_hist.pdf",width=2.5,height=3.6)
p=ggplot(data, aes(length)) +
  geom_histogram(aes(fill = "pink"))+theme_classic()
p
dev.off()


library(ComplexHeatmap)
library(preprocessCore)
data=read.table("sorted.merged.DEP.count")[,c(4:17)]
data=data[,c(1,2,3,4,11,12,13,14,5,6,7,8,9,10)]
#data=log(data+1)
pdf("Heatmap.pdf",width=2.5,height=3.6)
matrix=normalize.quantiles(as.matrix(data[4:17]))
#matrix=as.matrix(t(scale(t(matrix),scale=T,center=T)))
#matrix=as.matrix(t(scale(t(matrix),scale=T,center=T)))
pheatmap(data,scale='row',cluster_cols=F,cutree_rows=4,color = colorRampPalette(colors = c("blue","white","red"))(100))
dev.off()


matrix=as.matrix(t(scale(t(data),scale=T,center=T)))
matrix=normalize.quantiles(as.matrix(matrix))
Heatmap(matrix,cluster_columns=F,row_km=4)


dev.off()















pheatmap(matrix)

matrix=as.matrix(t(scale(t(data),scale=T,center=F)))
matrix=normalize.quantiles(matrix)
Heatmap(matrix,cluster_columns=F)

matrix=as.matrix(t(scale(t(data),scale=F,center=T)))
matrix=normalize.quantiles(matrix)
Heatmap(matrix,cluster_columns=F)




matrix=normalize.quantiles(matrix)
matrix=as.matrix(t(scale(t(matrix),scale=T,center=T)))
Heatmap(matrix,cluster_columns=F)



matrix=normalize.quantiles(matrix)
matrix=as.matrix(t(scale(t(matrix),scale=T,center=F)))
Heatmap(matrix,cluster_columns=F)



matrix=normalize.quantiles(matrix)
matrix=as.matrix(t(scale(t(matrix),scale=F,center=T)))
Heatmap(matrix,cluster_columns=F)


dev.off()




library("ComplexHeatmap")
library("pheatmap")                                                                                                                                                            
library("preprocessCore")
countsTable <- read.table( "/md01/changjy/data/quail/peak/all.sorted.merged.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
#length=c(32980464,36105232,25169816,44620304,34433734,28843022,50665778,31530834,34895780,45266950,37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*40000000
countsTable=normalize.quantiles(as.matrix(countsTable))
countsTable=countsTable[,c(1,2,3,5,6,7,8)]
colData <- data.frame(condition=factor(c("control","control","control","ld3","ld3","ld3","ld3")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds)
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"control_3d_DEP.txt", row.name=T,sep="\t",quote=F)


countsTable <- read.table( "/md01/changjy/data/quail/peak/all.sorted.merged.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*40000000
countsTable=countsTable[,c(1,2,3,4,9,10,11)]
colData <- data.frame(condition=factor(c("control","control","control","control","ld7","ld7","ld7")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds)
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"control_7d_DEP.txt", row.name=T,sep="\t",quote=F)


countsTable <- read.table( "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/Union_peak.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*40000000
countsTable=countsTable[,c(1,2,3,4,12,13,14)]
colData <- data.frame(condition=factor(c("control","control","control","control","ld28","ld28","ld28")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds)
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"control_28d_DEP.txt", row.name=T,sep="\t",quote=F)


countsTable <- read.table( "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/Union_peak.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*40000000
countsTable=countsTable[,c(5,6,7,8,9,10,11)]
colData <- data.frame(condition=factor(c("control","control","control","control","ld28","ld28","ld28")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds)
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"3d_7d_DEP.txt", row.name=T,sep="\t",quote=F)



countsTable <- read.table( "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/Union_peak.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*40000000
countsTable=countsTable[,c(5,6,7,8,12,13,14)]
colData <- data.frame(condition=factor(c("control","control","control","control","ld28","ld28","ld28")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds)
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"3d_28d_DEP.txt", row.name=T,sep="\t",quote=F)



countsTable <- read.table( "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/Union_peak.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*40000000
countsTable=countsTable[,c(9,10,11,12,13,14)]
colData <- data.frame(condition=factor(c("control","control","control","control","ld28","ld28","ld28")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds)
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1& (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"7d_28d_DEP.txt", row.name=T,sep="\t",quote=F)

countsTable <- read.table( "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/Union_peak.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*40000000
a=read.table("control_3d_DEP.txt")
b=read.table("control_7d_DEP.txt")
c=read.table("control_28d_DEP.txt")
d=read.table("3d_7d_DEP.txt")
e=read.table("3d_28d_DEP.txt")
f=read.table("7d_28d_DEP.txt")
aa=countsTable[rownames(a),]
bb=countsTable[rownames(a),]
cc=countsTable[rownames(a),]
dd=countsTable[rownames(a),]
ee=countsTable[rownames(a),]
ff=countsTable[rownames(a),]
peak=rbind(aa,bb)
peak=rbind(peak,cc)
peak=rbind(peak,dd)
peak=rbind(peak,ee)
peak=rbind(peak,ff)
pdf("quantile.Heatmap.pdf")
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
peak=t(scale(t(peak),center=T,scale=T))
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pheatmap(peak,cluster_cols=F,show_rownames=F)
peak=t(scale(t(peak),center=F,scale=T))
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pheatmap(peak,cluster_cols=F,show_rownames=F)
peak=t(scale(t(peak),center=T,scale=F))
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pheatmap(peak,cluster_cols=F,show_rownames=F)
peak=t(scale(t(peak),center=F,scale=F))
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pheatmap(peak,cluster_cols=F,show_rownames=F)
pheatmap(peak,cluster_cols=F,show_rownames=F)
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
peak=read.table("all.sorted.merged.DEP.count")
peak=peak[c(4,5,6,8,9,10,11,12,13,14,15,16,17)]/length*10000000
pdf("test.pdf")
colnames(peak)<-c("control_1","control_2","control_3","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

peak2=log2(peak+1)
colnames(peak2)<-c("control_1","control_2","control_3","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

pheatmap(peak2,cluster_cols=F,show_rownames=F)
rowmean=apply(peak,1,mean)
peak3=log2((peak+1)/(rowmean+1))
colnames(peak3)<-c("control_1","control_2","control_3","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

pheatmap(peak3,cluster_cols=F,show_rownames=F)
peak4=t(scale(t(peak)))
colnames(peak4)<-c("control_1","control_2","control_3","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

pheatmap(peak4,cluster_cols=F,show_rownames=F)

peak[peak>100]=100
pheatmap(peak,cluster_cols=F,show_rownames=F)
dev.off()





bedtools multicov -bed 
no_shift_DiffBind_Deseq2_peakAnno_all_sorted_merged.txt 
-bams /md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M1-1_L3_801D73.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-1_L2_801F10.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-2_L1_802F10.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M1-1_L1_807D75.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-1_L1_803F10.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-2_L1_806F10.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M5-2_L1_806D74.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M1-2_L1_808F10.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M3-2_L1_808D73.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M6-2_L1_808D74.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M1-1_L1_809D73.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M3-1_L1_809D74.unchrM.pe.q10.sort.rmdup.bam 
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M6-1_L2_801D75.unchrM.pe.q10.sort.rmdup.bam >no_shift_DiffBind_Deseq2_peakAnno_all_sorted_merged_count.txt

data=read.table("/md01/changjy/data/quail/next_1.0/diffbind/no_shift/no_shift_DiffBind_Deseq2_peakAnno_all_sorted_merged_count.txt")
data=data[4:17]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
data1=data/length*40000000
data2=log2(data1+1)
rowmean=apply(data1,1,mean)
data3=log2((data1+1/rowmean+1))
data4=t(scale(t(data1)))
colnames(data1)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

colnames(data2)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

colnames(data3)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

colnames(data4)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pdf("Heatmap.pdf")
Heatmap(data1,cluster_columns=F,cluster_rows=F,column_title="Raw count / length")
Heatmap(data1,cluster_columns=F,cluster_rows=T,column_title="Raw count / length")
Heatmap(data2,cluster_columns=F,cluster_rows=F,column_title="Raw count / length log2")
Heatmap(data2,cluster_columns=F,cluster_rows=T,column_title="Raw count / length log2")
Heatmap(data3,cluster_columns=F,cluster_rows=F,column_title="Raw count / length log2FC")
Heatmap(data3,cluster_columns=F,cluster_rows=T,column_title="Raw count / length log2FC")
Heatmap(data4,cluster_columns=F,cluster_rows=F,column_title="Raw count / length z-score")
Heatmap(data4,cluster_columns=F,cluster_rows=T,column_title="Raw count / length z-score")
dev.off()



countsTable <- read.table( "/md01/changjy/data/quail/next_1.0/diffbind/no_shift/Union_peak.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*40000000
a=read.table("control_3d_DEP.txt")
b=read.table("control_7d_DEP.txt")
c=read.table("control_28d_DEP.txt")
d=read.table("3d_7d_DEP.txt")
e=read.table("3d_28d_DEP.txt")
f=read.table("7d_28d_DEP.txt")
aa=countsTable[rownames(a),]
bb=countsTable[rownames(a),]
cc=countsTable[rownames(a),]
dd=countsTable[rownames(a),]
ee=countsTable[rownames(a),]
ff=countsTable[rownames(a),]
peak=rbind(aa,bb)
peak=rbind(peak,cc)
peak=rbind(peak,dd)
peak=rbind(peak,ee)
peak=rbind(peak,ff)
pdf("quantile.Heatmap.pdf")
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
peak=t(scale(t(peak),center=T,scale=T))
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pheatmap(peak,cluster_cols=F,show_rownames=F)
peak=t(scale(t(peak),center=F,scale=T))
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pheatmap(peak,cluster_cols=F,show_rownames=F)
peak=t(scale(t(peak),center=T,scale=F))
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pheatmap(peak,cluster_cols=F,show_rownames=F)
peak=t(scale(t(peak),center=F,scale=F))
colnames(peak)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pheatmap(peak,cluster_cols=F,show_rownames=F)
pheatmap(peak,cluster_cols=F,show_rownames=F)
pdf("test.pdf")
peak2=log2(peak+1)
colnames(peak2)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

pheatmap(peak2,cluster_cols=F,show_rownames=F)
rowmean=apply(peak,1,mean)
peak3=log2((peak/rowmean)+1)
colnames(peak3)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

pheatmap(peak3,cluster_cols=F,show_rownames=F)

dev.off()



##DiffBind找差异peak
library(DiffBind)
library(preprocessCore)
library(pheatmap)
library(ComplexHeatmap)
getwd()
setwd("/md01/changjy/data/quail/peak/peak/")
##导入基本信息
sample_info <-dba(sampleSheet="/md01/changjy/data/quail/peak/peak/sample_info.csv")
##count计数
dba_count=dba.count(sample_info,bUseSummarizeOverlaps=TRUE,bParallel=TRUE)

##构建差异组
dba_contrast=dba.contrast(dba_count, categories=DBA_CONDITION)
##差异分析
dba_diff=dba.analyze(dba_contrast)
##保存差异结果
dba_report_all_1 <- dba.report(dba_diff,contrast=1,method=DBA_DESEQ2)
dba_report_all_2 <- dba.report(dba_diff,contrast=2,method=DBA_DESEQ2)
dba_report_all_3 <- dba.report(dba_diff,contrast=3,method=DBA_DESEQ2)
dba_report_all_4 <- dba.report(dba_diff,contrast=4,method=DBA_DESEQ2)
dba_report_all_5 <- dba.report(dba_diff,contrast=5,method=DBA_DESEQ2)
dba_report_all_6 <- dba.report(dba_diff,contrast=6,method=DBA_DESEQ2)
write.table(dba_report_all_1,"control_3d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_2,"control_7d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_3,"control_28d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_4,"3d_7d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_5,"3d_28d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)
write.table(dba_report_all_6,"7d_28d_DEP.csv",row.names=T,sep="\t",quote=F,col.names=F)

##合并六组差异peak
cat *.csv|awk '{print $2"\t"$3"\t"$4}'|sort -k1,1 -k2,2n > DiffBind_DEP.sorted.csv
##peak计数
bedtools multicov -bed DiffBind_DEP.sorted.csv -bams \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M1-1_L3_801D73.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-1_L2_801F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-2_L1_802F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M1-1_L1_807D75.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-1_L1_803F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-2_L1_806F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M5-2_L1_806D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M1-2_L1_808F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M3-2_L1_808D73.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M6-2_L1_808D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M1-1_L1_809D73.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M3-1_L1_809D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M6-1_L2_801D75.unchrM.pe.q10.sort.rmdup.bam \
> DiffBind.sorted.DEP.count

##第二个bam有问题
bedtools multicov -bed DiffBind_DEP.sorted.csv \
-bams /md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bam \
>tmp

##merge
bedtools merge -i DiffBind_DEP.sorted.csv >DiffBind_DEP.sorted.merged.csv
##merge后计数
bedtools multicov -bed DiffBind_DEP.sorted.merged.csv -bams \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M1-1_L3_801D73.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-1_L2_801F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-2_L1_802F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M1-1_L1_807D75.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-1_L1_803F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-2_L1_806F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M5-2_L1_806D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M1-2_L1_808F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M3-2_L1_808D73.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M6-2_L1_808D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M1-1_L1_809D73.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M3-1_L1_809D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M6-1_L2_801D75.unchrM.pe.q10.sort.rmdup.bam \
> DiffBind.sorted.merged.DEP.count

##第二个bam有问题
bedtools multicov -bed DiffBind_DEP.sorted.merged.csv \
-bams /md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bam \
>tmp

##DEseq2找差异peak
dba_count$binding[dba_count$binding[,1]==1,][,1]="NC_029516.1"
dba_count$binding[dba_count$binding[,1]==2,][,1]="NC_029517.1"
dba_count$binding[dba_count$binding[,1]==3,][,1]="NC_029518.1"
dba_count$binding[dba_count$binding[,1]==4,][,1]="NC_029519.1"
dba_count$binding[dba_count$binding[,1]==5,][,1]="NC_029520.1"
dba_count$binding[dba_count$binding[,1]==6,][,1]="NC_029521.1"
dba_count$binding[dba_count$binding[,1]==7,][,1]="NC_029522.1"
dba_count$binding[dba_count$binding[,1]==8,][,1]="NC_029523.1"
dba_count$binding[dba_count$binding[,1]==9,][,1]="NC_029524.1"
dba_count$binding[dba_count$binding[,1]==10,][,1]="NC_029525.1"
dba_count$binding[dba_count$binding[,1]==11,][,1]="NC_029526.1"
dba_count$binding[dba_count$binding[,1]==12,][,1]="NC_029527.1"
dba_count$binding[dba_count$binding[,1]==13,][,1]="NC_029528.1"
dba_count$binding[dba_count$binding[,1]==14,][,1]="NC_029529.1"
dba_count$binding[dba_count$binding[,1]==15,][,1]="NC_029530.1"
dba_count$binding[dba_count$binding[,1]==16,][,1]="NC_029531.1"
dba_count$binding[dba_count$binding[,1]==17,][,1]="NC_029532.1"
dba_count$binding[dba_count$binding[,1]==18,][,1]="NC_029533.1"
dba_count$binding[dba_count$binding[,1]==19,][,1]="NC_029534.1"
dba_count$binding[dba_count$binding[,1]==20,][,1]="NC_029535.1"
dba_count$binding[dba_count$binding[,1]==21,][,1]="NC_029536.1"
dba_count$binding[dba_count$binding[,1]==22,][,1]="NC_029537.1"
dba_count$binding[dba_count$binding[,1]==23,][,1]="NC_029538.1"
dba_count$binding[dba_count$binding[,1]==24,][,1]="NC_029539.1"
dba_count$binding[dba_count$binding[,1]==25,][,1]="NC_029540.1"
dba_count$binding[dba_count$binding[,1]==26,][,1]="NC_029541.1"
dba_count$binding[dba_count$binding[,1]==27,][,1]="NC_029542.1"
dba_count$binding[dba_count$binding[,1]==28,][,1]="NC_029543.1"
dba_count$binding[dba_count$binding[,1]==29,][,1]="NC_029544.1"
dba_count$binding[dba_count$binding[,1]==30,][,1]="NC_029545.1"
dba_count$binding[dba_count$binding[,1]==31,][,1]="NC_029547.1"
##保存count矩阵
write.table(dba_count$binding,"dba_count.txt",sep="\t",quote=F,row.names=F,col.names=F)
##调整count矩阵格式
awk '{print $1"\t"$2"\t"$3}' dba_count.txt|sort -k1,1 -k2,2n> dba_count.bed

##矩阵计数
bedtools multicov -bed dba_count.bed -bams \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M1-1_L3_801D73.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-1_L2_801F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P0-M6-2_L1_802F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M1-1_L1_807D75.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-1_L1_803F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M2-2_L1_806F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P8-M5-2_L1_806D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M1-2_L1_808F10.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M3-2_L1_808D73.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P10-M6-2_L1_808D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M1-1_L1_809D73.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M3-1_L1_809D74.unchrM.pe.q10.sort.rmdup.bam \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/SD-P12-M6-1_L2_801D75.unchrM.pe.q10.sort.rmdup.bam \
> dba_count.bed.count

##第二个bam有问题
bedtools multicov -bed dba_count.bed \
-bams /md01/changjy/data/quail/mapping/sort.rmdup.bam/LD-P12-M2-1_L2_801D74.unchrM.pe.q10.sort.rmdup.bam \
>tmp

library("ComplexHeatmap")
library("pheatmap")                                                                                                                                                            
library("preprocessCore")
library(DESeq2)
###control vs ld_3 差异分析
countsTable <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
#length=c(32980464,36105232,25169816,44620304,34433734,28843022,50665778,31530834,34895780,45266950,37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*mean(length)
countsTable=normalize.quantiles(as.matrix(countsTable))
countsTable=countsTable[,c(1,2,3,4,5,6,7,8)]
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
colData <- data.frame(condition=factor(c("control","control","control","control","ld_3","ld_3","ld_3","ld_3")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds,contrast=c("condition","ld_3","control"))
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"test/control_3d_DEP.txt", row.name=T,sep="\t",quote=F)


###control vs ld_7 差异分析
countsTable <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
#length=c(32980464,36105232,25169816,44620304,34433734,28843022,50665778,31530834,34895780,45266950,37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*mean(length)
countsTable=normalize.quantiles(as.matrix(countsTable))
countsTable=countsTable[,c(1,2,3,4,9,10,11)]
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
colData <- data.frame(condition=factor(c("control","control","control","control","ld_7","ld_7","ld_7")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds,contrast=c("condition","ld_7","control"))
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"test/control_7d_DEP.txt", row.name=T,sep="\t",quote=F)


###control vs ld_28 差异分析
countsTable <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
#length=c(32980464,36105232,25169816,44620304,34433734,28843022,50665778,31530834,34895780,45266950,37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*mean(length)
countsTable=normalize.quantiles(as.matrix(countsTable))
countsTable=countsTable[,c(1,2,3,4,12,13,14)]
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
colData <- data.frame(condition=factor(c("control","control","control","control","ld_28","ld_28","ld_28")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds,contrast=c("condition","ld_28","control"))
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"test/control_28d_DEP.txt", row.name=T,sep="\t",quote=F)



###ld_3 vs ld_7 差异分析
countsTable <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
#length=c(32980464,36105232,25169816,44620304,34433734,28843022,50665778,31530834,34895780,45266950,37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*mean(length)
countsTable=normalize.quantiles(as.matrix(countsTable))
countsTable=countsTable[,c(5,6,7,8,9,10,11)]
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
colData <- data.frame(condition=factor(c("ld_3","ld_3","ld_3","ld_3","ld_7","ld_7","ld_7")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds,contrast=c("condition","ld_7","ld_3"))
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"test/ld3_7d_DEP.txt", row.name=T,sep="\t",quote=F)

###ld_3 vs ld_28 差异分析
countsTable <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
#length=c(32980464,36105232,25169816,44620304,34433734,28843022,50665778,31530834,34895780,45266950,37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*mean(length)
countsTable=normalize.quantiles(as.matrix(countsTable))
countsTable=countsTable[,c(5,6,7,8,12,13,14)]
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
colData <- data.frame(condition=factor(c("ld_3","ld_3","ld_3","ld_3","ld_28","ld_28","ld_28")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds,contrast=c("condition","ld_28","ld_3"))
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"test/ld3_ld28_DEP.txt", row.name=T,sep="\t",quote=F)

###ld_7 vs ld_28 差异分析
countsTable <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed.count", stringsAsFactors=TRUE )
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
#peak<-countsTable[,4]
#countsTable <- countsTable[ , 8:20]
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
#length=c(32980464,36105232,25169816,44620304,34433734,28843022,50665778,31530834,34895780,45266950,37775976,36153760,41528598)
countsTable=countsTable[4:17]/length*mean(length)
countsTable=normalize.quantiles(as.matrix(countsTable))
countsTable=countsTable[,c(9,10,11,12,13,14)]
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
colData <- data.frame(condition=factor(c("ld_7","ld_7","ld_7","ld_28","ld_28","ld_28")))
countsTable=round(countsTable)
dds<-DESeqDataSetFromMatrix(countsTable,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
res=results(dds,contrast=c("condition","ld_28","ld_7"))
summary(res)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.1 & (log2FoldChange > 1 | log2FoldChange < -1))
write.table(as.data.frame(diff_gene_deseq2),"test/ld7_ld28_DEP.txt", row.name=T,sep="\t",quote=F)

#======================================================================
#======================================================================
#======================================================================




###根据行名获得矩阵
countsTable <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed.count", stringsAsFactors=TRUE )
countsTable=countsTable[4:17]/length*mean(length)
countsTable=normalize.quantiles(as.matrix(countsTable))
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
###输入文库长度
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)

###读取数据
df1=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/control_3d_DEP.txt")
df2=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/control_7d_DEP.txt")
df3=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/control_28d_DEP.txt")
df4=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/ld3_7d_DEP.txt")
df5=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/ld3_ld28_DEP.txt")
df6=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/ld7_ld28_DEP.txt")
##按行名提取count
df1=countsTable[rownames(df1),]
df2=countsTable[rownames(df2),]
df3=countsTable[rownames(df3),]
df4=countsTable[rownames(df4),]
df5=countsTable[rownames(df5),]
df6=countsTable[rownames(df6),]

##合并矩阵
pic=rbind(df1,df2)
pic=rbind(pic,df3)
pic=rbind(pic,df4)
pic=rbind(pic,df5)
pic=rbind(pic,df6)

colnames(pic)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")

##画图
#pic=unique_pic
pic <- pic[!duplicated(rownames(pic)),]
pdf("Heatmap.pdf")
Heatmap(pic,,show_row_names=F,row_title="pic")
Heatmap(log2(pic+1),,show_row_names=F,row_title="log2(pic+1)")
Heatmap(log10(pic+1),,show_row_names=F,row_title="log10(pic+1)")
rowmean=apply(pic,1,mean)
Heatmap(log2((pic+1)/(rowmean+1)),show_row_names=F,row_title="log2((pic+1)/(rowmean+1))")
Heatmap(t(scale(t(pic))),show_row_names=F,row_title="t(scale(t(pic)))")
dev.off()

##画图
pdf("Heatmap_nocolumns.pdf")
Heatmap(pic,,show_row_names=F,row_title="pic",,cluster_columns=F)
Heatmap(log2(pic+1),,show_row_names=F,row_title="log2(pic+1)",cluster_columns=F)
Heatmap(log10(pic+1),,show_row_names=F,row_title="log10(pic+1)",cluster_columns=F)
rowmean=apply(pic,1,mean)
Heatmap(log2((pic+1)/(rowmean+1)),show_row_names=F,row_title="log2((pic+1)/(rowmean+1))",cluster_columns=F)
Heatmap(t(scale(t(pic))),show_row_names=F,row_title="t(scale(t(pic)))",cluster_columns=F)
dev.off()
=======================================================================================================================
##画图
library(circlize)
library("ComplexHeatmap")
library("pheatmap")                                                                                                                                                            
library("preprocessCore")
library(DESeq2)
###根据行名获得矩阵
countsTable <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed.count", stringsAsFactors=TRUE )
countsTable=countsTable[4:17]/length*mean(length)
countsTable=normalize.quantiles(as.matrix(countsTable))
rownames( countsTable ) <- paste(rep("peak_",dim(countsTable)[1]),c(1:dim(countsTable)[1]),sep="")
###输入文库长度
length=c(32980464,36105232,25169816,46231024,
    44620304,34433734,28843022,50665778,
    31530834,34895780,45266950,
    37775976,36153760,41528598)
###读取数据
df1=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/control_3d_DEP.txt")
df2=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/control_7d_DEP.txt")
df3=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/control_28d_DEP.txt")
df4=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/ld3_7d_DEP.txt")
df5=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/ld3_ld28_DEP.txt")
df6=read.table("/md01/changjy/data/quail/peak/DESeq2/test/Differential/ld7_ld28_DEP.txt")
##按行名提取count
df1=countsTable[rownames(df1),]
df2=countsTable[rownames(df2),]
df3=countsTable[rownames(df3),]
df4=countsTable[rownames(df4),]
df5=countsTable[rownames(df5),]
df6=countsTable[rownames(df6),]
##合并矩阵
pic=rbind(df1,df2)
pic=rbind(pic,df3)
pic=rbind(pic,df4)
pic=rbind(pic,df5)
pic=rbind(pic,df6)
colnames(pic)<-c("control_1","control_2","control_3","control_4","3d_1","3d_2","3d_3","3d_4","7d_1","7d_2","7d_3","28d_1","28d_2","28d_3")
pdf("Heatmap_nocolumns_cluster_9:16.pdf")
# Heatmap(pic,,show_row_names=F,row_title="pic",,cluster_columns=F,row_km=4)
# Heatmap(log2(pic+1),,show_row_names=F,row_title="log2(pic+1)",cluster_columns=F,row_km=4)
# Heatmap(log10(pic+1),,show_row_names=F,row_title="log10(pic+1)",cluster_columns=F,row_km=4)
rowmean=apply(pic,1,mean)
set.seed(222)
cluster_info=c("control","control","control","control","ld3","ld3","ld3","ld3","ld7","ld7","ld7","ld28","ld28","ld28")


col <- jdb_color_maps[1:4]
names(col) <- levels(cluster_info)

top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col), # 设置填充色
                       #labels = levels(cluster_info), 
                       labels_gp = gpar(cex = 0.5, col = "white")))

col_fun <- colorRamp2(
  c(-2, 0, 2), 
  c( "#FFA500","#AFEEEE","#FF4500")
  )
pdf("/md01/changjy/data/quail/peak/Heatmap.pdf",height=8,width=4.5)
set.seed(222)

Heatmap(log2((pic+1)/(rowmean+1)),
    #col=col_fun,
    show_row_names=F,
    #row_title="log2((pic+1)/(rowmean+1))",
    cluster_columns=F,
    column_names_rot = 60, 
    row_km=4,
    column_split = cluster_info,
    column_title = NULL ,
    top_annotation = top_anno,
    show_row_dend=F,
    cluster_rows=F,
    cluster_col=F,
    column_names_gp= gpar(
    col = "#1f78b4",
    fontsize = 12),
    heatmap_legend_param = list(title = "log2(FoldChange)",title_position = "leftcenter-rot"))

dev.off()
===============================================================================================================================


Heatmap(t(scale(t(pic))),show_row_names=F,row_title="t(scale(t(pic)))",cluster_columns=F,row_km=4)
dev.off()

###提取聚类信息
set.seed(222)
list<-Heatmap(log2((pic+1)/(rowmean+1)),
    #col=col_fun,
    show_row_names=F,
    #row_title="log2((pic+1)/(rowmean+1))",
    cluster_columns=F,
    column_names_rot = 60, 
    row_km=4,
    column_split = cluster_info,
    column_title = NULL ,
    top_annotation = top_anno,
    show_row_dend=F,
    cluster_rows=F,
    cluster_col=F,
    column_names_gp= gpar(
    col = "#1f78b4",
    fontsize = 12),
    heatmap_legend_param = list(title = "log2(FoldChange)",title_position = "leftcenter-rot"))
##提取聚类信息
row_order(list)
column_order(list)


##聚类信息赋值
cluster1 = pic[row_order(list)[[1]],column_order(list)] #提取第一簇的相关信息
cluster2 = pic[row_order(list)[[2]],column_order(list)] #同上
cluster3 = pic[row_order(list)[[3]],column_order(list)] #同上
cluster4 = pic[row_order(list)[[4]],column_order(list)] #同上

##查看各个类的数量
dim(cluster1)
dim(cluster2)
dim(cluster3)
dim(cluster4)

bed <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed.count", stringsAsFactors=TRUE )[,c(1,2,3)]
#bed <- read.table( "/md01/changjy/data/quail/peak/DESeq2/dba_count.bed", stringsAsFactors=TRUE )
rownames( bed ) <- paste(rep("peak_",dim( countsTable )[1]),c(1:dim( countsTable )[1]),sep="")
#data=merge(pic,bed,by='row.names')
cluster1=bed[rownames(bed)%in%rownames(cluster1),][,c(1,2,3)]
cluster2=bed[rownames(bed)%in%rownames(cluster2),][,c(1,2,3)]
cluster3=bed[rownames(bed)%in%rownames(cluster3),][,c(1,2,3)]
cluster4=bed[rownames(bed)%in%rownames(cluster4),][,c(1,2,3)]

# ##提取行名
# rownames(cluster1)<-row_order(list)[[1]]
# rownames(cluster2)<-row_order(list)[[2]]
# rownames(cluster3)<-row_order(list)[[3]]
# rownames(cluster4)<-row_order(list)[[4]]

###保存聚类结果
write.table(cluster1,"test/cluster1.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(cluster2,"test/cluster2.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(cluster3,"test/cluster3.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(cluster4,"test/cluster4.txt",sep="\t",quote=F,row.names=F,col.names=F)



####GO KEGG分析

# 导入包
library(AnnotationHub)
# 创建AnnotationHub实例
ah <- AnnotationHub()
# 获取物种信息，以鹌鹑为例
ah[ah$species=="Coturnix japonica"]
#ah[ah$species=="Japanese Quail"]
ah[ah$species=="Coturnix japonica" & ah$rdataclass=="OrgDb"]
##保存orgdb数据库
Coturnix_japonica_orgdb=ah[["AH86236"]]
saveDb(Coturnix_japonica_orgdb, file="./Coturnix_japonica.orgdb")


##GO
library("clusterProfiler")
# 读取基因名称数据
gene_symbol=read.table("/md01/changjy/data/quail/RNAseq/control_treat.DESeq2.select.txt",header=T)
gene=rownames(gene_symbol)

##基因名转换
gene_trans = bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="Coturnix_japonica_orgdb")

##基因kegg名转换
gene_trans_kegg = bitr(gene,fromType="SYMBOL",toType="GID",OrgDb="Coturnix_japonica_orgdb")

ego <- enrichGO(
          gene  = gene_trans$ENTREZID,
          keyType = "ENTREZID", 
          #universe = names(geneList), #背景基因集，可省
          OrgDb   = Coturnix_japonica_orgdb,
          ont     = "ALL",
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.01,
          qvalueCutoff  = 0.05,
          readable      = TRUE)

ego <- enrichKEGG(
          gene = gene_trans$SYMBOL,
          keyType = "kegg",
          organism  = 'hsa',
          pvalueCutoff  = 0.05,
          pAdjustMethod  = "BH",
          qvalueCutoff  = 0.05
)

##HINT-ATAC
## installation 
condaup
conda activate python3
pip install --user RGT # RGT only supports Python 3 

## Configuration of Genomic Data
cd /public/home/shipy3/DB/dm6/annotation
awk '$3=="gene"' edited_chrom_dmel-all-r6.33.gtf >edited_chrom_dmel-all-r6.33_genes.gtf
gff2bed < /public/home/shipy3/DB/dm6/annotation/edited_chrom_dmel-all-r6.33_genes.gtf > /public/home/shipy3/DB/dm6/annotation/edited_chrom_dmel-all-r6.33_genes.bed
cut -f 1,2,3,5,6 edited_chrom_dmel-all-r6.33_genes.bed > edited_chrom_dmel-all-r6.33_genes.bed.temp1
cut -f 10 edited_chrom_dmel-all-r6.33_genes.bed | cut -d ";" -f 2 | cut -d '"' -f 2 > edited_chrom_dmel-all-r6.33_genes.bed.temp2
paste -d "\t" edited_chrom_dmel-all-r6.33_genes.bed.temp1 edited_chrom_dmel-all-r6.33_genes.bed.temp2 > edited_chrom_dmel-all-r6.33_genes.bed
cut -f 9 edited_chrom_dmel-all-r6.33_genes.gtf | cut -d ";" -f 1 | cut -d '"' -f 2 > alias_dm6.txt.temp1
cut -f 9 edited_chrom_dmel-all-r6.33_genes.gtf | cut -d ";" -f 2 | cut -d '"' -f 2 > alias_dm6.txt.temp2
paste -d "&" alias_dm6.txt.temp1 alias_dm6.txt.temp2 > alias_dm6.txt.temp3
paste -d "\t" alias_dm6.txt.temp1 alias_dm6.txt.temp2 alias_dm6.txt.temp3 > alias_dm6.txt
cd ~/rgtdata
vim data.config.user # a customized genome for drosophila

###HINT-ATAC footprint
ls /md01/changjy/data/quail/mapping/sort.rmdup.bam/*.bam|while read bam
do
rgt-hint footprinting --organism cj22 \
--atac-seq \
$bam \
/md01/changjy/data/quail/peak/DESeq2/test/cluster/merge_cluster.bed \
--output-location=/md01/changjy/data/quail/HINT-ATAC/split \
--output-prefix=${bam%.bam*} \
--paired-end
done





rgt-hint footprinting --organism cj22 \
--atac-seq \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/merge/ld3.bam \
/md01/changjy/data/quail/HINT-ATAC/test/dba_sorted_merge_count_mpbs.bed \
--output-location=/md01/changjy/data/quail/HINT-ATAC \
--output-prefix=ld3 \
--pair-end

rgt-hint footprinting --organism cj22 \
--atac-seq \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/merge/ld7.bam \
/md01/changjy/data/quail/HINT-ATAC/test/dba_sorted_merge_count_mpbs.bed \
--output-location=/md01/changjy/data/quail/HINT-ATAC \
--output-prefix=ld7 \
--pair-end

rgt-hint footprinting --organism cj22 \
--atac-seq \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/merge/ld28.bam \
/md01/changjy/data/quail/HINT-ATAC/test/dba_sorted_merge_count_mpbs.bed \
--output-location=/md01/changjy/data/quail/HINT-ATAC \
--output-prefix=ld28 \
--pair-end


rgt-hint footprinting --organism cj22 \
--atac-seq \
/md01/changjy/data/quail/mapping/sort.rmdup.bam/merge/control.bam \
/md01/changjy/data/quail/HINT-ATAC/test/dba_sorted_merge_count_mpbs.bed \
--output-location=/md01/changjy/data/quail/HINT-ATAC \
--output-prefix=control \
--pair-end



TOBIAS BINDetect --motifs known10.motif \
--signals /md01/changjy/data/quail/TOBIAS/ScoreBigwig/control_footprints.bw \
/md01/changjy/data/quail/TOBIAS/ScoreBigwig/ld3_footprints.bw \
/md01/changjy/data/quail/TOBIAS/ScoreBigwig/ld7_footprints.bw \
/md01/changjy/data/quail/TOBIAS/ScoreBigwig/ld28_footprints.bw \
--genome /md01/changjy/data/quail/ref/ncbi_dataset/index/quail.fa \
--peaks cluster_1_known10.motif.position.motif.find \
--outdir RFX2 \
--cond_names control ld3 ld7 ld28 \
--cores 8



data=read.table("C:/Users/Jeff/OneDrive/桌面/quail/ref/gtf/quail.gtf",sep="\t",)
data=read.csv("C:/Users/Jeff/OneDrive/桌面/quail/motif/cluster1/knownResults.txt",sep="\t")
y=read.csv("C:/Users/Jeff/OneDrive/桌面/quail/motif.csv")
library(ggplot2)
ggplot(y,aes(cluster,toupper(name)))+
  geom_point(aes(size=enrichment,color=-log10(pvalue)))+
  scale_color_gradientn(colors = c(     "#006400", "#FC4E07")) +
  theme_minimal()+
  xlab(NULL)+
  ylab(NULL)+
  scale_size_continuous(range=c(0,8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=6,angle=-30),axis.text.y=element_text(size=6))
ggsave('C:/Users/Jeff/OneDrive/桌面/quail/plot/Motif_bubble.pdf', p, width = 5, height = 6)



library(ggplot2)
data<-read.csv("/md01/changjy/data/quail/RNAseq/control_vs_treatment.new.csv",header=T,row.names=1)



data[which(data$log2FoldChange >= 1 & data$padj < 0.05),'sig'] <- 'UP'
data[which(data$log2FoldChange <= -1 & data$padj< 0.05),'sig'] <- 'DOWN'
data[which(abs(data$log2FoldChange) < 1 | data$padj > 0.05),'sig'] <- 'NO DIFF'


data$label <-''
data[data$padj < 0.01 & data$log2FoldChange >= 2,]
data[rownames(data)%in%c("DIO3","GPR20","TSHB","CGA","SLC16A2","GHRH","DIO2","POMC","PER2","PER3"),]$label=c("DIO3","GPR20","TSHB","CGA","SLC16A2","GHRH","DIO2","POMC","PER2","PER3")


this_tile <- paste0("Control vs LD7",
    '\ncutoff for abs(logFC) and FDR is 1 and 0.05',
    '\nThe number of up gene is ',nrow(data[data$sig =='UP',]) ,
    '\nThe number of down gene is ',nrow(data[data$sig =='DOWN',]))
volcano <- ggplot(data, aes(log2FoldChange, -log(padj, 10))) +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=10,hjust = 0.5))+
  geom_point(aes(color = sig), alpha = 0.8, size = 2) +
  scale_color_manual(values = c('#000080', '#C0C0C0', '#8B0000')) +
  theme(panel.grid = element_blank(), 
      panel.background = element_rect(color = 'black', fill = 'transparent'), 
      ) +
  #geom_text(aes(label = label), size = 3,vjust=-0.5, alpha=0.8) +  
  theme(legend.title = element_blank(), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 padj', color = NA) +
  xlim(-8, 8)

 p=volcano+ geom_text_repel(data = data, aes(x = data$log2FoldChange, 
                                      y = -log10(data$padj), 
                                      label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black")


ggsave('/md01/changjy/data/quail/RNAseq/Volcano.pdf', p, width = 5, height = 6)


pdf('/md01/changjy/data/quail/RNAseq/Venn_ATAC_RNA.pdf')
venn.plot <- draw.pairwise.venn(  area1 = 931,  #区域1的数 
area2 =137 , #区域2的数 
cross.area = 7,  #重叠的个数 
category = c("ATAC-seq", "RNA-seq"),#分类命名
fill = c("#00AFBB", "#FC4E07"),#1 2 区域分别的填充颜色 
lty = "blank",  #1 2 区域的边框线类型 
cex = 2,        #1 2 区域内部数字的字体大小 
cat.cex = 2,    # 分类名称的字体大小 
cat.dist = 0.09,   #分类名称距离边的距离 实际调整 
cat.just = list(c(-1, -1), c(1, 1)),  #分类名称的位置  ，圈内或者圈外
ext.pos = 30,  #线的角度 默认是正上方12点位置 
ext.dist = -0.05,   #外部线的距离  跟根据圆圈的大小适当调整
ext.length = 0.85,  #外部线长度 
ext.line.lwd = 2,  #外部线的宽度 
ext.line.lty = "dashed" )  #外部线为虚线);
grid.draw(venn.plot)
dev.off()


data=read.csv("/md01/changjy/data/quail/peak/ATAC_KEGG.txt",sep="\t")
p <- ggplot(data[c(1:10),], aes(x=Term, y=Input.number, fill=-log10(P.Value))) + geom_bar(stat="identity") +
    theme_bw()+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +

    ylab("Gene number") +
    xlab("pathway name") +
     theme(panel.grid = element_blank())+
    #labs(fill=-log10(P.Value)) +
    scale_fill_gradient(low="#C1FFC1", high="#228B22")+
    coord_flip()
p
ggsave('/md01/changjy/data/quail/peak/DEP_KEGG.pdf', p)



library(ChIPseeker)
C1=readPeakFile("/md01/changjy/data/quail/supplementary document/Differential_peak/Differential_peak_C4.txt",sep="\t")
peakAnno=annotatePeak(C1,TxDb = txdb,tssRegion=c(-2000,2000),level = "gene")
write.table(as_tibble(peakAnno), file = "/md01/changjy/data/quail/supplementary document/Differential_peak/Differential_peak_C4.anno",sep = '\t', quote = FALSE, row.names = FALSE)




y=read.csv("C:/Users/Jeff/OneDrive/桌面/gProfiler_ggallus_2021-10-17 上午11-08-38__intersections.csv")
library(ggplot2)
ggplot(y,aes(cluster,Term))+
  geom_point(aes(size=Enrichment,color=-log(pvalue)))+
  scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_minimal()+
  xlab(NULL)+
  ylab(NULL)+
  scale_size_continuous(range=c(0,6))+
  theme(axis.text.x=element_text(size=6,angle=-30),axis.text.y=element_text(size=6))
ggsave('C:/Users/Jeff/OneDrive/桌面/quail/plot/ATAC_KEGG_bubble.pdf', p, width = 5, height = 6)
y$Term <- factor(y$term,levels = c(
"Propanoate metabolism",
"Cell adhesion molecules (CAMs)", 
"Peroxisome",
"Regulation of actin cytoskeleton", 
"Glycosaminoglycan degradation",


"Axon guidance", 
"beta-Alanine metabolism",
"Tyrosine metabolism",


"Calcium signaling pathway", 
"Nicotinate and nicotinamide metabolism",
"Intestinal immune network for IgA production",
"Adrenergic signaling in cardiomyocytes",

"TGF-beta signaling pathway",
"Fanconi anemia pathway",
"MAPK signaling pathway",
"One carbon pool by folate",
"Pentose and glucuronate interconversions",
))





species <- "Homo sapiens"
collection <- "CORE"
opts <- list()
opts["species"] <- species
opts["collection"] <- collection
out  <- TFBSTools::getMatrixSet(JASPAR2018, opts)

if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
  names(out) <- paste(names(out), TFBSTools::name(out), 
                      sep = "_")
motif <- out


data=read.csv("/md01/changjy/data/quail/footprint/cluster1/cluster_1_known10.motif.position.motif.find",header=F,sep="\t")
peaks <- GRanges(seqnames = data$V1,ranges=IRanges(start=data$V2,end=data$V3))
library(BSgenome.Quail.UCSC.taeGut2)
motif_ix <- matchMotifs(motif, peaks, genome = "BSgenome.Quail.UCSC.taeGut2")
motif_ix <- matchMotifs(motif, peaks, genome = "BSgenome.Quail.UCSC.taeGut2",
                        out = "positions")

library(TFBSTools)
library(motifmatchr)
library(tibble)
library(ChIPseeker)
library(GenomicFeatures)
txdb<-makeTxDbFromGFF("/md01/changjy/data/quail/ref/ncbi_dataset/GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff",organism="Coturnix japonica",dataSource="NCBI_ANNOTATION")
for (i in c(1:10)){
motif=paste("/md01/changjy/data/quail/motif/cluster4/knownResults/known",i,".motif",sep="")
data=read.csv(motif,sep="\t")[,c(1:4)]
num=as.vector(as.matrix(data*1000))
nam=colnames(data)[2]
## PFMatrix construction; Not all of the slots need to be initialised.
pfm <- PFMatrix(ID="MA0004.1", name=nam, 
                matrixClass="Zipper-Type",
                bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                tags=list(family="Helix-Loop-Helix", species="10090",
                          tax_group="vertebrates",medline="7592839", 
                          type="SELEX",ACC="P53762", pazar_tf_id="TF0000003",
                          TFBSshape_ID="11", TFencyclopedia_ID="580"),
                profileMatrix=matrix(num,
                                     byrow=TRUE, nrow=4,
                                     dimnames=list(c("A", "C", "G", "T"))
                                     )
                )
peak=read.table("/md01/changjy/data/quail/peak/DESeq2/test/cluster/cluster4.txt",sep="\t")
library(BSgenome.Quail.UCSC.taeGut2)
peaks <- GRanges(seqnames = peak[,1],
                 ranges = IRanges(start = peak[,2],
                                  end=peak[,3]))
motif_ix <- matchMotifs(  pfm, peaks, genome = "BSgenome.Quail.UCSC.taeGut2")

motif_ix <- matchMotifs(  pfm, peaks, genome = "BSgenome.Quail.UCSC.taeGut2",
                         out = "positions")

out=paste(motif,"_TFBSTools",".bed",sep="")
write.table(as_tibble(motif_ix)[,c(-1,-2)],out,sep = '\t', quote = FALSE,col.names=F, row.names = FALSE)
peak=readPeakFile(out)
peakAnno <- annotatePeak(peak,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
out_anno=paste(out,".anno",sep="")
write.table(cbind(as.data.frame(peak),as_tibble(peakAnno@anno)), file = out_anno,sep = '\t', quote = FALSE, row.names = FALSE)

}




library(JASPAR2018)
species <- "Homo sapiens"
collection <- "CORE"
opts <- list()
opts["species"] <- species
opts["collection"] <- collection
out  <- TFBSTools::getMatrixSet(JASPAR2018, opts)
 
if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
  names(out) <- paste(names(out), TFBSTools::name(out), 
                      sep = "_")
motif <- out



##R中合并find函数结果和peak对应关系生成总文件
for (i in c(1:10)){

position_file=paste("/md01/changjy/data/quail/motif/cluster4/knownResults/position/motif",i,".position",sep="")
position=read.table(position_file)
bed=read.table("/md01/changjy/data/quail/peak/DESeq2/test/cluster/cluster4.txt")
cbind(bed[position[,1],],position)
file=paste(position_file,".merge",sep="")
write.table(cbind(bed[position[,1],],position),file,sep="\t",row.names=F,col.names=F,quote=F)

}



##python扫描正负链序列打印motif位置
# _*_ coding: utf-8 _*_


import regex
import pysam
import pandas as pd
import sys

file_path=sys.argv[1]
print(file_path)
#ref=sys.argv[1]
def rev(seq):
    base_trans = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 'a':'t', 'c':'g', 't':'a', 'g':'c', 'n':'n'}
    rev_seq = list(reversed(seq))
    rev_seq_list = [base_trans[k] for k in rev_seq]
    rev_seq = ''.join(rev_seq_list)
    return(rev_seq)
def search_motif(str, substr):
    match = regex.finditer(substr, str, overlapped=True)
    position = [x.start() for x in match]
    return position
if __name__ == '__main__':
    file=pd.read_csv(file_path,header=None,sep="\t")
    ref = "/md01/changjy/data/quail/ref/ncbi_dataset/index/quail.fa"
    with pysam.FastaFile(ref) as fa_obj:
        res=[]
        for index,row in file.iterrows():
            seq = fa_obj.fetch(str(row[0]), row[1],row[2])
            subseq=str(row[5])
            #print(file)
            res1=search_motif(seq, subseq)
            res2=search_motif(rev(seq), subseq)
            #print(res1)
            #print(res2)
            #res1.extend(res2)
            #print(seq)
            #print(rev(seq))
            if res1:
                for pos in res1:
                    chr=str(row[0])
                    start=int(row[1])+int(pos)
                    end=start+len(subseq)
                    print(chr,start,end,subseq,str(row[6]),str(row[7]),sep="\t")
            if res2:
                for pos in res2:
                    chr=str(row[0])
                    start=int(row[1])+int(pos)
                    end=start+len(subseq)
                    print(chr,start,end,subseq,str(row[6]),str(row[7]),sep="\t")


ls /md01/changjy/data/quail/motif/cluster4/knownResults/position/*.merge|while read file
do
out=${file%}.bed
python motif.py  $file |sort -k1,1 -k2,2n > $out
done


library(tibble)
library(ChIPseeker)
library(GenomicFeatures)
txdb<-makeTxDbFromGFF("/md01/changjy/data/quail/ref/ncbi_dataset/GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff",organism="Coturnix japonica",dataSource="NCBI_ANNOTATION")
file_path="/md01/changjy/data/quail/motif/cluster4/knownResults/position/"
file_list <- list.files(path="/md01/changjy/data/quail/motif/cluster4/knownResults/position/", pattern="*.bed")
file_name=paste(file_path,file_list,sep="")
for (i in file_name){
    peak=readPeakFile(i)
    peakAnno <- annotatePeak(peak,TxDb = txdb,tssRegion=c(-3000,3000),level = "gene")
    out=paste(i,".anno",sep="")

    write.table(cbind(as.data.frame(peak),as_tibble(peakAnno@anno)), file = out,sep = '\t', quote = FALSE, row.names = FALSE)
}



