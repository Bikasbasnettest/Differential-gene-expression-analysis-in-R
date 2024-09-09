library(pheatmap)
require(RColorBrewer)
TA <- read.csv("D:/transccriptome analyis datasets/raw_counts.tsv", sep = '\t', row.names = 1)
SI<-read.csv("D:/transccriptome analyis datasets/design.tsv", sep = '\t')
###for setting the factor lvel
factors<-factor(SI$Group)
groups<-unique(SI$Group)
groups<-rev(groups)
groups
SI$Group<-factors
SI$Group
###creation of the object
dds<-DESeqDataSetFromMatrix(countData = TA, colData = SI, design = ~Group)
dds$Group<-relevel(dds$Group, ref = "control")
#creation of the group data
keep<-rowSums(counts(dds) >=10) >=min(table(SI$Group))
dds<-dds[keep,]
dds <- DESeq(dds, test = "Wald", sfType = 'poscount')
deseq_result<-results(dds)
deseq_result
deseq_result<-as.data.frame(deseq_result)
class(deseq_result)
head(deseq_result)
names(deseq_result)
deseq_result$GeneName<-row.names(deseq_result)
names(deseq_result)
head(deseq_result)
deseq_result <- subset(deseq_result, select = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "GeneName", "padj"))
names(deseq_result)
write.table(deseq_result, file = "deseq_result.all.tsv",row.names =F, sep = "\t")
deg<-subset(deseq_result, padj<0.05 & abs(log2FoldChange) >=1)
dim(deg)
dim(deseq_result)
deg<-deg[order(deg$padj),]
head(deg)
write.table(deg, file = "deseq_deg.tsv", row.names = FALSE, sep = "\t")
plotDispEsts(dds, main="GSE203159 Dispersion Estimates")
hist(deseq_result$padj, breaks = seq(0,1, length=21, col= "darkgreen", border = "red",
                                     xlab="", ylab="", main="GSE203159 frequency of padj-value"))
old.pal<-palette(c("green", "orange"))
par(mar=c(4,4,2,1),cex.main=1.5)
title=paste(groups[1], "vs", groups[2])
plot(deseq_result$log2FoldChange, -log10(deseq_result$padj), main="Bikas analysis result for the transcriptime analyis",
     xlab="log2fc", ylab="-log10(padj", pch=20, cex=0.5)
with(subset(deseq_result, padj < 0.05 & abs(log2FoldChange) >= 1), 
     points(log2FoldChange, -log10(padj), 
            pch = 20, 
            col = ifelse(sign(log2FoldChange) == 1, "darkgreen", "red"), 
            cex = 1))
normalized_count<-counts(dds, normalized=T)
head(normalized_count)
transfored_count<-log2(normalized_count+1)
tophits<-row.names(deg[1:15, ])
tophits
tophits<-transfored_count[tophits,]
pheatmap(tophits, cluster_rows = FALSE, cluster_cols = FALSE)
