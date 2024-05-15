library(DESeq2)
library(pheatmap)
library(RColorBrewer)

### Load the table of counts from the HT-seq files
sampleFiles = list.files("~/final_MSA/counts-files/", pattern="*.htseq.out")
sampleNames = gsub("\\.htseq\\.out", "", sampleFiles)
sampleStates = unlist(lapply(strsplit(sampleNames, split="_"), `[[`, 1))
sampleReps = unlist(lapply(strsplit(sampleNames, split="_"), `[[`, 2))
sampleTable = data.frame(sampleName = sampleNames,
                         fileName = sampleFiles,
                         condition = sampleStates,
                         rep = sampleReps)
### Create the DESeq2 object
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "~/final_MSA/counts-files", design = ~condition)

### Determine the size factors needed for normalization
dds = estimateSizeFactors(ddsHTSeq)

### Do Varaiance Stabilizing transformation
vsd = vst(dds, blind = T,nsub=600)

### Extract the vst matrix from the vsd object
vsd_mat = assay(vsd)

### Compute the pairwise correlation values and plot the heatmap with dendrogram
vsd_cor = cor(vsd_mat)
pheatmap(vsd_cor)

plotPCA(vsd, intgroup = "condition")


### Calculate the mean for each gene
readCounts = counts(ddsHTSeq)
mean_readCounts = apply(readCounts[,1:3], 1, mean)

### Calculate the variance for each gene
var_readCounts = apply(readCounts[,1:3], 1, var)

### Plot the mean versus variance in read count data
df = data.frame(mean_readCounts, var_readCounts)

library(ggplot2)

ggplot(df) +
  geom_point(aes(x=mean_readCounts, y= var_readCounts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  labs(title = "DESeq2 model - Dispersion")


dds = DESeq(dds)
res = results(dds)
write.table(res, file="DESeq_Results", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
