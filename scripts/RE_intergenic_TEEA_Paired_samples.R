pacman::p_load(pacman, TCGAbiolinks, IHW, apeglm, pheatmap, RColorBrewer, PCAtools, reshape2)
library(biomaRt)
library(psych, warn.conflicts=F, quietly=T)
library(dplyr, warn.conflicts=F, quietly=T)
library(magrittr)
library(SummarizedExperiment)
library(DESeq2)
pacman::p_load(GenomicDataCommons)
BiocManager::install("SummarizedExperiment", force = TRUE)
BiocManager::install("biomaRt", force = TRUE)
BiocManager::valid()
BiocManager::install()
BiocManager::install("DESeq2", force = TRUE) 
library(DESeq2)
setwd("~/Desktop/Intern/Trasposable element expression")

column_names <- colnames(Paired_counts_intergenic_final)
paired_clinical <- clinical[column_names,]
paired_intergenic_counts <- RE_intergenic_1_raw_counts[,column_names]
genomic_idx <- match(rownames(paired_clinical), colnames(paired_intergenic_counts))
genomic_idx

dataset_ordered <- paired_intergenic_counts[, genomic_idx]
all(rownames(paired_clinical) == colnames(dataset_ordered))


paired_clinical$definition <- gsub(" ", "_", paired_clinical$definition)

#Setting definition as factor as it will be the classification factor
paired_clinical$definition <- as.factor(paired_clinical$definition)

#Setting the levels
paired_clinical$definition <- relevel(paired_clinical$definition, ref = "Solid_Tissue_Normal")

#DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(dataset_ordered),
                              colData = paired_clinical,
                              design = ~definition)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

View(normalized_counts)
#write.csv(normalized_counts, "normalized_counts.csv")

vsd <- vst(dds, blind = TRUE)

#Extracting vst matrix
vsd_mat <- assay(vsd)

#Compute pairwise correlation analysis
vsd_cor <- cor(vsd_mat)
View(vsd_cor)
class(vsd_cor)
pheatmap(vsd_cor, annotation = select(paired_clinical, definition))

#Principal Component Analysis
plotPCA(vsd, intgroup="definition")

dds <- DESeq(dds)

plotDispEsts(dds)

dds_res <- results(dds,
                   contrast = c("definition", "Primary_solid_Tumor", "Solid_Tissue_Normal"),
                   alpha = 0.05,
                   altHypothesis = "greaterAbs", lfcThreshold = 1.5)

dds_res
plotMA(dds_res, ylim=c(-8,8))

#LFC Shrinkage
dds_res <- lfcShrink(dds,
                     coef=resultsNames(dds)[2], 
                     type="apeglm")
dds_res.Ordered <-  dds_res[with(dds_res, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]

plotMA(dds_res, ylim=c(-8,8)) #Shrinkage should allow for more accurate fold changes

mcols(dds_res)
head(dds_res, n=10)
summary(dds_res)

dds_res_sig <- subset(dds_res, padj < 0.05)
head(dds_res_sig)

rownames(dds_res_sig)
res_sig <- data.frame(normalized_counts[rownames(dds_res_sig), ])

heat_colors <- brewer.pal(6, "YlOrRd")

annotation_c <- dplyr::select(paired_clinical, definition)

rownames(annotation_c) <- colnames(res_sig)

# Run pheatmap
pheatmap(res_sig,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = T,
         annotation = annotation_c,
         scale = "row"
)


write.csv(dds_res_sig, file="significant_intergenic_paired.csv")
write.csv(dds_res, file="all_significant_intergenic_paired.csv")
