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

file_uuids = colnames(RE_intron_1_raw_counts)
head(file_uuids)

#write.csv(file_uuids, "RE_introns.csv", row.names = FALSE)

new_column_names <- (barcodes_RE_introns$`0`)
colnames(RE_intron_1_raw_counts) <- new_column_names
View(RE_intron_1_raw_counts)

genomic_idx <- match(rownames(clinical), colnames(RE_intron_1_raw_counts))
genomic_idx

dataset_ordered <- RE_intron_1_raw_counts[, genomic_idx]
all(rownames(clinical) == colnames(dataset_ordered))

#saveRDS(dataset_ordered, file="Dataset_RE_introns_TEEAnalysis.rds")

clinical$definition <- gsub(" ", "_", clinical$definition)

#Setting definition as factor as it will be the classification factor
clinical$definition <- as.factor(clinical$definition)

#Setting the levels
clinical$definition <- relevel(clinical$definition, ref = "Solid_Tissue_Normal")

dds <- DESeqDataSetFromMatrix(countData = round(dataset_ordered),
                              colData = clinical,
                              design = ~definition)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

View(normalized_counts)
nrow(dds)
mcols(dds)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

#Extracting vst matrix
vsd_mat <- assay(vsd)

#Compute pairwise correlation analysis
vsd_cor <- cor(vsd_mat)
View(vsd_cor)
class(vsd_cor)
pheatmap(vsd_cor, annotation = select(clinical, definition))

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

mcols(dds_res)$description
head(dds_res, n=10)
summary(dds_res)

dds_res_sig <- subset(dds_res, padj < 0.05)
head(dds_res_sig)
summary(dds_res_sig)

write.csv(dds_res_sig, file="significant_introns.csv")
write.csv(dds_res, file="all_significant_introns.csv")