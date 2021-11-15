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

file_uuids = colnames(RE_all_1_raw_counts)
head(file_uuids)

new_column_names <- (barcodes_RE_raw$`0`)
colnames(RE_all_1_raw_counts) <- new_column_names
View(RE_all_1_raw_counts)

genomic_idx <- match(rownames(clinical), colnames(RE_all_1_raw_counts))
genomic_idx

dataset_ordered <- RE_all_1_raw_counts[, genomic_idx]
all(rownames(clinical) == colnames(dataset_ordered))
head(dataset_ordered)

na_rows <- which(is.na(clinical$paper_Histologic.subtype))
na_rows
Histologic_subtype_clinical_df <- clinical[-(na_rows),]
View(Histologic_subtype_clinical_df)
Histologic_subtype_dataset_ordered <- subset(dataset_ordered, select=-(na_rows))
View(Histologic_subtype_dataset_ordered)

na_rows_histology <- which(Histologic_subtype_clinical_df$paper_Histologic.subtype == "ND")
na_rows_histology
Histologic_subtype_clinical_df <- Histologic_subtype_clinical_df[-(na_rows_histology),]
Histologic_subtype_dataset_ordered <- subset(Histologic_subtype_dataset_ordered, select=-(na_rows_histology))
all(rownames(Histologic_subtype_clinical_df) == colnames(Histologic_subtype_dataset_ordered))
Histologic_subtype_clinical_df$paper_Histologic.subtype <- gsub("-", "_", Histologic_subtype_clinical_df$paper_Histologic.subtype)

Histologic_subtype_clinical_df$paper_Histologic.subtype <- as.factor(Histologic_subtype_clinical_df$paper_Histologic.subtype)

#Setting the levels
Histologic_subtype_clinical_df$paper_Histologic.subtype <- relevel(Histologic_subtype_clinical_df$paper_Histologic.subtype, ref = "Papillary")

#write.csv(clinical, file="clinical.csv")
#DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(Histologic_subtype_dataset_ordered),
                              colData = Histologic_subtype_clinical_df,
                              design = ~paper_Histologic.subtype)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

View(normalized_counts)

vsd <- vst(dds, blind = TRUE)

#Extracting vst matrix
vsd_mat <- assay(vsd)

#Compute pairwise correlation analysis
vsd_cor <- cor(vsd_mat)
View(vsd_cor)
class(vsd_cor)
pheatmap(vsd_cor, annotation = select(Histologic_subtype_clinical_df, paper_Histologic.subtype))

#Principal Component Analysis
plotPCA(vsd, intgroup="paper_Histologic.subtype")

dds <- DESeq(dds)

plotDispEsts(dds)

dds_res <- results(dds,
                   contrast = c("paper_Histologic.subtype", "Papillary", "Non_Papillary"),
                   alpha = 0.05,
                   altHypothesis = "greaterAbs", lfcThreshold = 1.5)

dds_res
plotMA(dds_res, ylim=c(-6.5,6.5))

#LFC Shrinkage
dds_res <- lfcShrink(dds,
                     coef=resultsNames(dds)[2], 
                     type="apeglm")
dds_res.Ordered <-  dds_res[with(dds_res, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]

plotMA(dds_res, ylim=c(-6.5,6.5)) #Shrinkage should allow for more accurate fold changes

mcols(dds_res)
head(dds_res, n=10)
summary(dds_res)

dds_res_sig <- subset(dds_res, padj < 0.05)
head(dds_res_sig)

write.csv(dds_res_sig, file="significant_histologic_subtype.csv")
write.csv(dds_res, file="all_significant_raw_histologic_subtype.csv")
