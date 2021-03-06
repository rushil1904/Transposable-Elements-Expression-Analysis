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

which( colnames(clinical)=="paper_NMF.based.fraction.APOBEC.b.mutations")

na_rows <- which(is.na(clinical$paper_NMF.based.fraction.APOBEC.b.mutations))
na_rows
NMF.based.fraction.APOBEC.b.mutations_clinical_df <- clinical[-(na_rows),]
View(NMF.based.fraction.APOBEC.b.mutations_clinical_df)
NMF.based.fraction.APOBEC.b.mutations_dataset_ordered <- subset(dataset_ordered, select=-(na_rows))
View(NMF.based.fraction.APOBEC.b.mutations_dataset_ordered)
all(rownames(NMF.based.fraction.APOBEC.b.mutations_clinical_df) == colnames(NMF.based.fraction.APOBEC.b.mutations_dataset_ordered))

na_rows_NMF.based.fraction.APOBEC.b.mutations <- which(NMF.based.fraction.APOBEC.b.mutations_clinical_df$paper_NMF.based.fraction.APOBEC.b.mutations == "ND")
na_rows_NMF.based.fraction.APOBEC.b.mutations
NMF.based.fraction.APOBEC.b.mutations_clinical_df <- NMF.based.fraction.APOBEC.b.mutations_clinical_df[-(na_rows_NMF.based.fraction.APOBEC.b.mutations),]
NMF.based.fraction.APOBEC.b.mutations_dataset_ordered <- subset(NMF.based.fraction.APOBEC.b.mutations_dataset_ordered, select=-(na_rows_NMF.based.fraction.APOBEC.b.mutations))
all(rownames(NMF.based.fraction.APOBEC.b.mutations_clinical_df) == colnames(NMF.based.fraction.APOBEC.b.mutations_dataset_ordered))

NMF.based.fraction.APOBEC.b.mutations_clinical_df$paper_NMF.based.fraction.APOBEC.b.mutations <- as.factor(NMF.based.fraction.APOBEC.b.mutations_clinical_df$paper_NMF.based.fraction.APOBEC.b.mutations)

#Setting the levels
NMF.based.fraction.APOBEC.b.mutations_clinical_df$paper_NMF.based.fraction.APOBEC.b.mutations <- relevel(NMF.based.fraction.APOBEC.b.mutations_clinical_df$paper_NMF.based.fraction.APOBEC.b.mutations, ref = "")

#write.csv(clinical, file="clinical.csv")
#DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(NMF.based.fraction.APOBEC.b.mutations_dataset_ordered),
                              colData = NMF.based.fraction.APOBEC.b.mutations_clinical_df,
                              design = ~paper_NMF.based.fraction.APOBEC.b.mutations)

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
pheatmap(vsd_cor, annotation = select(NMF.based.fraction.APOBEC.b.mutations_clinical_df, paper_NMF.based.fraction.APOBEC.b.mutations))

#Principal Component Analysis
plotPCA(vsd, intgroup="paper_NMF.based.fraction.APOBEC.b.mutations")

dds <- DESeq(dds)

plotDispEsts(dds)

dds_res <- results(dds,
                   contrast = c("paper_NMF.based.fraction.APOBEC.b.mutations", "", ""),
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

write.csv(dds_res_sig, file="significant_raw_paper_NMF.based.fraction.APOBEC.b.mutations.csv")
write.csv(dds_res, file="all_significant_raw_paper_NMF.based.fraction.APOBEC.b.mutations.csv")
