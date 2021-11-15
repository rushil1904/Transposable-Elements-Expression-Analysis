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

head(dataset_ordered)

na_rows <- which(is.na(clinical$paper_AJCC.Tumor.category))
na_rows
tumor_category_clinical_df <- clinical[-(na_rows),]
View(tumor_category_clinical_df)
tumor_category_dataset_ordered <- subset(dataset_ordered, select=-(na_rows))
View(tumor_category_dataset_ordered)

na_rows_tumor_category <- which(tumor_category_clinical_df$paper_AJCC.Tumor.category == "ND")
na_rows_tumor_category
tumor_category_clinical_df <- tumor_category_clinical_df[-(na_rows_tumor_category),]
tumor_category_dataset_ordered <- subset(tumor_category_dataset_ordered, select=-(na_rows_tumor_category))

discard_rows_tumor_category <- which(tumor_category_clinical_df$paper_AJCC.Tumor.category == "TX")
discard_rows_tumor_category
tumor_category_clinical_df <- tumor_category_clinical_df[-(discard_rows_tumor_category),]
tumor_category_dataset_ordered <- subset(tumor_category_dataset_ordered, select=-(discard_rows_tumor_category))

discard_rows_tumor_category <- which(tumor_category_clinical_df$paper_AJCC.Tumor.category == "T1")
discard_rows_tumor_category
tumor_category_clinical_df <- tumor_category_clinical_df[-(discard_rows_tumor_category),]
tumor_category_dataset_ordered <- subset(tumor_category_dataset_ordered, select=-(discard_rows_tumor_category))

tumor_category_clinical_df$paper_AJCC.Tumor.category <- gsub("T4a", "T4", tumor_category_clinical_df$paper_AJCC.Tumor.category)
tumor_category_clinical_df$paper_AJCC.Tumor.category <- gsub("T3b", "T3", tumor_category_clinical_df$paper_AJCC.Tumor.category)

all(rownames(tumor_category_clinical_df) == colnames(tumor_category_dataset_ordered))

tumor_category_clinical_df$paper_AJCC.Tumor.category <- as.factor(tumor_category_clinical_df$paper_AJCC.Tumor.category)

#Setting the levels
tumor_category_clinical_df$paper_AJCC.Tumor.category <- relevel(tumor_category_clinical_df$paper_AJCC.Tumor.category, ref = "T2")

#write.csv(clinical, file="clinical.csv")
#DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(tumor_category_dataset_ordered),
                              colData = tumor_category_clinical_df,
                              design = ~paper_AJCC.Tumor.category)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

View(normalized_counts)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

#Extracting vst matrix
vsd_mat <- assay(vsd)

#Compute pairwise correlation analysis
vsd_cor <- cor(vsd_mat)
View(vsd_cor)
class(vsd_cor)
pheatmap(vsd_cor, annotation = select(tumor_category_clinical_df, paper_AJCC.Tumor.category))

#Principal Component Analysis
plotPCA(vsd, intgroup="paper_AJCC.Tumor.category")

dds <- DESeq(dds)

plotDispEsts(dds)

dds_res <- results(dds,
                   contrast = c("paper_AJCC.Tumor.category", "T2", "T4"),
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
rownames(dds_res_sig)
res_sig <- data.frame(normalized_counts[rownames(dds_res_sig), ])

heat_colors <- brewer.pal(6, "YlOrRd")

annotation_c <- dplyr::select(tumor_category_clinical_df, paper_AJCC.Tumor.category)

rownames(annotation_c) <- colnames(res_sig)

# Run pheatmap
pheatmap(res_sig,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = T,
         annotation = annotation_c,
         scale = "row"
)

write.csv(dds_res_sig, file="significant_introns_tumor_category_2.csv")
write.csv(dds_res, file="all_significant_introns_tumor_category_2.csv")
