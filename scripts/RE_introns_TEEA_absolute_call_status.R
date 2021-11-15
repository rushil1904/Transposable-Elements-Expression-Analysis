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

na_rows <- which(is.na(clinical$paper_ABSOLUTE.call.status))
na_rows
call_status_clinical_df <- clinical[-(na_rows),]
View(call_status_clinical_df)
call_status_dataset_ordered <- subset(dataset_ordered, select=-(na_rows))
View(call_status_dataset_ordered)
all(rownames(call_status_clinical_df) == colnames(call_status_dataset_ordered))

na_rows_call_status <- which(call_status_clinical_df$paper_ABSOLUTE.call.status == "NA")
na_rows_call_status
call_status_clinical_df <- call_status_clinical_df[-(na_rows_call_status),]
call_status_dataset_ordered <- subset(call_status_dataset_ordered, select=-(na_rows_call_status))
call_status_clinical_df$paper_ABSOLUTE.call.status <- gsub(" ", "_", call_status_clinical_df$paper_ABSOLUTE.call.status)
call_status_clinical_df$paper_ABSOLUTE.call.status <- gsub("-", "_", call_status_clinical_df$paper_ABSOLUTE.call.status)
all(rownames(call_status_clinical_df) == colnames(call_status_dataset_ordered))

call_status_clinical_df$paper_ABSOLUTE.call.status <- as.factor(call_status_clinical_df$paper_ABSOLUTE.call.status)

#Setting the levels
call_status_clinical_df$paper_ABSOLUTE.call.status <- relevel(call_status_clinical_df$paper_ABSOLUTE.call.status, ref = "called")

#write.csv(clinical, file="clinical.csv")
#DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(call_status_dataset_ordered),
                              colData = call_status_clinical_df,
                              design = ~paper_ABSOLUTE.call.status)

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
pheatmap(vsd_cor, annotation = select(call_status_clinical_df, paper_ABSOLUTE.call.status))

#Principal Component Analysis
plotPCA(vsd, intgroup="paper_ABSOLUTE.call.status")

dds <- DESeq(dds)

plotDispEsts(dds)

dds_res <- results(dds,
                   contrast = c("paper_ABSOLUTE.call.status", "called", "high_non_clonal"),
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

annotation_c <- dplyr::select(call_status_clinical_df, paper_ABSOLUTE.call.status)

rownames(annotation_c) <- colnames(res_sig)

# Run pheatmap
pheatmap(res_sig,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = T,
         annotation = annotation_c,
         scale = "row"
)

write.csv(dds_res_sig, file="significant_introns_call_status.csv")
write.csv(dds_res, file="all_significant_introns_call_status.csv")
