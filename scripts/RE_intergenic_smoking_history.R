pacman::p_load(pacman, TCGAbiolinks, IHW, apeglm, pheatmap, RColorBrewer, PCAtools, reshape2)
library(BiocManager)
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
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", force = TRUE) 
BiocManager::install("apeglm", force = TRUE)

setwd("~/Desktop/Intern/Trasposable element expression")

file_uuids = colnames(RE_intergenic_1_raw_counts)
head(file_uuids)

#write.csv(file_uuids, "RE_intergenic.csv", row.names = FALSE)

new_column_names <- (barcodes_RE_intergenic$`0`)
colnames(RE_intergenic_1_raw_counts) <- new_column_names
View(RE_intergenic_1_raw_counts)

genomic_idx <- match(rownames(clinical), colnames(RE_intergenic_1_raw_counts))
genomic_idx

dataset_ordered <- RE_intergenic_1_raw_counts[, genomic_idx]
all(rownames(clinical) == colnames(dataset_ordered))

na_rows <- which(is.na(clinical$paper_Tobacco.smoking.history))
na_rows
smoking_history_clinical_df <- clinical[-(na_rows),]
View(smoking_history_clinical_df)
smoking_history_dataset_ordered <- subset(dataset_ordered, select=-(na_rows))
View(smoking_history_dataset_ordered)
all(rownames(smoking_history_clinical_df) == colnames(smoking_history_dataset_ordered))

smoking_history_clinical_df$paper_Tobacco.smoking.history <- gsub(" ", "_", smoking_history_clinical_df$paper_Tobacco.smoking.history)
smoking_history_clinical_df$paper_Tobacco.smoking.history <- gsub("-", "_", smoking_history_clinical_df$paper_Tobacco.smoking.history)
smoking_history_clinical_df$paper_Tobacco.smoking.history <- gsub("reformed_smoker>15_years", "reformed_smoker_more_than_15_years", smoking_history_clinical_df$paper_Tobacco.smoking.history)
smoking_history_clinical_df$paper_Tobacco.smoking.history <- gsub("Current_reformed_smoker_for_<_or_=_15_years", "reformed_smoker_less_than_or_equal_to_15_years", smoking_history_clinical_df$paper_Tobacco.smoking.history)
smoking_history_clinical_df$paper_Tobacco.smoking.history <- gsub("Current_reformed_smoker_for_>_15_years", "Current_reformed_smoker_for_more_than_15_years", smoking_history_clinical_df$paper_Tobacco.smoking.history)
smoking_history_clinical_df$paper_Tobacco.smoking.history <- gsub("Current_Reformed_Smoker,_Duration_Not_Specified", "Current_Reformed_Smoker_Duration_Not_Specified", smoking_history_clinical_df$paper_Tobacco.smoking.history)
smoking_history_clinical_df$paper_Tobacco.smoking.history <- gsub("ND", "NA", smoking_history_clinical_df$paper_Tobacco.smoking.history)
smoking_history_clinical_df$paper_Tobacco.smoking.history <- gsub("Current_Reformed_Smoker_Duration_Not_Specified", "NA", smoking_history_clinical_df$paper_Tobacco.smoking.history)
na_rows_smoking <- which(smoking_history_clinical_df$paper_Tobacco.smoking.history == "NA")
na_rows_smoking
smoking_history_clinical_df <- smoking_history_clinical_df[-(na_rows_smoking),]
smoking_history_dataset_ordered <- subset(smoking_history_dataset_ordered, select=-(na_rows_smoking))

smoking_history_clinical_df$paper_Tobacco.smoking.history <- as.factor(smoking_history_clinical_df$paper_Tobacco.smoking.history)

#Setting the levels
smoking_history_clinical_df$paper_Tobacco.smoking.history <- relevel(smoking_history_clinical_df$paper_Tobacco.smoking.history, ref = "Lifelong_Non_smoker")


#DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(smoking_history_dataset_ordered),
                              colData = smoking_history_clinical_df,
                              design = ~paper_Tobacco.smoking.history)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

vignette("DESeq2")

View(normalized_counts)

vsd <- vst(dds, blind = TRUE)

#Extracting vst matrix
vsd_mat <- assay(vsd)

#Compute pairwise correlation analysis
vsd_cor <- cor(vsd_mat)
View(vsd_cor)
class(vsd_cor)
pheatmap(vsd_cor, annotation = select(smoking_history_clinical_df, paper_Tobacco.smoking.history))

#Principal Component Analysis
plotPCA(vsd, intgroup="paper_Tobacco.smoking.history")

dds <- DESeq(dds)

plotDispEsts(dds)

dds_res <- results(dds,
                   contrast = c("paper_Tobacco.smoking.history", "Lifelong_Non_smoker", "Current_smoker"),
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
View(res_sig)

heat_colors <- brewer.pal(6, "YlOrRd")

annotation_c <- dplyr::select(smoking_history_clinical_df, paper_Tobacco.smoking.history)

rownames(annotation_c) <- colnames(res_sig)

# Run pheatmap
pheatmap(res_sig,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = T,
         annotation = annotation_c,
         scale = "row"
)


write.csv(dds_res_sig, file="significant_intergenic_smoking_history.csv")
write.csv(dds_res, file="all_significant_intergenic_smoking_history.csv")