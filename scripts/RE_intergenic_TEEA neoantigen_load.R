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

na_rows <- which(is.na(clinical$paper_Neoantigen.load))
na_rows
neoantigen_load_clinical_df <- clinical[-(na_rows),]
View(neoantigen_load_clinical_df)
neoantigen_load_dataset_ordered <- subset(dataset_ordered, select=-(na_rows))
View(neoantigen_load_dataset_ordered)
all(rownames(neoantigen_load_clinical_df) == colnames(neoantigen_load_dataset_ordered))

class(neoantigen_load_clinical_df$paper_Neoantigen.load)
median(neoantigen_load_clinical_df$paper_Neoantigen.load)
mean(neoantigen_load_clinical_df$paper_Neoantigen.load)

neoantigen_load_clinical_df$paper_Neoantigen.load <- ifelse(neoantigen_load_clinical_df$paper_Neoantigen.load>=550, 'More_than_550', 'Less_than_550')
View(neoantigen_load_clinical_df)

neoantigen_load_clinical_df$paper_Neoantigen.load <- as.factor(neoantigen_load_clinical_df$paper_Neoantigen.load)

#Setting the levels
neoantigen_load_clinical_df$paper_Neoantigen.load <- relevel(neoantigen_load_clinical_df$paper_Neoantigen.load, ref = "Less_than_550")

#write.csv(clinical, file="clinical.csv")
#DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(neoantigen_load_dataset_ordered),
                              colData = neoantigen_load_clinical_df,
                              design = ~paper_Neoantigen.load)

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
pheatmap(vsd_cor, annotation = select(neoantigen_load_clinical_df, paper_Neoantigen.load))

#Principal Component Analysis
plotPCA(vsd, intgroup="paper_Neoantigen.load")

dds <- DESeq(dds)

plotDispEsts(dds)

dds_res <- results(dds,
                   contrast = c("paper_Neoantigen.load", "More_than_550", "Less_than_550"),
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

annotation_c <- dplyr::select(neoantigen_load_clinical_df, paper_Neoantigen.load)

rownames(annotation_c) <- colnames(res_sig)

# Run pheatmap
pheatmap(res_sig,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = T,
         annotation = annotation_c,
         scale = "row"
)

write.csv(dds_res_sig, file="significant_intergenic_neoantigen_load.csv")
write.csv(dds_res, file="all_significant_intergenic_neoantigen_load.csv")
