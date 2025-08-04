# This script performs Differential Gene Expression (DGE) analysis using unpaired sample 
# Baggerley’s test between Glucose and Cellobiose-treated yeast samples.
# Set working directory
setwd("D:/KCL2024/Courses/7BBG1002_Cloud_computing/Project")

# Load required packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

# Load sample info
sample_info <- read.csv("../metadata/sample_types.csv", sep = ",", header = TRUE)

# Load normalized RPKM data
glucose_data <- read.csv("../data/processed/Normalized_data/glucose_merged.csv", header = TRUE)
cellobiose_data <- read.csv("../data/processed/Normalized_data/cellobiose_merged.csv", header = TRUE)

# Set gene IDs as rownames and remove the Geneid column
gene_ids <- glucose_data$Geneid
row.names(glucose_data) <- glucose_data$Geneid
glucose_data <- glucose_data[, -1]
row.names(cellobiose_data) <- cellobiose_data$Geneid
cellobiose_data <- cellobiose_data[, -1]

# Define Baggerley’s test function
baggerley_test <- function(prop1, prop2, pseudo_count = 0.001) {
  prop1 <- as.numeric(prop1) + pseudo_count
  prop2 <- as.numeric(prop2) + pseudo_count
  n1 <- length(prop1)
  n2 <- length(prop2)
  mean1 <- mean(prop1)
  mean2 <- mean(prop2)
  var1 <- var(prop1)
  var2 <- var(prop2)
  pooled_var <- (((n1 - 1) * var1) + ((n2 - 1) * var2)) / (n1 + n2 - 2) + 1e-6
  z <- (mean1 - mean2) / sqrt(pooled_var * (1/n1 + 1/n2))
  p_value <- 2 * (1 - pnorm(abs(z)))
  return(p_value)
}

# Apply test to each gene
baggerley_results <- sapply(1:nrow(glucose_data), function(i) {
  baggerley_test(cellobiose_data[i, ], glucose_data[i, ], 0.001)
})

# Store results
result_df2 <- data.frame(
  glucose_norm_rpkm_mean = rowMeans(glucose_data),
  cellobiose_norm_rpkm_mean = rowMeans(cellobiose_data),
  p_value = baggerley_results
)

# Adjust p-values
result_df2$adjusted_p_value <- p.adjust(result_df2$p_value, method = "fdr")

# Calculate fold changes
result_df2$fold_change <- (result_df2$cellobiose_norm_rpkm_mean + 0.001) /
  (result_df2$glucose_norm_rpkm_mean + 0.001)
result_df2$log2_fold_change <- log2(abs(result_df2$fold_change))

# Identify significant DEGs
significant_genes <- result_df2[result_df2$adjusted_p_value <= 0.001 &
                                   abs(result_df2$log2_fold_change) >= 1.0, ]

# Annotate DEGs using GAF file
gaf_data <- as.data.frame(fread("../metadata/sgd.gaf/sgd_noheader.gaf",
                                sep = "\t", quote = "", fill = TRUE))
gaf_data <- gaf_data[gaf_data[,1] == "SGD", c(3,10,11)]
colnames(gaf_data) <- c("feature_id", "description", "synonym")
gaf_data$gene_id <- sub("\\|.*", "", gaf_data$synonym)
gaf_data <- distinct(gaf_data)

significant_genes$gene_id <- rownames(significant_genes)
annotated_degs <- right_join(gaf_data, significant_genes, by = "gene_id")
annotated_degs$feature_id <- coalesce(annotated_degs$feature_id, annotated_degs$gene_id)

# Export results
write.csv(annotated_degs, "../Output/results/Annotated_degs_Baggerley.csv", row.names = FALSE)
write(annotated_degs$feature_id, "../Output/results/FeatureId.txt", ncolumns = 1)
write(annotated_degs$gene_id, "../Output/results/GeneId.txt", ncolumns = 1)

# Transcription factor list
tf_list <- c('MET32', 'MET28', 'THI2', 'MIG2', 'UGA3', 'SIP4', 'MIG3', 'HMS2', 'KAR4',
             'MAL13', 'YAP5', 'DAL80','ADR1', 'USV1', 'CAT8', 'GSM1', 'XBP1', 'SUT1', 'HAP4')
tf_dge <- filter(annotated_degs, feature_id %in% tf_list)

# Volcano plot
result_df2$Significance <- ifelse(result_df2$adjusted_p_value <= 0.001 &
                                   abs(result_df2$log2_fold_change) >= 1.0,
                                   "Significant", "Non-significant")

volcano_plot <- ggplot(result_df2, aes(x = log2_fold_change, y = -log10(adjusted_p_value))) +
  geom_point(size = 0.7, aes(color = Significance)) +
  geom_hline(yintercept = -log10(0.001), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  scale_color_manual(values = c("Non-significant" = "red", "Significant" = "blue")) +
  labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-Value)", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(filename = "../Output/plots/VolcanoPlot_Baggerley.jpg", plot = volcano_plot,
       width = 6, height = 4, dpi = 300)

# Bar plot: Top 10 up/downregulated genes
deg_sorted <- annotated_degs[order(-annotated_degs$log2_fold_change), ]
top_upregulated <- head(deg_sorted, 10)
top_downregulated <- tail(deg_sorted, 10)
top_genes <- rbind(top_upregulated, top_downregulated)

top10 <- ggplot(top_genes, aes(x = reorder(feature_id, log2_fold_change),
                               y = log2_fold_change, fill = log2_fold_change < 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("skyblue", "pink"),
                    labels = c("Top Upregulated", "Top Downregulated")) +
  labs(x = "Genes", y = "Log2 Fold Change", fill = "Regulation Status") +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "top") +
  coord_flip()

ggsave(filename = "../Output/plots/Top_DEG_Baggerley.jpg", plot = top10,
       width = 6, height = 4, dpi = 300)

# Transcription factor plot
tf_fc_plot <- ggplot(tf_dge, aes(x = reorder(feature_id, fold_change),
                                 y = log2_fold_change,
                                 fill = log2_fold_change > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("pink", "skyblue"),
                    labels = c("Downregulated", "Upregulated")) +
  labs(x = "Transcription Factors", y = "Log2 Fold Change", fill = "Expression Levels") +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "top") +
  coord_flip()

ggsave(filename = "../Output/plots/Fold_change_tf_Baggerley.jpg", plot = tf_fc_plot,
       width = 6, height = 4, dpi = 300)

# Optional comparison with DESeq2 results
deseq <- read.csv("../Output/results/sig_DEG_DeSEQ.csv", header = TRUE)
deseq_annotated <- right_join(gaf_data, deseq, by = "gene_id")
deseq_annotated$log2FoldChange <- (-1) * deseq_annotated$log2FoldChange

deseq_sorted <- deseq_annotated[order(-deseq_annotated$log2FoldChange), ]
top_upregulated <- head(deseq_sorted, 10)
top_downregulated <- tail(deseq_sorted, 10)
top_genes <- rbind(top_upregulated, top_downregulated)

top10_deseq <- ggplot(top_genes, aes(x = reorder(feature_id, log2FoldChange),
                                     y = log2FoldChange, fill = log2FoldChange < 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("skyblue", "pink"),
                    labels = c("Top Upregulated", "Top Downregulated")) +
  labs(x = "Genes", y = "Log2 Fold Change", fill = "Regulation Status") +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "top") +
  coord_flip()

ggsave(filename = "../Output/plots/Top_DEG_DESeq.jpg", plot = top10_deseq,
       width = 6, height = 4, dpi = 300)

tf_deseq <- filter(deseq_annotated, feature_id %in% tf_list)

tf_fc_deseq <- ggplot(tf_deseq, aes(x = reorder(feature_id, log2FoldChange),
                                    y = log2FoldChange,
                                    fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("pink", "skyblue"),
                    labels = c("Downregulated", "Upregulated")) +
  labs(x = "Transcription Factors", y = "Log2 Fold Change", fill = "Expression Levels") +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "top") +
  coord_flip()

ggsave(filename = "../Output/plots/Fold_change_tf_DESeq.jpg", plot = tf_fc_deseq,
       width = 6, height = 4, dpi = 300)
