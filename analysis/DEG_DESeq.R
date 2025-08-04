# Load DESeq2 library
library(DESeq2)
library(ggplot2)
library(pheatmap)

setwd("D:/KCL2024/Courses/7BBG1002_Cloud_computing/Project/data/processed/counts_folder/")

# %%% Merge Data on Geneid %%%

# Specify Input Filepaths
filenames <- c("counts_alignment_SRR1166442.txt", 
               "counts_alignment_SRR1166443.txt", 
               "counts_alignment_SRR1166444.txt", 
               "counts_alignment_SRR1166445.txt", 
               "counts_alignment_SRR1166446.txt", 
               "counts_alignment_SRR1166447.txt")

# Initialize a list to store data
counts_data_list <- list()

# Read all files and store them in the list
# in [25]This line of code is using a regular expression to extract a specific part of the filename, which is the sample ID. 
#The 26 line renames the last column to the sample ID. 
#The 27 line appends the sample ID to any column names that contain "Length".
for (filename in filenames) {
  df <- read.delim(filename, header = TRUE, comment.char = "#")
  sample_id <- sub(".*_(SRR[0-9]+)\\.txt$", "\\1", filename)
  colnames(df)[ncol(df)] <- sample_id
  colnames(df)[grep("Length", colnames(df))] <- paste0("Length_", sample_id)
  
  # Store only necessary columns
  counts_data_list[[filename]] <- df[, c("Geneid", paste0("Length_", sample_id), sample_id)]
}

# Merge all data frames by "Geneid"
counts_df <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), counts_data_list)

# %%% RPKM Normalization %%%

# Extract length and count columns
length_columns <- grep("^Length", colnames(counts_df), value = TRUE)
count_columns <- grep("^SRR", colnames(counts_df), value = TRUE)

# Initialize a data frame to store RPKM normalized values
normalized_df <- counts_df[, "Geneid", drop = FALSE]

# Calculate total mapped reads for each sample
total_mapped_reads <- colSums(counts_df[, count_columns])

# Perform RPKM normalization
for (i in seq_along(count_columns)) {
  count_col <- count_columns[i]
  length_col <- length_columns[i]
  sample_id <- sub("_count$", "", count_col)
  
  # Compute RPKM
  rpkm_values <- (10^9 * counts_df[[count_col]]) / 
    (total_mapped_reads[i] * counts_df[[length_col]])
  
  # Add RPKM values to the RPKM data frame
  normalized_df[[paste0(sample_id, "_RPKM")]] <- rpkm_values
}

# Set Geneid as the index
row.names(normalized_df) <- normalized_df$Geneid
normalized_df <- normalized_df[, -1]

# Optionally, write to CSV
write.csv(normalized_df, "rpkm_normalized_counts.csv", row.names = TRUE)

setwd("D:/KCL2024/Courses/7BBG1002_Cloud_computing/Project/")

# %%% Perform Differential Gene Expression Analysis %%%

# Get dataframe and sample info (meta data)
df <- counts_df[, c("Geneid", count_columns)]
row.names(df) <- df$Geneid
df <- df[, -1]
info <- read.csv("metadata/sample_types.csv", header = TRUE)
head(df)

# Create a colData with sample names and conditions
colData <- DataFrame(
  Sample = info$Run,
  Condition = info$carbon_source)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = df, colData = colData, design = ~ Condition)

# remove lowly expressed entries
filtered_rows <- rowMeans(counts(dds)) >= 5
dds <- dds[filtered_rows, ]

# Run DESeq2 analysis
ddsDE <- DESeq(dds)
normCounts <- counts(ddsDE, normalized = TRUE)
threshold <- 0.001

# Extract differential expression results
results <- results(ddsDE, alpha = threshold)
summary(results)

# Filter for significantly differentially expressed genes
significant_results <- subset(results, padj < threshold & abs(log2FoldChange) >= 1.0)
summary(significant_results)
write.csv(significant_results, "Output/results/sig_DEG_DeSEQ.csv")

# write results to file
result_ordered = results[order(results$padj), ]
write.csv(result_ordered, "Condition_Glucose_vs_Cellobiose.csv")

# Convert DESeq2 results to a dataframe
result_ordered <- as.data.frame(result_ordered)

# Filter out rows with NA values
result_ordered <- result_ordered[complete.cases(result_ordered), ]

# Create a new column "Significant" in result_ordered based on padj

result_ordered$Significant <- ifelse(result_ordered$padj < threshold & abs(result_ordered$log2FoldChange) >= 1.0, "Significant", "Not Significant")



# Create a volcano plot with significance labeling
volcano_plot <- ggplot(result_ordered, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant)) +
  geom_hline(yintercept = -log10(threshold), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  labs(x = "log2(Fold Change)", y = "-log10(padj)") +
  ggtitle("Volcano Plot")

pdf("Output/plots/volcano_plot.pdf", width = 12, height = 8)
print(volcano_plot)
dev.off()

# Map gene expression data with significant results
allsig <- merge(normCounts, significant_results, by = 0)
sigCounts <- allsig[,2:7]
row.names(sigCounts) <- allsig$Row.names
write.csv(allsig, "Output/results/allSignificant_mapped_DESeq.csv")
head(allsig)

# Create a heatmap with significance labeling
heatmap_plot <- pheatmap(log2(sigCounts + 1e-6), scale = 'row', show_rownames = FALSE,
                         treeheight_row = 0, treeheight_col = 0,
                         angle_col = 45)

# Create a PDF device for saving the heatmap
pdf("Output/plots/heatmap_DESeq.pdf", width = 12, height = 8)
print(heatmap_plot)
# Close the PDF device
dev.off()

library(AnnotationDbi)

library(org.Sc.sgd.db)

# Map gene symbols to SGD IDs (example for yeast)
gene_symbols <- c("SUT1", "HAP4")

# Map gene symbols to Ensembl IDs (if applicable)
ensembl_ids <- mapIds(org.Sc.sgd.db, keys = gene_symbols, column = "ENSEMBL", keytype = "GENENAME", multiVals = "first")
print(ensembl_ids)
