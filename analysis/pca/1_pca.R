# === Set working directory ===

# === Load packages ===
library(ggplot2)
library(readr)
library(dplyr)
library(ggrepel)

# === Step 1: Read and clean metadata ===
md <- read_tsv("../../metadata/sample_information_from_prep_17425.tsv", show_col_types = FALSE)

skin_md <- md %>%
  filter(microbiome_type == "human_skin") %>%
  mutate(
    ad_status = tolower(ad_status),
    sample_name = gsub("\\.", "_", sample_name)
  ) %>%
  as.data.frame()

rownames(skin_md) <- skin_md$sample_name

# === Step 2: Read featureCounts output ===
counts <- read.delim("../../outputs/reference-map/featureCounts/counts.tsv", sep = "\t", header = TRUE, comment.char = "#", check.names = FALSE)
count_data <- counts[, 7:ncol(counts)]
rownames(count_data) <- counts$Geneid

# === Step 3: Filter, transpose, and scale ===
count_data <- count_data[rowSums(count_data) > 0, ]
count_data_t <- t(count_data)
count_data_scaled <- scale(log1p(count_data_t))

# === Step 4: PCA ===
pca <- prcomp(count_data_scaled, center = TRUE, scale. = TRUE)

# === Step 5: Compute % variance explained for axis labels ===
var_explained <- summary(pca)$importance[2, 1:2] * 100
x_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
y_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# === Step 6: Extract coordinates ===
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$sample <- rownames(pca_df)

# Clean sample names
pca_df <- pca_df %>%
  mutate(sample = gsub("_filt$", "", gsub("\\.bam$", "", basename(sample))))

# Set as rownames
rownames(pca_df) <- pca_df$sample
pca_df$sample <- NULL

# === Step 7: Merge metadata ===
skin_md <- skin_md %>%
  mutate(sample_name = gsub("\\.", "_", sample_name))

rownames(skin_md) <- skin_md$sample_name

# Subset metadata columns for merging
meta_subset <- skin_md[, c("ad_status", "lesion_status")]

# Debug overlap
cat("Overlap samples:", sum(rownames(pca_df) %in% rownames(meta_subset)), "\n")

# Merge by rownames
pca_df <- cbind(pca_df, meta_subset[rownames(pca_df), ])

# Save coordinates
write_tsv(pca_df %>% mutate(sample = rownames(pca_df)), "../../outputs/reference-map/pca/pca_coordinates.tsv")

# === Step 8: Filter for samples with ad_status ===
pca_df_filtered <- pca_df %>%
  filter(!is.na(ad_status))

# === Step 9: MANOVA test ===
manova_result <- manova(cbind(PC1, PC2) ~ ad_status, data = pca_df_filtered)
summary_manova <- summary(manova_result, test = "Pillai")
p_value <- summary_manova$stats[1, "Pr(>F)"]

cat("MANOVA p-value:", p_value, "\n")
print(summary_manova)

# === Step 10: PCA Plot with % variance explained ===
p <- ggplot(pca_df_filtered, aes(x = PC1, y = PC2, color = ad_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    data = subset(pca_df_filtered, ad_status == "ad_pos"),
    aes(label = lesion_status),
    size = 4,
    color = "gray50",
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("ad_pos" = "salmon", "healthy" = "lightblue")) +
  theme_minimal(base_size = 14) +
  labs(
    x = x_label,
    y = y_label,
    color = "AD Status"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("MANOVA p = ", signif(p_value, 3)),
    hjust = 1.1, vjust = 1.5,
    size = 4,
    color = "gray50"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

# Save plot as 600 DPI PNG
ggsave("../../figures/pca.png", p, width = 8, height = 6, dpi = 600)
cat("PCA plot saved as 'pca.png' (600 DPI)\n")

# === Step 11: Output all PCA feature loadings ===
# Extract loadings (feature contributions)
loadings_df <- as.data.frame(pca$rotation[, 1:2])
loadings_df$feature <- rownames(loadings_df)

# Add absolute values to rank contributions
loadings_ranked <- loadings_df %>%
  mutate(abs_PC1 = abs(PC1), abs_PC2 = abs(PC2))

# Save all features
write_tsv(loadings_ranked, "../../outputs/reference-map/pca/pca_feature_loadings.tsv")
cat("All PCA feature loadings saved to '../../outputs/reference-map/pca/pca_feature_loadings.tsv'\n")

# Save sorted by PC1 (descending)
sorted_PC1 <- loadings_ranked %>%
  arrange(desc(PC1))
write_tsv(sorted_PC1, "../../outputs/reference-map/pca/pca_feature_loadings_sorted_by_PC1.tsv")
cat("PCA feature loadings sorted by PC1 saved to '../../outputs/reference-map/pca/pca_feature_loadings_sorted_by_PC1.tsv'\n")

# Save sorted by PC2 (descending)
sorted_PC2 <- loadings_ranked %>%
  arrange(desc(PC2))
write_tsv(sorted_PC2, "../../outputs/reference-map/pca/pca_feature_loadings_sorted_by_PC2.tsv")
cat("PCA feature loadings sorted by PC2 saved to '../../outputs/reference-map/pca/pca_feature_loadings_sorted_by_PC2.tsv'\n")

