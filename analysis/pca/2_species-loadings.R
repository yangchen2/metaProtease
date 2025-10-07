# === Load species-level loadings and summarize ===
library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)

# Read PCA loadings
loadings <- read_tsv("pca_feature_loadings.tsv", show_col_types = FALSE)

# Extract species name (text before first "|")
loadings <- loadings %>%
  mutate(species = sub("\\|.*", "", feature))

# Summarize by species: take mean absolute loading per PC
species_loadings <- loadings %>%
  group_by(species) %>%
  summarize(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    abs_PC1 = mean(abs_PC1),
    abs_PC2 = mean(abs_PC2)
  ) %>%
  ungroup()

# Identify top separating species (strongest loadings along PC1 or PC2)
top_species <- species_loadings %>%
  mutate(abs_sum = abs_PC1 + abs_PC2) %>%
  arrange(desc(abs_sum)) %>%
  slice(1:15)  # adjust number of labels

# === Scatter plot of species-level contributions ===
p <- ggplot(species_loadings, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.6, color = "gray60") +
  geom_point(data = top_species, aes(x = PC1, y = PC2), color = "steelblue", size = 3) +
  geom_text_repel(
    data = top_species,
    aes(label = species),
    size = 4,
    box.padding = 0.3,
    point.padding = 0.2,
    color = "black",
    max.overlaps = Inf
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Species Contributions to PCA Axes",
    x = "PC1 Loadings",
    y = "PC2 Loadings"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# === Save plot as high-resolution PNG ===
ggsave("../../figures/species_loadings.png", plot = p, width = 8, height = 8, dpi = 600)
cat("Species loadings plot saved as 'species_loadings.png' (600 DPI)\n")

