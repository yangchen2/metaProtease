# === Load functional-level loadings and summarize ===
library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)

# Read PCA loadings
loadings <- read_tsv("pca_feature_loadings.tsv", show_col_types = FALSE)

# Extract gene function annotation (text after first "|")
loadings <- loadings %>%
  mutate(function_annotation = sub("^[^|]*\\|", "", feature))

# Initial classification into broad gene types
loadings <- loadings %>%
  mutate(
    gene_class = case_when(
      grepl("metalloprotease|zinc", function_annotation, ignore.case = TRUE) ~ "Metalloprotease",
      grepl("serine", function_annotation, ignore.case = TRUE) ~ "Serine protease",
      grepl("Clp|ATP-dependent", function_annotation, ignore.case = TRUE) ~ "ATP-dependent protease",
      grepl("peptidase", function_annotation, ignore.case = TRUE) ~ "Peptidase",
      grepl("RIP|RseP", function_annotation, ignore.case = TRUE) ~ "RIP protease",
      TRUE ~ "Other/Unknown"
    )
  )

# --- Step: Expand the 'Other/Unknown' category dynamically ---
# For "Other/Unknown", check which functional annotations appear >=3 times
other_map <- loadings %>%
  filter(gene_class == "Other/Unknown") %>%
  count(function_annotation, sort = TRUE) %>%
  mutate(dynamic_class = ifelse(n >= 3, function_annotation, "Other/Unknown"))

# Join back these refined labels
loadings <- loadings %>%
  left_join(other_map, by = "function_annotation") %>%
  mutate(gene_class_final = ifelse(!is.na(dynamic_class), dynamic_class, gene_class)) %>%
  select(-gene_class, -dynamic_class)

# Summarize mean loadings per functional category
function_loadings <- loadings %>%
  group_by(gene_class_final) %>%
  summarize(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    abs_PC1 = mean(abs_PC1),
    abs_PC2 = mean(abs_PC2),
    n_features = n()
  ) %>%
  ungroup()

# Identify top separating functions
top_functions <- function_loadings %>%
  mutate(abs_sum = abs_PC1 + abs_PC2) %>%
  arrange(desc(abs_sum)) %>%
  slice(1:15)

# === Plot ===
p <- ggplot(function_loadings, aes(x = PC1, y = PC2)) +
  geom_point(aes(size = n_features), alpha = 0.6, color = "gray60") +
  geom_point(data = top_functions, aes(x = PC1, y = PC2, color = gene_class_final, size = n_features)) +
  geom_text_repel(
    data = top_functions,
    aes(label = gene_class_final, color = gene_class_final),
    size = 4,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_size_continuous(name = "Number of Features") +
  theme_minimal(base_size = 14) +
  labs(
    x = "PC1 Loadings",
    y = "PC2 Loadings",
    color = "Gene Type"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# Save PNG
ggsave("../../figures/gene_type_loadings.png", p, width = 8, height = 8, dpi = 600)
cat("Functional gene type PCA loadings plot saved as 'gene_type_loadings.png' (600 DPI)\n")

