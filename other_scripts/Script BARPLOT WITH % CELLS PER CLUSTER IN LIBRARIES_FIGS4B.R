library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Load the Seurat object
seurat_object <- readRDS("Control-drought_last_30int-50pca-06_june 2024_tracking each REP.rds")

# Obtain the table with clusters and treatments
df <- seurat_object@meta.data %>%
  dplyr::group_by(seurat_clusters, treat) %>%
  dplyr::summarise(n = n(), .groups = "drop")

# Calculate the percentage of cells per cluster and treatment
df <- df %>%
  group_by(seurat_clusters) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

# Reorder the X-axis levels with the desired order
cluster_order <- c(5, 22, 6, 19, 0, 7, 3, 13, 16, 4, 10, 12, 18, 14, 8,
                   1, 2, 20, 17, 9, 11, 15, 23, 21)
df$seurat_clusters <- factor(df$seurat_clusters, levels = cluster_order)

# Create the stacked bar plot with custom order
ggplot(df, aes(x = seurat_clusters, y = percent, fill = treat)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  geom_text(
    aes(label = round(percent)),  # sin decimales y sin símbolo %
    position = position_stack(vjust = 0.5),
    size = 4,                    # tamaño más pequeño
    color = "black"
  ) +
  xlab("Cluster") +
  ylab("% of cells per treatment") +
  scale_fill_manual(
    values = c(
      "Control1" = "#1b9e77",
      "Control2" = "#a6d854",
      "Drought1" = "#8c510a",
      "Drought2" = "#e6ab02"
    ),
    name = "Treatment"
  ) +
  theme_minimal() +
  ggtitle("Percentage distribution of treatments within each cluster") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
