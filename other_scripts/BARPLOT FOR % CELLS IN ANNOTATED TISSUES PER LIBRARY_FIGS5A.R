library(Seurat)
library(dplyr)
library(ggplot2)

# Load the Seurat object
seurat_object <- readRDS("Control-drought_last_30int-50pca-06_june 2024_tracking each REP.rds")

# Create a table with the number of cells per treatment and cluster
df <- seurat_object@meta.data %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  mutate(tissue = case_when(
    cluster %in% c("5", "22", "6", "19", "0", "7", "3", "13") ~ "XY",
    cluster %in% c("16", "4") ~ "PC",
    cluster %in% c("10") ~ "VC",
    cluster %in% c("12", "18", "14") ~ "PH",
    cluster %in% c("8", "1", "2", "20") ~ "CP",
    cluster %in% c("17") ~ "CI",
    cluster %in% c("9", "11", "15") ~ "EP",
    cluster %in% c("23", "21") ~ "UNK",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(tissue)) %>%
  group_by(tissue, treat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

# Set the custom order for the X-axis
df$tissue <- factor(df$tissue, levels = c("XY", "PC", "VC", "PH", "CP", "CI", "EP", "UNK"))

# Colores para los tratamientos
colores_tratamientos <- c(
  "Control1" = "#228B22",
  "Control2" = "#66CD00",
  "Drought1" = "#E07B00",
  "Drought2" = "#FFAE42"
)


ggplot(df, aes(x = tissue, y = percent, fill = treat)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  geom_text(
    aes(label = round(percent)),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  xlab("Tissue") +
  ylab("% of cells per treatment") +
  scale_fill_manual(
    values = colores_tratamientos,
    name = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  ggtitle("Distribution of cells by tissue and treatment") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
