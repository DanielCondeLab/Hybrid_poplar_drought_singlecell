library(Seurat)
library(dplyr)
library(ggplot2)


seurat_object <- readRDS("Cluster X1_35_022_replicas.rds")


df <- seurat_object@meta.data %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  mutate(tissue = case_when(
    cluster == "0" ~ "X1-0",
    cluster == "1" ~ "X1-1",
    cluster == "2" ~ "X1-2",
    cluster == "3" ~ "X1-3"
  )) %>%
  filter(!is.na(tissue)) %>%
  group_by(treat, tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(treat) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()


df$treat <- factor(df$treat, levels = c("Control1", "Control2", "Drought1", "Drought2"))
df$tissue <- factor(df$tissue, levels = c("X1-0", "X1-1", "X1-2", "X1-3"))


ggplot(df, aes(x = treat, y = percent, fill = tissue)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  geom_text(
    aes(label = round(percent)),  
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  ylab("% of cells per cluster") +
  xlab("Treatment") +
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  theme_minimal(base_size = 14) +
  ggtitle("Distribution of cells by treatment and cluster") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
