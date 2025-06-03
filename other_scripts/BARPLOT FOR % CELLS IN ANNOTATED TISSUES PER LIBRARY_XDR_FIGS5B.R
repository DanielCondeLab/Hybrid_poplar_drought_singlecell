library(Seurat)
library(dplyr)
library(ggplot2)


seurat_object <- readRDS("XYLEM CELL DIV PARENCHYMA 40 012_tracking each REP.rds")


df <- seurat_object@meta.data %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  mutate(tissue = case_when(
    cluster %in% c("5") ~ "VC",
    cluster %in% c("2") ~ "PC",
    cluster %in% c("0") ~ "XMC-0",
    cluster %in% c("1") ~ "XMC-1",
    cluster %in% c("4") ~ "XP",
    cluster %in% c("6") ~ "SCW/PCD",
    cluster %in% c("3") ~ "SCW",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(tissue)) %>%
  group_by(tissue, treat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()


df$tissue <- factor(df$tissue, levels = c("VC", "PC", "XMC-0", "XMC-1", "XP", "SCW", "SCW/PCD"))


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
  scale_fill_manual(values = colores_tratamientos, name = "Treatment") +
  theme_minimal(base_size = 14) +
  ggtitle("Distribution of cells by tissue and treatment") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
