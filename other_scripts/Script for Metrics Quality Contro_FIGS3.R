library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gt)



seurat_object <- readRDS("Control-drought_last_30int-50pca-06_june 2024_tracking each REP.rds")



# Table to represent the metrics

meta <- seurat_object@meta.data
counts <- GetAssayData(seurat_object, slot = "counts", assay = "RNA")

# Add columns: expressed genes and number of UMIs per cell
meta$nFeature <- Matrix::colSums(counts > 0)
meta$nUMI <- Matrix::colSums(counts)

# be sure that "treat" is a factor
meta$treat <- as.factor(meta$treat)

# number of detected genes per treatment
genes_by_treat <- sapply(levels(meta$treat), function(tr) {
  cells <- rownames(meta)[meta$treat == tr]
  detected_genes <- rowSums(counts[, cells] > 0)
  sum(detected_genes > 0)
})

# Summary per treatment
summary_table <- meta %>%
  group_by(treat) %>%
  summarise(
    `Number of cells` = n(),
    `Genes per cell (mean)` = round(mean(nFeature), 1),
    `UMIs per cell (mean)` = round(mean(nUMI), 1)
  ) %>%
  mutate(`Detected genes` = genes_by_treat[treat]) %>%
  arrange(treat)

# Create a table
summary_table %>%
  gt(rowname_col = "treat") %>%
  tab_header(
    title = "Summary of metrics by library"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 1
  ) %>%
  cols_label(
    treat = "Library"
  ) %>%
  
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(everything())
  ) %>%

  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body()
  ) %>%
 
  tab_style(
    style = cell_borders(sides = "right", color = "gray", weight = px(1)),
    locations = cells_column_labels(columns = -last_col())
  ) %>%
 
  tab_style(
    style = cell_borders(sides = "right", color = "gray", weight = px(1)),
    locations = cells_body(columns = -last_col())
  ) %>%
 
  tab_options(
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    column_labels.border.bottom.style = "solid",
    heading.align = "center"
  )

#Violin plot with UMIs and Genes per cell per library
library(Seurat)
seurat_object <- readRDS("Control-drought_last_30int-50pca-06_june 2024_tracking each REP.rds")
seurat_object$treat <- as.factor(seurat_object$treat)
VlnPlot(
  object = seurat_object,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "treat",
  pt.size = 0,
  ncol = 2  # Número de columnas para organizar los plots
)
