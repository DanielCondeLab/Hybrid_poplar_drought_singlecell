library(Seurat)
library(dplyr)


seurat_object <- readRDS("XYLEM CELL DIV PARENCHYMA 40 012_tracking each REP.rds")
gene_of_interest <- "PtXaAlbH.08G129200"


cluster1_cells <- WhichCells(seurat_object, ident = 1)

# Filter by treatment within cluster 1
control1_cells <- colnames(seurat_object)[seurat_object@meta.data$treat == "Control1"]
control2_cells <- colnames(seurat_object)[seurat_object@meta.data$treat == "Control2"]
drought1_cells <- colnames(seurat_object)[seurat_object@meta.data$treat == "Drought1"]
drought2_cells <- colnames(seurat_object)[seurat_object@meta.data$treat == "Drought2"]

# Intersection: treatment + cluster 1
expr_control1 <- GetAssayData(seurat_object, slot = "data")[gene_of_interest, intersect(cluster1_cells, control1_cells)]
expr_control2 <- GetAssayData(seurat_object, slot = "data")[gene_of_interest, intersect(cluster1_cells, control2_cells)]
expr_drought1 <- GetAssayData(seurat_object, slot = "data")[gene_of_interest, intersect(cluster1_cells, drought1_cells)]
expr_drought2 <- GetAssayData(seurat_object, slot = "data")[gene_of_interest, intersect(cluster1_cells, drought2_cells)]

# Mean expression
mean_expr_control1 <- mean(expr_control1, na.rm = TRUE)
mean_expr_control2 <- mean(expr_control2, na.rm = TRUE)
mean_expr_drought1 <- mean(expr_drought1, na.rm = TRUE)
mean_expr_drought2 <- mean(expr_drought2, na.rm = TRUE)

# Show results
print(paste("Expresión media en Control1:", mean_expr_control1))
print(paste("Expresión media en Control2:", mean_expr_control2))
print(paste("Expresión media en Drought1:", mean_expr_drought1))
print(paste("Expresión media en Drought2:", mean_expr_drought2))

# Percentage of cells with expression > mean
percent_high_control1 <- sum(expr_control1 > mean_expr_control1) / length(expr_control1) * 100
percent_high_control2 <- sum(expr_control2 > mean_expr_control2) / length(expr_control2) * 100
percent_high_drought1 <- sum(expr_drought1 > mean_expr_drought1) / length(expr_drought1) * 100
percent_high_drought2 <- sum(expr_drought2 > mean_expr_drought2) / length(expr_drought2) * 100

cat("Percentage of cells with expression > mean (cluster 1 only)")
cat(paste("Control1:", round(percent_high_control1, 2), "%\n"))
cat(paste("Control2:", round(percent_high_control2, 2), "%\n"))
cat(paste("Drought1:", round(percent_high_drought1, 2), "%\n"))
cat(paste("Drought2:", round(percent_high_drought2, 2), "%\n"))

df_expr <- data.frame(
  treat = c("Control1", "Control2", "Drought1", "Drought2"),
  percent = round(c(
    percent_high_control1,
    percent_high_control2,
    percent_high_drought1,
    percent_high_drought2
  )),
  tissue = "Cluster 1"
)


library(ggplot2)
ggplot(df_expr, aes(x = treat, y = percent, fill = tissue)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  geom_text(
    aes(label = round(percent)),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  ylab("% of cells with expression > mean") +
  xlab("Treatment") +
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
