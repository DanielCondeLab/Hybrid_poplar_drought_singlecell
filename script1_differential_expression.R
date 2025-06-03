library (Seurat)
seurat_obj <- readRDS("XYLEM CELL DIV PARENCHYMA 40 012.rds")
Idents(seurat_obj) <- seurat_obj$seurat_clusters
clusters <- unique(seurat_obj$seurat_clusters)
deg_results <- list()

for (cluster in clusters) {
  # Subset the Seurat object to include only the current cluster
  cluster_subset <- subset(seurat_obj, idents = cluster)
  
  # Perform differential expression analysis
  degs <- FindMarkers(
    cluster_subset, 
    ident.1 = "Drought", 
    ident.2 = "Control", assay = "RNA", slot = "data",
    test.use = "bimod", group.by = "treat",
    min.pct = 0.03,         # Include all genes regardless of percentage
    logfc.threshold = 0.35, # Include all genes regardless of fold-change
    verbose = TRUE
  )
  
  # Add cluster information to the results
  degs$Cluster <- cluster
  
  # Store results in the list
  deg_results[[paste0("Cluster_", cluster)]] <- degs
}
for (i in seq_along(deg_results)) {
  # Example of adding a gene ID column if it’s missing
  if (!"gene" %in% colnames(deg_results[[i]])) {
    deg_results[[i]]$gene <- rownames(deg_results[[i]])
  }
}
deg_combined <- do.call(rbind, deg_results)
rownames(deg_combined) <- NULL
print(colnames(deg_combined))
filtered_degs <- deg_combined[deg_combined$p_val_adj < 0.05, ]
output_file <- "xylem_BIMOD_degs_003PC_035FC_normalized log_R script.csv"
write.csv(filtered_degs, file = output_file, row.names = FALSE)
