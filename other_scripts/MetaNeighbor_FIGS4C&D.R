library(Seurat)            
library(MetaNeighbor)       
library(SummarizedExperiment)  
library(pheatmap)           

# -----------------------------------------------
# 📂 2. Load Seurat object and filter treatments
# -----------------------------------------------
# Load the .rds file
seurat_object <- readRDS("Control-drought_last_30int-50pca-06_june 2024_tracking each REP.rds")

# Filter only the control treatment replicates
seurat_sub <- subset(seurat_object, subset = treat %in% c("Control1", "Control2"))

# -----------------------------------------------
# 📊 3. Obtain the normalized expression matrix
# -----------------------------------------------
# Extract the log-normalized matrix (genes × cells)
expr_matrix <- GetAssayData(seurat_sub, assay = "RNA", slot = "data")

# -----------------------------------------------
# 🧬 4. Prepare the SummarizedExperiment object
# -----------------------------------------------
# Extract relevant metadata
meta <- seurat_sub@meta.data

# Create SummarizedExperiment object required by MetaNeighbor
se <- SummarizedExperiment(
  assays = list(exprs = as.matrix(expr_matrix)),
  colData = data.frame(
    study_id = meta$treat,                  
    cell_type = as.character(meta$seurat_clusters),  
    row.names = colnames(expr_matrix)
  )
)

# -----------------------------------------------
# 📈 5. Identify highly variable genes between treatments
# -----------------------------------------------
# This improves the sensitivity of MetaNeighbor
var_genes <- variableGenes(dat = se, exp_labels = se$study_id)

# -----------------------------------------------
# 🤖 6. Run MetaNeighborUS for unsupervised reproducibility
# -----------------------------------------------
# Calculate the AUROC matrix between each pair of clusters across treatments
results <- MetaNeighborUS(
  var_genes = var_genes,
  dat = se,
  study_id = se$study_id,
  cell_type = se$cell_type
)

# -----------------------------------------------
# 🔥 7. Visualize AUROC matrix (reproducibility by cell type)
# -----------------------------------------------
cluster_order <- c(5, 22, 6, 19, 0, 7, 3, 13, 16, 4, 10, 12, 18, 14, 8,
                   1, 2, 20, 17, 9, 11, 15, 23, 21)

head(rownames(results))


rows_custom <- paste0("Control1|",cluster_order)
cols_custom <- paste0("Control2|",cluster_order)

existing_rows <- intersect(rows_custom, rownames(results))
existing_cols <- intersect(cols_custom, colnames(results))

results_final <- results[rows_custom, cols_custom]

pheatmap(results_final,
         main = "",  
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize = 10,
         angle_col = 315)

