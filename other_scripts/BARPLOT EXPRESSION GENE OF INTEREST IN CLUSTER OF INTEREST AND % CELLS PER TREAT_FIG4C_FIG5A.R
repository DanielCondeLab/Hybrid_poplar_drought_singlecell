
library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)

seurat_object <- readRDS("XYLEM CELL DIV PARENCHYMA 40 012.rds")
DefaultAssay(seurat_object) <- "RNA"
# Define Inputs
gene_name <- "PtXaAlbH.05G186400"  # Replace with your gene
cluster_id <- "1"                  # Replace with your cluster

# Extract Expression Data from Seurat Object
cluster_cells <- WhichCells(seurat_object, ident = cluster_id)  # Get cells in the cluster

# Fetch expression values for the gene in selected cluster
expr_data <- FetchData(seurat_object, vars = gene_name, cells = cluster_cells)
expr_data$Treatment <- seurat_object$treat[rownames(expr_data)]  # Add treatment info
colnames(expr_data) <- c("Expression", "Treatment")

# Compute Summary Statistics for Average Expression
avg_expression <- expr_data %>%
  group_by(Treatment) %>%
  summarise(
    Average = mean(Expression),
    SE = sd(Expression) / sqrt(n())  # Standard Error (SE)
  )

# Compute % of Cells Expressing the Gene
percent_expressing <- expr_data %>%
  group_by(Treatment) %>%
  summarise(Percent = mean(Expression > 0) * 100)  # Convert to percentage

# Plot 1: Average Expression Bar Plot with SE
p1 <- ggplot(avg_expression, aes(x = Treatment, y = Average, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Average - SE, ymax = Average + SE), width = 0.2) +
  scale_fill_manual(values = c("lightblue", "lightcoral")) +
  theme_minimal() +
  ggtitle(paste("Average Expression of", gene_name, "in Cluster", cluster_id)) +
  ylab("Average Expression (±SE)") +
  theme(legend.position = "none")

# Plot 2: % of Cells Expressing Bar Plot
p2 <- ggplot(percent_expressing, aes(x = Treatment, y = Percent, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_manual(values = c("lightblue", "lightcoral")) +
  theme_minimal() +
  ggtitle(paste("% of Cells Expressing", gene_name, "in Cluster", cluster_id)) +
  ylab("Percentage of Cells (%)") +
  theme(legend.position = "none")

# Combine Both Plots Side by Side
grid.arrange(p1, p2, ncol = 2)
