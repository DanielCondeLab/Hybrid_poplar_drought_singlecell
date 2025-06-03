library(Seurat)
library(ggplot2) 
integrated_data <- readRDS("Cluster X1_35_022_replicas.rds")
integrated_data <- NormalizeData(integrated_data, normalization.method = "LogNormalize", scale.factor = 10000)
genes_list <- read.csv("Vessels figure early vessels reviewer2 ordered.csv", header = TRUE)
gene_ids <- genes_list$Gene
desired_order <- c("0", "1", "2", "3")

                      

dot_plot <- DotPlot(integrated_data, features = gene_ids, , assay = "RNA", dot.scale = 5, cols = c("aliceblue", "red")) + ggtitle("X") + scale_y_discrete(limits = desired_order) + scale_color_gradientn(limits = c(-1.25, 2), colors = c("azure2", 'azure3', "green4")) + coord_flip() + scale_x_discrete (position = "top") +
  theme(
    axis.text.x = element_text(size = 8),  # Adjust x-axis text size
    axis.text.y = element_text(size = 8)   # Adjust y-axis text size
    ,  # Adjust y-axis text size
    axis.title.x = element_blank(),         # Remove x-axis title
    axis.title.y = element_blank(),          # Remove y-axis title , 
    legend.text = element_text(size = 8),  # Adjust legend text size
legend.title = element_text(size = 8))  # Adjust legend title size

                      
dot_plot2 <- dot_plot +
scale_x_discrete(position="top", limits = rev(levels(dot_plot$data$feature)))

dot_plot2
                      