library(UpSetR)
library(tidyr)
library(dplyr)
deg_data <- read.csv("DEGS FOR UPSET PLOT.csv")
head(deg_data)
deg_data <- deg_data %>%filter(!is.na(Cluster))

deg_wide <- deg_data %>%
  mutate(present = 1) %>%  # Add a column to indicate presence of the gene in the cluster
  tidyr::spread(key = Cluster, value = present, fill = 0)  

colnames(deg_wide)[-1] <- paste0("Cluster", colnames(deg_wide)[-1])

cluster_colors <- c(
  Cluster1 = "dodgerblue",
  Cluster2 = "darkorange",
  Cluster3 = "forestgreen",
  Cluster5 = "purple"
)

svg("upset_plot2.svg", width = 10, height = 8)

desired_order <- c("Cluster3", "Cluster1", "Cluster2", "Cluster5")


upset(deg_wide, 
      sets = desired_order[desired_order %in% colnames(deg_wide)],  # Only include clusters that exist in the data
      order.by = "freq", 
      text.scale = c(2, 1.5, 2, 1.5, 2, 2),
      main.bar.color = "black",  
      sets.bar.color = cluster_colors[1:length(desired_order)],  # Adjust colors according to the new order
      matrix.color = "red",  
      number.angles = 30,
      keep.order = TRUE)  # Keep order of sets as in the input

dev.off()
