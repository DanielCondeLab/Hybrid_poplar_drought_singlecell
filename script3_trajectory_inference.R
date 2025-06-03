suppressMessages(require(docopt))
suppressMessages(require(Seurat))
suppressMessages(require(tidyverse))
suppressMessages(require(grDevices))
suppressMessages(require(slingshot))
suppressMessages(require(RColorBrewer))
suppressMessages(require(ComplexHeatmap))
suppressMessages(require(grDevices))
suppressMessages(require(SingleCellExperiment))
suppressMessages(require(ggthemes))
suppressMessages(require(tradeSeq))
suppressMessages(require(UpSetR))
suppressMessages(require(ggplot2))
suppressMessages(require(ggridges))

set.seed(1407)

## Load the dataset
sample_name <- "XYLEM_CELL_DIV_35_01.rds"
cds <- readRDS( sample_name)

## Check the datasets
( p <- Seurat::DimPlot(cds,
                       label.size = 5,
                       label = T) )

( p_treat <- Seurat::DimPlot(cds,
                             split.by = "treat",
                             label.size = 5,
                             label = T) )

system("mkdir -p results")
ggsave("results/XYLEM_CELL_DIV_35_01_with_all_clusters.svg",
       p,
       width = 12,
       height = 10)

ggsave("results/XYLEM_CELL_DIV_35_01_with_all_clusters_by_treat.svg",
       p_treat,
       width = 24,
       height = 10)

## Removing cells from C0

### Get the cell IDs for cluster 0
cells_to_remove <- WhichCells(cds, idents = 0)

### Subset the Seurat object to exclude these cells
sing_cell_data <- subset(cds,
                         cells = setdiff(colnames(cds),
                                         cells_to_remove))

( p1 <- Seurat::DimPlot(sing_cell_data,
                        label.size = 5,
                        label = T) )

( p1_treat <- Seurat::DimPlot(sing_cell_data,
                              split.by = "treat",
                              label.size = 5,
                              label = T) )

ggsave("results/XYLEM_CELL_DIV_35_01_without_C0.svg",
       p1,
       width = 12,
       height = 10)

ggsave("results/XYLEM_CELL_DIV_35_01_without_C0_by_treat.svg",
       p1_treat,
       width = 24,
       height = 10)

saveRDS(sing_cell_data, file = "results/XYLEM_CELL_DIV_35_01_without_C0.rds")

########################################################### 
### Trajectory inference and DEGs within the trajectory ###
###########################################################

## Prepare the data for trajectory
COUNTS <- as.matrix(sing_cell_data@assays$RNA@counts)
COUNTS <- COUNTS[rowSums(COUNTS) > 0, ]
length(rownames(COUNTS))

sce <- SingleCellExperiment(assays = List(counts = COUNTS ) )
reducedDims(sce) <- SimpleList(UMAP = sing_cell_data@reductions$umap@cell.embeddings)
colData(sce)$seurat_clusters <- sing_cell_data$seurat_clusters
colData(sce)$treat <- sing_cell_data$treat

## Trajectory inference using slingshot
sds <- slingshot::slingshot(
    sce,
    clusterLabels = "seurat_clusters",
    start.clus = 4,
    reducedDim = 'UMAP',
    shrink = 1L,
    reweight = TRUE,
    reassign = TRUE,
    maxit = 10L,
    smoother = "smooth.spline",
)

# Visualize the inferred lineages
slingshot::slingLineages(sds)

### Plots the cells dynamics on each lineage by treatment ###

## Lineage 1
pseudotime_control_lin1 <- slingPseudotime(sds)[colData(sds)$treat == "Control", 1]
pseudotime_control_lin1 <- pseudotime_control_lin1[!is.na(pseudotime_control_lin1)]

pseudotime_drought_lin1 <- slingPseudotime(sds)[colData(sds)$treat == "Drought", 1]
pseudotime_drought_lin1 <- pseudotime_drought_lin1[!is.na(pseudotime_drought_lin1)]

den_control_lin1 <- density( pseudotime_control_lin1 )
den_drought_lin1 <- density( pseudotime_drought_lin1 )

svg(filename = "results/density_plot_pseudotime_lineage1.svg", width = 10, height = 8)
plot( den_control_lin1, 
      xlab = "Pseudotime - 'lineage' 1", main = "")
lines( den_drought_lin1)
legend("topright",
       legend=c( paste("Control:", length(pseudotime_control_lin1), "cells"),
                 paste("Drought:", length(pseudotime_drought_lin1), "cells") ),
       fill=c(rgb(2/255, 99/255, 227/255, alpha = 0.6),
              rgb(227/255, 62/255, 2/255, alpha = 0.6)))
polygon(den_control_lin1, col = rgb(2/255, 99/255, 227/255, alpha = 0.6))
polygon(den_drought_lin1, col = rgb(227/255, 62/255, 2/255, alpha = 0.6))
dev.off()

# Perform a statistical test to check if the differences in the distribution are significant 
ks.test(slingPseudotime(sds)[colData(sds)$treat == "Control", 1],
        slingPseudotime(sds)[colData(sds)$treat == "Drought", 1])

## Lineage 2
pseudotime_control_lin2 <- slingPseudotime(sds)[colData(sds)$treat == "Control", 2]
pseudotime_control_lin2 <- pseudotime_control_lin2[!is.na(pseudotime_control_lin2)]

pseudotime_drought_lin2 <- slingPseudotime(sds)[colData(sds)$treat == "Drought", 2]
pseudotime_drought_lin2 <- pseudotime_drought_lin2[!is.na(pseudotime_drought_lin2)]

den_control_lin2 <- density( pseudotime_control_lin2 )
den_drought_lin2 <- density( pseudotime_drought_lin2 )

svg(filename = "results/density_plot_pseudotime_lineage2.svg", width = 10, height = 8)
plot( den_control_lin2, 
      xlab = "Pseudotime - 'lineage' 2", main = "")
lines( den_drought_lin2)
legend("topright",
       legend=c( paste("Control:", length(pseudotime_control_lin2), "cells"),
                 paste("Drought:", length(pseudotime_drought_lin2), "cells") ),
       fill=c(rgb(2/255, 99/255, 227/255, alpha = 0.6),
              rgb(227/255, 62/255, 2/255, alpha = 0.6)))
polygon(den_control_lin2, col = rgb(2/255, 99/255, 227/255, alpha = 0.6))
polygon(den_drought_lin2, col = rgb(227/255, 62/255, 2/255, alpha = 0.6))
dev.off()

ks.test(slingPseudotime(sds)[colData(sds)$treat == "Control", 2],
        slingPseudotime(sds)[colData(sds)$treat == "Drought", 2])

### Split the ITs by treatment ###

### Filtering by the treat
to_filter_spec1 <- sing_cell_data@meta.data %>%
    dplyr::filter( treat ==  "Control") %>%
    rownames()

to_filter_spec2 <- sing_cell_data@meta.data %>%
    dplyr::filter(treat == "Drought") %>%
    rownames()

sc_data_traj_spec1 <- base::subset(sing_cell_data,
                                   cells = to_filter_spec1)
sc_data_traj_spec2 <- base::subset(sing_cell_data,
                                   cells = to_filter_spec2)

species_data <- list(Control = sc_data_traj_spec1,
                     Drought = sc_data_traj_spec2)

for (t in 1:length(species_data) ) {
    
    sing_cell_data_sub <- species_data[[t]]
    
    COUNTS <- as.matrix(sing_cell_data_sub@assays$RNA@counts)
    COUNTS <- COUNTS[rowSums(COUNTS) > 0, ]
    
    sce <- SingleCellExperiment(assays = List(counts = COUNTS ) )
    reducedDims(sce) <- SimpleList(UMAP = sing_cell_data_sub@reductions$umap@cell.embeddings)
    colData(sce)$seurat_clusters <- sing_cell_data_sub$seurat_clusters
    
    sds <- slingshot::slingshot(
        sce,
        clusterLabels = "seurat_clusters",
        start.clus = 4,
        reducedDim = 'UMAP',
        shrink = 1L,
        reweight = TRUE,
        reassign = TRUE,
        maxit = 10L,
        smoother = "smooth.spline",
    )
    
    slingshot::slingLineages(sds)
    
    sds2 <- SlingshotDataSet(sds)
    counts <- as.matrix(assays(sds)$counts)
    
    sce_fitted <- fitGAM(counts = counts,
                         sds = sds2,
                         nknots = 6,
                         verbose = T)
    
    saveRDS(sce_fitted,
            file = paste0("TRADEseq_",
                          names(species_data)[t],
                          ".rds")
    )
    
}

sce_fitted_control <- readRDS("TRADEseq_Control.rds")
sce_fitted_drought <- readRDS("TRADEseq_Drought.rds")

# test for dynamic expression
feat_importances_control <- associationTest(sce_fitted_control, lineages = TRUE,
                                            contrastType="end", inverse="Chol")
feat_importances_drought <- associationTest(sce_fitted_drought, lineages = TRUE,
                                            contrastType="end", inverse="Chol")

feat_importances_control$fdr_lin1 <- stats::p.adjust(feat_importances_control$pvalue_1,
                                                     method = "fdr",
                                                     n = length(feat_importances_control$pvalue_1))

feat_importances_control$fdr_lin2 <- stats::p.adjust(feat_importances_control$pvalue_2,
                                                     method = "fdr",
                                                     n = length(feat_importances_control$pvalue_2))

feat_importances_drought$fdr_lin1 <- stats::p.adjust(feat_importances_drought$pvalue_1,
                                                     method = "fdr",
                                                     n = length(feat_importances_drought$pvalue_1))

feat_importances_drought$fdr_lin2 <- stats::p.adjust(feat_importances_drought$pvalue_2,
                                                     method = "fdr",
                                                     n = length(feat_importances_drought$pvalue_2))

feat_importances_control$genes <- rownames(feat_importances_control)
feat_importances_drought$genes <- rownames(feat_importances_drought)

write.table(feat_importances_control,
            "results/association_test_output_control.tsv",
            quote = F, sep = "\t", col.names = T, row.names = F)

write.table(feat_importances_drought,
            "results/association_test_output_drought.tsv",
            quote = F, sep = "\t", col.names = T, row.names = F)

## Adjust the p-values
ControlGenes_lin1 <- rownames(feat_importances_control)[
    which( feat_importances_control$fdr_lin1 <= 0.05)
]
ControlGenes_lin2 <- rownames(feat_importances_control)[
    which( feat_importances_control$fdr_lin2 <= 0.05)
]
DroughtGenes_lin1 <-  rownames(feat_importances_drought)[
    which( feat_importances_drought$fdr_lin1 <= 0.05)
]
DroughtGenes_lin2 <-  rownames(feat_importances_drought)[
    which( feat_importances_drought$fdr_lin2 <= 0.05)
]

length(ControlGenes_lin1)
length(DroughtGenes_lin1)
length(ControlGenes_lin2)
length(DroughtGenes_lin2)

png("results/Upset_plot_filtered.png",
    width = 26, height = 14, units = "cm", res = 300)
UpSetR::upset(
    fromList( list( Control_lin1 = ControlGenes_lin1,
                    Drought_lin1 = DroughtGenes_lin1,
                    Control_lin2 = ControlGenes_lin2,
                    Drought_lin2 = DroughtGenes_lin2) )
)
dev.off()

png("results/Upset_plot_all_intersections.png",
    width = 26, height = 14, units = "cm", res = 300)
UpSetR::upset(
    fromList( list(
        Control_lin1 = ControlGenes_lin1,
        Drought_lin1 = DroughtGenes_lin1,
        Control_lin2 = ControlGenes_lin2,
        Drought_lin2 = DroughtGenes_lin2)
    ),empty.intersections = TRUE
)
dev.off()

Control_lin1_df <- tibble(DEG = ControlGenes_lin1, Condition = "Control Lineage 1")
Control_lin2_df <- tibble(DEG = ControlGenes_lin2, Condition = "Control Lineage 2")
Drought_lin1_df <- tibble(DEG = DroughtGenes_lin1, Condition = "Drought Lineage 1")
Drought_lin2_df <- tibble(DEG = DroughtGenes_lin2, Condition = "Drought Lineage 2")

all_DEGs_tradeseq <- rbind(Control_lin1_df,
                           Control_lin2_df,
                           Drought_lin1_df,
                           Drought_lin2_df)

write.table(all_DEGs_tradeseq,
            "results/all_DEGs_tradeseq.tsv",
            quote = F, sep = "\t", col.names = T, row.names = F)

list_filter <- list("Control lin 1" = ControlGenes_lin1,
                    "Drought lin 1" = DroughtGenes_lin1,
                    "Control lin 2" = ControlGenes_lin2,
                    "Drought lin 2" = DroughtGenes_lin2)

df2 <- data.frame(gene=unique(unlist(list_filter)))

df1 <- lapply(list_filter,function(x){
    data.frame(gene = x)
}) %>% 
    bind_rows(.id = "path")

df_int <- lapply(df2$gene,function(x){
    
    intersection <- df1 %>% 
        dplyr::filter(gene==x) %>% 
        arrange(path) %>% 
        pull("path") %>% 
        paste0(collapse = "|")
    
    data.frame(gene = x,int = intersection)
}) %>% 
    bind_rows()

df_int %>% 
    group_by(int) %>% 
    summarise(n=n()) %>% 
    arrange(desc(n))

write.table(df_int,
            "results/all_DEGs_tradeseq_by_intersection.tsv",
            quote = F, sep = "\t", col.names = T, row.names = F)