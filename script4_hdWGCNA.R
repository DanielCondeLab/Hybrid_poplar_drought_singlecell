set.seed(1407)

# single-cell analysis package
require(Seurat)

# plotting and data science packages
require(tidyverse)
require(cowplot)
require(patchwork)

# co-expression network analysis packages:
require(WGCNA)
require(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

## Help functions
my_ggsave <- function( name = name,
                       plot = plot,
                       height = 12,
                       width = 14) {
  
  ggplot2::ggsave(
    filename = name,
    plot = plot,
    height = height,
    width = width,
    units = "in",
    bg = "#FFFFFF",
    dpi = 300,
    limitsize = F)
  
}

#############################
## Setting main parameters ##
#############################

## Number of cores
ncores = 36

## RDS file containing the Seurat (clustered) dataset
INPUT = "data/Control-drought_last_30int-50pca-06_june 2024.rds"

## Name of the output folder
sample_name = "whole_dataset_pearson"

## what variable to use for batch correction (in the modules)
## >> Create a variable here telling with you should normalize or not
use_harmony = NULL
harmony_group = "treat"

## What variables to consider when building the metacells.
metacells_grouping = c("cell_type", harmony_group)

## Number of hub genes to be reported
n_of_hub_genes = 100

## What dimension reduction to use during the analysis (this if for the network construction, not data representation. Plots will always use UMAP).
my_reduction = "pca"

# optionally enable multithreading
enableWGCNAThreads(nThreads = ncores)

## Creates dir to save all outputs
dir_name <- paste0("mkdir -p ", sample_name)
system( dir_name)

# load the snRNA-seq dataset
seurat_obj <- readRDS(INPUT)
DefaultAssay(seurat_obj) <- "RNA"

(p1 <- DimPlot(seurat_obj, group.by='seurat_clusters', label=TRUE) +
    umap_theme() + ggtitle('') + NoLegend())

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_input_dataset.png"),
          plot = p1)

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_input_dataset.svg"),
          plot = p1)

(p1_treat <- DimPlot(seurat_obj, group.by='seurat_clusters', label=TRUE, split.by = "treat") +
    umap_theme() + ggtitle('') + NoLegend())

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_input_dataset_by_treat.png"),
          plot = p1_treat)

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_input_dataset_by_treat.svg"),
          plot = p1_treat)

seurat_obj <- Seurat::FindVariableFeatures(object = seurat_obj,
                                           assay = "RNA")
# Scaling the object
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

if ( is.null(seurat_obj@meta.data$cell_type) ) {
  
  seurat_obj@meta.data$cell_type <- seurat_obj@meta.data$seurat_clusters
  
}

# Creating a group variable to be used to focus the analysis
seurat_obj@meta.data$cell_type <- ifelse(
  seurat_obj@meta.data$seurat_clusters == 0,
  "vasculature",
  ifelse( seurat_obj@meta.data$seurat_clusters == 3,
          "vasculature",
          ifelse( seurat_obj@meta.data$seurat_clusters == 4,
                  "vasculature",
                  ifelse( seurat_obj@meta.data$seurat_clusters == 5,
                          "vasculature",
                          ifelse( seurat_obj@meta.data$seurat_clusters == 6,
                                  "vasculature",
                                  ifelse( seurat_obj@meta.data$seurat_clusters == 7,
                                          "vasculature",
                                          ifelse( seurat_obj@meta.data$seurat_clusters == 10,
                                                  "vasculature",
                                                  ifelse( seurat_obj@meta.data$seurat_clusters == 13,
                                                          "vasculature",
                                                          ifelse( seurat_obj@meta.data$seurat_clusters == 16,
                                                                  "vasculature",
                                                                  ifelse( seurat_obj@meta.data$seurat_clusters == 4,
                                                                          "vasculature",
                                                                          ifelse( seurat_obj@meta.data$seurat_clusters == 10,
                                                                                  "vasculature",
                                                                                  ifelse( seurat_obj@meta.data$seurat_clusters == 0,
                                                                                          "vasculature",
                                                                                          "others") ) ) ) ) ) ) ) ) ) ) )

## What cell type (or cluster) is of interest
group_of_interest = "vasculature"

## Set up Seurat object for WGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.001, # fraction of cells that a gene needs to be expressed in order to be included in the network.
  wgcna_name = sample_name # the name of the hdWGCNA experiment
)

# Construct metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = metacells_grouping, # specify the columns in seurat_obj@meta.data to group by
  reduction = my_reduction, # select the dimensionality reduction to perform KNN on
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type' # set the Idents of the metacell seurat object
)

seurat_obj <- NormalizeMetacells(seurat_obj)

# ####################################
# ## Co-expression network analysis ##
# ####################################

## Set up the expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = 'vasculature', # the name of the group of interest in the group.by column
  group.by = 'cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups#
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Select soft-power threshold
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
( p6 <- wrap_plots(plot_list, ncol=2) )
my_ggsave(name = paste0(sample_name, "/soft_power_decision.png"),
          plot = p6)

my_ggsave(name = paste0(sample_name, "/soft_power_decision.svg"),
          plot = p6)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  setDatExpr=FALSE,
  tom_name = sample_name, # name of the topological overlap matrix written to disk
  overwrite_tom = TRUE
)

saveRDS(seurat_obj, "ConstructNetwork_output.rds")

png(filename = paste0(sample_name, "/hdWGCNA_Dendrogram.png"),
    height = 10, width = 20,units = "cm",res = 300, bg = "white")
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

svg(filename = paste0(sample_name, "/hdWGCNA_Dendrogram.svg"),
    height = 10, width = 20, bg = "white")
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')

##########################################
### Module Eigengenes and Connectivity ###
##########################################

# Compute harmonized module eigengenes
# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars=harmony_group ###>>> We may remove this normalization, if it is affecting the biological variation that is expected
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
write.csv(hMEs,
          file = paste0(sample_name, "/harmonized_module_eigengenes.csv"),
          row.names = TRUE)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
write.csv(MEs, file = paste0(sample_name, "/module_eigengenes.csv"),
          row.names = TRUE)

###################################
### Compute module connectivity ###

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type',
  group_name = group_of_interest,
  corFnc = 'cor',
  corOptions = "use = 'p',
  method = 'pearson'")

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(group_of_interest, "-M")
)

# plot genes ranked by kME for each module
n_whole_datasets <- length(unique(seurat_obj@misc$whole_dataset_pearson$wgcna_modules$module))
p7 <- PlotKMEs(seurat_obj, ncol = round(n_whole_datasets/3) )

my_ggsave(name = paste0(sample_name, "/genes_ranked_by_kME_per_module.png"),
          plot = p7, width = round(n_whole_datasets/3)*5, height = 15)

my_ggsave(name = paste0(sample_name, "/genes_ranked_by_kME_per_module.svg"),
          plot = p7, width = round(n_whole_datasets/3)*5, height = 15)

## Getting the module assignment table
# get the module assignment table:
modules <- GetModules(seurat_obj)
write.csv(modules, 
          file = paste0(sample_name, "/module_assignment_table.csv"),
          row.names = T)

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = n_of_hub_genes)
write_csv(hub_df, file = paste0(sample_name, "/hub_genes_of_modules.csv") )

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = n_of_hub_genes,
  method='UCell')

###########################
### Basic Visualization ###
###########################

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "umap",
  features='hMEs', # plot the hMEs
  order=TRUE, # order so the points with highest hMEs are on top
  ucell = T, # depending on Seurat vs UCell for gene scoring
  raster = F
)

# stitch together with patchwork
p8 <- wrap_plots(plot_list, ncol=round(n_whole_datasets/3)) 
my_ggsave(name = paste0(sample_name, "/UMAP_hMEs.png"),
          plot = p8, width = round(n_whole_datasets/3)*3, 
          height = 15)

my_ggsave(name = paste0(sample_name, "/UMAP_hMEs.svg"),
          plot = p8, width = round(n_whole_datasets/3)*3, 
          height = 15)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "umap",
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE, # depending on Seurat vs UCell for gene scoring
  raster = F
)

# stitch together with patchwork
p9 <- wrap_plots(plot_list, ncol=round(n_whole_datasets/3)) 

my_ggsave(name = paste0(sample_name, "/UMAP_hub_scores_provided_by_UCell.png"),
          plot = p9, width = round(n_whole_datasets/3)*3, height = 15)

my_ggsave(name = paste0(sample_name, "/UMAP_hub_scores_provided_by_UCell.svg"),
          plot = p9, width = round(n_whole_datasets/3)*3, height = 15)

# plot module correlagram
png(filename = paste0(sample_name, "/module_correlagram.png"),
    height = round(n_whole_datasets/3)*3,
    width = round(n_whole_datasets/3)*3,
    units = "cm",res = 300, bg = "white")
ModuleCorrelogram(seurat_obj)
dev.off()

svg(filename = paste0(sample_name, "/module_correlagram.svg"),
    height = round(n_whole_datasets/3)*3,
    width = round(n_whole_datasets/3)*3, bg = "white")
ModuleCorrelogram(seurat_obj)
dev.off()

#################################
### Seurat plotting functions ###
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs)
mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p12 <- DotPlot(seurat_obj, features=mods, group.by = 'cell_type') +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
my_ggsave(name = paste0(sample_name, "/dot_plot_of_module_expression.png"),
          plot = p12, width = round(n_whole_datasets/3)*3, height = 15)
my_ggsave(name = paste0(sample_name, "/dot_plot_of_module_expression.svg"),
          plot = p12, width = round(n_whole_datasets/3)*3, height = 15)

#Plot violin oplots per module
for (i in 1:length(mods) ) {
  
  md_name <- paste0(group_of_interest, "-M", i)
  
  md_plot <- VlnPlot(
    seurat_obj,
    features = md_name,
    pt.size = 0,
    group.by = 'cell_type') +
    # add box-and-whisker plots on top:
    geom_boxplot(width=.25, fill='white')+
    # change axis labels and remove legend:
    xlab('') + ylab('hME') + NoLegend()
  
  my_ggsave(name = paste0(sample_name, 
                          "/Violin_plot_", md_name, ".png"),
            plot = md_plot, height = 8, width = 16)
  my_ggsave(name = paste0(sample_name, 
                          "/Violin_plot_", md_name, ".svg"),
            plot = md_plot, height = 8, width = 16)
  
}

#############################
### Network Visualization ###
require(igraph)

ModuleNetworkPlot(seurat_obj, outdir = sample_name)

## Applying UMAP to co-expression networks
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 25, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
p13 <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors.png"),
          plot = p13, width = 12, height = 12)
my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors.svg"),
          plot = p13, width = 12, height = 12)

options(future.globals.maxSize = 10000 * 1024^2)
        
p14 <- ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE,
)

my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors_with_labeled_top_3.png"),
          plot = p14, width = 12, height = 12)
my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors_with_labeled_top_3.svg"),
          plot = p14, width = 12, height = 12)

####################################################
### Differential module eigengene (DME) analysis ###

## Comparing all possibilities for all time-points
md_comparisons <- data.frame( Var1 = c("Control", "Drought"),
                             Var2 = c("Drought", "Control") )
                             
DMEs_all <- data.frame()
for (c in 1:nrow(md_comparisons) ) {
  
  g1 <- md_comparisons[c, 1]
  g2 <- md_comparisons[c, 2]
  
  group1 <- seurat_obj@meta.data %>% 
    subset(cell_type == group_of_interest & treat == g1) %>%
    rownames
  
  group2 <- seurat_obj@meta.data %>% 
    subset(cell_type == group_of_interest & treat == g2) %>%
    rownames
  
  DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    wgcna_name=NULL, 
    verbose = T ) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(comparison = paste0(g1, "_vs_", g2) ) %>%
    dplyr::mutate(p_val_adj = round(p_val_adj, 5) )# %>%
  #tibble::rownames_to_column("Module")
  
  DMEs_all <- rbind(DMEs_all, DMEs)
  
}

write_csv(DMEs_all, file = paste0(sample_name, "/Sig_Diff_Modules_between_time_points.csv") )

## This removed the scaled dataset from the seurat object, to reduce it size before writting it to the disk.
seurat_obj@assays$RNA$scale.data <-  NULL

saveRDS(seurat_obj,
        file = paste0(sample_name, "/processed_seurat_obg.rds") )

## save the network 
net <- GetTOM(seurat_obj)
saveRDS(net, "scRNAseq_populus_drought_network_TOM_file.rds")