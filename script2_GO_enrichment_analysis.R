library(topGO)
genes_of_interest <- read.csv("UP1.csv")$Gene
#Load the list of poplar IDs with the GOs of their Arabidopsis homologs
all_genes_with_go <- read.csv("Palba all genome GOs_no ARA.csv")
geneUniverse <- all_genes_with_go$Gene
#Identify the genes in the list that have assigned GOs
geneList <- factor(as.integer(geneUniverse %in% genes_of_interest))
names(geneList) <- geneUniverse
GOdata <- new("topGOdata",
              ontology = "BP",  # 'BP' for Biological Process, 'MF' for Molecular Function, 'CC' for Cellular Component
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = split(all_genes_with_go$GO, all_genes_with_go$Gene))
result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
total_GO_terms <- length(usedGO(GOdata))
GOresult <- GenTable(GOdata, weight01Fisher = result,
                     orderBy = "weight01Fisher", 
                     ranksOf = "weight01Fisher", 
                     topNodes = total_GO_terms)
GOresult$weight01Fisher <- as.numeric(GOresult$weight01Fisher)
#Sometimes the p-value is NA when it is less than 1e-30. We fix it with this
GOresult$weight01Fisher[is.na(GOresult$weight01Fisher)] <- 1e-30
print(GOresult)

#Calculate the adjusted p-values
GOresult$adj_p_value <- p.adjust(GOresult$weight01Fisher, method = "BH")
print(GOresult)
library(ggplot2)
GOresult_filtered <- GOresult[GOresult$adj_p_value < 0.05, ]
GOresult_filtered$GeneRatio <- as.numeric(GOresult_filtered$Significant) / as.numeric(GOresult_filtered$Annotated)
print(GOresult_filtered)
write.csv(GOresult_filtered,"GO enrichment UP1.csv")
GOresults_plot <- GOresult_filtered[1:50,]
ggplot(GOresults_plot, aes(x = GeneRatio, y = reorder(Term, GeneRatio), 
                              size = GeneRatio, 
                              color = adj_p_value)) +
  geom_point() +
  scale_color_gradient(limits=c(0, 0.05), low = "purple", high = "orange", 
                       name = "p.adjusted") + 
  scale_size_continuous(name = "GeneRatio", range = c(2, 6)) +
  labs(x = "GeneRatio", y = NULL,
       title = "GO Enrichment DEGS UP cluster 0 TOP 50") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),  # Adjust font size for x-axis
    axis.text.y = element_text(size = 7, hjust = 0, vjust = 0.5),  # Adjust font size, and alignment for y-axis
    axis.title.x = element_text(size = 14),  # Adjust font size for x-axis label
    axis.title.y = element_text(size = 14),  # Adjust font size for y-axis label
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title size and style
    legend.title = element_text(size = 12),  # Adjust font size for legend title
    legend.text = element_text(size = 10),  # Adjust font size for legend text
    plot.margin = margin(1, 1, 1.5, 1.5, "cm")  # Increase plot margins to prevent clipping
  )

#Identify the genes of interest within each enriched GO 
resultTable <- read.csv("GO enrichment UP1.csv")
genesOfInterest <- read.csv("UP1.csv")
enrichedGOs <- resultTable$GO.ID
genesOfInterest <- genesOfInterest$Gene
allGenesInTerms <- genesInTerm(GOdata)
genesPerEnrichedGO <- list()

for(goTerm in enrichedGOs) {
  # Get genes associated with this GO term from the genesInTerm() result
  genesInGO <- allGenesInTerms[[goTerm]]
  
  # Filter to keep only genes of interest
  genesOfInterestInGO <- genesInGO[genesInGO %in% genesOfInterest]
  
  # Store the result if there are any genes of interest in this GO term
  if(length(genesOfInterestInGO) > 0) {
    genesPerEnrichedGO[[goTerm]] <- genesOfInterestInGO
  }
}

# Convert the result to a data frame for easier viewing/saving
resultDF <- data.frame(
  GO.Term = rep(names(genesPerEnrichedGO), lengths(genesPerEnrichedGO)),
  Genes = unlist(genesPerEnrichedGO)
)

# Save the results to a CSV file if needed
write.csv(resultDF, "genes_in_enriched_GOs_COMMON_WAT1-2_200coexpressed.csv", row.names = FALSE)

# Load necessary packages
library(dplyr)
library(stringr)

# Read your CSV file
data <- read.csv("Genes_GOs.csv")

# Group genes by GO term and collapse them into a single comma-separated string
grouped_data <- data %>%
  group_by(GO.Term) %>%
  summarise(Genes = str_c(Genes, collapse = ","))

# Write the output to a new CSV file
write.csv(grouped_data, "grouped_genes_by_GO_UP1.csv", row.names = FALSE)
