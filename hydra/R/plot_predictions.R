#!/usr/bin/env Rscript

###############################################

# Manoj M Wagle (USydney; CMRI)

###############################################


suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(data.table)
    library(SingleCellExperiment)
    library(scater)
})

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
modality <- args[2]
cell_type_predicted <- args[3]

# Load the dataset
dataset <- readRDS(dataset)

if (inherits(dataset, "SingleCellExperiment")) {
    dataset <- as.Seurat(dataset)
}

# Read the predicted cell types from the CSV file
predicted_labels <- read.csv(cell_type_predicted)

# Make sure the row names of the Seurat object match the order of the predictions
predicted_labels <- predicted_labels[order(as.numeric(rownames(predicted_labels))),]

# Add the predicted cell type labels to the Seurat object
dataset$predicted_cell_type <- predicted_labels$x

# Normalize the data
dataset <- NormalizeData(dataset)

# Find variable features
dataset <- FindVariableFeatures(dataset)

# Scale the data
dataset <- ScaleData(dataset)

# Run PCA
dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset))

# Run t-SNE
dataset <- RunTSNE(dataset, dims = 1:10)

# Common theme for the plots
common_theme <- theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
    )

# Generate t-SNE plot colored by predicted cell types
tsne_plot <- DimPlot(dataset, reduction = "tsne", group.by = "predicted_cell_type", pt.size = 0.6, alpha = 1, raster = FALSE) +
    ggtitle("t-SNE Plot of Predicted Cell Types") +
    common_theme

# Create the output directory if it doesn't exist
output_dir <- "Results/Plots"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Save the t-SNE plot
ggsave(filename = file.path(output_dir, paste0("Hydra_predicted_cell_types_", modality, ".pdf")), plot = tsne_plot, width = 15, height = 10, units = "in")
