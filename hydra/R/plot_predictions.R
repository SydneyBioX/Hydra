#!/usr/bin/env Rscript

###############################################

# Manoj M Wagle (USydney)

###############################################


suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
dataset_path <- args[1]
modality <- args[2]
cell_type_predicted <- args[3]

# Determine file type based on extension
file_ext <- tools::file_ext(dataset_path)

# Initialize dataset
dataset <- NULL

# Load dataset based on file type
if (tolower(file_ext) == "rds") {
    # Load Seurat-specific packages
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(scater)
    })
    
    # Load Seurat object
    dataset <- readRDS(dataset_path)
    if (inherits(dataset, "SingleCellExperiment")) {
        dataset <- as.Seurat(dataset)
    }
    
} else if (tolower(file_ext) %in% c("h5ad")) {
    # Load AnnData-specific packages
    suppressPackageStartupMessages({
        library(anndata)
        library(zellkonverter)
    })
    
    # Read AnnData object
    adata <- anndata::read_h5ad(dataset_path)
    
    # Convert AnnData to Seurat
    dataset <- zellkonverter::as.Seurat(adata)
    
} else {
    stop("Unsupported file format. Please provide a .rds or .h5ad file.")
}

# Read the predicted cell types from the CSV file
predicted_labels <- fread(cell_type_predicted)

# Ensure row names match and add predicted labels
predicted_labels <- predicted_labels[order(as.numeric(rownames(predicted_labels))),]

# Add predicted cell type labels to the Seurat object
if (!"x" %in% colnames(predicted_labels)) {
    stop("Predicted labels CSV must contain a column named 'x' with cell type predictions.")
}
dataset$predicted_cell_type <- predicted_labels$x

# Normalize the data
dataset <- NormalizeData(dataset)

# Find variable features
dataset <- FindVariableFeatures(dataset)

# Scale the data
dataset <- ScaleData(dataset)

# Run PCA
dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset))

# Run UMAP instead of t-SNE
dataset <- RunUMAP(dataset, dims = 1:10)

# Define a common theme for plots without gridlines
common_theme <- theme_bw() +
    theme(
        panel.grid.major = element_blank(),    # Remove major gridlines
        panel.grid.minor = element_blank(),    # Remove minor gridlines
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
    )

# Generate UMAP plot colored by predicted cell types
umap_plot <- DimPlot(
    dataset, 
    reduction = "umap", 
    group.by = "predicted_cell_type", 
    pt.size = 0.6, 
    alpha = 1, 
    raster = FALSE
) +
    ggtitle("UMAP Plot of Predicted Cell Types") +
    common_theme

# Create the output directory if it doesn't exist
output_dir <- "Results/Plots"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Save the UMAP plot at 300 DPI
ggsave(
    filename = file.path(output_dir, paste0("Hydra_predicted_cell_types_", modality, ".png")),
    plot = umap_plot,
    width = 15,
    height = 10,
    units = "in",
    dpi = 300
)
