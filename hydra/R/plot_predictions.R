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

# Load dataset based on file type
if (tolower(file_ext) == "rds") {
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(scater)
    })
    
    dataset <- readRDS(dataset_path)
    if (inherits(dataset, "SingleCellExperiment")) {
        dataset <- as.Seurat(dataset)
    }
    
} else if (tolower(file_ext) %in% c("h5ad")) {
    suppressPackageStartupMessages({
        library(anndata)
        library(zellkonverter)
    })
    
    adata <- anndata::read_h5ad(dataset_path)
    
    # Convert AnnData to Seurat
    dataset <- zellkonverter::as.Seurat(adata)
    
} else {
    stop("Unsupported file format. Please provide a .rds or .h5ad file.")
}

predicted_labels <- fread(cell_type_predicted)

# Ensure row names match and add predicted labels
predicted_labels <- predicted_labels[order(as.numeric(rownames(predicted_labels))),]

# Add predicted cell type labels to the Seurat object
if (!"x" %in% colnames(predicted_labels)) {
    stop("Predicted labels CSV must contain a column named 'x' with cell type predictions.")
    quit()
}
dataset$predicted_cell_type <- predicted_labels$x

dataset <- NormalizeData(dataset)
dataset <- FindVariableFeatures(dataset)
dataset <- ScaleData(dataset)
dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset))
dataset <- RunUMAP(dataset, dims = 1:10)

common_theme <- theme_bw() +
    theme(
        panel.grid.major = element_blank(),    
        panel.grid.minor = element_blank(),    
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

ggsave(
    filename = file.path(output_dir, paste0("Hydra_predicted_cell_types_", modality, ".png")),
    plot = umap_plot,
    width = 15,
    height = 10,
    units = "in",
    dpi = 300
)
