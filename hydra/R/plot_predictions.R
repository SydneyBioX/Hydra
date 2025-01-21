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
        library(scater)
    })
    
    dataset <- readRDS(dataset_path)
    if (inherits(dataset, "SingleCellExperiment")) {
        library(Seurat)
        library(SingleCellExperiment)
        dataset <- as.Seurat(dataset)
    }
    
} else if (tolower(file_ext) %in% c("h5ad")) {
    suppressPackageStartupMessages({
        library(reticulate)
        library(Seurat)
    })
    anndata <- import("anndata", convert = FALSE)
    
    adata <- anndata$read_h5ad(dataset_path)

    if (py_has_attr(adata$X, "toarray")) {
        assay_data <- t(py_to_r(adata$X$toarray()))  
    } else {
        assay_data <- t(py_to_r(adata$X))  
    }
    
    gene_names <- py_to_r(adata$var$index$to_list())
    rownames(assay_data) <- gene_names

    cell_names <- py_to_r(adata$obs$index$to_list())
    colnames(assay_data) <- cell_names

    dataset <- CreateSeuratObject(counts = assay_data, project = "SeuratProject", assay = modality)
    
} else {
    stop("Unsupported file format. Please provide a .rds or .h5ad file.")
}

predicted_labels <- read.csv(cell_type_predicted)
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
    ggtitle("UMAP plot of Hydra predicted cell types") +
    common_theme

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
