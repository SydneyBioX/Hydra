#!/usr/bin/env Rscript

###############################################

# Manoj M Wagle (USydney)

###############################################


suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(ggridges)
    library(rlang)
    library(SingleCellExperiment)
    library(scater)
})

args <- commandArgs(trailingOnly = TRUE)
dataset_path <- args[1]              
modality <- args[2]                   
cell_type_label <- args[3]           
gene_name <- ifelse(args[4] == "None", NA, args[4])   # Gene name for expression plotting
cell_type_of_interest <- ifelse(length(args) >= 5 && args[5] != "None", args[5], NA)  # Specific cell type

# Determine file type based on extension
file_ext <- tools::file_ext(dataset_path)

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
    
} else if (tolower(file_ext) %in% c("h5ad", "h5")) {
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
    stop("Unsupported file format. Please provide a .rds (Seurat) or .h5ad (AnnData) file.")
}

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

# Generate UMAP plot based on cell type label
if (cell_type_label %in% colnames(dataset@meta.data)) {
    umap_plot <- DimPlot(
        dataset, 
        reduction = "umap", 
        group.by = cell_type_label, 
        pt.size = 0.6, 
        alpha = 1, 
        raster = FALSE
    ) +
        ggtitle("UMAP Plot Colored by Cell Types") +
        common_theme
} else {
    # If cell type label does not exist, perform clustering
    suppressPackageStartupMessages({
        library(Seurat)
    })
    dataset <- FindNeighbors(dataset, dims = 1:10)
    dataset <- FindClusters(dataset, resolution = 0.5)
    umap_plot <- DimPlot(
        dataset, 
        reduction = "umap", 
        group.by = "seurat_clusters", 
        pt.size = 0.6, 
        alpha = 1, 
        raster = FALSE
    ) +
        ggtitle("UMAP Plot Colored by Clusters") +
        common_theme
}

# Create the output directory if it doesn't exist
output_dir <- "Results/Plots"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Save the UMAP plot at 300 DPI
ggsave(
    filename = file.path(output_dir, paste0("umap_cell_types_", modality, ".png")),
    plot = umap_plot,
    width = 15,
    height = 10,
    units = "in",
    dpi = 300
)

# If gene_name is provided, plot UMAP with gene expression highlighted
if (!is.na(gene_name)) {
    if (!(gene_name %in% rownames(dataset))) {
        warning(paste("Gene", gene_name, "not found in the dataset. Skipping gene expression plot."))
    } else {
        gene_plot <- FeaturePlot(
            dataset, 
            features = gene_name, 
            reduction = "umap", 
            pt.size = 0.6, 
            alpha = 1, 
            raster = FALSE
        ) +
            scale_color_gradient(low = "lightgrey", high = "red", name = "Log-normalized\nexpression") +
            ggtitle(paste("UMAP Plot of", gene_name, "Expression")) +
            common_theme
        # Save the gene expression UMAP plot at 300 DPI
        ggsave(
            filename = file.path(output_dir, paste0("umap_gene_expression_", gene_name, ".png")),
            plot = gene_plot,
            width = 15,
            height = 10,
            units = "in",
            dpi = 300
        )
    }
}

# If cell_type_of_interest is provided, create a ridgeline plot
if (!is.na(cell_type_of_interest) && !is.na(gene_name)) {
    if (!(gene_name %in% rownames(dataset))) {
        warning(paste("Gene", gene_name, "not found in the dataset. Skipping ridgeline plot."))
    } else if (!(cell_type_label %in% colnames(dataset@meta.data))) {
        warning(paste("Cell type label", cell_type_label, "not found in the dataset metadata. Skipping ridgeline plot."))
    } else {
        expression_data <- FetchData(dataset, vars = c(gene_name, cell_type_label), layer = "data")
        expression_data$Group <- ifelse(
            tolower(expression_data[[cell_type_label]]) == tolower(cell_type_of_interest), 
            cell_type_of_interest, 
            "Other cell types"
        )
        expression_data$Group <- factor(expression_data$Group, levels = c("Other cell types", cell_type_of_interest))
        
        ridgeline_plot <- ggplot(expression_data, aes(x = .data[[gene_name]], y = Group, fill = Group)) +
            geom_density_ridges() +
            scale_fill_manual(values = c("lightgray", "red")) +
            labs(
                x = "Log-normalized expression",
                y = "Cell type",
                title = paste("Expression Distribution of", gene_name)
            ) +
            ggridges::theme_ridges(grid = FALSE) +
            theme(
                legend.position = "none",
                axis.title.x = element_text(hjust = 0.5, vjust = -0.5),
                axis.title.y = element_text(hjust = 0.5)
            )
        
        # Save the ridgeline plot at 300 DPI
        ggsave(
            filename = file.path(output_dir, paste0("ridgeline_expression_plot_", cell_type_of_interest, "_", gene_name, ".png")),
            plot = ridgeline_plot,
            width = 8,
            height = 6,
            units = "in",
            dpi = 300
        )
    }
}
