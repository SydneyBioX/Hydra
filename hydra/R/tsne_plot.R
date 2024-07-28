#!/usr/bin/env Rscript

###############################################
# Manoj M Wagle (USydney; CMRI)
###############################################

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(ggridges)
    library(rlang)
})


args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
modality <- args[2]
cell_type_label <- args[3]
gene_name <- ifelse(args[4] == "None", NA, args[4])
cell_type_of_interest <- ifelse(length(args) >= 5 && args[5] != "None", args[5], NA)

# Load the dataset
dataset <- readRDS(dataset)

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

# Check if cell type label exists
if (cell_type_label %in% colnames(dataset@meta.data)) {
    tsne_plot <- DimPlot(dataset, reduction = "tsne", group.by = cell_type_label, pt.size = 0.6, alpha = 1, raster = FALSE) +
        ggtitle("t-SNE Plot colored by cell types") +
        common_theme
} else {
    # If cell type label does not exist, use clustering
    dataset <- FindNeighbors(dataset, dims = 1:10)
    dataset <- FindClusters(dataset, resolution = 0.5)
    tsne_plot <- DimPlot(dataset, reduction = "tsne", group.by = "seurat_clusters", pt.size = 0.6, alpha = 1, raster = FALSE) +
        ggtitle("t-SNE Plot colored by clusters") +
        common_theme
}

# Create the output directory if it doesn't exist
output_dir <- "Results/Plots"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Save the t-SNE plot
ggsave(filename = file.path(output_dir, paste0("tsne_cell_types_", modality, ".pdf")), plot = tsne_plot, width = 15, height = 10, units = "in")

# If gene_name is provided, plot t-SNE with gene expression highlighted
if (!is.na(gene_name)) {
    gene_plot <- FeaturePlot(dataset, features = gene_name, reduction = "tsne", pt.size = 0.6, alpha = 1, raster = FALSE) +
        scale_color_gradient(low = "lightgrey", high = "red", name = "Log-normalized\nexpression") +
        ggtitle(paste("t-SNE Plot of", gene_name, "Expression")) +
        common_theme
    # Save the gene expression t-SNE plot
    ggsave(filename = file.path(output_dir, paste0("tsne_gene_expression_", gene_name, ".pdf")), plot = gene_plot, width = 15, height = 10, units = "in")
}


# If cell_type_of_interest is provided, create a ridgeline plot
if (!is.na(cell_type_of_interest) && !is.na(gene_name)) {
    expression_data <- FetchData(dataset, vars = c(gene_name, cell_type_label), layer = "data")
    expression_data$Group <- ifelse(tolower(expression_data[[cell_type_label]]) == tolower(cell_type_of_interest), cell_type_of_interest, "Other cell types")
    expression_data$Group <- factor(expression_data$Group, levels = c("Other cell types", cell_type_of_interest))
    
    ridgeline_plot <- ggplot(expression_data, aes(x = .data[[gene_name]], y = Group, fill = Group)) +
        geom_density_ridges() +
        scale_fill_manual(values = c("lightgray", "red")) +
        labs(x = "Log-normalized expression",
             y = "Cell type",
             title = paste("Expression Distribution of", gene_name)) +
        theme_ridges(grid = FALSE) +
        theme(legend.position = "none",
              axis.title.x = element_text(hjust = 0.5, vjust = -0.5),
              axis.title.y = element_text(hjust = 0.5))
    
    ggsave(filename = file.path(output_dir, paste0("ridgeline_expression_plot_", cell_type_of_interest, "_", gene_name, ".pdf")), plot = ridgeline_plot, width = 8, height = 6, units = "in")
}