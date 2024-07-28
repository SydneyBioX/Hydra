#!/usr/bin/env Rscript

##############################################

# Manoj M Wagle (USydney; CMRI)

##############################################


args <- commandArgs(trailingOnly = TRUE)
train_file <- args[1]
test_file <- args[2]
cell_type_label <- args[3]
data_type <- args[4]


##############################################

suppressPackageStartupMessages({
  library(rhdf5)
  library(HDF5Array)
  library(Seurat)
  library(caret)
  library(data.table)
  library(SingleCellExperiment)
  library(Signac)
  library(glue)
  library(scater)
  library(reticulate)
})


##############################################

# Writing H5 files

write_h5_scJoint <- function(exprs_list, h5file_list) {
  
  if (length(unique(lapply(exprs_list, rownames))) != 1) {
    stop("rownames of exprs_list are not identical.")
  }
  
  for (i in seq_along(exprs_list)) {
    if (file.exists(h5file_list[i])) {
      warning("h5file exists! will rewrite it.")
      system(paste("rm", h5file_list[i]))
    }
    
    h5createFile(h5file_list[i])
    h5createGroup(h5file_list[i], "matrix")
    writeHDF5Array(t((exprs_list[[i]])), h5file_list[i], name = "matrix/data")
    h5write(rownames(exprs_list[[i]]), h5file_list[i], name = "matrix/features")
    h5write(colnames(exprs_list[[i]]), h5file_list[i], name = "matrix/barcodes")
    print(h5ls(h5file_list[i]))
    
  } 
}

write_csv_scJoint <- function(cellType_list, csv_list) {
  
  for (i in seq_along(cellType_list)) {
    
    if (file.exists(csv_list[i])) {
      warning("csv_list exists! will rewrite it.")
      system(paste("rm", csv_list[i]))
    }
    
    names(cellType_list[[i]]) <- NULL
    write.csv(cellType_list[[i]], file = csv_list[i])
    
  }
}


##############################################

# Processing datasets
preprocess_dataset_train <- function(dataset_file, cell_type_label) {
  file_extension <- tools::file_ext(dataset_file)
  
  if (file_extension == "rds") {
    dataset <- readRDS(dataset_file)
  } else if (file_extension == "h5ad") {
    anndata <- import("anndata", convert = FALSE)
    dataset <- anndata$read_h5ad(dataset_file)
  } else {
    stop("Unsupported file format")
  }

  if ("Seurat" %in% class(dataset)) 
  {
    if (!cell_type_label %in% colnames(dataset@meta.data)) {
      stop(glue("The specified cell type label column '{cell_type_label}' does not exist in the Seurat object. Please specify the correct column that corresponds to cell type labels in your train dataset."))
    }
    assay_data <- Seurat::GetAssayData(dataset, layer = "counts")
    cell_type_vector <- dataset@meta.data[[cell_type_label]]
  } 
  else if ("SingleCellExperiment" %in% class(dataset)) 
  {
    if (!cell_type_label %in% colnames(colData(dataset))) {
      stop(glue("The specified cell type label column '{cell_type_label}' does not exist in the SingleCellExperiment object. Please specify the correct column that corresponds to cell type labels in your train dataset."))
    }
    assay_data <- counts(dataset)
    cell_type_vector <- colData(dataset)[[cell_type_label]]
  } 
  else if (inherits(dataset, "Anndata")) {
    if (!cell_type_label %in% colnames(dataset$obs)) {
      stop(glue("The specified cell type label column '{cell_type_label}' does not exist in the AnnData object. Please specify the correct column that corresponds to cell type labels in your train dataset."))
    }
    assay_data <- dataset$X$toarray()
    cell_type_vector <- dataset$obs[[cell_type_label]]$to_list()
  } 
  else {
    stop("Unsupported object type")
  }

  sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
  dataset <- dataset[sel,]

  if ("Seurat" %in% class(dataset)) {
    modality.filt <- as.matrix(Seurat::GetAssayData(dataset, layer = "counts"))
  }
  if ("SingleCellExperiment" %in% class(dataset)) {
    modality.filt <- as.matrix(counts(dataset))
  }

  cty <- factor(cell_type_vector)
  cty <- as.character(cty)

  modality.filt <- log2(modality.filt + 1)
  modality.filt_scaled <- scale(modality.filt)

  num_samples_per_cty <- table(cty)

  return(list(modality = modality.filt_scaled, modality_noscale = modality.filt, cty = cty, num_samples_per_cty = num_samples_per_cty))
}

preprocess_dataset_test <- function(dataset_file) {  
  file_extension <- tools::file_ext(dataset_file)
  
  if (file_extension == "rds") {
    dataset <- readRDS(dataset_file)
  } else if (file_extension == "h5ad") {
    anndata <- import("anndata", convert = FALSE)
    dataset <- anndata$read_h5ad(dataset_file)
  } else {
    stop("Unsupported file format")
  }

  if ("Seurat" %in% class(dataset)) {
    assay_data <- Seurat::GetAssayData(dataset, layer = "counts")
  } else if ("SingleCellExperiment" %in% class(dataset)) {
    assay_data <- counts(dataset)
  } else if (inherits(dataset, "Anndata")) {
    assay_data <- dataset$X$toarray()
  } else {
    stop("Unsupported object type")
  }

  sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
  dataset <- dataset[sel,]

  if ("Seurat" %in% class(dataset)) {
    modality.filt <- as.matrix(Seurat::GetAssayData(dataset, layer = "counts"))
  } else if ("SingleCellExperiment" %in% class(dataset)) {
    modality.filt <- as.matrix(counts(dataset))
  }

  modality.filt <- log2(modality.filt + 1)
  modality.filt_scaled <- scale(modality.filt)

  return(list(modality = modality.filt_scaled, modality_noscale = modality.filt))
}

dataset_files <- c(train_file)
dataset_files1 <- c(test_file)

preprocessed_datasets <- lapply(dataset_files, function(x) preprocess_dataset_train(x, cell_type_label))
preprocessed_datasets1 <- lapply(dataset_files1, preprocess_dataset_test)

# Ensure train and test daatset have same features
common_features <- sort(intersect(rownames(preprocessed_datasets[[1]]$modality), rownames(preprocessed_datasets1[[1]]$modality)))
preprocessed_datasets[[1]]$modality <- preprocessed_datasets[[1]]$modality[common_features,]
preprocessed_datasets1[[1]]$modality <- preprocessed_datasets1[[1]]$modality[common_features,]


common_features <- sort(intersect(rownames(preprocessed_datasets[[1]]$modality_noscale), rownames(preprocessed_datasets1[[1]]$modality_noscale)))
preprocessed_datasets[[1]]$modality_noscale <- preprocessed_datasets[[1]]$modality_noscale[common_features,]
preprocessed_datasets1[[1]]$modality_noscale <- preprocessed_datasets1[[1]]$modality_noscale[common_features,]

# Ensuring Feature order is same in both train and test datasets
if ((identical(rownames(preprocessed_datasets[[1]]$modality), rownames(preprocessed_datasets1[[1]]$modality))) & identical(rownames(preprocessed_datasets[[1]]$modality_noscale), rownames(preprocessed_datasets1[[1]]$modality_noscale)))
{
  print("Feature order is same in train and test...")
} else 
{
  print("Error: Train and test do not have the same features/feature order...Exiting!")
  quit(status = 1)
}



##############################################

# Save processed data
for (dataset_idx in seq_along(preprocessed_datasets)) {
  print(glue("Processing train data..."))

  dataset <- preprocessed_datasets[[dataset_idx]]
  modality.filt <- dataset$modality
  modality.filt_noscale <- dataset$modality_noscale
  cty <- dataset$cty

  dataset_folder <- glue("Input_Processed")
  dir.create(dataset_folder, showWarnings = FALSE, recursive = TRUE)
  split_folder <- glue("{dataset_folder}/split_1")
  dir.create(split_folder, showWarnings = FALSE)

  # Write train data
  write_h5_scJoint(exprs_list = list(modality = modality.filt),
                   h5file_list = c(glue("{split_folder}/{data_type}_train.h5")))

  write_h5_scJoint(exprs_list = list(modality = modality.filt_noscale),
                   h5file_list = c(glue("{split_folder}/{data_type}_train_noscale.h5")))

  write_csv_scJoint(cellType_list = list(ct = cty),
                    csv_list = c(glue("{split_folder}/ct_train.csv")))

  # Handle test data
  print(glue("Processing test data...}"))

  dataset_test <- preprocessed_datasets1[[dataset_idx]]
  modality_test_filt <- dataset_test$modality

  dataset_test_folder <- glue("Input_Processed")
  dir.create(dataset_test_folder, showWarnings = FALSE, recursive = TRUE)
  split_test_folder <- glue("{dataset_test_folder}/split_1")
  dir.create(split_test_folder, showWarnings = FALSE)

  write_h5_scJoint(exprs_list = list(modality = modality_test_filt),
                   h5file_list = c(glue("{split_test_folder}/{data_type}_test.h5")))
}