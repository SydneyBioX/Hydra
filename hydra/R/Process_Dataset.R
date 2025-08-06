#!/usr/bin/env Rscript

##############################################

# Manoj M Wagle (USydney, MIT CSAIL)

##############################################


args <- commandArgs(trailingOnly = TRUE)

train_file <- args[1]
test_file <- if (length(args) >= 2 && args[2] != "None") args[2] else NULL
cell_type_label <- args[3]
# peak <- if (length(args) >= 5 && args[5] %in% c("TRUE", "True", "true")) TRUE else FALSE
data_type <- args[4]
batch_size <- args[5]


##############################################

suppressPackageStartupMessages({
  library(rhdf5)
  library(HDF5Array)
  library(glue)
  library(reticulate)
  library(Matrix)
})


##############################################

# Writing H5 files
write_h5 <- function(exprs_list, h5file_list) {
  
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

write_csv <- function(cellType_list, csv_list) {
  
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
preprocess_dataset_train <- function(dataset_file, cell_type_label, batch_size = 1000) {
  file_extension <- tools::file_ext(dataset_file)
  
  # Load the dataset based on file extension
  if (file_extension == "rds") {
    dataset <- readRDS(dataset_file)
  } else if (file_extension == "h5ad") {
    anndata <- import("anndata", convert = FALSE)
    dataset <- anndata$read_h5ad(dataset_file)
  } else {
    stop("Unsupported file format")
  }
  
  # Process based on dataset type
  if ("Seurat" %in% class(dataset)) {
    library(Seurat)
    if (!cell_type_label %in% colnames(dataset@meta.data)) {
      stop(glue("The specified cell type label column '{cell_type_label}' does not exist in the Seurat object. Please specify the correct column that corresponds to cell type labels in your reference dataset."))
    }
    assay_data <- Seurat::GetAssayData(dataset, layer = "counts")
    cell_type_vector <- dataset@meta.data[[cell_type_label]]
    
    # # Filtering logic if peak is FALSE
    # if (!peak) {
    #   sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    #   dataset <- dataset[sel, ]
    # }
    sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    dataset <- dataset[sel, ]
    modality.filt <- as(Seurat::GetAssayData(dataset, layer = "counts"), "sparseMatrix")
    rm(assay_data)
    rm(dataset)
    gc()
  } else if ("SingleCellExperiment" %in% class(dataset)) {
    library(SingleCellExperiment)
    if (!cell_type_label %in% colnames(colData(dataset))) {
      stop(glue("The specified cell type label column '{cell_type_label}' does not exist in the SingleCellExperiment object. Please specify the correct column that corresponds to cell type labels in your reference dataset."))
    }
    assay_data <- counts(dataset)
    cell_type_vector <- colData(dataset)[[cell_type_label]]
    
    # # Filtering logic if peak is FALSE
    # if (!peak) {
    #   sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    #   dataset <- dataset[sel, ]
    # }
    sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    dataset <- dataset[sel, ]
    modality.filt <- as(counts(dataset), "sparseMatrix")
    rm(assay_data)
    rm(dataset)
    gc()
  } else if (any(grepl("AnnData", class(dataset)))) {
    column_exists <-  py_to_r(dataset$obs$columns$`__contains__`(cell_type_label))
    if (!column_exists) {
      stop(glue("The specified cell type label column '{cell_type_label}' does not exist in the AnnData object. Please specify the correct column that corresponds to cell type labels in your reference dataset."))
    }
    if (py_has_attr(dataset$X, "toarray")) {
      assay_data <- t(py_to_r(dataset$X$toarray()))  # Convert sparse matrix to dense array
    } else {
      assay_data <- t(py_to_r(dataset$X))  # Already a dense array
    }
    gene_names <- py_to_r(dataset$var$index$to_list())
    rownames(assay_data) <- gene_names
    cell_type_vector <- py_to_r(dataset$obs[[cell_type_label]]$to_list())

    # Assign column names (cells or samples)
    cell_names <- py_to_r(dataset$obs$index$to_list())
    colnames(assay_data) <- cell_names
    
    # # Filtering logic if peak is FALSE
    # if (!peak) {
    #   sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    #   assay_data <- assay_data[sel, ]
    # }
    sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    assay_data <- assay_data[sel, ]
    modality.filt <- as(assay_data, "sparseMatrix")
    rm(assay_data)
    rm(dataset)
    gc()
  } else {
    stop("Unsupported object type")
  }

  # Initialize containers for batched processing
  num_batches <- ceiling(ncol(modality.filt) / batch_size)
  modality.filt_list <- vector("list", num_batches)
  modality.filt_scaled_list <- vector("list", num_batches)

  # Process data in batches
  for (i in seq_len(num_batches)) {
    batch_indices <- ((i - 1) * batch_size + 1):min(i * batch_size, ncol(modality.filt))
    batch_data <- modality.filt[, batch_indices, drop = FALSE]

    # Apply log2 transformation and scaling to the batch
    batch_data <- log2(batch_data + 1)
    batch_data_dense <- as.matrix(batch_data)
    batch_data_scaled <- scale(batch_data_dense)

    # Store the processed batch
    modality.filt_list[[i]] <- batch_data_dense
    modality.filt_scaled_list[[i]] <- batch_data_scaled
  }

  # Combine all batches into final datasets
  modality.filt <- do.call(cbind, modality.filt_list)
  rm(modality.filt_list)
  modality.filt_scaled <- do.call(cbind, modality.filt_scaled_list)
  rm(modality.filt_scaled_list)

  num_samples_per_cty <- table(factor(cell_type_vector))

  return(list(modality = modality.filt_scaled, modality_noscale = modality.filt, cty = as.character(factor(cell_type_vector)), num_samples_per_cty = num_samples_per_cty))
}

preprocess_dataset_test <- function(dataset_file, batch_size = 1000) {
  file_extension <- tools::file_ext(dataset_file)
  
  # Load the dataset based on file extension
  if (file_extension == "rds") {
    dataset <- readRDS(dataset_file)
  } else if (file_extension == "h5ad") {
    anndata <- import("anndata", convert = FALSE)
    dataset <- anndata$read_h5ad(dataset_file)
  } else {
    stop("Unsupported file format")
  }
  
  # Process based on dataset type
  if ("Seurat" %in% class(dataset)) {
    library(Seurat)
    assay_data <- Seurat::GetAssayData(dataset, layer = "counts")

    # # Filtering logic if peak is FALSE
    # if (!peak) {
    #   sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    #   dataset <- dataset[sel, ]
    # }
    sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    dataset <- dataset[sel, ]
    modality.filt <- as(Seurat::GetAssayData(dataset, layer = "counts"), "sparseMatrix")
    rm(assay_data)
    rm(dataset)
    gc()
  } else if ("SingleCellExperiment" %in% class(dataset)) {
    library(SingleCellExperiment)
    assay_data <- counts(dataset)

    # # Filtering logic if peak is FALSE
    # if (!peak) {
    #   sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    #   dataset <- dataset[sel, ]
    # }
    sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    dataset <- dataset[sel, ]
    modality.filt <- as(counts(dataset), "sparseMatrix")
    rm(assay_data)
    rm(dataset)
    gc()
  } else if (any(grepl("AnnData", class(dataset)))) {
    if (py_has_attr(dataset$X, "toarray")) {
      assay_data <- t(py_to_r(dataset$X$toarray()))  # Convert sparse matrix to dense array
    } else {
      assay_data <- t(py_to_r(dataset$X))  # Already a dense array
    }
    gene_names <- py_to_r(dataset$var$index$to_list())
    rownames(assay_data) <- gene_names

    # Assign column names (cells or samples)
    cell_names <- py_to_r(dataset$obs$index$to_list())
    colnames(assay_data) <- cell_names

    # # Filtering logic if peak is FALSE
    # if (!peak) {
    #   sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    #   assay_data <- assay_data[sel, ]
    # }
    sel <- names(which(rowSums(assay_data == 0) / ncol(assay_data) < 0.99))
    assay_data <- assay_data[sel, ]
    modality.filt <- as(assay_data, "sparseMatrix")
    rm(assay_data)
    rm(dataset)
    gc()
  } else {
    stop("Unsupported object type")
  }

  # Initialize containers for batched processing
  num_batches <- ceiling(ncol(modality.filt) / batch_size)
  modality.filt_list <- vector("list", num_batches)
  modality.filt_scaled_list <- vector("list", num_batches)

  # Process data in batches
  for (i in seq_len(num_batches)) {
    batch_indices <- ((i - 1) * batch_size + 1):min(i * batch_size, ncol(modality.filt))
    batch_data <- modality.filt[, batch_indices, drop = FALSE]

    # Apply log2 transformation and scaling to the batch
    batch_data <- log2(batch_data + 1)
    batch_data_dense <- as.matrix(batch_data)
    batch_data_scaled <- scale(batch_data_dense)

    # Store the processed batch
    modality.filt_list[[i]] <- batch_data_dense
    modality.filt_scaled_list[[i]] <- batch_data_scaled
  }

  # Combine all batches into final datasets
  modality.filt <- do.call(cbind, modality.filt_list)
  rm(modality.filt_list)
  modality.filt_scaled <- do.call(cbind, modality.filt_scaled_list)
  rm(modality.filt_scaled_list)

  return(list(modality = modality.filt_scaled, modality_noscale = modality.filt))
}

dataset_files <- c(train_file)
print("Now processing train dataset...")
preprocessed_datasets <- lapply(dataset_files, function(x) preprocess_dataset_train(x, cell_type_label))

if (!is.null(test_file)) {
  dataset_files1 <- c(test_file)
  print("Now processing test dataset...")
  preprocessed_datasets1 <- lapply(dataset_files1, preprocess_dataset_test)

  # Ensure train and test dataset have same features
  common_features <- sort(intersect(rownames(preprocessed_datasets[[1]]$modality), rownames(preprocessed_datasets1[[1]]$modality)))
  preprocessed_datasets[[1]]$modality <- preprocessed_datasets[[1]]$modality[common_features,]
  preprocessed_datasets1[[1]]$modality <- preprocessed_datasets1[[1]]$modality[common_features,]

  common_features <- sort(intersect(rownames(preprocessed_datasets[[1]]$modality_noscale), rownames(preprocessed_datasets1[[1]]$modality_noscale)))
  preprocessed_datasets[[1]]$modality_noscale <- preprocessed_datasets[[1]]$modality_noscale[common_features,]
  preprocessed_datasets1[[1]]$modality_noscale <- preprocessed_datasets1[[1]]$modality_noscale[common_features,]

  # Ensuring Feature order is same in both train and test datasets
  if ((identical(rownames(preprocessed_datasets[[1]]$modality), rownames(preprocessed_datasets1[[1]]$modality))) & identical(rownames(preprocessed_datasets[[1]]$modality_noscale), rownames(preprocessed_datasets1[[1]]$modality_noscale))) {
      print("Feature order is same in train and test...")
    } else {
      print("Error: Train and test do not have the same features/feature order...Exiting!")
      quit(status = 1)
    }
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
  write_h5(exprs_list = list(modality = modality.filt),
                   h5file_list = c(glue("{split_folder}/{data_type}_train.h5")))

  write_h5(exprs_list = list(modality = modality.filt_noscale),
                   h5file_list = c(glue("{split_folder}/{data_type}_train_noscale.h5")))

  write_csv(cellType_list = list(ct = cty),
                    csv_list = c(glue("{split_folder}/ct_train.csv")))

  if (!is.null(test_file)) {                  
    # Handle test data
    print(glue("Processing test data...}"))

    dataset_test <- preprocessed_datasets1[[dataset_idx]]
    modality_test_filt <- dataset_test$modality

    dataset_test_folder <- glue("Input_Processed")
    dir.create(dataset_test_folder, showWarnings = FALSE, recursive = TRUE)
    split_test_folder <- glue("{dataset_test_folder}/split_1")
    dir.create(split_test_folder, showWarnings = FALSE)

    write_h5(exprs_list = list(modality = modality_test_filt),
                    h5file_list = c(glue("{split_test_folder}/{data_type}_test.h5")))
    }                    
  }
