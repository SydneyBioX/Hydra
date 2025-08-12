# 🚀 Installation guide

[![PyPI version](https://img.shields.io/pypi/v/Hydra-tools?color=orange)](https://pypi.org/project/Hydra-tools/)
[![Downloads](https://img.shields.io/pypi/dm/Hydra-tools?color=blue)](https://pypi.org/project/Hydra-tools/)
[![Docs](https://img.shields.io/badge/docs-passing-brightgreen)](https://sydneybiox.github.io/Hydra/)
![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.0-blueviolet)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/SydneyBioX/Hydra?tab=MIT-1-ov-file#readme)


### **Pre-requisites**

- Python >= 3.8
- R >= 4.0

<span style="display: block; height: 1px;"></span>

!!! note
    We recommend creating a separate environment such as **[Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html#)** to avoid package conflicts.

### **R dependencies**

Before installing <code><span style="color: red;">Hydra</span></code>, please make sure you install the following packages:

```r
mamba install -c conda-forge -c bioconda \
  bioconductor-hdf5array \
  bioconductor-singlecellexperiment \
  bioconductor-rhdf5 \
  r-seurat \
  r-glue \
  r-reticulate \
  r-matrix \
  r-ggplot2 \
  r-rlang \
  r-ggridges \
  r-anndata \
  bioconductor-zellkonverter
```



### **Installing Hydra**

Install <code><span style="color: red;">Hydra</span></code> via pip:

```bash
pip3 install hydra-tools
```

### **Verifying installation**

To check the <code><span style="color: red;">Hydra</span></code> installation, please run:

```bash
hydra --help
```
<span style="display: block; height: 1px;"></span>

You should see an output like this:
<div style="border-left: 1px solid purple; padding-left: 10px; overflow: auto; font-size: 14px;">
<pre>
Thank you for using Hydra 😄, an interpretable deep generative tool for single-cell omics. Please refer to the full documentation available at [xxx] for detailed usage instructions. If you encounter any issues running the tool - Please open an issue on Github, and we will get back to you as soon as possible!!


📍 NOTE 📍: You need to run feature selection (`fs`) on the train datatset before annotating the cell types in the query dataset. If you have already run feature selection on the train & want to annotate (`annotation`) a different related query dataset, please process the data (`processdata`) first and then provide the path to the directory containing this processed data.

usage: Hydra [-h] [--seed SEED] [--train TRAIN] [--test TEST] [--celltypecol CELLTYPECOL] [--modality {rna,adt,atac}] [--base_dir DIR] [--gene GENE] [--ctofinterest CTOFINTEREST]
             [--predictions PREDICTIONS] [--ctpredictions CTPREDICTIONS] [--processdata_batch_size PROCESSDATA_BATCH_SIZE] [--batch_size BATCH_SIZE] [--attr_batch_size ATTR_BATCH_SIZE]
             [--epochs EPOCHS] [--lr LR] [--gpu GPU] [--z_dim Z_DIM] [--hidden_rna HIDDEN_RNA] [--hidden_adt HIDDEN_ADT] [--hidden_atac HIDDEN_ATAC] [--num_models NUM_MODELS] --setting
             {processdata,fs,plot,annotation}
             ...

positional arguments:
  annotation_args       Additional arguments for annotation script

options:
  -h, --help                  show this help message and exit
  --seed SEED                 seed
  --train TRAIN               Path to the training dataset (Seurat, SCE or Anndata object)
  --test TEST                 Path to the test dataset (Seurat or SCE object)
  --celltypecol CELLTYPECOL
                              Cell type label column in your input dataset (Seurat, SCE or Anndata object). Default: `cell_type`
  --modality {rna,adt,atac}
                              Input data modality. Default: `rna`
  --base_dir DIR              Path to the directory containing processed data directory. Default: Current working directory
  --gene GENE                 Name of the gene whose expression is to be highlighted in the plot
  --ctofinterest CTOFINTEREST
                              Name of the cell type for which a ridgeline plot of gene expression should be generated
  --predictions PREDICTIONS
                              Generate UMAP plot for Hydra predicted cell types
  --ctpredictions CTPREDICTIONS
                              Path to the csv file containing cell types predicted by Hydra
  --processdata_batch_size PROCESSDATA_BATCH_SIZE
                              batch size for processing reference and query datasets
  --batch_size BATCH_SIZE
                              batch size for processing data during training
  --attr_batch_size ATTR_BATCH_SIZE
                              batch size for feature atrribution. Please adjust this based on your GPU memory
  --epochs EPOCHS             num of training epochs
  --lr LR                     learning rate
  --gpu GPU                   Please specify the GPU to use
  --z_dim Z_DIM               Number of neurons in latent space
  --hidden_rna HIDDEN_RNA
                              Number of neurons for RNA layer
  --hidden_adt HIDDEN_ADT
                              Number of neurons for ADT layer
  --hidden_atac HIDDEN_ATAC
                              Number of neurons for ATAC layer
  --num_models NUM_MODELS
                              Number of models for Ensemble Learning
  --setting {processdata,fs,plot,annotation}
                            `processdata` for processing input train and test Seurat, SCE or Anndata objects;
                            `fs` for feature selection to obtain cell-identity genes;
                            `plot` for generating UMAP plot of the dataset (Additionally, highlights gene expression when called with the `--gene` argument; Generates a ridgeline plot of expression of the specified gene in cell type of interest vs all other cell types when called with `--ctofinterest` argument; Generates a UMAP plot of Hydra predicted labels when called with `--predictions` argument);
                            `annotation` for automated annotation of the query dataset
</pre>
</div>

<br>

---
<p style="text-align: left; font-size: 15px">
  Documentation by <a href="http://manojmw.github.io" target="_blank">Manoj M Wagle</a>
</p>