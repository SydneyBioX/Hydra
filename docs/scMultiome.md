# Single-cell multiomics

[![PyPI version](https://img.shields.io/pypi/v/Hydra-tools?color=orange)](https://pypi.org/project/Hydra-tools/)
[![Docs](https://img.shields.io/badge/docs-passing-brightgreen)](https://sydneybiox.github.io/Hydra/)
![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.0-blueviolet)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/SydneyBioX/Hydra?tab=MIT-1-ov-file#readme)


In this tutorial, we will use the example Skin SHARE-seq dataset (Ma et al.) introduced in the data processing step earlier. All analysis steps remain exactly the same as those used in the scRNA-seq analysis. For simplicity, we will demonstrate how to use <code><span style="color: red;">Hydra</span></code> to predict cell types in the query scMultiome dataset.

### **Running feature selection**

As described in the scRNA-seq analysis, the first step is to run feature selection. Since we are running <code><span style="color: blue;">fs</span></code> from the same directory containing the <code><span style="color: brown;">Input_Processed</span></code> directory, we can skip the <code><span style="color: blue;">base_dir</span></code> argument.

```bash
hydra --setting fs
```
<span style="display: block; height: 1px;"></span>

<div style="border-left: 1px solid purple; padding-left: 10px; overflow: auto; font-size: 14px; max-height: 700px">
<pre>
Thank you for using Hydra 😄, an interpretable deep generative tool for single-cell omics. Please refer to the full documentation available at https://sydneybiox.github.io/Hydra/ for detailed usage instructions. If you encounter any issues running the tool - Please open an issue on Github, and we will get back to you as soon as possible!!


===============================

Device to be used: CUDA 

===============================

INFO - 2025-08-12 05:38:39,834 - Starting to run 

INFO - 2025-08-12 05:38:39,834 - Training model... 

INFO - 2025-08-12 05:38:56,926 - The Dataset is: scRNA+scATAC 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 40/40 [00:13<00:00,  2.95it/s]
INFO - 2025-08-12 05:39:11,652 - 

Refining Model: 1 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 30/30 [00:05<00:00,  5.41it/s]
INFO - 2025-08-12 05:39:18,494 - 

Refining Model: 2 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 48/48 [00:05<00:00,  8.97it/s]
INFO - 2025-08-12 05:39:24,651 - 

Refining Model: 3 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 49/49 [00:06<00:00,  8.09it/s]
INFO - 2025-08-12 05:39:31,526 - 

Refining Model: 4 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 48/48 [00:06<00:00,  7.75it/s]
INFO - 2025-08-12 05:39:38,465 - 

Refining Model: 5 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 41/41 [00:04<00:00,  9.48it/s]
INFO - 2025-08-12 05:39:43,544 - 

Refining Model: 6 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 46/46 [00:05<00:00,  8.89it/s]
INFO - 2025-08-12 05:39:49,483 - 

Refining Model: 7 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 33/33 [00:04<00:00,  7.34it/s]
INFO - 2025-08-12 05:39:54,774 - 

Refining Model: 8 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 37/37 [00:04<00:00,  8.76it/s]
INFO - 2025-08-12 05:39:59,777 - 

Refining Model: 9 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 45/45 [00:04<00:00,  9.61it/s]
INFO - 2025-08-12 05:40:05,233 - 

Refining Model: 10 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 50/50 [00:04<00:00, 10.58it/s]
INFO - 2025-08-12 05:40:10,719 - 

Refining Model: 11 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 42/42 [00:04<00:00,  8.49it/s]
INFO - 2025-08-12 05:40:16,439 - 

Refining Model: 12 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 37/37 [00:04<00:00,  8.71it/s]
INFO - 2025-08-12 05:40:21,446 - 

Refining Model: 13 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 35/35 [00:04<00:00,  8.22it/s]
INFO - 2025-08-12 05:40:26,477 - 

Refining Model: 14 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 33/33 [00:04<00:00,  7.98it/s]
INFO - 2025-08-12 05:40:31,382 - 

Refining Model: 15 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 49/49 [00:05<00:00,  9.23it/s]
INFO - 2025-08-12 05:40:37,455 - 

Refining Model: 16 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 37/37 [00:04<00:00,  8.76it/s]
INFO - 2025-08-12 05:40:42,486 - 

Refining Model: 17 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 36/36 [00:04<00:00,  8.18it/s]
INFO - 2025-08-12 05:40:47,631 - 

Refining Model: 18 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 44/44 [00:04<00:00,  9.05it/s]
INFO - 2025-08-12 05:40:53,323 - 

Refining Model: 19 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 32/32 [00:03<00:00,  8.33it/s]
INFO - 2025-08-12 05:40:57,953 - 

Refining Model: 20 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 50/50 [00:05<00:00,  9.19it/s]
INFO - 2025-08-12 05:41:04,263 - 

Refining Model: 21 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 38/38 [00:04<00:00,  8.02it/s]
INFO - 2025-08-12 05:41:09,805 - 

Refining Model: 22 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 41/41 [00:04<00:00,  9.22it/s]
INFO - 2025-08-12 05:41:15,028 - 

Refining Model: 23 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 44/44 [00:04<00:00,  9.39it/s]
INFO - 2025-08-12 05:41:20,492 - 

Refining Model: 24 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 43/43 [00:05<00:00,  8.58it/s]
INFO - 2025-08-12 05:41:26,249 - 

Refining Model: 25 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 32/32 [00:04<00:00,  7.82it/s]
INFO - 2025-08-12 05:41:31,118 - 

Running feature selection... 

100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [01:33<00:00,  4.69s/it]
INFO - 2025-08-12 05:43:18,247 - Completed successfully! 
</pre>
</div>

<span style="display: block; height: 1px;"></span>

The results are stored in the <code><span style="color: brown;">Feature_Selection</span></code> sub-directory located within the <code><span style="color: brown;">Results</span></code> directory.

<span style="display: block; height: 1px;"></span>
<div style="border-left: 1px solid purple; padding-left: 10px; overflow: auto; font-size: 14px;">
<pre>
Results/
└── Feature_Selection
    └── Hydra-25
        ├── fs.ahighCD34+ bulge_Hydra.csv
        ├── fs.alowCD34+ bulge_Hydra.csv
        ├── fs.Basal_Hydra.csv
        ├── fs.Dermal Fibrobalst_Hydra.csv
        ├── fs.Dermal Papilla_Hydra.csv
        ├── fs.Dermal Sheath_Hydra.csv
        ├── fs.Endothelial_Hydra.csv
        ├── fs.Granular_Hydra.csv
        ├── fs.Hair Shaft-cuticle.cortex_Hydra.csv
        ├── fs.Infundibulum_Hydra.csv
        ├── fs.IRS_Hydra.csv
        ├── fs.Isthmus_Hydra.csv
        ├── fs.K6+ Bulge Companion Layer_Hydra.csv
        ├── fs.Macrophage DC_Hydra.csv
        ├── fs.Medulla_Hydra.csv
        ├── fs.Melanocyte_Hydra.csv
        ├── fs.ORS_Hydra.csv
        ├── fs.Spinous_Hydra.csv
        ├── fs.TAC-1_Hydra.csv
        └── fs.TAC-2_Hydra.csv
</pre>
</div>

<span style="display: block; height: 1px;"></span>


### **Automated cell type annotation**

```bash
hydra --setting annotation
```

<span style="display: block; height: 1px;"></span>

<div style="border-left: 1px solid purple; padding-left: 10px; overflow: auto; font-size: 14px; max-height: 700px">
<pre>
Thank you for using Hydra 😄, an interpretable deep generative tool for single-cell omics. Please refer to the full documentation available at https://sydneybiox.github.io/Hydra/ for detailed usage instructions. If you encounter any issues running the tool - Please open an issue on Github, and we will get back to you as soon as possible!!


===============================

Device to be used: CUDA 

===============================

INFO - 2025-08-12 05:44:16,051 - Starting to run 

Device to be used:  cuda 

INFO - 2025-08-12 05:44:18,360 - Training classifier: 1 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 13.90it/s]
INFO - 2025-08-12 05:44:37,639 - Training classifier: 2 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 25.91it/s]
INFO - 2025-08-12 05:44:38,162 - Training classifier: 3 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.12it/s]
INFO - 2025-08-12 05:44:38,685 - Training classifier: 4 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.18it/s]
INFO - 2025-08-12 05:44:39,206 - Training classifier: 5 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 24.45it/s]
INFO - 2025-08-12 05:44:39,729 - Training classifier: 6 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.03it/s]
INFO - 2025-08-12 05:44:40,234 - Training classifier: 7 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.10it/s]
INFO - 2025-08-12 05:44:51,072 - Training classifier: 8 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 24.89it/s]
INFO - 2025-08-12 05:44:51,584 - Training classifier: 9 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.42it/s]
INFO - 2025-08-12 05:44:52,104 - Training classifier: 10 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 24.46it/s]
INFO - 2025-08-12 05:44:52,620 - Training classifier: 11 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 24.38it/s]
INFO - 2025-08-12 05:44:53,137 - Training classifier: 12 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.58it/s]
INFO - 2025-08-12 05:45:04,071 - Training classifier: 13 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.78it/s]
INFO - 2025-08-12 05:45:04,580 - Training classifier: 14 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 25.97it/s]
INFO - 2025-08-12 05:45:15,420 - Training classifier: 15 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 25.88it/s]
INFO - 2025-08-12 05:45:15,951 - Training classifier: 16 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.10it/s]
INFO - 2025-08-12 05:45:16,476 - Training classifier: 17 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 24.83it/s]
INFO - 2025-08-12 05:45:17,026 - Training classifier: 18 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 24.84it/s]
INFO - 2025-08-12 05:45:17,566 - Training classifier: 19 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.09it/s]
INFO - 2025-08-12 05:45:18,094 - Training classifier: 20 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.16it/s]
INFO - 2025-08-12 05:45:18,624 - Training classifier: 21 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.25it/s]
INFO - 2025-08-12 05:45:19,151 - Training classifier: 22 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 25.95it/s]
INFO - 2025-08-12 05:45:19,684 - Training classifier: 23 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 25.72it/s]
INFO - 2025-08-12 05:45:20,216 - Training classifier: 24 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 25.83it/s]
INFO - 2025-08-12 05:45:20,756 - Training classifier: 25 

100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 26.37it/s]
INFO - 2025-08-12 05:45:21,286 - Annotating the query dataset... 

INFO - 2025-08-12 05:46:32,062 - Completed successfully! 
</pre>
</div>

<span style="display: block; height: 1px;"></span>

The annotation results generated will be stored in the directory <code><span style="color: brown;">Results/Annotation</span></code>.

<span style="display: block; height: 1px;"></span>

### **Plotting cell type predictions**

We will now use <span style="color: red;">Hydra</span></code> to plot predicted labels for the test dataset:

```bash
hydra --setting plot --predictions True --test scMultiome/Ma_Skin_scRNA_test.h5ad --modality rna --ctpredictions Results/Annotation/Hydra-25/cell_type_predicted_Hydra-25.csv
```
<span style="display: block; height: 1px;"></span>

<div style="border-left: 1px solid purple; padding-left: 10px; overflow: auto; font-size: 14px; max-height: 700px">
<pre>
Thank you for using Hydra 😄, an interpretable deep generative tool for single-cell omics. Please refer to the full documentation available at https://sydneybiox.github.io/Hydra/ for detailed usage instructions. If you encounter any issues running the tool - Please open an issue on Github, and we will get back to you as soon as possible!!


===============================

Device to be used: CUDA 

===============================

INFO - 2025-08-12 05:48:54,487 - Starting to run 

INFO - 2025-08-12 05:48:54,488 - Generating plot for Hydra predicted cell types... 

Warning: Data is of class matrix. Coercing to dgCMatrix.
Normalizing layer: counts
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Finding variable features for layer counts
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |======================================================================| 100%
PC_ 1 
Positive:  EBF1, SPARC, FBN1, LAMA2, COL6A3, ZEB2, COL3A1, FBXL7, COL5A2, FSTL1 
           PRRX1, COL4A1, MEG3, COL1A2, VIM, ZEB1, COL5A1, LDB2, PREX2, GRB10 
           COL1A1, COL15A1, ZFP521, COL4A2, NID1, PDE3A, IGFBP7, VCAN, THSD7A, PTPRM 
Negative:  CLCA3A2, CHL1, DSG1A, SERPINB2, ABCA12, IL20RA, SKINT6, SKINT5, SLC24A3, SOX6 
           PKIB, KRT77, IPMK, SPTLC3, EPGN, ANO9, FOSL1, RGS20, FCGBP, KRT10 
           ITGA6, KRT79, SKINT9, TPRG, LRP4, DMKN, KRT1, P3H2, ADGRL3, TNFAIP3 
PC_ 2 
Positive:  LAMA2, RSPO2, COL6A3, VCAN, COL5A2, COL5A1, HHIP, PRRX1, COL3A1, FREM1 
           CORIN, SVEP1, RSPO4, DCC, COL1A2, TBX15, CDH11, TMEM132C, PRR16, NRG2 
           MEG3, NCAM1, FBN1, KCNQ3, ABCA8A, COL14A1, FBN2, FGF10, NAV3, LINGO2 
Negative:  PTPRB, CD36, PECAM1, ADGRF5, FLT1, MECOM, ADGRL4, CDH5, FABP4, CYYR1 
           SHANK3, DYSF, CD93, KDR, AQP1, RASGRF2, SCARB1, FLI1, PTPRM, ETS1 
           EMCN, SDPR, ESAM, NRP1, MEOX2, FGD5, CCDC85A, PREX2, RASGRP3, S1PR1 
PC_ 3 
Positive:  MEG3, COL1A2, CCDC80, COL5A1, COL1A1, SPARC, GFPT2, ADGRD1, COL14A1, FBN1 
           IFI205, TNXB, NID1, COL3A1, CLEC3B, COL6A2, OPCML, ANXA1, ZEB1, ACKR3 
           TNFAIP6, EBF2, COL6A1, NOVA1, DPT, COL5A2, TMEFF2, DLC1, ESR1, BICC1 
Negative:  LEF1, ROBO1, CUX1, BNC2, NHS, GJA1, MKI67, TRPS1, TOP2A, HIST1H1B 
           NCL, CSGALNACT1, GLI2, NRP2, HIST1H1D, HIST1H1E, SLC30A1, GTF2IRD1, MYO5A, ROR2 
           KIF11, HIVEP3, VAV3, NEO1, HIST1H1C, KIF15, BRAF, JAG1, NOTCH1, CENPE 
PC_ 4 
Positive:  RSPO2, CORIN, RSPO4, DCC, HHIP, KCNQ3, FREM1, NRG2, ENPP2, RSPO3 
           MASP1, TMEM132C, LINGO2, MIR155HG, NCAM1, SHISA9, SNAP91, PTPRZ1, INHBA, FGF10 
           PLSCR4, PRR16, CHODL, PAPPA, NDNF, CNTN1, LEPR, GRID1, WIF1, PAPPA2 
Negative:  MEG3, DCN, FN1, COL1A1, ADGRD1, GFPT2, CELF2, CCDC80, COL5A1, COL1A2 
           BNC2, COL14A1, IFI205, OPCML, MKI67, SOX5, CLEC3B, NOVA1, ACKR3, TMEFF2 
           TOP2A, HIST1H1B, GJA1, CUX1, DPT, PI16, EBF2, FNDC1, COL6A2, ADAMTSL1 
PC_ 5 
Positive:  MKI67, ROBO1, TOP2A, VCAN, RSPO4, KIF11, RSPO2, KIF15, HIST1H1B, FREM1 
           WIF1, ANLN, HHIP, CENPE, LEF1, CENPF, MT2, KCNQ3, PTPRZ1, CORIN 
           COL23A1, TMEM132C, THSD7A, PLSCR4, PREX2, MASP1, CASC5, EBF1, HIST1H1A, FLI1 
Negative:  KRT79, NEBL, ADGRL3, PLXDC2, KRT17, SKINT3, ABCA12, ATP6V1C2, ANO9, DEFB6 
           CST6, SPINK5, TPRG, LDLRAD4, SLIT3, SBSN, CSNK2A2, AGPAT4, MCCC2, RGS20 
           KRT77, KLK7, ABI3BP, GAREM, EPHA4, KALRN, EPS8L1, DMKN, MYO5B, CHL1 
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
05:49:23 UMAP embedding parameters a = 0.9922 b = 1.112
05:49:23 Read 23597 rows and found 10 numeric columns
05:49:23 Using Annoy for neighbor search, n_neighbors = 30
05:49:23 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
05:49:24 Writing NN index file to temp file /tmp/RtmpDpQLyc/file1788b42fa8a747
05:49:24 Searching Annoy index using 1 thread, search_k = 3000
05:49:33 Annoy recall = 100%
05:49:33 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
05:49:34 Initializing from normalized Laplacian + noise (using RSpectra)
05:49:34 Commencing optimization for 200 epochs, with 1017522 positive edges
05:49:34 Using rng type: pcg
Using method 'umap'
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
05:49:43 Optimization finished
INFO - 2025-08-12 05:49:44,996 - Completed successfully!
</pre>
</div>

<br>

UMAP visualization of <code><span style="color: red;">Hydra</span></code> predicted cell types for scRNA-seq Skin dataset (Ma <i>et al</i>):


<span style="display: block; height: 1px;"></span>

<div style="text-align: center;">
    <img src="../images/Hydra_predicted_cell_types_rna_Ma.png" alt="Hydra Predicted cell types"; style="width: 600px; height: 400px";>
</div>

<span style="display: block; height: 1px;"></span>

```bash
hydra --setting plot --predictions True --test scMultiome/Ma_Skin_scATAC_test.h5ad --modality atac --ctpredictions Results/Annotation/Hydra-25/cell_type_predicted_Hydra-25.csv
```
<span style="display: block; height: 1px;"></span>

<div style="border-left: 1px solid purple; padding-left: 10px; overflow: auto; font-size: 14px; max-height: 700px">
<pre>
Thank you for using Hydra 😄, an interpretable deep generative tool for single-cell omics. Please refer to the full documentation available at https://sydneybiox.github.io/Hydra/ for detailed usage instructions. If you encounter any issues running the tool - Please open an issue on Github, and we will get back to you as soon as possible!!


===============================

Device to be used: CUDA 

===============================

INFO - 2025-08-12 05:50:04,545 - Starting to run 

INFO - 2025-08-12 05:50:04,545 - Generating plot for Hydra predicted cell types... 

Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
Warning: Data is of class matrix. Coercing to dgCMatrix.
Normalizing layer: counts
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Finding variable features for layer counts
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |======================================================================| 100%
PC_ 1 
Positive:  PLEKHG1, CELF2, KIF6, FANCC, STK32B, SHC4, CD82, FOXP1, GPR156, COL1A2 
           RGL1, LEF1, SGMS2, KAT7, MROH8, COL5A1, GM15802, USP46, GTF2IRD1, BASP1 
           ZFP521, GM40383, CDH6, SPSB4, HAPLN1, GM29345, ST7L, ADAMTS2, GM26820, PSMB7 
Negative:  GM34047, GM19585, SKINT4, FABP7, GM48236, BC064078, OLFR420, CLEC2G, PHYKPL, CLEC2F 
           PPP2R2D, GM47964, GM49776, 4930544I03RIK, OLFR1341, MALRD1, 2510002D24RIK, 4930511M06RIK, SIRT1, LGALS3 
           GM47465, GM35037, GM43123, USP38, SLC7A12, GM44148, DEFA-PS8, GM28500, GM35240, ASPHD2 
PC_ 2 
Positive:  SHC4, COL1A2, BICC1, ARHGAP15, ZFP521, GM40438, GM40383, PAPPA2, THSD7A, MAPK9 
           ADAMTS2, PDE10A, COL5A1, TH, COL4A1, GM48616, RGL1, FNDC3B, ADCY2, GM45338 
           TBX15, HAGH, GM2415, DCLK1, CEP162, VIM, GM12153, SH2D1B1, KIF6, KIFAP3 
Negative:  KAT7, SGMS2, GTF2IRD1, GM34022, LEF1, GM49891, 9530036O11RIK, ROBO1, HBS1L, NFIL3 
           TTC9, 4930534D22RIK, FOXP1, GM35208, PLEKHG1, RAP1GAP2, AFP, STK32B, SNX32, TNFAIP8L3 
           GM26820, BNC2, RXYLT1, G630018N14RIK, SPSB4, GM8520, 4930553J12RIK, 4931406G06RIK, RNF180, TMEFF1 
PC_ 3 
Positive:  SGMS2, KAT7, CD82, STK32B, LEF1, PLEKHG1, GM34022, CELF2, CDH6, 4930534D22RIK 
           COL1A2, NFIL3, GM29345, HBS1L, GM35208, BASP1, RNF180, RXYLT1, SPSB4, 9530036O11RIK 
           SLC25A36, KIF6, SHC4, MYL12B, FAR2OS2, BCAT1, 4930553J12RIK, METTL21A, GM40383, GPR156 
Negative:  4930511M06RIK, GM34047, SKINT4, MCF2L, NFXL1, ANK, MALRD1, GM47465, GM48236, ABCG1 
           PDIA6, OLFR1341, GM49868, 1700018A04RIK, VTI1B, ANAPC2, CLEC2G, KLK1B4, FABP7, CLEC2F 
           OSGEP, GAL3ST3, 2510002D24RIK, GPATCH4, BC064078, RDH13, FXYD3, KRT18, GM47732, GM44442 
PC_ 4 
Positive:  DPH5, SLC30A7, KCNMB4, GM48616, SEMA3C, GM28053, SCEL, SDR42E1, HSD17B13, PHOX2A 
           HAGH, SPAAR, TTC16, DNAJB8, GM13391, DHX35, GNG11, LTA4H, GM27246, SIAE 
           9330175M20RIK, TBC1D20, SLC25A45, TMPRSS6, COL4A1, GM15544, TOR2A, CEP95, KCNK15, GM15802 
Negative:  COL1A2, GM12153, BICC1, GM40438, GM40383, ADAMTS2, MYL12B, MAPK9, GM45338, PAPPA2 
           CEP162, COL5A1, TBX15, ARHGAP15, TH, NOVA1, ADCY2, GM45883, IQSEC3, STARD3NL 
           SHC4, CDO1, GM35554, RDH1, CELF2, LRRC58, SH2D1B1, DCLK1, CFAP299, PROKR1 
PC_ 5 
Positive:  NCALD, NR3C2, ANK, 1700018A04RIK, GM47465, COMMD10, GM29040, DCLK1, ADCY1, PHF11C 
           GM36633, ACAD9, DDX4, GM47732, BNC2, SLC9C1, GM30954, C130026L21RIK, BTG3, TTC27 
           PAKAP, MALRD1, GM20744, CDK5, GM17771, SYCP1, GM42932, PLG, GM11680, SNX5 
Negative:  CFAP299, ADORA2B, SIRT1, CNOT10, GM19585, GM44058, GM49539, GM26008, OGT, 4930544I03RIK 
           TYSND1, BC064078, ROBO1, IRF1, PHYKPL, GM35037, PIGV, GM38207, HS3ST1, GM43123 
           TAF15, GM48662, TNFAIP8L3, FGD6, GM28500, ODR4, GM49776, GPR176, TRDC, OLFR420 
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
05:50:43 UMAP embedding parameters a = 0.9922 b = 1.112
05:50:43 Read 23597 rows and found 10 numeric columns
05:50:43 Using Annoy for neighbor search, n_neighbors = 30
05:50:43 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
05:50:44 Writing NN index file to temp file /tmp/RtmpytiF6Z/file1791f8a3aec57
05:50:44 Searching Annoy index using 1 thread, search_k = 3000
05:50:53 Annoy recall = 100%
05:50:53 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
05:50:54 Initializing from normalized Laplacian + noise (using RSpectra)
05:50:55 Commencing optimization for 200 epochs, with 993682 positive edges
05:50:55 Using rng type: pcg
Using method 'umap'
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
05:51:03 Optimization finished
INFO - 2025-08-12 05:51:04,433 - Completed successfully!
</pre>
</div>

<br>

UMAP visualization of <code><span style="color: red;">Hydra</span></code> predicted cell types for the scATAC-seq skin dataset (Ma <i>et al.</i>).

<span style="display: block; height: 1px;"></span>

<div style="text-align: center;">
  <img src="../images/Hydra_predicted_cell_types_atac_Ma.png" alt="Hydra predicted cell types" style="width: 600px; height: 400px;"> 
</div>

<span style="display: block; height: 1px;"></span>

Below, we include an external figure showing the agreement between Hydra predicted cell type labels and the authors’ labels

<div style="text-align: center;">
  <img src="../images/heatmap_CT_Ma.png" 
       alt="Heatmap comparing authors’ labels and Hydra-predicted labels" 
       style="width: 700px; height: 600px; margin-left: -100px;">

<span style="display: block; height: 1px;"></span>

<br>

---
<p style="text-align: left; font-size: 15px">
  Documentation by <a href="http://manojmw.github.io" target="_blank">Manoj M Wagle</a>
</p>
