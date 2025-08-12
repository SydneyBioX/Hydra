# 🛠️ Data processing

[![PyPI version](https://img.shields.io/pypi/v/Hydra-tools?color=orange)](https://pypi.org/project/Hydra-tools/)
[![Downloads](https://img.shields.io/pypi/dm/Hydra-tools?color=blue)](https://pypi.org/project/Hydra-tools/)
[![Docs](https://img.shields.io/badge/docs-passing-brightgreen)](https://sydneybiox.github.io/Hydra/)
![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.0-blueviolet)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/SydneyBioX/Hydra?tab=MIT-1-ov-file#readme)


### **Overview**

Before using <code><span style="color: red;">Hydra</span></code>, the datasets need to be preprocessed. The required input format is **RDS/H5AD files** (Reference and Query) containing raw single-cell count data in one of the following formats:

- Seurat
- Single-cell Experiment (SCE)
- Anndata

!!! note
    If you are using an AnnData object, please make sure the default data matrix (`.X`) contains raw counts.

### **Input data**

<code><span style="color: red;">Hydra</span></code> supports both unimodal and multimodal single-cell omcis data. The supported modalities are:

- Single-cell RNA sequencing (scRNA-seq)
- Single-cell ATAC sequencing (scATAC-seq)
- Single-cell ADT sequencing (scADT-seq)
- scMultiome (Joint profiling of any of the above modalities - e.g., CITE-seq, TEA-seq, SHARE-seq, sci-CAR, and other similar technologies)

<span style="display: block; height: 1px;"></span>

### **Processing scRNA-seq datasets**

Below, we demonstrate the processing of two Lung scRNA-seq datasets. For simplicity, we will use only a subset of the data.

<span style="display: block; height: 1px;"></span>

#### Example dataset

##### Lung (<a href="https://doi.org/10.1186/s13059-019-1906-x" target="_blank">Madissoon <i>et al</i></a>)


```bash
wget "XXX"
```

<span style="display: block; height: 1px;"></span>

##### Lung (<a href="https://doi.org/10.1016/j.cell.2022.11.005" target="_blank">He <i>et al</i></a>)


```bash
wget "XXX"
```
<span style="display: block; height: 1px;"></span>

#### Process with Hydra

Datasets can be easily processed using the <code><span style="color: blue;">processdata</span></code> setting of <code><span style="color: red;">Hydra</span></code>.

```bash
# Example command
hydra --setting processdata --modality [rna, adt or atac] --train [Path to the reference dataset] --test [Path to the query dataset]

```

<span style="display: block; height: 1px;"></span>

!!! note
    By default, Hydra assumes the label column in the reference dataset as `cell_type`. If this is not the case, please specify using the <code><span style="color: blue;">celltypecol</span></code> argument. If you are only running feature selection, you can skip the query dataset.


<span style="display: block; height: 1px;"></span>

##### scRNA-seq
```bash
hydra --setting processdata --modality rna --train scRNA/Lung_Madissoon.h5ad --test scRNA/Lung_He.h5ad
```

<span style="display: block; height: 1px;"></span>

<div style="border-left: 1px solid purple; padding-left: 10px; overflow: auto; font-size: 14px; max-height: 700px">
<pre>
Thank you for using Hydra 😄, an interpretable deep generative tool for single-cell omics. Please refer to the full documentation available at https://sydneybiox.github.io/Hydra/ for detailed usage instructions. If you encounter any issues running the tool - Please open an issue on Github, and we will get back to you as soon as possible!!


===============================

Device to be used: CUDA 

===============================

INFO - 2025-08-10 07:13:14,402 - Starting to run 

INFO - 2025-08-10 07:13:14,402 - Processing datasets... 

[1] "Now processing train dataset..."
[1] "Now processing test dataset..."
[1] "Feature order is same in train and test..."
Processing train data...
                   group           name       otype dclass          dim
0                      /         matrix   H5I_GROUP                    
1                /matrix .data_dimnames   H5I_GROUP                    
2 /matrix/.data_dimnames              1 H5I_DATASET STRING         5950
3 /matrix/.data_dimnames              2 H5I_DATASET STRING        10728
4                /matrix       barcodes H5I_DATASET STRING         5950
5                /matrix           data H5I_DATASET  FLOAT 5950 x 10728
6                /matrix       features H5I_DATASET STRING        10728
                   group           name       otype dclass          dim
0                      /         matrix   H5I_GROUP                    
1                /matrix .data_dimnames   H5I_GROUP                    
2 /matrix/.data_dimnames              1 H5I_DATASET STRING         5950
3 /matrix/.data_dimnames              2 H5I_DATASET STRING        10728
4                /matrix       barcodes H5I_DATASET STRING         5950
5                /matrix           data H5I_DATASET  FLOAT 5950 x 10728
6                /matrix       features H5I_DATASET STRING        10728
Processing test data...
                   group           name       otype dclass          dim
0                      /         matrix   H5I_GROUP                    
1                /matrix .data_dimnames   H5I_GROUP                    
2 /matrix/.data_dimnames              1 H5I_DATASET STRING         6483
3 /matrix/.data_dimnames              2 H5I_DATASET STRING        10728
4                /matrix       barcodes H5I_DATASET STRING         6483
5                /matrix           data H5I_DATASET  FLOAT 6483 x 10728
6                /matrix       features H5I_DATASET STRING        10728
INFO - 2024-07-25 16:54:40,967 - Completed successfully! 
</pre>
</div>

<span style="display: block; height: 1px;"></span>

The processed files will be stored in the directory - <code><span style="color: brown;">Input_Processed</span></code>.

<br>


### **Processing scMultiome datasets**

Below, we demonstrate the processing of a SHARE-seq dataset by <a href="https://doi.org/10.1016/j.cell.2020.09.056">Ma <i>et al</i></a> that jointly profiles single-cell transcriptomic (RNA) and chomatin accessibility (ATAC) profiles. 

!!! note
    If using scMultiome, the Reference–Query files of each modality need to be processed separately (Refer below)

<span style="display: block; height: 1px;"></span>

#### Example dataset

##### Skin (<a href="https://doi.org/10.1016/j.cell.2020.09.056">Ma <i>et al</i></a>)


```bash
wget "XXX"
```

<span style="display: block; height: 1px;"></span>

#### Process with Hydra

##### scRNA-seq
```bash
hydra --setting processdata --modality rna --train scMultiome/Ma_Skin_scRNA_train.h5ad --test scMultiome/Ma_Skin_scRNA_test.h5ad
```

<span style="display: block; height: 1px;"></span>


<div style="border-left: 1px solid purple; padding-left: 10px; overflow: auto; font-size: 14px; max-height: 700px">
<pre>
Thank you for using Hydra 😄, an interpretable deep generative tool for single-cell omics. Please refer to the full documentation available at https://sydneybiox.github.io/Hydra/ for detailed usage instructions. If you encounter any issues running the tool - Please open an issue on Github, and we will get back to you as soon as possible!!


===============================

Device to be used: CUDA 

===============================

INFO - 2025-08-12 05:24:58,529 - Starting to run 

INFO - 2025-08-12 05:24:58,529 - Processing datasets... 

[1] "Now processing train dataset..."
[1] "Now processing test dataset..."
[1] "Feature order is same in train and test..."
Processing train data...
                   group           name       otype dclass         dim
0                      /         matrix   H5I_GROUP                   
1                /matrix .data_dimnames   H5I_GROUP                   
2 /matrix/.data_dimnames              1 H5I_DATASET STRING        4463
3 /matrix/.data_dimnames              2 H5I_DATASET STRING        8628
4                /matrix       barcodes H5I_DATASET STRING        4463
5                /matrix           data H5I_DATASET  FLOAT 4463 x 8628
6                /matrix       features H5I_DATASET STRING        8628
                   group           name       otype dclass         dim
0                      /         matrix   H5I_GROUP                   
1                /matrix .data_dimnames   H5I_GROUP                   
2 /matrix/.data_dimnames              1 H5I_DATASET STRING        4463
3 /matrix/.data_dimnames              2 H5I_DATASET STRING        8628
4                /matrix       barcodes H5I_DATASET STRING        4463
5                /matrix           data H5I_DATASET  FLOAT 4463 x 8628
6                /matrix       features H5I_DATASET STRING        8628
Processing test data...
                   group           name       otype dclass          dim
0                      /         matrix   H5I_GROUP                    
1                /matrix .data_dimnames   H5I_GROUP                    
2 /matrix/.data_dimnames              1 H5I_DATASET STRING        23597
3 /matrix/.data_dimnames              2 H5I_DATASET STRING         8628
4                /matrix       barcodes H5I_DATASET STRING        23597
5                /matrix           data H5I_DATASET  FLOAT 23597 x 8628
6                /matrix       features H5I_DATASET STRING         8628
INFO - 2025-08-12 05:26:28,556 - Completed successfully!
</pre>
</div>

<span style="display: block; height: 1px;"></span>

##### scATAC-seq

!!! note
    <code><span style="color: red;">Hydra</span></code> requires gene activity scores for scATAC data

<span style="display: block; height: 1px;"></span>

```bash
hydra --setting processdata --modality atac --train scMultiome/Ma_Skin_scATAC_train.h5ad --test scMultiome/Ma_Skin_scATAC_test.h5ad
```

<span style="display: block; height: 1px;"></span>


<div style="border-left: 1px solid purple; padding-left: 10px; overflow: auto; font-size: 14px; max-height: 700px">
<pre>
Thank you for using Hydra 😄, an interpretable deep generative tool for single-cell omics. Please refer to the full documentation available at https://sydneybiox.github.io/Hydra/ for detailed usage instructions. If you encounter any issues running the tool - Please open an issue on Github, and we will get back to you as soon as possible!!


===============================

Device to be used: CUDA 

===============================

INFO - 2025-08-12 05:29:34,072 - Starting to run 

INFO - 2025-08-12 05:29:34,072 - Processing datasets... 

[1] "Now processing train dataset..."
[1] "Now processing test dataset..."
[1] "Feature order is same in train and test..."
Processing train data...
                   group           name       otype dclass          dim
0                      /         matrix   H5I_GROUP                    
1                /matrix .data_dimnames   H5I_GROUP                    
2 /matrix/.data_dimnames              1 H5I_DATASET STRING         4463
3 /matrix/.data_dimnames              2 H5I_DATASET STRING        13755
4                /matrix       barcodes H5I_DATASET STRING         4463
5                /matrix           data H5I_DATASET  FLOAT 4463 x 13755
6                /matrix       features H5I_DATASET STRING        13755
                   group           name       otype dclass          dim
0                      /         matrix   H5I_GROUP                    
1                /matrix .data_dimnames   H5I_GROUP                    
2 /matrix/.data_dimnames              1 H5I_DATASET STRING         4463
3 /matrix/.data_dimnames              2 H5I_DATASET STRING        13755
4                /matrix       barcodes H5I_DATASET STRING         4463
5                /matrix           data H5I_DATASET  FLOAT 4463 x 13755
6                /matrix       features H5I_DATASET STRING        13755
Processing test data...
                   group           name       otype dclass           dim
0                      /         matrix   H5I_GROUP                     
1                /matrix .data_dimnames   H5I_GROUP                     
2 /matrix/.data_dimnames              1 H5I_DATASET STRING         23597
3 /matrix/.data_dimnames              2 H5I_DATASET STRING         13755
4                /matrix       barcodes H5I_DATASET STRING         23597
5                /matrix           data H5I_DATASET  FLOAT 23597 x 13755
6                /matrix       features H5I_DATASET STRING         13755
Warning message:
In write_csv(cellType_list = list(ct = cty), csv_list = c(glue("{split_folder}/ct_train.csv"))) :
  csv_list exists! will rewrite it.
INFO - 2025-08-12 05:32:17,150 - Completed successfully!
</pre>
</div>

<span style="display: block; height: 1px;"></span>

The processed files will be stored in the directory - <code><span style="color: brown;">Input_Processed</span></code>

<br>

---
<p style="text-align: left; font-size: 15px">
  Documentation by <a href="http://manojmw.github.io" target="_blank">Manoj M Wagle</a>
</p>