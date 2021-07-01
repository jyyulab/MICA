# MICA
[![Build Status](https://travis-ci.com/jyyulab/MICA.svg?token=HDr9KWz2yFUbD2psHJxJ&branch=master)](https://travis-ci.com/jyyulab/MICA)

MICA is a mutual information-based non-linear clustering algorithm for single-cell RNA-seq data 
that consists of following steps:
1. Mutual information-based k-nearest neighbor graph construction
2. Graph embedding for dimension reduction
3. Clustering on dimension reduced space
4. UMAP or tSNE visualization


## Prerequisites
* [python>=3.7.0](https://www.python.org/downloads/) (developed and tested on python 3.7.4)
    * See [requirements.txt](https://github.com/jyyulab/MICA/blob/million/requirements.txt) file for other python library 
    dependencies


## Installation
#### Using conda to create a virtual environment (recommended)
The recommended method of setting up the required Python environment and dependencies 
is to use the [conda](https://conda.io/docs/) dependency manager:
```
$ conda create -n py374 python=3.7.4          # Create a python3.7 virtual environment
$ source activate py37                        # Activate the virtual environment
$ conda install --file requirements.txt       # Install dependencies
```

#### Install using pip
```
$ pip install MICA
```


#### Install from source
```
$ git clone https://github.com/jyyulab/MICA     # Clone the repo
$ cd MICA                                       # Switch to the MICA root directory
$ python setup.py install                       # Install MICA from source
$ mica -h                                       # Check if mica works correctly
```


## Usage
```
$ mica -h
usage: mica [-h] {auto,ge,mds} ...

MICA - Mutual Information-based Clustering Analysis tool.

optional arguments:
  -h, --help     show this help message and exit

subcommands:
  {auto,ge,mds}  versions
    auto         automatic version
    ge           graph embedding version
    mds          MDS version
```

#### Inputs
The main input for MICA is a tab-separated cells/samples by genes/proteins (rows are cells/samples) expression 
matrix or an [adata](https://anndata.readthedocs.io/en/latest/index.html) file after preprocessing.


#### Outputs
After the completion of the pipeline, `mica` will generate the following outputs:
* Clustering results plot with clustering label mapped to each cluster
* Clustering results txt file with visualization coordinates and clustering label


## Examples
#### Running MICA MDS version
`mica mds -i /research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/GoldernStd/Yan/Yan_MICA_input.txt -o 
/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/outputs/GoldernStd/Yan/mds_1_2 -pn Yan -nc 8 
-dk 12 13 14 15 16 17 18 19`

#### Running MICA GE version
`mica ge -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs/cwl_lsf 
-ar 10.0`

## Reference
To be added
