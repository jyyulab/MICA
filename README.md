# MICA
[![Build Status](https://travis-ci.com/jyyulab/MICA.svg?token=HDr9KWz2yFUbD2psHJxJ&branch=master)](https://travis-ci.com/jyyulab/MICA)

MICA is a mutual information-based non-linear clustering algorithm for single-cell RNA-seq data 
that consists of following steps:
1. Mutual information-based k-nearest neighbor graph construction
2. Graph embedding for dimension reduction
3. Clustering on dimension reduced space
4. UMAP or tSNE visualization


## Prerequisites
* [python>=3.6.1](https://www.python.org/downloads/) (developed and tested on python 3.6.1)
    * See [requirements.txt](https://github.com/jyyulab/MICA/blob/million/requirements.txt) file for other python library 
    dependencies
* [node2vec](https://github.com/snap-stanford/snap/tree/master/examples/node2vec) (default) or 
[deepwalk](https://github.com/phanein/deepwalk) (Executable much be available in your environment)


## Installation
#### Using conda to create a virtual environment (recommended)
The recommended method of setting up the required Python environment and dependencies 
is to use the [conda](https://conda.io/docs/) dependency manager:
```
$ conda create -n py36 python=3.6.1           # Create a python3.6 virtual environment
$ source activate py36                        # Activate the virtual environment
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
usage: mica [-h] -i FILE -o DIR [-m STR] [-d INT] [-e FLOAT] [-v STR]

MICA is a Mutual Information-based nonlinear Clustering Analysis tool designed for scRNA-seq data. This version uses a graph embedding method for dimension reduction on MI-kNN graph.

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --input-file FILE
                        Path to an input file (h5ad file or tab-delimited text
                        file)
  -o DIR, --output-dir DIR
                        Path to final output directory
  -m STR, --dr-method STR
                        Dimension reduction method [node2vec | deepwalk]
                        (default: node2vec)
  -d INT, --dr-dim INT  Number of dimensions to reduce to (default: 12)
  -e FLOAT, --resolution FLOAT
                        Determines size of the communities. (default: 1.0)
  -v STR, --visual-method STR
                        Visualization embedding method [umap | tsne] (default:
                        umap)
```

#### Inputs
The main input for MICA is a tab-separated cells/samples by genes/proteins (rows are cells/samples) expression 
matrix or an [adata](https://anndata.readthedocs.io/en/latest/index.html) file after preprocessing.


#### Outputs
After the completion of the pipeline, `mica` will generate the following outputs:
* Clustering results plot with clustering label mapped to each cluster
* Clustering results txt file with visualization coordinates and clustering label


## Examples
#### Running MICA graph embedding version
`mica -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad 
-o ./test_data/outputs -d 12 -e 1.0`


## Reference
To be added
