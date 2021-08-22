# MICA
[![Build Status](https://travis-ci.com/jyyulab/MICA.svg?token=HDr9KWz2yFUbD2psHJxJ&branch=master)](https://travis-ci.com/jyyulab/MICA)

MICA is a clustering tool for single-cell RNA-seq data. MICA takes a preprocessed gene expression matrix as input and
efficiently cluster the cells.
MICA consists of the following main components:
1. Mutual information estimation for cell-cell distance quantification
2. Dimension reduction on the non-linear mutual information-based distance space
3. Consensus clustering on dimension reduced spaces
4. Clustering visualization and cell type annotation

MICA workflow:

<img src="images/MICA_workflow.png" width="500">


## Prerequisites
* [python>=3.7.0](https://www.python.org/downloads/) (developed and tested on python 3.7.4)
    * See [requirements.txt](https://github.com/jyyulab/MICA/blob/million/requirements.txt) file for other dependencies


## Installation
#### Using conda to create a virtual environment (recommended)
The recommended method of setting up the required Python environment and dependencies 
is to use the [conda](https://conda.io/docs/) dependency manager:
```
$ conda create -n py374 python=3.7.4          # Create a python3.7 virtual environment
$ source activate py37                        # Activate the virtual environment
$ conda install --file requirements.txt       # Install dependencies
```

<!--- 
#### Install using pip
```
$ pip install MICA
```
-->

#### Install from source
```
$ git clone https://github.com/jyyulab/MICA     # Clone the repo
$ cd MICA                                       # Switch to the MICA root directory
$ python setup.py install                       # Install MICA from source
$ mica -h                                       # Check if mica works correctly
```


## Usage
MICA workflow has two build-in dimension reduction methods. The auto mode (```mica``` or ```mica auto```) 
selects a dimension reduction method automatically based on the cell count of the preprocessed matrix. 
Users can select graph embedding method (```mica ge```)  or MDS (```mica mds```) method manually using subcommand 
```ge``` or ```mds``` respectively. 
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
matrix or an [anndata](https://anndata.readthedocs.io/en/latest/index.html) file after preprocessing.


#### Outputs
After the completion of the pipeline, `mica` will generate the following outputs:
* Clustering results plot with clustering label mapped to each cluster
* Clustering results txt file with visualization coordinates and clustering label


## Examples
#### Running MICA auto mode
MICA auto mode reduces the dimensionality using either the multidimensional scaling method (<= 3,000 cells) or 
the graph embedding method (> 3,000 cells), where the number of cells cutoff was chosen based on performance
evaluation of datasets of various sizes. 

`mica -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs`

#### Running MICA GE mode
MICA GE mode reduces the dimensionality using the graph embedding method. It sweeps a range of resolutions
of Louvain clustering algorithm. ```-ar``` parameter sets the upper bound of the range.

`mica ge -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs
-ar 10.0`

#### Running MICA MDS mode
MICA MDS mode reduces the dimensionality using the multidimensional scaling method. ```-pn``` parameter sets the
project name; ```-nc``` specifies the numbers of clusters (k in k-mean clustering algorithm); ```-dk``` is the
the numbers of dimensions used in performing k-mean clusterings in the dimension reduced matrix.

`mica mds -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o 
./test_data/outputs -pn Yan -nc 8 9 10 -dk 12 13 14 15 16 17 18 19`

## Reference
To be added
