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
usage: mica [-h] -i FILE -o DIR -pn STR -nc INT [INT ...] [-pl STR] [-sr]
            [-cj FILE] [-rr STR] [-bs INT] [-df STR] [-dm STR] [-sn INT]
            [-tn INT] [-dk INT [INT ...]] [-dp INT] [-nn INT]
            [-dg INT [INT ...]] [-vm STR] [-md FLOAT] [-pp INT]

MICA - Mutual Information-based Clustering Analysis tool. This version uses a 
    multidimensional scaling method for dimension reduction.

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -i FILE, --input-file FILE
                        Path to an input file (h5ad file or tab-delimited text
                        file)
  -o DIR, --output-dir DIR
                        Path to final output directory
  -pn STR, --project-name STR
                        Project name/ID.
  -nc INT [INT ...], --num-clusters INT [INT ...]
                        Number of cluster to be specified in kmeans

platform arguments:
  -pl STR, --platform STR
                        Platform local (cwltool) or lsf (cwlexec) (default:
                        local)
  -sr, --serial-run     run cwltool in serial mode
  -cj FILE, --config-json FILE
                        [Required if use lsf platform] LSF-specific
                        configuration file in JSON format to be used for
                        workflow execution
  -rr STR, --rerun STR  Rerun an exited workflow with the given cwlexec
                        workflow ID.

additional arguments:
  -bs INT, --bootstrap INT
                        Maximum number of iterations per dimension (default:
                        10)
  -df STR, --dist-func STR
                        Method for distance matrix calculation [mi | euclidean
                        | spearman | pearson](default:mi)
  -dm STR, --dr-method STR
                        Transformation method used for dimension reduction
                        [MDS | PCA | LPL | LPCA] (default: MDS)
  -sn INT, --slice-size INT
                        Number of cells in each MI sub-matrix (default: 1000)
  -tn INT, --thread-number INT
                        Number of poolings used for multiple kmeans
                        iterations,usually equals to iterations_km (default:
                        10)
  -dk INT [INT ...], --dims-km INT [INT ...]
                        Dimensions used in k-mean clustering, array inputs are
                        supported (default: [19])
  -dp INT, --dims-plot INT
                        Number of dimensions used in visualization (default:
                        19)

graph clustering arguments:
  -nn INT, --num-neighbor INT
                        Number of neighbors of a cell for building k-nearest
                        neighbor graph
  -dg INT [INT ...], --dims-graph INT [INT ...]
                        Dimensions used in graph clustering, array inputs are
                        supported (default: 19)

visualization arguments:
  -vm STR, --visual-method STR
                        Visualization method UMAP or t-SNE (default: UMAP)
  -md FLOAT, --min-dist FLOAT
                        min_dist parameter in UMAP, minimum distance of points
                        in the embedded space (default: 0.6)
  -pp INT, --perplexity INT
                        Perplexity parameter for tSNE visualization
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
`mica -pl lsf -i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad -o ./test_data/outputs/cwl_lsf
-pn "cwl_lsf" -nc 3 4 -cj ./MICA/config/config_cwlexec.json`


## Reference
To be added
