# MICA
[![Build Status](https://travis-ci.com/jyyulab/MICA.svg?token=HDr9KWz2yFUbD2psHJxJ&branch=master)](https://travis-ci.com/jyyulab/MICA)

MICA is a mutual information-based clustering algorithm that consists of following steps:
1. Mutual information calculation
2. Dimension Reduction
3. Ensemble Clustering
4. Consensus Clustering


## Prerequisites
* [python==3.6.1](https://www.python.org/downloads/)
    * [pandas>=0.22.0](https://pandas.pydata.org/)
    * [numpy>=1.17.0](https://www.scipy.org/scipylib/download.html)
    * [scikit-learn>=0.19.1](http://scikit-learn.org/stable/install.html#)
    * [matplotlib==2.2.2](https://matplotlib.org/users/installing.html)
    * [scipy>=1.0.1](https://www.scipy.org/install.html)
    * [tables>=3.5.1](https://github.com/PyTables/PyTables)
    * [h5py>=2.10.0](https://www.h5py.org/)
    * [anndata>=0.7.4](https://anndata.readthedocs.io/en/latest/index.html#)
    * [scanpy>=1.6.0](https://scanpy-tutorials.readthedocs.io/en/latest/index.html)
    * [python-louvain>=0.14](https://github.com/taynaud/python-louvain)
    * [toil[cwl]@9af87a094f4a8b31b3c8a265e6dbcb2cc7595c91](https://toil.readthedocs.io/en/latest/running/cwl.html)

Note: the pipeline is written in [common workflow language](https://www.commonwl.org/) for multi-platform compatibility. Toil is used as the cwl-runner, which has batch scheduler implementation natively in python. It is installed at a specific git commit due to new implementation which is not yet in pypi.

## Installation
#### Using conda to create a virtual environment (recommended)
The recommended method of setting up the required Python environment and dependencies 
is to use the [conda](https://conda.io/docs/) dependency manager:
```
$ conda create -n py36 python=3.6.1           # Create a python3.6 virtual environment
$ source activate py36                        # Activate the virtual environment
$ conda install --file requirements.txt       # Install dependencies
$ pip install cwlref-runner                   # cwltool is not available in conda, install with pip
```

#### Install using pip
Note: installation requires pip>20.0. You can easily upgrade to the latest version of pip using `pip install --upgrade pip`.
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
usage: mica [-h] {local,batch} ...

MICA is a scalable tool to perform unsupervised scRNA-seq clustering analysis.

optional arguments:
  -h, --help   show this help message and exit

Subcommands:
  {local,batch}  platforms
    local        run cwltool in a local workstation
    batch        runs toil-cwl-runner for a given batch system
```
`mica` workflow is implemented with CWL. It supports multiple computing platforms. 
We have tested it locally using cwltool and on an IBM LSF cluster using cwlexec / toil. 
For the convenience, a python wrapper is developed for you to choose computing platform 
using subcommand.

The local mode (sjaracne local) runs in parallel by default using cwltool's --parallel option. 
To run it in serial, use --serial option.

To use MICA with a batch scheduler in an HPC environment, MICA implements toil with a variety of schedulers. The currently supported
schedulers are lsf, grid_engine, htcondor, torque, and slurm. The memory requirements are baked into the cwl workflows. By passing
the `-dm` flag, MICA will be run with memory doubling functionality on the systems using lsf as the scheduler. This will enable
toil to automatically resubmit jobs with double memory which fail due to hitting memory limits.
For all other schedulers, the cwl files `./cwl/mica.cwl` and `./cwl/mica_g.cwl` will need to be modified until toil supports the
`--overrides` flag.

#### Inputs
The main input for MICA is a tab-separated cells/samples by genes/proteins (rows are cells/samples) expression 
matrix and a list of cluster sizes to be specified in multiple runs of K-means.


#### Outputs
After the completion of the pipeline, `mica` will generate the following outputs:
1. Cell-cell mutual information matrix 
2. Dimension reduced distance matrix 
3. Clustering results plot with clustering label mapped to each cluster
4. Clustering results txt file with visualization coordinates and clustering label


## Examples
#### Running on a single machine using k-mean clustering (Linux/OSX)
`mica local 
-i ./test_data/inputs/PBMC_Demo_MICA_input_mini.txt 
-p "cwl_local" 
-k 3 4 
-o ./test_data/outputs/cwl_local/ 
--dist "spearman"`


#### Running on an IBM LSF cluster using k-mean clustering
```
mica batch \
-i ./test_data/inputs/PBMC_Demo_MICA_input_mini.txt \
-p "cwl_lsf" \
-k 3 4 \
-o ./test_data/outputs/cwl_lsf/ \
-s lsf \
-w wordir \
-js MICA_jobstore 
```

#### Running on an IBM LSF clustering using graph-based clustering
```
mica batch \
-i ./test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad \
-p "cwl_lsf_graph" \
-n 10 \
-o ./test_data/outputs/cwl_lsf/ \
-s lsf \
-w workdir \
-js MICA_jobstore
```

#### Rerun a failed workflow on an IBM LSF cluster
**THIS IS CURRENTLY BEING REIMPLEMENTED.** For a failed workflow, e.g., due to memory limit, MICA supports rerunning the workflow starting from the failed step 
with an input of the workflow ID using option `-r`. The default settings in the config file `config_cwlexec.json` may 
be updated to increase the memory limit for the failed step. 

```
mica lsf \
-r c62fb0be-cdb0-4bf6-b17f-2758ac0b51d7 \
-j ./MICA/config/config_cwlexec.json
```

## Reference
To be added
