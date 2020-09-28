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
    * [cwltool>=1.0.2](https://github.com/common-workflow-language/cwltool)
    * [h5py>=2.10.0](https://www.h5py.org/)
    * [anndata>=0.7.4](https://anndata.readthedocs.io/en/latest/index.html#)
    * [scanpy>=1.6.0](https://scanpy-tutorials.readthedocs.io/en/latest/index.html)
* [cwlexec>=0.2.2](https://github.com/IBMSpectrumComputing/cwlexec) (required for running on IBM LSF)

Note: the pipeline is written in [common workflow language](https://www.commonwl.org/) for multi-platform compatibility.


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

Note: cwlexec requires manual installation if needed.


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
usage: mica [-h] {local,lsf} ...

MICA is a scalable tool to perform unsupervised scRNA-seq clustering analysis.

optional arguments:
  -h, --help   show this help message and exit

Subcommands:
  {local,lsf}  platforms
    local      run cwltool in a local workstation
    lsf        run cwlexec in a IBM LSF interface
```
`mica` workflow is implemented with CWL. It supports multiple computing platforms. 
We have tested it locally using cwltool and on an IBM LSF cluster using cwlexec. 
For the convenience, a python wrapper is developed for you to choose computing platform 
using subcommand.

The local mode (sjaracne local) runs in parallel by default using cwltool's --parallel option. 
To run it in serial, use --serial option.

To use LSF mode, editing the LSF-specific configuration file (a copy of the file is MICA/config/config_cwlexec.json)
to change the default queue and adjust memory reservation for each step is necessary. Consider 
increasing memory reservation for bootstrap step and consensus step if the dimension of your expression 
matrix file is large.


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
#### Running on a single machine (Linux/OSX)
`mica local 
-i ./test_data/inputs/PBMC_Demo_MICA_input_mini.txt 
-p "cwl_local" 
-k 3 4 
-o ./test_data/outputs/cwl_local/ 
--dist "spearman"`


#### Running on an IBM LSF cluster
`mica lsf 
-i ./test_data/inputs/PBMC_Demo_MICA_input_mini.txt 
-p "cwl_lsf" 
-k 3 4 
-o ./test_data/outputs/cwl_lsf/ 
-c ./MICA/config/config_cwlexec.json`


#### Rerun a failed workflow on an IBM LSF cluster
For a failed workflow, e.g., due to memory limit, MICA supports rerunning the workflow starting from the failed step 
with an input of the workflow ID using option `-r`. The default settings in the config file can be updated to 
increase the memory limit for the failed step. 

`mica lsf 
-r c62fb0be-cdb0-4bf6-b17f-2758ac0b51d7 
-j ./MICA/config/config_cwlexec.json`


## Reference
To be added
