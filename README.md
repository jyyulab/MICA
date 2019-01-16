# MICA

MICA is an ensemble and mutual information-based clustering algorithm that consists of three steps:
1. Decomposition and Dimension Reduction
2. Ensemble Clustering
3. Consensus Clustering

## Prerequisites
* [python==3.6.1](https://www.python.org/downloads/)
    * [pandas>=0.22.0](https://pandas.pydata.org/)
    * [numpy>=1.14.2](https://www.scipy.org/scipylib/download.html)
    * [scikit-learn>=0.19.1](http://scikit-learn.org/stable/install.html#)
    * [matplotlib>=2.2.2](https://matplotlib.org/users/installing.html)
    * [scipy>=1.0.1](https://www.scipy.org/install.html)
* [R==3.4.0](https://www.r-project.org/)
* [Phoenix] (Required only for parallel processing on a LSF cluster)

## Download
```
git clone https://github.com/jyyulab/MICA   # Clone the repo
```

## Installation
MICA can be run without an installation. So this is optional.

Setup a virtual python environment using [conda](https://conda.io/docs/) (not required but highly recommended). 
```
conda create -n py36 python=3.6.1           # Create a python3.6 virtual environment
source activate py36                        # Activate the virtual environment
```
Install MICA in the virtual environment.
```
cd MICA                     # Switch to the MICA root directory
python setup.py install     # Install MICA from source
mica -h                      # Check if mica works correctly
```

## Demo
This is an example of running MICA using the input file in `test_data` without an installation on
St. Jude Research LSF Cluster (assume you have the access to CompBio environment).
```
ssh hpc                             # ssh to a head node
hpcf_interactive                    # login an interactive node
setcbenv prod                       # set CompBio environment to prod
cbload phoenix                      # load CompBio modules
cbload util-python
module load python/3.6.1            # load python
module load R/3.4.0                 # load R
./MICA/mica.py Clust test_no_install_LSF ./test_data/inputs/test_local.whole.h5 ./test_data/inputs/test_local_mi.h5 
./test_data/outputs/test_no_install_LSF/ test_no_install_LSF --k 3 4 5 6 --perplexity 30 --retransformation False 
--host LSF
```
After the completion of the pipeline, the outputs are saved in `./test_data/outputs/MICA_test_no_intall_LSF`.
