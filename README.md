# MICA

MICA is a mutual information-based clustering algorithm that consists of following steps:
1. Mutual information calculation
2. Decomposition and Dimension Reduction
3. Ensemble Clustering
4. Consensus Clustering
In thie pipeline, we utilized common workflow languages as a wrapper for multi-platform compatibility.

## Prerequisites
* [python==3.6.1](https://www.python.org/downloads/)
    * [pandas>=0.22.0](https://pandas.pydata.org/)
    * [numpy>=1.14.2](https://www.scipy.org/scipylib/download.html)
    * [scikit-learn>=0.19.1](http://scikit-learn.org/stable/install.html#)
    * [matplotlib>=2.2.2](https://matplotlib.org/users/installing.html)
    * [scipy>=1.0.1](https://www.scipy.org/install.html)
    * [cwltool>=1.0.2](https://github.com/common-workflow-language/cwltool)   
		or
    * [cwlexec>=0.2.1](https://github.com/IBMSpectrumComputing/cwlexec)
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

mica lsf \
-i [path_to_input_txt] \
-p [your_project_name] \
-k 3 4 \ # specify number of cluster
-o [path_to_outputs] \
-q standard # which queue

```

After the completion of the pipeline, MICA will generate following outputs:
1. Cell-cell Mutual information matrix 
2. Dimension reduced distance matrix 
3. Clustering results plot with clustering label mapped to each cluster
4. Clustering results txt file with visualization coordinates and clustering label

