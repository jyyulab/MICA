# How to install MICA hnsw-ann version

Step 1:

download MICA from Jiayu branch

```bash
module load git
git clone -b mica031 https://github.com/jyyulab/MICA.git
```

Step 2:

change current directory

```bash
cd MICA
```

Step 3:

create conda env

```bash
module load conda3
conda create -y -n mica_031 python=3.7.6
source activate mica_031
```

Step 4:

install customized hnswlib first:

```bash
cd hnswlib
pip install setuptools==65.6.3 anndata==0.8.0 pandas==1.2.3 scipy==1.7.3 typing==3.7.4.3 typing-extensions==4.7.1 numba==0.53.1 networkx==2.6.3
python -m pip install . numpy==1.21.6
```

To check hnswlib:

```bash
Python
```

```python
import hnswlib
# no error means ok
quit()
```

if any errors occur, try to switch conda in module

Step 5:

install other requirements

```bash
cd ../
pip install -r mica_hnsw.txt
```

Step 6:

install MICA

```bash
cd MICA
python setup.py install
```

Step 7:

check MICA:

```bash
mica -h
mica ge -h
mica mds -h
```

In mica ge, you will see 3 new options from version 0.2.3:

-nnt, —nn-type: select type of the nearest neighbor algo: knn or ann

-annm: set m value of ann, bigger the m is, more accuract the algo is, but more time will be cost

- `M` - the number of bi-directional links created for every new element during construction. Reasonable range for `M` is 2-100. Higher `M` work better on datasets with high intrinsic dimensionality and/or high recall, while low `M` work better for datasets with low intrinsic dimensionality and/or low recalls. The parameter also determines the algorithm's memory consumption, which is roughly `M * 8-10` bytes per stored element.
    
    As an example for `dim`=4 random vectors optimal `M` for search is somewhere around 6, while for high dimensional datasets (word embeddings, good face descriptors), higher `M` are required (e.g. `M`=48-64) for optimal performance at high recall. The range `M`=12-48 is ok for the most of the use cases. When `M` is changed one has to update the other parameters. Nonetheless, ef and ef_construction parameters can be roughly estimated by assuming that `M`*`ef_{construction}` is a constant.
    

-annef: set ef value,

 `ef` - the size of the dynamic list for the nearest neighbors (used during the search). Higher `ef` leads to more accurate but slower search. `ef` cannot be set lower than the number of queried nearest neighbors `k`. The value `ef` of can be anything between `k` and the size of the dataset.

In mica get, you will see 1 new option from version 0.3.1:

-cldis,--clustering-distance: select the distance base of the clustering process, cosine is considered to be more efficient and robust when the data is high dimonsional.

Also new in mica 0.3.1:

The package requirement range has been broaden, and it will be more convenient in the future.



Example:

mica ge -i /input_root -o /output_root -nnt ann -annm 8 -annef 1000 

## Common questions/issues:

1. Try to install MICA 0.3.1 but when checking, it’s still old version: 

Solve: make sure to create a new folder and change direction into there before you start Step 1, if you’ve have a like ‘MICA’ at the current direction.

```bash
mkdir MICA031
cd MICA031
#start step1
```

1. Illegal instruction(core dumped) or other problem when verifying hnswlib

Often because of the conda version or gcc version

try:

```bash
module load conda3/202311
module load gcc/10.2.0
#others may also works like gcc 9.5.0
```

1. Web package
Use pip install instead of conda install

reference: 

[https://github.com/nmslib/hnswlib](https://github.com/nmslib/hnswlib)
