# How to install MICA hnsw-ann version

Step 1:

download MICA from Jiayu branch
```bash
module load git
git clone -b mica_hnsw_test https://github.com/jyyulab/MICA.git
```

Step 2:

unzip file
```bash
unzip MICA-mica_hnsw_test
cd MICA-mica_hnsw_test
```

Step 3:

create conda env

```bash
module load conda3
conda create -n mica_hnsw python=3.7.6
#press y if needed
source activate mica_hnsw
```

Step 4:

install customized hnswlib first:

```bash
cd hnswlib
pip install setuptools==57.5.0 anndata==0.7.5 pandas==1.2.3 scipy==1.6.1 typing==3.7.4.3 typing-extensions==3.7.4.3 numba==0.53.1 networkx==2.5
python -m pip install . numpy==1.20.1
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

in mica ge, you will see 3 new options

-nnt, —nn-type: select type of the nearest neighbor algo: knn or ann

-annm: set m value of ann, bigger the m is, more accuract the algo is, but more time will be cost

- `M` - the number of bi-directional links created for every new element during construction. Reasonable range for `M` is 2-100. Higher `M` work better on datasets with high intrinsic dimensionality and/or high recall, while low `M` work better for datasets with low intrinsic dimensionality and/or low recalls. The parameter also determines the algorithm's memory consumption, which is roughly `M * 8-10` bytes per stored element.
    
    As an example for `dim`=4 random vectors optimal `M` for search is somewhere around 6, while for high dimensional datasets (word embeddings, good face descriptors), higher `M` are required (e.g. `M`=48-64) for optimal performance at high recall. The range `M`=12-48 is ok for the most of the use cases. When `M` is changed one has to update the other parameters. Nonetheless, ef and ef_construction parameters can be roughly estimated by assuming that `M`*`ef_{construction}` is a constant.
    

-annef: set ef value,

 `ef` - the size of the dynamic list for the nearest neighbors (used during the search). Higher `ef` leads to more accurate but slower search. `ef` cannot be set lower than the number of queried nearest neighbors `k`. The value `ef` of can be anything between `k` and the size of the dataset.

Example:

mica ge -i /input_root -o /output_root -nnt ann -annm 8 -annef 1000

reference: 

[https://github.com/nmslib/hnswlib](https://github.com/nmslib/hnswlib)