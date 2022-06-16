from setuptools import setup, find_packages

version = {}
with open("./version.py") as fp:
    exec(fp.read(), version)

setup(
    name="MICA",
    version=version["__version__"],
    description="Mutual Information-based Clustering Analysis tool designed for scRNA-seq data",
    url="https://github.com/jyyulab/MICA",
    author="Liang Ding, Hao Shi",
    author_email="liang.ding@stjude.org, hao.shi@stjude.org",
    license="See LICENSE.md",
    install_requires=[
        "pandas==1.2.3",
        "numpy==1.20.1",
        "scikit-learn==0.24.1",
        "matplotlib==3.4.0rc3",
        "scipy==1.6.1",
        "tables==3.6.1",
        "cwltool==3.0.20210124104916",
        "h5py==3.2.1",
        "scanpy==1.7.1",
        "python-louvain==0.15",
        "fast_histogram==0.9",
        "pecanpy==1.0.1",
        "numba==0.53.1",
        "anndata==0.7.5",
        "umap-learn==0.5.1",
        "gensim==4.1.2",
        "pynndescent==0.5.2",
        "llvmlite==0.36.0",
        "setuptools==57.5.0",
        "rdflib_jsonld==0.5.0",
        ],
    python_requires="==3.7.6",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    test_suite="tests",
    entry_points={"console_scripts": ["mica=MICA.mica:main"]},
)
