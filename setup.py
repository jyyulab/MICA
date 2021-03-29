from setuptools import setup, find_packages

version = {}
with open("./version.py") as fp:
    exec(fp.read(), version)

setup(
    name="MICA",
    version=version["__version__"],
    description="Mutual Information-based Clustering Analysis tool designed for scRNA-seq data",
    url="https://github.com/jyyulab/MICA",
    author="Liang Ding, Chenxi Qian, Joao P. Veloso",
    author_email="liang.ding@stjude.org",
    license="See LICENSE.md",
    install_requires=[
        "pandas>=1.1.3",
        "numpy>=1.19.2",
        "scikit-learn>=0.23.2",
        "matplotlib>=3.3.2",
        "scipy>=1.5.2",
        "tables>=3.5.2",
        "cwltool>=1.0.2",
        "h5py>=2.10.0",
        "scanpy>=1.6.0",
        "python-louvain>=0.14",
        "fast_histogram>=0.9",
        "node2vec>=0.4.1",
        "numba>=0.52.0",
        "anndata>=0.7.4",
        "umap-learn>=0.5.1",
        ],
    python_requires=">=3.6.1",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    test_suite="tests",
    entry_points={"console_scripts": ["mica=MICA.mica_ge:main"]},
)
