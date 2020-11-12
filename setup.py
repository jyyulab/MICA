from setuptools import setup, find_packages

version = {}
with open("./version.py") as fp:
    exec(fp.read(), version)

setup(
    name="MICA",
    version=version["__version__"],
    description="Mutual Information based Clustering Analysis",
    url="https://github.com/jyyulab/MICA",
    author="Liang Ding, Chenxi Qian, Joao P. Veloso",
    author_email="liang.ding@stjude.org",
    license="See LICENSE.md",
    install_requires=[
        "pandas>=1.0",
        "numpy>=1.17.0",
        "scikit-learn>=0.19.1",
        "matplotlib>=3.1.2",
        "scipy>=1.0.1",
        "tables>=3.5.2",
        "h5py>=2.10.0",
        "anndata>=0.7.4",
        "scanpy>=1.6.0",
        "python-louvain>=0.14",
        "toil[cwl] @ git+https://github.com/DataBiosphere/toil.git@9af87a094f4a8b31b3c8a265e6dbcb2cc7595c91",
        ],
    python_requires=">=3.6.1",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    test_suite="tests",
    entry_points={"console_scripts": ["mica=MICA.mica:main"]},
)
