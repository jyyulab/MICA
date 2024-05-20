import os
import sys
import subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install

version = {}
with open("./version.py") as fp:
    exec(fp.read(), version)


class PreInstall(install):
    def run(self):
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-r', 'requirements.txt'])
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', './mihnsw'])
        install.run(self)
        

setup(
    name="MICA",
    version=version["__version__"],
    description="Mutual Information-based Clustering Analysis tool designed for scRNA-seq data",
    url="https://github.com/jyyulab/MICA",
    author="Liang Ding, Hao Shi, Jiayu Zhou",
    author_email="liang.ding@stjude.org, hao.shi@stjude.org, jiayu.zhou@stjude.org",
    license="See LICENSE.md",
    # install_requires=[
    #     "mihnsw==0.8.6",
    #     "numpy>=1.20.1,<=1.21.6",
    #     "scipy>=1.6.1,<=1.7.3",
    #     "h5py>=3.2.1, <=3.8.0",
    #     "scanpy==1.7.1",
    #     "python-louvain==0.15",
    #     "fast_histogram==0.9",
    #     "pecanpy==1.0.1",
    #     "numba==0.53.1",
    #     "anndata>=0.7.5, <=0.8.0",
    #     "networkx>=2.5.0, <=2.6.3"
    #     ],
    python_requires="==3.7.6",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    # test_suite="tests",
    entry_points={"console_scripts": ["mica=MICA.mica:main"]},
    cmdclass={
        'install': PreInstall
    }
)
