from setuptools import setup, find_packages

version = {}
with open("./version.py") as fp:
    exec(fp.read(), version)

setup(
    name="MICA",
    version=version["__version__"],
    description="Mutual Information based Clustering analysis",
    url="https://github.com/jyyulab/MICA",
    author="Chenxi Qian, Liang Ding, Alireza Khatamian",
    author_email="chenxi.qian@stjude.org, liang.ding@stjude.org",
    license="See LICENSE.md",
    install_requires=[
        "pandas>=0.22.0",
        "numpy>=1.14.2",
        "scikit-learn>=0.19.1",
        "matplotlib>=2.2.2",
        "scipy>=1.0.1",
    ],
    python_requires="==3.6.1",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    test_suite="tests",
    entry_points={"console_scripts": ["mica=MICA.MICA:main"]},
)
