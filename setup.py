from setuptools import setup, find_packages

version = {}
with open("./version.py") as fp:
    exec(fp.read(), version)

setup(
    name="MIE",
    version=version["__version__"],
    description="Mutual Information-based Clustering Analysis",
    url="https://github.com/jyyulab/MICA",
    author="Liang Ding, Chenxi Qian, Alireza Khatamian",
    author_email="liang.ding@stjude.org",
    license="Apache License 2.0",
    install_requires=[
        "pandas >= 0.22.0",
        "numpy >= 1.14.2",
        "scikit-learn >= 0.19.1",
        "matplotlib >= 2.2.2",
        "scipy >= 1.0.1",
    ],
    python_requires="==3.6.1",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    package_data={
        "MIE": ["test_data/*"],
    },
    include_package_data=True,
    test_suite="tests",
    entry_points={
        "console_scripts": [
            "mica=MICA.mica:main",
        ]
    },
)
