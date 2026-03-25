from setuptools import setup, find_packages

setup(
    name="plantvelo",
    version="0.1.0",
    description="Plant-specific RNA velocity with intron-retention aware three-state kinetic model",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "velocyto>=0.17",
        "numpy>=1.17",
        "loompy>=3.0",
        "click>=7.0",
        "pysam>=0.16",
    ],
    entry_points={
        "console_scripts": [
            "plantvelo=plantvelo.commands:cli",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
)
