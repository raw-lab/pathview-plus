"""
setup.py
Installation script for pathview-plus
"""

from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="pathview-plus",
    version="2.0.2",
    author="Richard Allen White III",
    description="Complete pathway visualization: KEGG + SBGN + highlighting + splines",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/raw-lab/pathview-plus",
    scripts=["bin/pathview-cli.py"],
    packages=["pathview"],
    package_dir=dict(pathview='lib'),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.10",
    install_requires=[
        "polars>=0.19.0",
        "numpy>=1.24.0",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        "Pillow>=10.0.0",
        "networkx>=3.1",
        "requests>=2.31.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "black>=23.0",
            "mypy>=1.0",
        ],
        "fast": [
            "lxml>=4.9.0",  # Faster XML parsing
        ],
    },
)
