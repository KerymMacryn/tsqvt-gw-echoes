"""
Setup script for TSQVT Gravitational Wave Echoes Analysis

Install with:
    pip install -e .

Author: Mohamed H.M. Makraini
"""

from setuptools import setup, find_packages
import os

# Read README
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="tsqvt-gw-echoes",
    version="0.1.0",
    author="Mohamed H.M. Makraini",
    author_email="mhamed34@alumno.uned.es",
    description="Search for gravitational wave echoes predicted by TSQVT in LIGO/Virgo data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/tsqvt-gw-echoes",
    project_urls={
        "Bug Tracker": "https://github.com/yourusername/tsqvt-gw-echoes/issues",
        "Documentation": "https://github.com/yourusername/tsqvt-gw-echoes/blob/main/README.md",
        "arXiv": "https://arxiv.org/abs/2501.xxxxx",
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "mypy>=0.950",
        ],
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
            "nbsphinx>=0.8.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "tsqvt-analyze=scripts.run_complete_analysis:main",
        ],
    },
    keywords="gravitational-waves ligo tsqvt quantum-gravity echoes black-holes",
    include_package_data=True,
    zip_safe=False,
)
