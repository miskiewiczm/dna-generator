"""Setup configuration for dna-generator package."""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

setup(
    name="dna-generator",
    version="1.0.0",
    author="",  # TODO: Add author name
    author_email="",  # TODO: Add author email
    description="Deterministic DNA sequence generator with biochemical quality control",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/miskiewiczm/dna-generator",
    packages=find_packages(exclude=["tests", "tests.*", "examples", "benchmarks"]),
    package_data={
        'dna_generator': ['*.json'],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    python_requires=">=3.8",
    install_requires=[
        "dna-commons>=0.1.0",  # Our base library
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
        ],
        "vis": [
            "matplotlib>=3.5.0",  # For plot_validation if needed
        ],
    },
    entry_points={
        "console_scripts": [
            "dna-generator=dna_generator.__main__:main",
        ],
    },
    include_package_data=True,
)
