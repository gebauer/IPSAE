from setuptools import setup, find_packages

setup(
    name="ipsae",
    version="0.1.8",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "biopython>=1.79",
        "scipy>=1.7.0",
    ],
    author="Jan Gebauer",
    author_email="jan.gebauer@tum.de",
    description="A library for calculating ipSAE scores for protein structures",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/gebauerj/IPSAE",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
)
