import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GWAS_Miner",  # Replace with your own username
    version="0.0.1",
    author="Thomas Rowlands",
    author_email="thomas.s.rowlands@gmail.com",
    description="A python package for extracting key information from GWAS publications.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Thomas-Rowlands/GWAS-Miner/tree/Dev1",
    packages=setuptools.find_packages(),
    install_requires=[
       "python-dateutil>=2.8.1",
       "rdflib>=5.0.0",
       "owlready2>=0.24",
       "lxml>=4.5.2",
       "jsonschema>=3.2.0",
       "rtgo>=0.0.7",
       "networkx>=2.5",
       "spacy>=2.3.2,<2.4",
       "en_core_sci_md @ https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.2.5/en_core_sci_md-0.2.5.tar.gz",
       "svglib>=1.0.1",
       "reportlab>=3.5.49",
       "PyQt5>=5.15.0",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
