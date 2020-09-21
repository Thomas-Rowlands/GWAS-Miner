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
        "python-dateutil",
        "rdflib",
        "owlready2",
        "lxml",
        "jsonschema",
        "rtgo",
        "networkx",
        "spacy",
        "en_core_sci_md @ https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.2.5/en_core_sci_md-0.2.5.tar.gz",
        "svglib",
        "reportlab",
        "PyQt5"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
