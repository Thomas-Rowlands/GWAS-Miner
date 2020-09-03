import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GWAS_Miner-QuantumJava", # Replace with your own username
    version="0.0.1",
    author="Thomas Rowlands",
    author_email="thomas.s.rowlands@gmail.com",
    description="A python package for extracting key information from GWAS publications.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Thomas-Rowlands/GWAS-Miner/tree/Dev1",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)