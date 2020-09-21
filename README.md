# GWAS Miner

## Overview
GWAS Miner was created as part of my PhD project at the University of Leicester, tackling the problem of extracting
 meaningful data from GWAS publication text.
 
## Features
- Extraction of genotype to phenotype associations, including genetic marker, disease and significance score (p-value).
- Visualisation of entity recognition and sentence structure within GWAS publication text.
- View publication statistics such as number of ontology disease term occurrences.

## Running GWAS Miner
GWAS Miner can be utilised both with a graphical user interface and through passing command line parameters.
Launching the graphical user interface can be done by passing the `-g` parameter to GWASMiner.py, allowing quick
 and easy access to all of it's features.
 
### Graphical User Interface 
```
python GWASMiner.py -g
```

### Command line
The following subset of features are available without launching the graphical user interface.

##### Process files within a directory
```
python GWASMiner.py -d <path_to_directory>
```

##### Update ontology cache
```
python GWASMiner.py -u
```

##### Visualise entities identified within a document
```
python GWASMiner.py -d <path_to_directory> -g "ents"
```

##### Visualise sentence dependencies within a document
```
python GWASMiner.py -d <path_to_directory> -g "sents"
```

## Dependencies
The following Python packages are required to run GWAS Miner, using at least python3.5 or later:
- [python-dateutil](https://pypi.org/project/python-dateutil/)
- [rdflib](https://pypi.org/project/rdflib/)
- [owlready2](https://pypi.org/project/Owlready2/)
- [lxml](https://pypi.org/project/lxml/)
- [jsonschema](https://pypi.org/project/jsonschema/)
- [rtgo](https://pypi.org/project/rtgo/)
- [networkx](https://pypi.org/project/networkx/)
- [spacy](https://pypi.org/project/spacy/)
- [SciSpaCy](https://allenai.github.io/scispacy/) pre-trained data model: `pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.2.5/en_core_sci_md-0.2.5.tar.gz`
- [svglib](https://pypi.org/project/svglib/)
- [reportlab](https://pypi.org/project/reportlab/)
- [PyQt5](https://pypi.org/project/PyQt5/)

## Contact
For issue reporting and feedback/recommendations please email Thomas Rowlands at [tr142@leicester.ac.uk](mailto:tr142@leicester.ac.uk).
