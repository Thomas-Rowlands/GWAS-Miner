# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 08:47:22 2019

@author: Thomas Rowlands

Based on the tutorial by Sowmya Vivek
Found here: https://medium.com/analytics-vidhya/automated-keyword-extraction-from-articles-using-nlp-bfd864f41b34

"""
import os
import DB
from Ontology import HPO, Mesh
from DataPreparation import PreProcessing
from DataStructures import Study
import numpy as np
import config
import sys
import pprint

connection = DB.Connection(config.user, config.password, config.database, config.server, config.port)
meshTerms = None
hpoTerms = None
hpoSyns = None
hpo2Mesh = None

def get_ontology_terms():
    global mesh_data, hpo_data, hpo_syns, hpo2Mesh
    # Update ontology file & extracted data.
    # Mesh.update_mesh_file()

    # Retrieve MeSH terms
    mesh_data = Mesh.get_mesh_phenotypes()

    # Retrieve HPO terms
    hpo_data = HPO.get_hpo()  # Array of lists (id, label string literal)

    # Retrieve HPO synonyms
    hpo_syns = HPO.get_hpo_synonyms()

    # Retrieve HPO to Mesh mappings
    query_string = """SELECT hpoID, meshID
                FROM hpo2mesh
                            """

    # hpo2Mesh = ontology.get_hpo_2_mesh(connection.query(query_string))


def process_studies(directory):
    #  Retrieve ontology terms for tagging
    get_ontology_terms()
    tagging_data = {"HPO": [], "MeSH": []}
    for (id, term) in hpo_data:
        tagging_data["HPO"].append(term)
    #for (id, synonym) in hpo_syns:
    #    tagging_data["HPO_Syn"].append(synonym)
    for (id, label) in mesh_data:
        tagging_data["MeSH"].append(label)
    print(tagging_data)
    # Load study data
    file_data = []
    for filename in os.listdir(directory):
        with open(directory + "/" + filename, 'r') as file:
            file_data.append([filename, file.read()])

    studies = []
    for (pmid, xmlText) in file_data:
        print("-------------------------\nProcessing study " + str(pmid) + "\n-------------------------")
        studies.append(PreProcessing.strip_xml(pmid, xmlText))



    # for (id, term) in hpoSyns:
    #    tagging_data[term] = "HPS"

    # Create array of text from main body of study
    from NLP import PhenoNLP
    pheno_nlp = PhenoNLP(tagging_data)
    for study in studies:
        corpus = study.abstract
        for section in study.sections:
            corpus += " " + section[:][1]
        pheno_nlp.process_corpus(corpus)
      # pre_processor = PreProcessing(np.array(data), tagging_data)

        sys.exit("Stopping after 1st study")


def main():
    args = sys.argv[1:]
    cores = 1
    docs = "study_data"
    for i in range(len(args)):
        if args[i] == "-cores":
            args[i + 1] = cores
        if args[i] == "-docs":
            args[i + 1] = docs
    process_studies(docs)


if __name__ == '__main__':
    main()
