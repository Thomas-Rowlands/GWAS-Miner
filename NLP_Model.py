# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 08:47:22 2019

@author: Thomas Rowlands

Based on the tutorial by Sowmya Vivek
Found here: https://medium.com/analytics-vidhya/automated-keyword-extraction-from-articles-using-nlp-bfd864f41b34

"""
import os

import DB
from Ontology import OntologyMap
from DataPreparation import PreProcessing
from Study import Study
import numpy as np
import config
import sys

ontology = OntologyMap
connection = DB.Connection(config.user, config.password, config.database, config.server, config.port)
meshTerms = None
hpoTerms = None
hpoSyns = None
hpo2Mesh = None


def old_get_ontology_terms():
    global meshTerms, hpoTerms, hpoSyns, hpo2Mesh
    # Retrieve MeSH terms
    query_string = """SELECT meshID, termID, termName
        FROM mesh_term
        WHERE meshID in 
            (SELECT meshID
                FROM mesh_termtree
                WHERE LEFT(treeID,1) = 'C');
                            """
    meshTerms = ontology.get_mesh(connection.query(query_string))

    # Retrieve HPO terms
    query_string = """SELECT hpoID, name
                FROM hpo_term"""

    hpoTerms = ontology.get_hpo(connection.query(query_string))

    # Retrieve HPO synonyms
    query_string = """SELECT hpoID, synonymText
                FROM hpo_synonym
                            """
    hpoSyns = ontology.get_hpo_synonyms(connection.query(query_string))

    # Retrieve HPO to Mesh mappings
    query_string = """SELECT hpoID, meshID
                FROM hpo2mesh
                            """
    hpo2Mesh = ontology.get_hpo_2_mesh(connection.query(query_string))


def get_ontology_terms():
    global meshTerms, hpoTerms, hpoSyns, hpo2Mesh
    # Retrieve MeSH terms

    # meshTerms = ontology.get_mesh()

    # Retrieve HPO terms
    hpoTerms = ontology.get_hpo()  # Array of lists (id, label string literal)

    # Retrieve HPO synonyms
    # hpoSyns = ontology.get_hpo_synonyms()

    # Retrieve HPO to Mesh mappings
    query_string = """SELECT hpoID, meshID
                FROM hpo2mesh
                            """

    # hpo2Mesh = ontology.get_hpo_2_mesh(connection.query(query_string))


def process_studies(directory):
    # Load study data
    file_data = []
    for filename in os.listdir(directory):
        with open(directory + "/" + filename, 'r') as file:
            file_data.append([filename, file.read()])

    studies = []
    for (pmid, xmlText) in file_data:
        studies.append(PreProcessing.strip_xml(None, pmid, xmlText))

    get_ontology_terms()

    #  Retrieve ontology terms for tagging

    tagging_data = {}
    for (id, term) in hpoTerms:
        tagging_data[term] = "HP"
    print(tagging_data)
    # for (id, term) in hpoSyns:
    #    tagging_data[term] = "HPS"

    #  Create array of text from main body of study
    for study in studies:
        data = []
        for section in study.sections:
            data.append(section[1])
        pre_processor = PreProcessing(np.array(data), tagging_data)
        print(pre_processor.tagged_text)
        for (word, tag) in pre_processor.tagged_text:
            if 'HP' in tag:
                print((word, tag))
    # dataset = pre_processor.study_data


def main(directory="study_data"):
    process_studies(directory)


if __name__ == '__main__':
    main()
