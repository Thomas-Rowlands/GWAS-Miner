# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 08:47:22 2019

@author: Thomas Rowlands

Based on the tutorial by Sowmya Vivek
Found here: https://medium.com/analytics-vidhya/automated-keyword-extraction-from-articles-using-nlp-bfd864f41b34

"""
import DB
from Ontology import OntologyMap
import pandas
from DataPreparation import PreProcessing
from Study import Study
import numpy as np
import config
import sys


ontology = OntologyMap
connection = DB.Connection(config.server, config.user, config.password, config.database)
meshTerms = None
hpoTerms = None
hpoSyns = None
hpo2Mesh = None


def get_ontology_terms():
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


get_ontology_terms()

# Get Abstracts from data set
queryString = """SELECT identifier, Title, StudyAbstract,  
	FROM Study
    WHERE StudyAbstract NOT IN ("", "Not supplied")"""
resultCursor = connection.query(queryString)

# add header row
tmp = []
Studies = []
for (identifier, Title, studyAbstract) in resultCursor:
    tmp.append([identifier, Title, studyAbstract])
    study = Study(gwas_id=identifier, title=Title, abstract=studyAbstract)
    Studies.append(study)
resultCursor.close()

pre_processor = PreProcessing(np.array(tmp))
dataset = pre_processor.study_data


