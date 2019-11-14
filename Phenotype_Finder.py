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
import pprint

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
    #  Retrieve ontology terms for tagging
    get_ontology_terms()
    tagging_data = {}
    for (id, term) in hpoTerms:
        tagging_data[term] = "HP"
    # print(tagging_data)

    # Load study data
    file_data = []
    for filename in os.listdir(directory):
        with open(directory + "/" + filename, 'r') as file:
            file_data.append([filename, file.read()])

    studies = []
    for (pmid, xmlText) in file_data:
        studies.append(PreProcessing.strip_xml(None, pmid, xmlText))



    # for (id, term) in hpoSyns:
    #    tagging_data[term] = "HPS"

    #  Create array of text from main body of study
    for study in studies:
        data = []
        #test = ""
        for section in study.sections:
            data.append(section[:][1])
            #test = test + str(study.sections[:][1])
        #import SpaceJam
        #SpaceJam.process_text(test, tagging_data)
        #sys.exit()
        pre_processor = PreProcessing(np.array(data), tagging_data)
        tmp = []
        snps = []
        snp_count = 0
        for (word, tag) in pre_processor.tagged_text:
            if tag == 'RSID':
                snps.append(word)
            if tag in ['HP', 'TEST']:
                tmp.append((word, tag))
        snps = list(dict.fromkeys(snps))
        snp_count = len(snps)
        print("PMID " + study.pmid + " contains " + str(snp_count) + " unique SNP identifiers")
        print(tmp)
        pprint.pprint(pre_processor.chunked_text)
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
