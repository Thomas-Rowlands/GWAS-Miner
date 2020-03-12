# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 08:47:22 2019

@author: Thomas Rowlands
"""
import os
import DB
import Utility_Functions
from Ontology import HPO, Mesh, EFO
from DataPreparation import PreProcessing
from DataStructures import Study, MasterLexicon
import config
import sys
from pprint import pprint

connection = DB.Connection(config.user, config.password, config.database, config.server, config.port)
meshTerms = None
hpoTerms = None
hpoSyns = None
hpo2Mesh = None


def get_ontology_terms():
    global mesh_data, hpo_data, hpo_syns, hpo2Mesh, efo_terms, efo_syns
    # Update ontology file & extracted data.
    # Mesh.update_mesh_file() #  Could be dangerous to use if structural changes are made.

    # Retrieve MeSH & HPO terms/descriptors/concepts
    mesh_data = Mesh.get_descriptors()
    hpo_data = HPO.get_hpo_from_cache()
    efo_terms, efo_syns = EFO.get_efo_from_cache()

    # Retrieve HPO synonyms
    hpo_syns = HPO.get_syns()

    # Retrieve HPO to Mesh mappings
    query_string = """SELECT hpoID, meshID
                FROM hpo2mesh
                            """

    # hpo2Mesh = ontology.get_hpo_2_mesh(connection.query(query_string))


def update_ontology_cache():
    Mesh.set_descriptors()
    HPO.set_terms()
    HPO.set_hpo_synonyms()
    EFO.set_terms()


def process_studies(directory):
    # Retrieve ontology terms for tagging
    get_ontology_terms()
    tagging_data = {"HPO": [], "MeSH": [], "EFO": []}
    for (id, label) in hpo_data:
        tagging_data["HPO"].append(label)
    for (id, label) in hpo_syns:
        tagging_data["HPO"].append(label)
    for (id, label) in mesh_data:
        tagging_data["MeSH"].append(label)
    for (id, label) in efo_terms:
        tagging_data["EFO"].append(label)
    for (id, label) in efo_syns:
        tagging_data["EFO"].append(label)

    # Load study data
    file_data = []
    for filename in os.listdir(directory):
        with open(directory + "/" + filename, 'r') as file:
            file_data.append([filename, file.read()])

    from NLP import Interpreter
    lexicon = MasterLexicon().parse(tagging_data)
    nlp = Interpreter(lexicon)

    for (pmid, xmlText) in file_data[0:1]:
        print("-------------------------\nProcessing study " + str(pmid) + "\n-------------------------")
        study = PreProcessing.strip_xml(pmid, xmlText)
        doc = nlp.process_corpus(nlp.replace_all_abbreviations(study.get_fulltext()))
        nlp.display_ents(doc)
        nlp.display_structure([sent for sent in doc.sents])
        test = nlp.extract_phenotypes(doc)
        print("\nPhenotypes matched from study text:\n")
        pprint(test)
        print("\nPhenotypes matched from results tables:\n")
        marker_count = 0
        ontology_matches = 0
        temp = []
        for snp in study.get_snps():
            if (not snp.rs_identifier) or (not snp.gee_p_val) or (snp.phenotype is None):
                continue
            output = snp.rs_identifier
            temp.append(snp.rs_identifier)
            if snp.gee_p_val:
                output += F" | GEE => {snp.gee_p_val}"
            if snp.fbat_p_val:
                output += F" | FBAT => {snp.fbat_p_val}"
            if snp.misc_p_val:
                output += F" | MISC => {snp.misc_p_val}"
            output += " | Phenotype => "
            snp.phenotype = nlp.replace_abbreviations(snp.phenotype, study.original)
            if snp.phenotype is None:
                continue
            marker_count += 1
            matches = nlp.onto_match(snp.phenotype)
            if matches:
                for match in matches:
                    output += F"{match[0]}<{match[1]}>"
                ontology_matches += 1
            else:
                output += F"{snp.phenotype}<NO MATCH>"
            print(output)
        print(Utility_Functions.Utility.remove_duplicates(temp))
        print("Total markers /w P-vals: " + str(marker_count))
        print(F"Total markers matched to an ontology: {ontology_matches}")
        # corpus = study.abstract
        # for section in study.sections:
        #    corpus += " " + section[:][1]
        # nlp.process_corpus(corpus)
        # nlp.process_table_corpus(study.tables)
        # sys.exit("Stopping after 1st study")


def main():
    args = sys.argv[1:]
    cores = 1
    docs = ""
    for i in range(len(args)):
        if args[i] == "-cores":
            cores = args[i + 1]
        if args[i] == "-docs":
            docs = args[i + 1]
        if args[i] == "-update-ont":
            update_ontology_cache()
            sys.exit("Ontology cache updated.")
    if docs == "":
        sys.exit("Document directory must be passed via -docs <path>")
    process_studies(docs)


if __name__ == '__main__':
    main()
