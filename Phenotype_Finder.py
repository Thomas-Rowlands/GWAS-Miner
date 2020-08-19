# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 08:47:22 2019

Main access point, 

@author: Thomas Rowlands
"""
import os
import DB
import Utility_Functions
from DataPreparation import PreProcessing
from Ontology import HPO, Mesh, EFO
from DataStructures import MasterLexicon
import config
import sys
from pprint import pformat
import logging
import time
import Experimental
from rtgo import ReadyThready

connection = DB.Connection(
    config.user, config.password, config.database, config.server, config.port)
meshTerms = None
hpoTerms = None
hpoSyns = None
hpo2Mesh = None

test = time.strftime("%Y%m%d-%H%M%S")

logger = logging.getLogger("Phenotype Finder")


def get_ontology_terms():
    """[Assigns ontology data from cache files ready for processing]
    """

    global mesh_data, hpo_data, hpo_syns, hpo2Mesh, efo_terms, efo_syns

    # Update ontology file & extracted data.
    # Mesh.update_mesh_file() #  Could be dangerous to use if structural changes are made.

    # Retrieve MeSH & HPO terms/descriptors/concepts

    # Retrieve HPO synonyms
    logger.info("Starting ontology data loading.")
    ont_data = ReadyThready.go_cluster([Mesh.get_descriptors, HPO.get_hpo_from_cache,
                                        EFO.get_efo_from_cache, HPO.get_syns])
    mesh_data = ont_data[0]
    hpo_data = ont_data[1]
    efo_terms, efo_syns = ont_data[2]
    hpo_syns = ont_data[3]
    logger.info("Finished ontology data loading.")

    # Retrieve HPO to Mesh mappings
    # query_string = """SELECT hpoID, meshID
    #             FROM hpo2mesh
    #                         """

    # hpo2Mesh = ontology.get_hpo_2_mesh(connection.query(query_string))


def update_ontology_cache():
    """[Updates the ontology cache files with data from the source ontology files.]
    """
    logger.info("Updating ontology cache files.")
    funcs = [Mesh.set_descriptors,
             HPO.set_terms,
             HPO.set_hpo_synonyms,
             EFO.set_terms]
    ReadyThready.go_cluster(funcs)
    logger.info("Finished updating ontology cache files.")


def process_studies(directory, visualise=None):
    """[Processes each file within the provided directory for GWAS information extraction.]

    Args:
        directory ([string]): [directory containing publication files.]
        visualise ([string], optional): [ents = entity visualisation, sents = dependency parsing visualisation]. Defaults to None.
    """
    # Retrieve ontology terms for tagging
    get_ontology_terms()
    tagging_data = {"HPO": [], "MeSH": [], "EFO": []}
    for (id, label) in hpo_data:
        tagging_data["HPO"].append(label)
    for (id, label) in hpo_syns:
        tagging_data["HPO"].append(label)
    for (id, label) in mesh_data:
        tagging_data["MeSH"].append(label)
    # for (id, label) in efo_terms:
    #     tagging_data["EFO"].append(label)
    # for (id, label) in efo_syns:
    #     tagging_data["EFO"].append(label)

    # Structure ontology data ready for NLP tagging
    from NLP import Interpreter
    lexicon = MasterLexicon().parse(tagging_data)

    # Initialise Interpreter module for NLP processing using master lexicon of ontology data.
    nlp = Interpreter(lexicon)

    # Process each publication file in turn
    for file_name in os.listdir(F"{directory}/full_text"):
        study = Experimental.load_study(directory, file_name)
        if not study:
            continue
        logger.info(F"Processing PMC {study.pmid}")
        # study = PreProcessing.strip_xml(pmid, xmlText) # unused while testing ICL files
        doc = nlp.process_corpus(
            nlp.replace_all_abbreviations(study.get_fulltext()))
        #print(study.get_fulltext())
       # nlp.display_ents(doc)
        if visualise:
            if visualise == "ents":
                nlp.display_ents(doc)
            elif visualise == "sents":
                nlp.display_structure([sent for sent in doc.sents])
        blah = []
        for sent in doc.sents:
            blah.append(sent)
        test = nlp.extract_phenotypes(doc)
        if test:
            logger.info(F"Phenotypes matched from study text: {pformat(test)}")
        marker_count = 0
        ontology_matches = 0
        temp = []
        if study.get_snps():
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
                # snp.phenotype = nlp.replace_abbreviations( May not be needed now, testing required.
                #     snp.phenotype, study.get_fulltext)
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
                logger.info(output)
            logger.info(F"Phenotypes matched from results tables: {Utility_Functions.Utility.remove_duplicates(temp)}")
            logger.info(F"Total markers /w P-vals: {str(marker_count)}")
            logger.info(F"Total markers matched to an ontology: {ontology_matches}")
        else:
            logger.info(F"No markers found for PMC{study.pmid}")
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
    skip = False
    visualise = None

    # Setup folder for log files.
    if not os.path.isdir("logs"):
        try:
            os.makedirs("logs")
        except:
            sys.exit("Unable to create logs folder")

    logging.basicConfig(filename=F"logs/{test}.log", format="%(asctime)s %(levelname)-8s %(message)s",
                        level=logging.INFO)

    for i in range(len(args)):
        if skip:
            skip = False
            continue
        if args[i] == "-cores":
            cores = args[i + 1]
            skip = True
        elif args[i] == "-docs":
            docs = args[i + 1]
            skip = True
        elif args[i] == "-update-ont":
            update_ontology_cache()
            sys.exit("Ontology cache updated.")
        elif args[i] == "-visualise":
            visualise = args[i + 1]
            skip = True
        else:
            sys.exit(F"Argument not recognised: {args[i]}")

    if docs == "":
        sys.exit("Document directory must be passed via -docs <path>")
    process_studies(docs, visualise)


if __name__ == '__main__':
    main()
