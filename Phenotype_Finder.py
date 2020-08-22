# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 08:47:22 2019

Main access point, 

@author: Thomas Rowlands
"""
import os
import Utility_Functions
from DataStructures import MasterLexicon
import sys
from pprint import pformat
import logging
import time
import Experimental
import Ontology
import argparse

test = time.strftime("%Y%m%d-%H%M%S")
logger = logging.getLogger("Phenotype Finder")
lexicon = None
gui = None


def prepare_ontology_data():
    global lexicon
    lexicon = MasterLexicon().parse(Ontology.get_tagging_data())


def process_studies(directory, visualise=None):
    """[Processes each file within the provided directory for GWAS information extraction.]

    Args: directory ([string]): [directory containing publication files.] visualise ([string], optional): [ents =
    entity visualisation, sents = dependency parsing visualisation]. Defaults to None.
    """
    # Structure ontology data ready for NLP tagging
    from NLP import Interpreter
    global lexicon
    prepare_ontology_data()
    # Initialise Interpreter module for NLP processing using master lexicon of ontology data.
    nlp = Interpreter(lexicon)

    # Process each publication file in turn
    for file_name in os.listdir(directory):
        if not file_name.endswith("maintext.json"):
            continue
        study = Experimental.load_study(directory, file_name)
        if not study:
            continue
        logger.info(F"Processing PMC {study.pmid}")
        if gui:
            gui.update_progress_text(F"Processing {file_name}")
        # study = PreProcessing.strip_xml(pmid, xmlText) # unused while testing ICL files
        doc = nlp.process_corpus(
            nlp.replace_all_abbreviations(study.get_fulltext()))
        # print(study.get_fulltext())
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
    # Configure command line arguments
    parser = argparse.ArgumentParser(description='GWAS Information Extraction')
    parser.add_argument('-c', '--cores', type=int, default=1, help='Number of CPU cores to utilise. 0 = Max '
                                                                   'available. Default = 1.')
    parser.add_argument('-d', '--docs', type=str, default="study_data", help='Directory containing the study JSON '
                                                                             'files.')
    parser.add_argument('-u', '--update_ont', action='store_true', help='Update ontology cache files from the source '
                                                                        'ontology files.')
    parser.add_argument('-v', '--visualise', type=str, help='Start displacy visualisation server for entities or '
                                                            'dependencies by specifying ents or sents respectively.')
    parser.add_argument('-g', '--interface', action='store_true', help='Launch using the graphical user interface.')

    # Parse input arguments
    args = parser.parse_args()
    cores = args.cores
    docs = args.docs
    visualise = args.visualise
    using_gui = args.interface
    update_ont = args.update_ont

    # Setup folder for log files.
    if not os.path.isdir("logs"):
        try:
            os.makedirs("logs")
        except:
            sys.exit("Unable to create logs folder")

    # Configure logger
    logging.basicConfig(filename=F"logs/{test}.log", format="%(asctime)s %(levelname)-8s %(message)s",
                        level=logging.INFO)

    # Update ontology cache files if requested.
    if update_ont:
        Ontology.update_ontology_cache()

    # Begin running either the GUI or processing studies immediately.
    if using_gui:
        from GUI import MainForm
        global gui
        gui = MainForm()
        gui.open()
    else:
        process_studies(docs, visualise)


if __name__ == '__main__':
    main()
