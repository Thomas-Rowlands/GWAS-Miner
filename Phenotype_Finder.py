# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 08:47:22 2019

Main access point, 

@author: Thomas Rowlands
"""
import json
import os
from DataStructures import MasterLexicon
import sys
import logging
import time
import Experimental
import Ontology
import argparse
from NLP import Interpreter

current_time = time.strftime("%Y%m%d-%H%M%S")
logger = logging.getLogger("Phenotype Finder")
lexicon = None
gui = None


def __prepare_ontology_data(gui_ptr=None):
    global lexicon
    lexicon = MasterLexicon().parse(Ontology.get_tagging_data())


def load_nlp_object():
    global lexicon
    if not lexicon:
        __prepare_ontology_data()
    nlp = Interpreter(lexicon)
    return nlp


def update_gui_progress(qt_progress_signal, text):
    if not qt_progress_signal:
        return
    qt_progress_signal.emit(text)


def output_study_results(study, qt_finished_signal=None):
    if len(study.get_snps()) == 0 and qt_finished_signal:
        from GUI import QtFinishedResponse
        response = QtFinishedResponse(False, F"PMC{study.pmid}")
        qt_finished_signal.emit(response)
        return
    result_json = {"PMCID": study.pmid, "Mutations": [snp.__dict__ for snp in study.get_snps()]}
    with open(F"output/PMC{study.pmid}_result.json", "w", encoding="utf-8") as out_file:
        out_file.write(json.dumps(result_json, indent=4))
    if qt_finished_signal:
        from GUI import QtFinishedResponse
        response = QtFinishedResponse(True, F"PMC{study.pmid}")
        qt_finished_signal.emit(response)


def process_study(nlp, study, visualise=None, qt_progress_signal=None, qt_finished_signal=None):
    if not study:
        return False
    update_gui_progress(qt_progress_signal, F"Processing study {study.pmid}...")
    doc = nlp.process_corpus(
        nlp.replace_all_abbreviations(study.get_fulltext()))
    if visualise:
        update_gui_progress(qt_progress_signal, F"Launching visualisation of study {study.pmid}...")
        if visualise == "ents":
            nlp.display_ents(doc)
        elif visualise == "sents":
            nlp.display_structure([sent for sent in doc.sents])
    blah = []  # for debugging only
    for sent in doc.sents:
        blah.append(sent)
    update_gui_progress(qt_progress_signal, F"Identifying data from study {study.pmid}...")
    study.set_snps(nlp.extract_phenotypes(doc))
    marker_count = 0
    ontology_matches = 0
    temp = []
    output_study_results(study, qt_finished_signal)
    return True


def prepare_study(directory, file_name):
    study = Experimental.load_study(directory, file_name)
    return study


def process_studies(directory, visualise=None, qt_progress_signal=None, qt_finished_signal=None):
    """[Processes each file within the provided directory for GWAS information extraction.]

    Args: directory ([string]): [directory containing publication files.] visualise ([string], optional): [ents =
    entity visualisation, sents = dependency parsing visualisation]. Defaults to None.
    """
    # Structure ontology data ready for NLP tagging
    update_gui_progress(qt_progress_signal, "Gathering ontology data...")
    __prepare_ontology_data()
    # Initialise Interpreter module for NLP processing using master lexicon of ontology data.
    update_gui_progress(qt_progress_signal, "Loading NLP pipeline...")
    nlp = load_nlp_object()

    # Process each publication file in turn
    for file_name in os.listdir(directory):
        if not file_name.endswith("maintext.json"):
            continue
        study = prepare_study(directory, file_name)
        if not study:
            from GUI import QtFinishedResponse
            response = QtFinishedResponse(False, file_name)
            qt_finished_signal.emit(response)
            continue
        logger.info(F"Processing PMC {study.pmid}")
        result = process_study(nlp, study, visualise, qt_progress_signal, qt_finished_signal)
        if not result:
            update_gui_progress(qt_progress_signal, F"Unable to process study {file_name}. Skipping...")


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
    cores = args.cores  # not yet implemented
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

    # Setup folder for output files.
    if not os.path.isdir("output"):
        try:
            os.makedirs("output")
        except:
            sys.exit("Unable to create output folder")

    # Configure logger
    logging.basicConfig(filename=F"logs/{current_time}.log", format="%(asctime)s %(levelname)-8s %(message)s",
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
