# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 08:47:22 2019

Main access point, 

@author: Thomas Rowlands
"""
import logging
from datetime import datetime

# import bioc
# from bioc import BioCFileType

import BioC
import json
# import OutputConverter
from GWAS_Miner import OutputConverter
from TableExtractor import parse_tables, output_tables


def __load_config():
    """
    Retrieve configurable settings from the config.ini file.
    """
    from configparser import ConfigParser
    import sys
    config_settings = ConfigParser()
    try:
        config_settings.read("settings//config.ini")
    except FileNotFoundError:
        sys.exit(F"Could not locate the config file, has it been moved/deleted?")
    except IOError as e:
        sys.exit(F"Error reading config file: {e}")
    return config_settings


def __prepare_ontology_data():
    import Ontology
    return Ontology.get_master_lexicon()


logger = logging.getLogger("GWAS Miner")
config = __load_config()
lexicon = None
nlp = None
gui = None
is_cancelled = False
global output_xml
processed_files = []


def theme():
    return config.get(section="preferences", option="theme").lower()


def save_config():
    with open("settings//config.ini", "w") as f_out:
        config.write(f_out)


def load_nlp_object(qt_progress_signal=None, qt_finished_signal=None):
    from NLP import Interpreter
    global lexicon, nlp
    if not lexicon:
        update_gui_progress(qt_progress_signal, "Loading Ontology Data...")
        lexicon = __prepare_ontology_data()
    if not nlp:
        update_gui_progress(qt_progress_signal, "Loading NLP Pipeline...")
        nlp = Interpreter(lexicon)
        if qt_finished_signal:
            qt_finished_signal.emit(True)
            return
    return nlp


def update_gui_progress(qt_progress_signal, text):
    if not qt_progress_signal:
        return
    qt_progress_signal.emit(text)


def output_study_results(study, qt_study_finished_signal=None):
    # if qt_study_finished_signal:
    #     from GUI import QtFinishedResponse
    #     response = QtFinishedResponse(False, F"PMC{study['documents'][0]['id']}")
    #     qt_study_finished_signal.emit(response)
    #     return
    global output_xml
    # if output_xml:
    # Convert study to json, then convert to BioC XML document for output.
    # xml_study = json.dumps(study, default=BioC.ComplexHandler)
    # xml_study = bioc.loads(xml_study, BioCFileType.BIOC_JSON)
    # with open(F"output/PMC{study['documents'][0]['id']}_result.xml", "w", encoding="utf-8") as out_file:
    #     bioc.dump(xml_study, out_file)
    # # else:
    # with open(F"output/PMC{study['documents'][0]['id']}_result.json", "w", encoding="utf-8") as out_file:
    #     json.dump(study, out_file, default=BioC.ComplexHandler)
    OutputConverter.output_xml(json.dumps(study, default=BioC.ComplexHandler),
                               F"output/xml/{study['documents'][0]['id']}_result.xml")
    OutputConverter.output_json(study, F"output/json/{study['documents'][0]['id']}_result.json")

    if qt_study_finished_signal:
        from GUI import QtFinishedResponse
        response = QtFinishedResponse(True, F"{study['documents'][0]['id']}")
        qt_study_finished_signal.emit(response)


def get_study_visualisations(study, qt_progress_signal=None, qt_finished_signal=None):
    nlp_object = load_nlp_object(qt_progress_signal)
    if study is None:
        from GUI import QtFinishedResponse
        response = QtFinishedResponse(False, F"Failed to process document")
        qt_finished_signal.emit(response)
        return
    update_gui_progress(qt_progress_signal, F"Processing study {study['documents'][0]['id']}...")
    study_fulltext = "<!break!>".join([x['text'] for x in study['documents'][0]['passages']])
    doc = nlp_object.process_corpus(nlp_object.replace_all_abbreviations(study_fulltext))
    phenotype_stats = nlp_object.get_phenotype_stats(doc, lexicon)
    update_gui_progress(qt_progress_signal, "Generating HTML")
    entities = nlp_object.display_ents(doc, True, theme())
    dependencies = nlp_object.display_structure([sent for sent in doc.sents], True, theme())
    if qt_finished_signal:
        from GUI import QtFinishedResponse
        response = QtFinishedResponse(True, F"PMC{study['documents'][0]['id']}",
                                      [entities, dependencies, phenotype_stats])
        qt_finished_signal.emit(response)


def process_study(nlp, study, qt_progress_signal=None, qt_study_finished_signal=None):
    global is_cancelled
    if not study or is_cancelled:
        return False

    update_gui_progress(qt_progress_signal, F"Identifying data from study {study['documents'][0]['id']}...")
    t, m, p = 0, 0, 0
    # abbreviations = nlp.get_all_abbreviations(study_fulltext)
    for passage in study['documents'][0]['passages']:
        if "section_title_1" in passage["infons"].keys() and passage["infons"]["section_title_1"].lower() != "results":
            continue
        passage_text = passage['text']
        # if abbreviations:
        #     for abbrev in abbreviations:
        #         passage_text = passage_text.replace(abbrev[0], abbrev[1])

        doc = nlp.process_corpus(passage_text, )
        # for sent in doc.sents:
        #     training_sent = [x.label_ for x in sent.ents]
        #     if training_sent:
        #         if [x for x in training_sent if x[0] == "D"] and "RSID" in training_sent and "PVAL" in training_sent:
        #             training_string = sent.text_with_ws
        #             for ent in sent.ents:
        #                 ent_string = F"<!TRAIT:{ent.text_with_ws}!>" if ent.label_[
        #                                                                     0] == "D" else F"<!{ent.label_}:{ent.text_with_ws}!>"
        #                 training_string = training_string.replace(ent.text_with_ws, ent_string)
        #             with open("training_input/training_input.txt", "a+", encoding="utf-8") as f_in:
        #                 f_in.write(training_string + "\n")

        annotations = nlp.get_entities(doc)
        used_annots = []
        if annotations:
            for annot in annotations:
                loc = BioC.BioCLocation(offset=annot["offset"] + passage["offset"], length=annot["length"])
                current_datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
                if annot["text"] in used_annots:
                    for old_annot in passage['annotations']:
                        if old_annot.text == annot["text"]:
                            old_annot.locations.append(loc)
                if "RSID" not in annot["entity_type"] and "PVAL" not in annot["entity_type"]:
                    genomic_trait = BioC.BioCAnnotation(id=F"T{t}", infons={"type": "trait", "identifier": annot["id"],
                                                                            "annotator": "tr142@le.ac.uk",
                                                                            "updated_at": current_datetime},
                                                        locations=[loc], text=annot["text"])
                    passage['annotations'].append(genomic_trait)
                    t += 1
                elif "RSID" in annot["entity_type"]:
                    marker_identifier = BioC.BioCAnnotation(id=F"M{m}",
                                                            infons={"type": "genomic_marker", "identifier": annot["id"],
                                                                    "annotator": "tr142@le.ac.uk",
                                                                    "updated_at": current_datetime},
                                                            locations=[loc], text=annot["text"])
                    passage['annotations'].append(marker_identifier)
                    m += 1
                elif "PVAL" in annot["entity_type"]:
                    p_value = BioC.BioCAnnotation(id=F"P{p}", infons={"type": "significance", "identifier": annot["id"],
                                                                      "annotator": "tr142@le.ac.uk",
                                                                      "updated_at": current_datetime},
                                                  locations=[loc], text=annot["text"])
                    passage['annotations'].append(p_value)
                    p += 1
                used_annots.append(annot["text"])

        # relations = nlp.extract_phenotypes(doc)
        relations = None
        if relations:
            print(relations)
    output_study_results(study, qt_study_finished_signal)
    return True


def prepare_study(directory, file_name):
    import Experimental
    # study = Experimental.load_study(directory, file_name)
    study = Experimental.load_bioc_study(directory, file_name)
    return study


def process_studies(directory, visualise=None, shortlist=None, qt_progress_signal=None, qt_study_finished_signal=None):
    """[Processes each file within the provided directory for GWAS information extraction.]

    Args: directory ([string]): [directory containing publication files.] visualise ([string], optional): [ents =
    entity visualisation, sents = dependency parsing visualisation]. Defaults to None.
    """
    import os
    cancel_response = None

    if qt_progress_signal:
        from GUI import QtFinishedResponse
        cancel_response = QtFinishedResponse(True, "Cancelled.", 1)
    # Structure ontology data ready for NLP tagging
    update_gui_progress(qt_progress_signal, "Gathering ontology data...")

    if is_cancelled:
        qt_study_finished_signal.emit(cancel_response)
        return

    # Initialise Interpreter module for NLP processing using master lexicon of ontology data.
    update_gui_progress(qt_progress_signal, "Loading NLP pipeline...")
    nlp_object = load_nlp_object()

    if is_cancelled:
        qt_study_finished_signal.emit(cancel_response)
        return

    # Process each publication file in turn
    for file_name in os.listdir(directory):
        if shortlist:
            if file_name not in shortlist:
                continue
        if is_cancelled:
            qt_study_finished_signal.emit(cancel_response)
            return

        logger.info(F"Extracting data for file: {file_name}")
        update_gui_progress(qt_progress_signal, F"Extracting data for file: {file_name}")

        if file_name.endswith("tables.json"):
            tables, contains_annotations = parse_tables(os.path.join(directory, file_name), nlp_object)
            if tables and contains_annotations:
                output_tables(F"output/json/{file_name}", tables)
            elif tables and not contains_annotations:
                update_gui_progress(qt_progress_signal, F"No annotations found for {file_name}...")
            else:
                update_gui_progress(qt_progress_signal, F"Unable to process study {file_name}. Skipping...")
            continue
        study = prepare_study(directory, file_name)

        if not study:
            response = QtFinishedResponse(False, file_name)
            qt_study_finished_signal.emit(response)
            continue

        logger.info(F"Processing PMC {study['documents'][0]['id']}")
        result = process_study(nlp_object, study, qt_progress_signal, qt_study_finished_signal)

        if not result:
            update_gui_progress(qt_progress_signal, F"Unable to process study {file_name}. Skipping...")

    if qt_study_finished_signal:
        response = QtFinishedResponse(True, "Finished processing.", 1)
        qt_study_finished_signal.emit(response)


def visualise_study(file, visualisation_type):
    import os
    study = prepare_study(directory=os.path.dirname(file), file_name=os.path.basename(file))
    nlp = load_nlp_object()
    doc = nlp.process_corpus(nlp.replace_all_abbreviations(study.get_fulltext()))
    if visualisation_type.lower() == "ents":
        nlp.display_ents(doc, theme=theme())
    elif visualisation_type.lower() == "sents":
        nlp.display_structure([sent for sent in doc.sents], theme=theme())


def main():
    import argparse
    import os
    import sys
    # Configure command line arguments
    parser = argparse.ArgumentParser(description='GWAS Information Extraction')
    parser.add_argument('-c', '--cores', type=int, default=1, help='Number of CPU cores to utilise. 0 = Max '
                                                                   'available. Default = 1.')
    parser.add_argument('-d', '--docs', type=str, help='Directory containing the study JSON files.')
    parser.add_argument('-u', '--update_ont', action='store_true', help='Update ontology cache files from the source '
                                                                        'ontology files.')
    parser.add_argument('-v', '--visualise', type=str, help='Start displacy visualisation server for entities or '
                                                            'dependencies by specifying ents or sents respectively.')
    parser.add_argument('-g', '--interface', action='store_true', help='Launch using the graphical user interface.')
    parser.add_argument('-x', '--xml', action='store_true', help='Output results in BioC XML format rather than JSON.')

    # Parse input arguments
    args = parser.parse_args()
    # cores = args.cores  # not yet implemented
    docs = args.docs
    visualise = args.visualise
    using_gui = args.interface
    update_ont = args.update_ont
    global output_xml
    output_xml = args.xml

    # Setup folder for log files.
    if not os.path.isdir("logs"):
        try:
            os.makedirs("logs")
        except IOError as e:
            sys.exit(F"Unable to create logs folder: {e}")

    # Setup folder for output files.
    if not os.path.isdir("output"):
        try:
            os.makedirs("output")
        except IOError as e:
            sys.exit(F"Unable to create output folder: {e}")

    # Configure logger
    import time
    current_time = time.strftime("%Y%m%d-%H%M%S")
    logging.basicConfig(filename=F"logs/{current_time}.log", format="%(asctime)s %(levelname)-8s %(message)s",
                        level=logging.INFO)

    # Update ontology cache files if requested.
    if update_ont:
        import Ontology
        Ontology.update_ontology_cache()

    if not using_gui:
        global lexicon
        lexicon = __prepare_ontology_data()
    # Begin running either the GUI or processing studies immediately.
    if using_gui:
        os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
        from GUI import MainForm
        global gui
        gui = MainForm()
        gui.open()
    else:
        if visualise and docs:
            visualise_study(docs, visualise)
        elif docs:
            process_studies(docs)
        else:
            print("No valid parameters given.")


if __name__ == '__main__':
    main()
