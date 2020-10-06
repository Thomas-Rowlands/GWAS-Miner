import codecs
import csv
import json
import logging
import pickle
from collections import namedtuple

import config
import owlready2
import rdflib
from lxml import etree
from rtgo import ReadyThready

from GWAS_Miner.DataStructures import Lexicon, LexiconEntry, MasterLexicon

logger = logging.getLogger("GWAS Miner")


def validate_data(ont_data):
    if ont_data:
        for i in ont_data:
            if not i:
                return False
    else:
        return False
    if len(ont_data[2]) != 2:
        return False
    return True


def set_master_lexicon():
    mesh_lexicon = Mesh.get_lexicon()
    hpo_lexicon = HPO.get_lexicon()
    master = MasterLexicon()
    master.add_lexicon(mesh_lexicon)
    master.add_lexicon(hpo_lexicon)
    master.set_priority_order({mesh_lexicon.name: 1, hpo_lexicon.name: 2})
    try:
        with open("../ontology_data/lexicon.lexi", "wb") as file:
            pickle.dump(master, file)
    except IOError as io:
        logger.error(F"Unable to create lexicon cache: {io}")
    except Exception as ex:
        logger.error(F"An unexpected error occurred creating lexicon cache: {ex}")
    return master


def get_master_lexicon():
    master = None
    try:
        with open("../ontology_data/lexicon.lexi", "rb") as file:
            master = pickle.load(file)
    except FileNotFoundError:
        logger.info(F"Cache missing, creating new cache...")
        return set_master_lexicon()
    except IOError as io:
        logger.error(F"Unable to read lexicon cache: {io}")
    except Exception as ex:
        logger.error(F"An unexpected error occurred reading lexicon cache: {ex}")
    return master


def __get_tagging_data():
    """[Assigns ontology data from cache files ready for processing]
    """
    # Update ontology file & extracted data.
    # Mesh.update_mesh_file() #  Could be dangerous to use if structural changes are made.

    # Retrieve MeSH & HPO terms/descriptors/concepts

    # Retrieve HPO synonyms
    logger.info("Starting ontology data loading.")
    ont_data = ReadyThready.go_cluster([Mesh.get_lexicon, HPO.get_lexicon,
                                        EFO.get_efo_from_cache, HPO.get_syns])
    if not validate_data(ont_data):
        update_ontology_cache()
        ont_data = ReadyThready.go_cluster([Mesh.get_lexicon, HPO.get_lexicon,
                                            EFO.get_efo_from_cache, HPO.get_syns])
        if not validate_data(ont_data):
            logger.error("Critical failure to extract ontology data.")
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
    tagging_data = {"HPO": [], "MeSH": [], "EFO": []}
    for (id, label) in hpo_data:
        tagging_data["HPO"].append(label)
    for (id, label) in hpo_syns:
        tagging_data["HPO"].append(label)
    tagging_data["HPO"] = list(set(tagging_data["HPO"]))
    for (id, label) in mesh_data:
        tagging_data["MeSH"].append(label)
    tagging_data["MeSH"] = list(set(tagging_data["MeSH"]))
    # for (id, label) in efo_terms:
    #     tagging_data["EFO"].append(label)
    # for (id, label) in efo_syns:
    #     tagging_data["EFO"].append(label)
    return tagging_data


def update_ontology_cache(qt_progress_signal=None, qt_finished_signal=None):
    """[Updates the ontology cache files with data from the source ontology files.]
    """
    logger.info("Updating ontology cache files.")
    funcs = [Mesh.extract_mesh_data,
             HPO.set_terms,
             HPO.set_hpo_synonyms,
             EFO.set_terms]
    ReadyThready.go_cluster(funcs)
    logger.info("Finished updating ontology cache files.")
    if qt_finished_signal:
        from GUI import QtFinishedResponse
        response = QtFinishedResponse(status=True, text="Updated ontology data.")
        qt_finished_signal.emit(response)


class EFO:
    efo_namespace = rdflib.namespace.Namespace("http://www.ebi.ac.uk/efo/efo.owl")
    efo_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS}

    @staticmethod
    def set_terms():
        g = rdflib.Graph()
        try:
            g.parse(config.efo_file)
            efo_ontology_terms = g.query(config.efo_terms_statement, initNs=EFO.efo_namespaces)
            with open("../ontology_data/efo_terms.json", "w") as out_file:
                results = {}
                for (id, label, exactSyn) in efo_ontology_terms:
                    if id not in results:
                        results[id] = [label, exactSyn]
                    else:
                        logger.warning("EFO ID was not found...")
                json.dump(results, out_file)
        except IOError as io:
            logger.error(F"The following IO error occurred querying EFO terms: {io.strerror}")
        except BaseException as ex:
            logger.error(F"The following error occurred querying EFO terms: {ex}")
        g.close()

    @staticmethod
    def get_efo_from_cache():
        results, terms, syns = None, None, None
        try:
            results = json.load(codecs.open('../ontology_data/efo_terms.json', 'r', 'utf-8-sig'))
            terms = [[x["a.id"], x["a.FSN"]] for x in results]
            temp_syns = [[x["a.id"], x["synonym"]] for x in results]
            syns = []
            for i in range(len(temp_syns)):
                for syn in temp_syns[i][1]:
                    syns.append([temp_syns[i][0], syn])
        except IOError as io:
            logger.error(F"The following IO error occurred retrieving cached EFO terms: {io.strerror}")
        except BaseException as ex:
            logger.error(F"The following IO error occurred retrieving cached EFO terms: {ex}")

        return terms, syns


class Mesh:
    mesh_namespace = rdflib.namespace.Namespace("http://id.nlm.nih.gov/mesh/vocab#")
    mesh_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS,
                       "meshv": mesh_namespace}

    @staticmethod
    def get_mesh_phenotypes():
        data = None
        with open("../ontology_data/MeSH_Phenotypes.csv", "r") as input_file:
            reader = csv.reader(input_file)
            data = list(reader)[1:]
        for i in range(len(data)):
            data[i][0], data[i][1] = data[i][1], data[i][0].replace("\"", "")  # Swap order and remove quotes
        return data

    @staticmethod
    def get_lexicon():
        new_lexicon = Lexicon(name="MESH")
        data = []
        with open("../ontology_data/mesh.json") as file:
            terms = json.load(file)
        for id, label in terms:
            if new_lexicon.identifier_used(id):
                new_lexicon.assign_synonym(id, label)
                continue
            entry = LexiconEntry(identifer=id, name=label)
            new_lexicon.add_entry(entry)
        return new_lexicon

    @staticmethod
    def validate_branch(tree_num):
        desired_branches = ["A", "C", "G", "F", "D", "E01", "N06.850"]
        branch_exceptions = ["G17", "G05", "F02.463.425.069", "F04.754.720.346", "F01.829.263",
                             "I01.880.853.150"]
        result = False
        for desired in desired_branches:
            if desired in tree_num:
                result = True
                for undesired in branch_exceptions:
                    if undesired in tree_num:
                        return False
        return result

    @staticmethod
    def extract_mesh_data():
        filtered_descriptors = []
        try:
            parser = etree.XMLParser(encoding='utf-8')
            tree = etree.parse("../ontology_data/desc2020.xml")
            unwanted_descriptors = []
            all_descriptors = [x for x in tree.xpath("//DescriptorRecord", smart_string=False)]
            for desc in all_descriptors:
                for tree_num in desc.xpath(".//TreeNumberList//TreeNumber//text()"):
                    result = Mesh.validate_branch(tree_num)
                    if not result:
                        unwanted_descriptors.append(tree_num)
                        break
                    else:
                        filtered_descriptors.append(desc)
            for desc in unwanted_descriptors:
                if desc in filtered_descriptors:
                    filtered_descriptors.remove(desc)

        except IOError as io:
            logger.error(F"IO error parsing MeSH descriptors XML file: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            logger.error(F"An unexpected error occurred whilst parsing MeSH descriptors XML file: {ex}")

        results = []
        try:
            for desc in filtered_descriptors:
                descriptors = [str(x) for x in desc.xpath("./DescriptorUI/text() | "
                                                          "./DescriptorName/String/text()", smart_string=False)]
                concepts = [str(x) for x in desc.xpath("./ConceptList/Concept/ConceptUI/text() | "
                                                       "./ConceptList/Concept/ConceptName/String/text()",
                                                       smart_string=False)]
                terms = [str(x) for x in desc.xpath("./ConceptList/Concept/TermList/Term/TermUI/text() | "
                                                    "./ConceptList/Concept/TermList/Term/String/text()",
                                                    smart_string=False)]
                descriptors = [descriptors[i:i + 2] for i in range(0, len(descriptors), 2)]
                concepts = [concepts[i:i + 2] for i in range(0, len(concepts), 2)]
                terms = [terms[i:i + 2] for i in range(0, len(terms), 2)]

                for (id, descriptor) in descriptors:
                    results.append([id, descriptor])
                for (id, concept) in concepts:
                    results.append([id, concept])
                for (id, term) in terms:
                    results.append([id, term])

            output = json.dumps(results)
            file = open("../ontology_data/mesh.json", "w")
            file.write(output)
            file.close()
        except IOError as io:
            logger.error(F"IO error saving MeSH descriptors JSON file: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            logger.error(F"An unexpected error occurred whilst parsing MeSH descriptors XML file: {ex}")
        return results

    @staticmethod
    def update_mesh_file():
        from ftplib import FTP
        import gzip
        import os
        import datetime
        from dateutil import parser
        import shutil

        is_updated = False
        is_file_missing = False
        current_file_time = datetime.datetime.now()
        if os.path.exists("mesh.nt"):
            current_file_time = datetime.datetime.fromtimestamp(os.path.getmtime("mesh.nt"))
        else:
            is_file_missing = True
        ftp_address = "ftp.nlm.nih.gov"
        ftp = FTP(ftp_address)
        ftp.login()
        ftp.cwd("online/mesh/rdf/")
        time = ftp.voidcmd("MDTM mesh.nt")[4:].strip()
        time = parser.parse(time)
        if time > current_file_time or is_file_missing:
            with open("mesh.nt.gz", "wb") as file:
                ftp.retrbinary("RETR mesh.nt.gz", file.write)
            is_updated = True
        ftp.close()
        if is_updated:
            # Unzip mesh file
            with gzip.open('mesh.nt.gz', 'rb') as input_file:
                with open('mesh.nt', 'wb') as output_file:
                    shutil.copyfileobj(input_file, output_file)
            logger.info("MeSH file updated")
            Mesh.extract_mesh_data()
            logger.info("New MeSH data extracted.")
        else:
            logger.info("Already using the latest MeSH file.")
        return True


class HPO:
    hp_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS}

    @staticmethod
    def set_terms():
        results = None
        try:
            ont = owlready2.get_ontology(config.hpo_file)
            ont.load()
            g = owlready2.default_world.as_rdflib_graph()
            hpo_ontology_terms = g.query(config.hpo_statement, initNs=HPO.hp_namespaces)
            g.close()
            results = [[id[id.rfind("/") + 1:], term.toPython().lower()] for (id, term) in hpo_ontology_terms]
            output = json.dumps(results)
            file = open("../ontology_data/hp.json", "w")
            file.write(output)
            file.close()
        except IOError as io:
            logger.error(F"IO error storing new HPO terms: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            logger.error(F"An unexpected error occurred whilst storing new HPO terms: {ex}")
        return results

    @staticmethod
    def get_lexicon():
        terms = None
        syns = None
        new_lexicon = Lexicon(name="HPO")
        try:
            with open("../ontology_data/hp.json", "r") as file:
                terms = json.load(file)
            with open("../ontology_data/hp_syns.json", "r") as file:
                syns = json.load(file)
        except IOError as io:
            logger.error(F"IO error retrieving cached HPO terms: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            logger.error(F"An unexpected error occurred whilst retrieving cached HPO terms: {ex}")
        for id, label in terms:
            new_entry = LexiconEntry(identifer=id, name=label)
            new_lexicon.add_entry(new_entry)
        for id, label in syns:
            new_lexicon.assign_synonym(id, label)

        return new_lexicon

    @staticmethod
    def set_hpo_synonyms():
        ont = owlready2.get_ontology(config.hpo_file)
        ont.load()
        g = owlready2.default_world.as_rdflib_graph()
        hpo_ontology_syns = g.query(config.hpo_syns_statement, initNs=HPO.hp_namespaces)
        g.close()
        results = [[id, synonym.toPython().lower()] for (id, synonym) in hpo_ontology_syns]
        output = json.dumps(results)
        file = open("../ontology_data/hp_syns.json", "w")
        file.write(output)
        file.close()
        return results

    @staticmethod
    def get_syns():
        file = open("../ontology_data/hp_syns.json", "r")
        results = json.load(file)
        return results

    @staticmethod
    def get_hpo_2_mesh(cursor):
        hpo_to_mesh = {"hpoID": [], "meshID": []}
        for (hpoID, meshID) in cursor:
            hpo_to_mesh["hpoID"].append(hpoID)
            hpo_to_mesh["meshID"].append(meshID)
        cursor.close()
        return hpo_to_mesh
