import codecs
import json
import logging
import pickle

import config
import owlready2
import rdflib
from rtgo import ReadyThready
from py2neo import Graph
from DataStructures import Lexicon, LexiconEntry, MasterLexicon, MeshDescriptor, MeshTerm, MeshConcept

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
    master = get_graph_ontology_data()
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
    tagging_data = {"HPO": [], "MESH": [], "EFO": []}
    for (id, label) in hpo_data:
        tagging_data["HPO"].append(label)
    for (id, label) in hpo_syns:
        tagging_data["HPO"].append(label)
    tagging_data["HPO"] = list(set(tagging_data["HPO"]))
    for (id, label) in mesh_data:
        tagging_data["MESH"].append(label)
    tagging_data["MESH"] = list(set(tagging_data["MESH"]))
    # for (id, label) in efo_terms:
    #     tagging_data["EFO"].append(label)
    # for (id, label) in efo_syns:
    #     tagging_data["EFO"].append(label)
    return tagging_data


def update_ontology_cache(qt_progress_signal=None, qt_finished_signal=None):
    """[Updates the ontology cache files with data from the source ontology files.]
    """
    logger.info("Updating ontology cache.")
    # set_master_lexicon()
    get_graph_ontology_data()
    logger.info("Finished updating ontology cache.")
    if qt_finished_signal:
        from GUI import QtFinishedResponse
        response = QtFinishedResponse(status=True, text="Updated ontology data.")
        qt_finished_signal.emit(response)


def filter_mesh_lexicon(lexi):
    included = ["A", "C", "G", "F", "D", "E01", "N06.850"]
    excluded = ["G05", "F02.463.425.069", "F04.754.720.346", "F01.829.263", "I01.880.853.150"]
    for entry in lexi.get_entries():
        in_included = False
        in_excluded = False
        for branch in included:
            for tree_id in entry.tree_id():
                if tree_id and tree_id.startswith(branch):
                    in_included = True
                    break
            if in_included:
                break
        if not in_included:
            for branch in excluded:
                for tree_id in entry.tree_id():
                    if tree_id and tree_id.startswith(branch):
                        in_excluded = True
                        break
                if in_excluded:
                    break
        if not in_included:
            lexi.remove_entry(entry)
    return lexi


def get_graph_ontology_data():
    ontologies = ["MESH", "HPO"]
    master_lexi = MasterLexicon()
    for ont in ontologies:
        lexi = __retrieve_ont_lexicon(ont)
        if lexi.name == "MESH":
            lexi = filter_mesh_lexicon(lexi)
            output_lexicon_terms(lexi)
        master_lexi.add_lexicon(lexi)
    return master_lexi


def output_lexicon_terms(lexi):
    with open("lexi_dump.tsv", "w", encoding="utf-8") as fout:
        for entry in lexi.get_entries():
            fout.write(", ".join(entry.tree_id()) + "\t" + entry.identifier + "\t" + entry.name())


def __retrieve_ont_lexicon(ontology_name):
    try:
        graph = Graph(scheme="bolt", host="localhost", password="12345")
        mesh_query = F"""
        MATCH (n:{ontology_name})-[:HAS_SYNONYM*0..]->(m)
        WHERE n:Term
        WITH n, COLLECT(m.FSN) AS syns
        RETURN n.id, n.FSN, n.treeid, syns
        """
        cursor = graph.run(mesh_query)
        lexi = Lexicon(name=ontology_name)
        while cursor.forward():
            old_entry = lexi.get_entry_by_term(str(cursor.current[1]))
            if old_entry:
                if ontology_name == "MESH":
                    old_entry.add_tree_id(str(cursor.current[2]))
                continue
            entry = LexiconEntry(str(cursor.current[0]), name=str(cursor.current[1]), tree_id=str(cursor.current[2]))
            for syn in cursor.current[3]:
                entry.add_synonym(id=entry.identifier, name=syn)
            lexi.add_entry(entry)
        return lexi
    except ConnectionRefusedError as cre:
        raise cre


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

    @staticmethod
    def get_lexicon():
        new_lexicon = Lexicon(name="MESH")
        data = Mesh.extract_mesh_data()
        for descriptor in data:
            entry = LexiconEntry(identifer=descriptor.ui, name=descriptor.name, mesh_descriptor=descriptor)
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
        descriptors = []
        try:
            import xml.etree.cElementTree as ET
            context = ET.iterparse("../ontology_data/desc2020.xml", events=("start", "end"))
            current_descriptor = None
            current_concept = None
            current_term = None
            path = []
            for event, elem in context:
                if event == "start":
                    path.append(elem.tag)
                    if elem.tag == "Concept":
                        current_concept = MeshConcept()
                        current_concept.is_preferred = elem.attrib["PreferredConceptYN"]
                    elif elem.tag == "DescriptorRecord":
                        current_descriptor = MeshDescriptor()
                    elif elem.tag == "Term":
                        current_term = MeshTerm()
                        current_term.is_preferred = elem.attrib["ConceptPreferredTermYN"]
                elif event == "end":
                    if elem.tag == "DescriptorUI":
                        current_descriptor.ui = elem.text
                    elif elem.tag == "String":
                        if "DescriptorName" in path:
                            current_descriptor.name = elem.text
                        elif "ConceptName" in path:
                            current_concept.name = elem.text
                        elif "Term" in path:
                            current_term.name = elem.text
                    elif elem.tag == "TreeNumber":
                        if "TreeNumberList" in path:
                            current_descriptor.tree_nums.append(elem.text)
                    elif elem.tag == "ConceptUI":
                        if "Concept" in path:
                            current_concept.ui = elem.text
                    elif elem.tag == "TermUI":
                        current_term.ui = elem.text
                    elif elem.tag == "Term":
                        current_concept.terms.append(current_term)
                    elif elem.tag == "Concept":
                        current_descriptor.concepts.append(current_concept)
                    elif elem.tag == "DescriptorRecord":
                        descriptors.append(current_descriptor)
                    path.pop()
        except IOError as io:
            logger.error(F"IO error parsing MeSH descriptors XML file: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            logger.error(F"An unexpected error occurred whilst parsing MeSH descriptors XML file: {ex}")
        descriptors = [x for x in descriptors if False not in [Mesh.validate_branch(y) for y in x.tree_nums]]
        return descriptors

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
    def get_terms():
        results = None
        try:
            ont = owlready2.get_ontology(config.hpo_file)
            ont.load()
            g = owlready2.default_world.as_rdflib_graph()
            hpo_ontology_terms = g.query(config.hpo_statement, initNs=HPO.hp_namespaces)
            g.close()
            results = [[id[id.rfind("/") + 1:], term.toPython().lower()] for (id, term) in hpo_ontology_terms]
        except IOError as io:
            logger.error(F"IO error reading new HPO terms: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            logger.error(F"An unexpected error occurred whilst storing new HPO terms: {ex}")
        return results

    @staticmethod
    def get_lexicon():
        terms = None
        syns = None
        new_lexicon = Lexicon(name="HPO")
        try:
            terms = HPO.get_terms()
            syns = HPO.get_hpo_synonyms()
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
    def get_hpo_synonyms():
        ont = owlready2.get_ontology(config.hpo_file)
        ont.load()
        g = owlready2.default_world.as_rdflib_graph()
        hpo_ontology_syns = g.query(config.hpo_syns_statement, initNs=HPO.hp_namespaces)
        g.close()
        results = [[id, synonym.toPython().lower()] for (id, synonym) in hpo_ontology_syns]
        return results

    @staticmethod
    def get_hpo_2_mesh(cursor):
        hpo_to_mesh = {"hpoID": [], "meshID": []}
        for (hpoID, meshID) in cursor:
            hpo_to_mesh["hpoID"].append(hpoID)
            hpo_to_mesh["meshID"].append(meshID)
        cursor.close()
        return hpo_to_mesh
