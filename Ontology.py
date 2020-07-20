import sys

import rdflib
import owlready2
from lxml import etree
import json
import config
import csv
import codecs
import logging

logger = logging.getLogger("Phenotype Finder")

class EFO:
    efo_namespace = rdflib.namespace.Namespace("http://www.ebi.ac.uk/efo/efo.owl")
    efo_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS}

    @staticmethod
    def set_terms():
        g = rdflib.Graph()
        try:
            g.parse(config.efo_file)
            efo_ontology_terms = g.query(config.efo_terms_statement, initNs=EFO.efo_namespaces)
            with open("ontology_data/efo_terms.json", "a") as out_file:
                results = {}
                for (id, label, exactSyn) in efo_ontology_terms:
                    if id not in results:
                        results[id] = [label, exactSyn]
                    else:
                        logger.warning("EFO ID was not found...")
                json.dump(results, out_file)
        except IOError as io:
            sys.exit(F"The following IO error occurred querying EFO terms: {io.strerror}")
        except BaseException as ex:
            sys.exit(F"The following error occurred querying EFO terms: {ex}")
        g.close()

    @staticmethod
    def get_efo_from_cache():
        results, terms, syns = None, None, None
        try:
            results = json.load(codecs.open('ontology_data/efo_terms.json', 'r', 'utf-8-sig'))
            terms = [[x["a.id"], x["a.FSN"]] for x in results]
            temp_syns = [[x["a.id"], x["synonym"]] for x in results]
            syns = []
            for i in range(len(temp_syns)):
                for syn in temp_syns[i][1]:
                    syns.append([temp_syns[i][0], syn])
        except IOError as io:
            sys.exit(F"The following IO error occurred retrieving cached EFO terms: {io.strerror}")
        except BaseException as ex:
            sys.exit(F"The following IO error occurred retrieving cached EFO terms: {ex}")

        return terms, syns


class Mesh:
    mesh_namespace = rdflib.namespace.Namespace("http://id.nlm.nih.gov/mesh/vocab#")
    mesh_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS,
                       "meshv": mesh_namespace}

    @staticmethod
    def get_mesh_phenotypes():
        data = None
        with open("ontology_data/MeSH_Phenotypes.csv", "r") as input_file:
            reader = csv.reader(input_file)
            data = list(reader)[1:]
        for i in range(len(data)):
            data[i][0], data[i][1] = data[i][1], data[i][0].replace("\"", "")  # Swap order and remove quotes
        return data

    @staticmethod
    def get_descriptors():
        data = []
        with open("ontology_data/mesh.json") as file:
            terms = json.load(file)
        for id, label in terms:
            data.append([id, label.lower()])
        return data

    @staticmethod
    def set_descriptors():
        filtered_descriptors = []
        try:
            parser = etree.XMLParser(encoding='utf-8')
            tree = etree.parse("ontology_data/desc2020.xml")
            desired_branches = ["A", "C", "G", "F", "D", "E01", "N06.850"]
            branch_exceptions = ["G17", "G05", "F02.463.425.069", "F04.754.720.346", "F01.829.263",
                                 "I01.880.853.150"]

            all_descriptors = [x for x in tree.xpath("//DescriptorRecord", smart_string=False)]
            for desc in all_descriptors:
                is_valid = False
                for num in desc.xpath(".//TreeNumberList//TreeNumber//text()"):
                    for desired in desired_branches:
                        if desired in num:
                            is_valid = True
                    for undesired in branch_exceptions:
                        if undesired in num:
                            is_valid = False
                            break
                if is_valid:
                    filtered_descriptors.append(desc)
        except IOError as io:
            sys.exit(F"IO error parsing MeSH descriptors XML file: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            sys.exit(F"An unexpected error occurred whilst parsing MeSH descriptors XML file: {ex}")

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
            file = open("ontology_data/mesh.json", "w")
            file.write(output)
            file.close()
        except IOError as io:
            sys.exit(F"IO error saving MeSH descriptors JSON file: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            sys.exit(F"An unexpected error occurred whilst parsing MeSH descriptors XML file: {ex}")
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
            Mesh.set_descriptors()
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
            file = open("ontology_data/hp.json", "w")
            file.write(output)
            file.close()
        except IOError as io:
            sys.exit(F"IO error storing new HPO terms: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            sys.exit(F"An unexpected error occurred whilst storing new HPO terms: {ex}")
        return results

    @staticmethod
    def get_hpo_from_cache():
        results = None
        try:
            file = open("ontology_data/hp.json", "r")
            results = json.load(file)
        except IOError as io:
            sys.exit(F"IO error retrieving cached HPO terms: {io.errno} -> {io.strerror}")
        except BaseException as ex:
            sys.exit(F"An unexpected error occurred whilst retrieving cached HPO terms: {ex}")
        return results

    @staticmethod
    def set_hpo_synonyms():
        ont = owlready2.get_ontology(config.hpo_file)
        ont.load()
        g = owlready2.default_world.as_rdflib_graph()
        hpo_ontology_syns = g.query(config.hpo_syns_statement, initNs=HPO.hp_namespaces)
        g.close()
        results = [[id, synonym.toPython().lower()] for (id, synonym) in hpo_ontology_syns]
        output = json.dumps(results)
        file = open("ontology_data/hp_syns.json", "w")
        file.write(output)
        file.close()
        return results

    @staticmethod
    def get_syns():
        file = open("ontology_data/hp_syns.json", "r")
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
