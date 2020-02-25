import json

import rdflib
import owlready2
import requests
from lxml import etree
import json
import config
import csv
import codecs

class EFO:
    efo_namespace = rdflib.namespace.Namespace("http://www.ebi.ac.uk/efo/efo.owl")
    efo_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS}

    @staticmethod
    def get_efo_data():
        query_string = """
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            PREFIX go: <http://www.geneontology.org/formats/oboInOwl#>
            SELECT ?term ?termLabel ?exactSynLabel
            WHERE {
                ?term rdfs:subClassOf ?parent .
                ?term rdfs:label ?termLabel .
                ?term <http://www.geneontology.org/formats/oboInOwl#hasExactSynonym> ?exactSyn .
                ?exactSyn rdfs:label ?exactSynLabel .
            } 
            ORDER BY ?term
                """
        g = rdflib.Graph()
        g.parse(config.efo_file)
        mesh_ontology_terms = g.query(query_string, initNs=EFO.efo_namespaces)
        with open("ontology_data/efo_terms.json", "a") as out_file:
            results = {}
            for (id, label, exactSyn) in mesh_ontology_terms:
                if id not in results:
                    results[id] = [label, exactSyn]
                else:
                    print("Mistakes were made...")
            json.dump(results, out_file)
        g.close()

    @staticmethod
    def get_efo_from_cache():
        results = json.load(codecs.open('ontology_data/efo_terms.json', 'r', 'utf-8-sig'))
        terms = [[x["a.id"], x["a.FSN"]] for x in results]
        temp_syns = [[x["a.id"], x["synonym"]] for x in results]
        syns = []
        for i in range(len(temp_syns)):
            for syn in temp_syns[i][1]:
                syns.append([temp_syns[i][0], syn])
        return terms, syns


class Mesh:
    mesh_namespace = rdflib.namespace.Namespace("http://id.nlm.nih.gov/mesh/vocab#")
    mesh_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS,
                       "meshv": mesh_namespace}

    @staticmethod
    def get_mesh_data(use_xml=False):
        if use_xml:
            return Mesh.__get_descriptors_xml()
        else:
            terms = Mesh.__get_terms()
            descriptors = Mesh.__get_descriptors()
            data = []
            for id, term in terms:
                data.append([id, term])
            for id, descriptor in descriptors:
                data.append([id, descriptor])
        return data

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
    def __get_descriptors():
        data = []
        with open("mesh_descriptors.json") as file:
            terms = json.load(file)
        for id in terms.keys():
            for label in terms[id]:
                data.append([id, label.lower()])
        return data

    @staticmethod
    def __get_terms():
        data = []
        with open("mesh_terms.json") as file:
            terms = json.load(file)
        for id in terms.keys():
            for label in terms[id]:
                data.append([id, label.lower()])
        return data

    @staticmethod
    def __set_terms():
        try:
            query_string = """
                PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
                PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
                PREFIX owl: <http://www.w3.org/2002/07/owl#>
                PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
                PREFIX mesh: <http://id.nlm.nih.gov/mesh/>

                SELECT ?t ?tName ?tAltName
                WHERE {
                  ?t a/(rdfs:subClassOf|owl:equivalentClass)* meshv:Term .
                  ?t meshv:prefLabel ?tName .
                  ?t meshv:altLabel ?tAltName .
                } 
                ORDER BY ?t
                    """
            g = rdflib.Graph()
            g.parse(config.mesh_file, format="nt")
            mesh_ontology_terms = g.query(query_string, initNs=Mesh.mesh_namespaces)
            with open("mesh_terms.json", "a") as out_file:
                results = {}
                for (id, label, altLabel) in mesh_ontology_terms:
                    if id not in results:
                        results[id] = [label, altLabel]
                    else:
                        results[id].append(altLabel)
                json.dump(results, out_file)
            g.close()
        except IOError as e:
            print("IO error storing new mesh terms: {0} -> {1}".format(e.errno, e.strerror))
            g.close()
            return False
        except Exception as e:
            print("An unexpected error occurred whilst retrieving/storing new MeSH terms: {0}".format(str(e)))
            g.close()
            return False
        return True

    @staticmethod
    def __get_descriptors_xml():
        parser = etree.XMLParser(encoding='utf-8')
        tree = etree.parse("desc2020.xml")
        desired_branches = ["A", "C", "G", "F", "D", "E01", "N06.850"]
        branch_exceptions = ["G17", "G05", "F02.463.425.069", "F04.754.720.346", "F01.829.263",
                             "I01.880.853.150"]
        # results = [str(x) for x in tree.xpath("//DescriptorRecord/DescriptorUI/text() | "
        #                                       "//DescriptorRecord/DescriptorName/String/text()", smart_string=False)]
        # concepts = [str(x) for x in tree.xpath("//DescriptorRecord/ConceptList/Concept/ConceptUI/text() | "
        #                                        "//DescriptorRecord/ConceptList/Concept/ConceptName/String/text()",
        #                                        smart_string=False)]
        # terms = [str(x) for x in tree.xpath("//DescriptorRecord/ConceptList/Concept/TermList/Term/TermUI/text() | "
        #                                     "//DescriptorRecord/ConceptList/Concept/TermList/Term/String/text()",
        #                                     smart_string=False)]

        all_descriptors = [x for x in tree.xpath("//DescriptorRecord", smart_string=False)]
        filtered_descriptors = []
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
        results = []
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

        # sys.exit()
        # results = []  # delete me
        # results = [results[i:i + 2] for i in range(0, len(results), 2)]
        # concepts = [concepts[i:i + 2] for i in range(0, len(concepts), 2)]
        # terms = [terms[i:i + 2] for i in range(0, len(terms), 2)]
        #
        # for (id, concept) in concepts:
        #     results.append([id, concept])
        # for (id, term) in terms:
        #     results.append([id, term])

        output = json.dumps(results)
        file = open("ontology_data/mesh.json", "w")
        file.write(output)
        file.close()
        return results

    @staticmethod
    def get_mesh_from_cache():
        file = open("ontology_data/mesh.json", "r")
        results = json.load(file)
        return results

    @staticmethod
    def __set_descriptors():
        try:
            query_string = """
                PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
                PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
                PREFIX owl: <http://www.w3.org/2002/07/owl#>
                PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
                PREFIX mesh: <http://id.nlm.nih.gov/mesh/>
                
                SELECT DISTINCT ?d ?dName
                WHERE {
                  ?a meshv:hasDescriptor ?d.
                  ?d rdfs:label ?dName .
                } 
                ORDER BY ?d 
            """
            g = rdflib.Graph()
            g.parse(config.mesh_file, format="nt")
            mesh_descriptors = g.query(query_string, initNs=Mesh.mesh_namespaces)
            with open("mesh_descriptors.json", "a") as out_file:
                results = {}
                for (id, label) in mesh_descriptors:
                    if id not in results:
                        results[id] = [label]
                    else:
                        results[id].append(label)
                json.dump(results, out_file)
            g.close()
        except IOError as e:
            print("IO error storing new mesh terms: {0} -> {1}".format(e.errno, e.strerror))
            g.close()
            return False
        except Exception as e:
            print("An unexpected error occurred whilst retrieving/storing new MeSH descriptors: {0}".format(str(e)))
            g.close()
            return False
        return True

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
            print("MeSH file updated")
            Mesh.__set_descriptors()
            Mesh.__set_terms()
            print("New MeSH data extracted.")
        else:
            print("Already using the latest MeSH file.")
        return True


class HPO:
    hp_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS}

    @staticmethod
    def get_hpo():
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
        return results

    @staticmethod
    def get_hpo_from_cache():
        file = open("ontology_data/hp.json", "r")
        results = json.load(file)
        return results

    @staticmethod
    def get_hpo_synonyms():
        ont = owlready2.get_ontology(config.hpo_file)
        ont.load()
        g = owlready2.default_world.as_rdflib_graph()
        hpo_ontology_syns = g.query(config.hpo_syns_statement, initNs=HPO.hp_namespaces)
        g.close()
        results = [[id, synonym.toPython().lower()] for (id, term, synonym) in hpo_ontology_syns]
        output = json.dumps(results)
        file = open("ontology_data/hp_syns.json", "w")
        file.write(output)
        file.close()
        return results

    @staticmethod
    def get_hpo_syns_from_cache():
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
