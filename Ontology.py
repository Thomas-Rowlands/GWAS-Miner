import json

import rdflib
import owlready2
import requests

import config
import csv


class Mesh:
    mesh_namespace = rdflib.namespace.Namespace("http://id.nlm.nih.gov/mesh/vocab#")
    mesh_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS,
                       "meshv": mesh_namespace}

    @staticmethod
    def get_mesh_data():
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
        results = [[id, term.toPython().lower()] for (id, term) in hpo_ontology_terms]
        return results

    @staticmethod
    def get_hpo_synonyms():
        ont = owlready2.get_ontology(config.hpo_file)
        ont.load()
        g = owlready2.default_world.as_rdflib_graph()
        hpo_ontology_syns = g.query(config.hpo_syns_statement, initNs=HPO.hp_namespaces)
        g.close()
        results = [[id, synonym.toPython().lower()] for (id, term, synonym) in hpo_ontology_syns]
        return results

    @staticmethod
    def get_hpo_2_mesh(cursor):
        hpo_to_mesh = {"hpoID": [], "meshID": []}
        for (hpoID, meshID) in cursor:
            hpo_to_mesh["hpoID"].append(hpoID)
            hpo_to_mesh["meshID"].append(meshID)
        cursor.close()
        return hpo_to_mesh
