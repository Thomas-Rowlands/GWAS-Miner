import json

import rdflib
import owlready2
import requests

import config


class OntologyMap:
    mesh_namespace = rdflib.namespace.Namespace("http://id.nlm.nih.gov/mesh/vocab#")
    hp_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS}
    mesh_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS,
                       "meshv": mesh_namespace}

    @staticmethod
    def get_mesh_file():
        result = []
        with open("mesh_descriptors", "r") as in_file:
            line = in_file.readline()
            while line:
                id = line[:line.index(" <|> ")].lower()
                term = line[(line.index(" <|> ") + len(" <|> ")):].lower().replace("\n", "")
                result.append((id, term))
                line = in_file.readline()
        results = [[id, term.lower()] for (id, term) in result]
        return results

    @staticmethod
    def get_mesh():
        g = rdflib.Graph()
        g.parse(config.mesh_file, format="nt")
        mesh_ontology_terms = g.query(config.mesh_statement, initNs=OntologyMap.mesh_namespaces)
        g.close()
        return mesh_ontology_terms

    @staticmethod
    def get_mesh_descriptors():
        url = "http://id.nlm.nih.gov/mesh/sparql?format=JSON&inference=true&query="
        url_offset = "http://id.nlm.nih.gov/mesh/sparql?format=JSON&inference=true&offset="
        query = requests.utils.quote(config.mesh_statement)
        response = requests.get(url + query)
        content = json.loads(response.content)
        offset = 1000
        result_count = len(content["results"]["bindings"])
        while result_count > 0:
            with open("mesh_descriptors", "a") as out_file:
                for result in content["results"]["bindings"]:
                    out_file.write(str(result["d"]["value"]) + " <|> " + result["dName"]["value"] + "\n")
            if result_count < 1000:
                break
            response = requests.get(url_offset + str(offset) + "&query=" + query)
            if not response.content:
                break
            content = json.loads(response.content)
            result_count = len(content["results"]["bindings"])
            offset += 1000

    @staticmethod
    def get_hpo():
        ont = owlready2.get_ontology(config.hpo_file)
        ont.load()
        g = owlready2.default_world.as_rdflib_graph()
        hpo_ontology_terms = g.query(config.hpo_statement, initNs=OntologyMap.hp_namespaces)
        g.close()
        results = [[id, term.toPython().lower()] for (id, term) in hpo_ontology_terms]
        return results

    @staticmethod
    def get_hpo_synonyms():
        ont = owlready2.get_ontology(config.hpo_file)
        ont.load()
        g = owlready2.default_world.as_rdflib_graph()
        hpo_ontology_syns = g.query(config.hpo_syns_statement, initNs=OntologyMap.hp_namespaces)
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
