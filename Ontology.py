import rdflib
import owlready2
import config


class OntologyMap:
    mesh_namespace = rdflib.namespace.Namespace("http://id.nlm.nih.gov/mesh/vocab#")
    hp_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS}
    mesh_namespaces = {"owl": rdflib.namespace.OWL, "rdf": rdflib.namespace.RDF, "rdfs": rdflib.namespace.RDFS,
                       "meshv": mesh_namespace}

    @staticmethod
    def get_mesh():
        g = rdflib.Graph()
        g.parse(config.mesh_file, format="nt")
        mesh_ontology_terms = g.query(config.mesh_statement, initNs=OntologyMap.mesh_namespaces)
        g.close()
        return mesh_ontology_terms

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
        g = rdflib.Graph()
        g.parse(config.hpo_file, format="owl")
        hpo_ontology_syns = g.query(config.hpo_syns_statement, initNs=OntologyMap.namespaces)
        g.close()
        return hpo_ontology_syns

    @staticmethod
    def get_hpo_2_mesh(cursor):
        hpo_to_mesh = {"hpoID": [], "meshID": []}
        for (hpoID, meshID) in cursor:
            hpo_to_mesh["hpoID"].append(hpoID)
            hpo_to_mesh["meshID"].append(meshID)
        cursor.close()
        return hpo_to_mesh
