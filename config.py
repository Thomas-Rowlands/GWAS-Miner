# connection details
server = "127.0.0.1"
database = "ontology_partial"
study_db = "study_0219_partial"
user = "root"
password = "Maggie7803GB!"
port = 3306

# files
hpo_file = "ontology_data/hp.owl"
mesh_file = "mesh.nt"
efo_file = "ontology_data/efo.owl"

# Ontology SPARQL Queries

efo_terms_statement = """
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

hpo_statement = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>

SELECT DISTINCT ?id ?name
WHERE { ?temp rdf:type owl:Class .
        ?temp rdfs:label "Phenotypic abnormality" . 
        ?id rdfs:subClassOf* ?temp . 
        ?id rdfs:label ?name
        }
"""

hpo_syns_statement = """PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX hp2: <http://purl.obolibrary.org/obo/hp.owl#>
PREFIX obo: <http://www.geneontology.org/formats/oboInOwl#>

SELECT DISTINCT ?id ?synonym
WHERE { ?temp rdf:type owl:Class .
        ?temp rdfs:label "Phenotypic abnormality" . 
        ?id rdfs:subClassOf* ?temp . 
        ?id <http://www.geneontology.org/formats/oboInOwl#hasExactSynonym> ?synonym .
        }
"""

study_query = """SELECT identifier, Title, StudyAbstract
	FROM Study
    WHERE StudyAbstract NOT IN ("", "Not supplied");"""