# Ontology files
hpo_file = "ontology_data/hp.owl"
mesh_file = "ontology_data/desc2020.xml"
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

hpo_syns_statement = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
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

# NLP Variables
regex_entity_patterns = {
    "PVAL": [
        r"((\(?\b[pP][  =<-]{1,}(val{1,}[ue]{0,})?[  <≥=×xX-]{0,}[  \(]?\d+[\.]?[\d]{0,}[-−_^*()  \d×xX]{0,}))",
        r"(\d?\..?\d[ ]?[*×xX]{1}[ ]?\d{1,}[ (]?[-−]\d{1,}[ )]?)",
        r"((\(?\b[pP][  =<-]{1,}(val{1,}[ue]{0,})?[  <≥=×xX-]{0,}[  \(]?\d+[\.]?[\d]{0,}[-^*()  \d×xX]{0,"
        r"})|(\d?\.?\d[  ]?[*×xX]{1}[  ]?\d{1,}[  ]?-\d{1,}))",
        "([pP][- ]{1,2}[val]{0,3}[ue]{0,}[^0-9]{0,9})([0-9]+)([0-9.e-]+)",
        r"([pP][VAL]{0,3}[ =]+[xX× _\-−]+[0-9]+)"
    ],
    "PTYPE": r"(\(?GEE\)?)|(\(?FBAT\)?)",
    "Table Ref": r"(table[- ]{0,}\d{1,})"

}
