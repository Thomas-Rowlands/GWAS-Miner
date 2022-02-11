import logging
import pickle

from py2neo import Graph
from DataStructures import Lexicon, LexiconEntry, MasterLexicon

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


def update_ontology_cache(qt_progress_signal=None, qt_finished_signal=None):
    """[Updates the ontology cache files with data from the source ontology files.]
    """
    logger.info("Updating ontology cache.")
    set_master_lexicon()
    # get_graph_ontology_data()
    logger.info("Finished updating ontology cache.")
    if qt_finished_signal:
        from GUI import QtFinishedResponse
        response = QtFinishedResponse(status=True, text="Updated ontology data.")
        qt_finished_signal.emit(response)


def get_graph_ontology_data():
    ontologies = ["MESH", "HPO"]
    master_lexi = MasterLexicon()
    for ont in ontologies:
        lexi = __retrieve_ont_lexicon(ont)
        master_lexi.add_lexicon(lexi)
    return master_lexi


def output_lexicon_terms(lexi):
    with open("included_only.tsv", "w", encoding="utf-8") as fout:
        for entry in lexi.get_entries():
            fout.write(", ".join(entry.tree_id()) + "\t" + entry.identifier + "\t" + entry.name() + "\n")


def __retrieve_ont_lexicon(ontology_name):
    try:
        graph = Graph(scheme="bolt", host="localhost", password="12345")
        ont_query = ""
        if ontology_name == "MESH":
            ont_query = """
                MATCH (n:MESH)
                WHERE n.miner_included = true
                WITH COLLECT(n.id) AS included
                
                MATCH (m:MESH)
                WHERE m.miner_included = false
                WITH COLLECT(m.id) AS excluded, included
                
                MATCH (n:MESH)-[:HAS_SYNONYM*0..]->(syn{miner_included: True})
                WHERE n.id IN included AND NOT n.id in excluded AND n.miner_included = TRUE
                WITH n, COLLECT(syn.FSN) AS syns
                RETURN DISTINCT n.id, n.FSN, n.treeid, syns
            """
        else:
            ont_query = F"""
            MATCH (n:{ontology_name})-[:HAS_SYNONYM*0..]->(m)
            WHERE n:Term
            WITH n, COLLECT(m.FSN) AS syns
            RETURN n.id, n.FSN, n.treeid, syns
            """
        cursor = graph.run(ont_query)
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
