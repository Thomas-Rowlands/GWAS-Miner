class OntologyMap:
    @staticmethod
    def get_mesh(cursor):
        test_mesh = {"meshID": [], "termID": [], "termName": []}
        for (meshID, termID, termName) in cursor:
            test_mesh["meshID"].append(meshID)
            test_mesh["termName"].append(termName)
            test_mesh["termID"].append(termID)
        cursor.close()
        return test_mesh

    @staticmethod
    def get_hpo(cursor):
        test_hpo = {"hpoID": [], "name": []}
        for (hpoID, name) in cursor:
            test_hpo["hpoID"].append(hpoID)
            test_hpo["name"].append(name)
        cursor.close()
        return test_hpo

    @staticmethod
    def get_hpo_synonyms(cursor):
        test_hpo_syns = {"hpoID": [], "synText": []}
        for (hpoID, synText) in cursor:
            test_hpo_syns["hpoID"].append(hpoID)
            test_hpo_syns["synText"].append(synText)
        cursor.close()
        return test_hpo_syns

    @staticmethod
    def get_hpo_2_mesh(cursor):
        hpo_to_mesh = {"hpoID": [], "meshID": []}
        for (hpoID, meshID) in cursor:
            hpo_to_mesh["hpoID"].append(hpoID)
            hpo_to_mesh["meshID"].append(meshID)
        cursor.close()
        return hpo_to_mesh
