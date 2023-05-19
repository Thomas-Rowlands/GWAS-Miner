import json
import logging

from DataStructures import Study

logger = logging.getLogger("GWAS Miner")


def validate_json_maintext(json_data):
    title = list(json_data.keys())[0]
    if not json_data[title]:
        return False
    else:
        return True


def load_bioc_study(directory, file_name):
    bioc_study = None
    try:
        with open(F"{directory}/{file_name}", "r", encoding="utf-8") as fin:
            bioc_study = json.loads(fin.read())
    except IOError:
        print(F"Unable to locate/open file: {file_name}")
    return bioc_study


def load_study(directory, file_name):
    json_table_data = None
    pmc_id = None
    try:
        with open(F"{directory}/{file_name}", 'r', encoding="utf-8") as file:
            json_text_data = json.load(file)
    except IOError as io:
        logger.error(F"IO Error: {io}")
        return None
    except Exception as e:
        logger.error(
            F"An error occurred attempting to read publication file {file_name}:\n {e}")
        return None
    if json_text_data:
        pmc_id = file_name[3:-15]
        try:
            with open(F"{directory}/{file_name.replace('._maintext', '')}", 'r', encoding="utf-8") as file:
                json_table_data = json.load(file)
        except FileNotFoundError:
            logger.error(F"Tables file was not found for PMCID: {pmc_id}, continuing without tables.")
            if not validate_json_maintext(json_text_data):
                return None
            else:
                study = Study(json_text_data)
                study.pmid = pmc_id
                return study
        except IOError as io:
            logger.error(F"IO Error: {io}")
            return None
        except Exception as e:
            logger.error(F"An error occurred attempting to read publication table file PMC{pmc_id}_tables.json:\n {e}")
            return None
    if not validate_json_maintext(json_text_data) or not json_table_data:
        logger.error(F"Study text or table text missing: {file_name}")
        return None
    study = Study(json_text_data)
    study.pmid = pmc_id

    return study
