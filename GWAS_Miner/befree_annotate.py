import html
import json
import re
from datetime import datetime
from difflib import SequenceMatcher
from os import listdir

from os.path import isfile, join

from GWAS_Miner import Ontology, Experimental, BioC, OutputConverter
from GWAS_Miner.BioC import BioCLocation, BioCAnnotation


def get_gc_content():
    gc_data = {}
    pmid_pmcid_dict = {}
    headers_skipped = False
    bioc_pmcids = [x.replace(".json", "").replace("_abbreviations", "") for x in listdir("BioC_Studies") if
                   isfile(join("BioC_Studies", x))]
    # retrieve matching data.
    with open("GC_content.tsv", "r", encoding="utf-8") as f_in:
        for line in f_in.readlines():
            if not headers_skipped:
                headers_skipped = not headers_skipped
                continue
            line = line.split("\t")
            if line[0] not in bioc_pmcids:
                continue
            if line[0] not in gc_data.keys():
                gc_data[line[0]] = []
            gc_data[line[0]].append([line[5], line[6], line[8]])
            if line[1] not in pmid_pmcid_dict.keys():
                pmid_pmcid_dict[line[1]] = line[0]
    return gc_data, pmid_pmcid_dict


def get_befree_data(pmid_pmcid_dict):
    befree_data = {}
    headers_skipped = False
    # retrieve befree data
    with open("BeFree_Data.tsv", "r", encoding="utf-8") as f_in:
        for line in f_in.readlines():
            if not headers_skipped:
                headers_skipped = not headers_skipped
                continue
            line = line.split("\t")
            if line[0] not in pmid_pmcid_dict.keys():
                continue
            if line[0] not in befree_data.keys():
                befree_data[line[0]] = []
            befree_data[line[0]].append({"sentence_number": line[3], "variantid": line[4], "variant_offset": line[6],
                                         "diseaseid": line[7], "disease_text": line[8], "disease_offset": line[9],
                                         "sentence": html.unescape(line[10].lstrip('"').rstrip('"')),
                                         "meshid": line[11],
                                         "mapping_source": line[13]})
    return befree_data


def get_bioc_annotations(annotations, used_annots, offset, t, m, p, g, r):
    current_datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    for annot in annotations:
        loc = BioCLocation(offset=annot["offset"] + offset, length=annot["length"])
        if annot["text"] in [x.text for x in used_annots]:
            for old_annot in used_annots:
                if old_annot.text == annot["text"] and loc not in old_annot.locations:
                    old_annot.locations.append(loc)
        if "RSID" not in annot["entity_type"] and "PVAL" not in annot["entity_type"]:
            genomic_trait = BioCAnnotation(id=F"T{t}", infons={"type": "trait", "identifier": annot["id"],
                                                               "annotator": "BeFree@example.com",
                                                               "updated_at": current_datetime},
                                           locations=[loc], text=annot["text"])
            used_annots.append(genomic_trait)
            t += 1
        elif "RSID" in annot["entity_type"]:
            marker_identifier = BioCAnnotation(id=F"M{m}",
                                               infons={"type": "genomic_marker", "identifier": annot["id"],
                                                       "annotator": "BeFree@example.com",
                                                       "updated_at": current_datetime},
                                               locations=[loc], text=annot["text"])
            used_annots.append(marker_identifier)
            m += 1
        elif "PVAL" in annot["entity_type"]:
            p_value = BioCAnnotation(id=F"P{p}", infons={"type": "significance", "identifier": annot["id"],
                                                         "annotator": "BeFree@example.com",
                                                         "updated_at": current_datetime},
                                     locations=[loc], text=annot["text"])
            used_annots.append(p_value)
            p += 1
    return used_annots, t, m, p, g, r


def get_closest_index(text, search, target):
    closest_index = None
    closest_difference = None
    for match in re.finditer(search, text, flags=re.MULTILINE):
        if closest_index is None:
            closest_index = match.start()
            closest_difference = abs(target - closest_index)
        elif abs(target - match.start()) < closest_difference:
            closest_index = match.start()
            closest_difference = abs(target - closest_index)
    return closest_index


def get_befree_annotations():

    gc_data, pmid_pmcid_dict = get_gc_content()
    befree_data = get_befree_data(pmid_pmcid_dict)

    for pmid in befree_data.keys():
        pmc_id = pmid_pmcid_dict[pmid]
        study = Experimental.load_bioc_study("BioC_Studies", F"{pmc_id}.json")
        study_befree_data = befree_data[pmid]
        used_annots = []
        relations = []
        t, m, p, g, r = 0, 0, 0, 0, 0
        if not study:
            print(F"{pmc_id} study was not found.")
            continue
        for passage in study['documents'][0]['passages']:
            if passage['infons']['section_type'] == "ABSTRACT":
                text = passage["text"]
                for entry in study_befree_data:
                    sent = sorted([(SequenceMatcher(None, x, entry['sentence']).ratio(), x) for x in text.split(". ")],
                                  key=lambda y: y[0], reverse=True)[0]
                    sent_score = sent[0]
                    sent = sent[1]
                    if not sent or sent_score < 0.6:
                        continue
                    # weird issues with P(combined)
                    annotations = []
                    sentence_offset = text.index(sent)
                    annotations.append({
                        "entity_type": "DISEASE",
                        "id": entry["meshid"],
                        "text": entry["disease_text"],
                        "offset": get_closest_index(sent, entry["disease_text"], int(
                            entry["disease_offset"][:entry["disease_offset"].index("#")])) + sentence_offset,
                        "length": len(entry["disease_text"])
                    })
                    annotations.append({
                        "entity_type": "RSID",
                        "id": entry["variantid"],
                        "text": entry["variantid"],
                        "offset": get_closest_index(sent, entry["variantid"], int(
                            entry["variant_offset"][:entry["variant_offset"].index("#")])) + sentence_offset,
                        "length": len(entry["variantid"])
                    })
                    disease_node_id = t
                    marker_node_id = m
                    used_annots, t, m, p, g, r = get_bioc_annotations(annotations, used_annots, passage["offset"], t, m, p,
                                                                      g, r)
                    phenotype_node = BioC.BioCNode(refid=F"T{disease_node_id}", role="")
                    marker_node = BioC.BioCNode(refid=F"M{marker_node_id}", role="")
                    bioc_relation = BioC.BioCRelation(id=F"R{r}",
                                                      infons={"type": "BeFree_Association",
                                                              "annotator": "tr142@le.ac.uk",
                                                              "updated_at": datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")},
                                                      nodes=[phenotype_node, marker_node])
                    r += 1
                    relations.append(bioc_relation)
                for annot in used_annots:
                    passage['annotations'].append(annot)
        study['documents'][0]['relations'] = relations
        OutputConverter.output_xml(json.dumps(study, default=BioC.ComplexHandler),
                                   F"output/PMC{study['documents'][0]['id']}_result.xml")
