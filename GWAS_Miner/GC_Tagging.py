import json
import re
from datetime import datetime
from os import listdir
from os.path import isfile, join

import networkx as nx
from spacy.matcher import PhraseMatcher
from spacy.tokens import Span, Token

import Ontology
from GWAS_Miner import BioC, OutputConverter, Experimental, befree_annotate, GCTableExtractor, TableExtractor
from GWAS_Miner.DataStructures import Marker, Significance, Phenotype, Association
from GWAS_Miner.PostProcessing import clean_output_annotations
from NLP import Interpreter


class GCInterpreter(Interpreter):

    def __init__(self, lexicon, ontology_only=False):
        super().__init__(lexicon, ontology_only)
        self.gc_relations = []

    def set_ontology_terms(self, term_ids):
        term_ids = list(set(term_ids))
        new_matcher = PhraseMatcher(self.nlp.vocab, attr="LOWER")
        for lexicon in self.lexicon.get_ordered_lexicons():
            if lexicon.name == "HPO":
                continue
            for term_id in term_ids:
                entry = lexicon.get_entry_by_id(term_id)
                if not entry:
                    print(term_id)
                    continue
                patterns = Interpreter.get_term_variations(entry)
                patterns = self.nlp.tokenizer.pipe(patterns)
                new_matcher.add(entry.identifier, patterns, on_match=self._Interpreter__on_match)
        self.__phrase_matcher = new_matcher

    def process_corpus(self, corpus, **kwargs):
        """[Applies tokenization, entity recognition and dependency parsing to the supplied corpus.]

        Args:
            corpus ([string]): [corpus text for information extraction]

        Returns: [SpaCy doc object]: [Parsed SpaCy doc object containing the processed input text with entities,
        tokens and dependencies.]
        """

        doc = self.nlp(corpus)
        doc.user_data["relations"] = {"PHENO_ASSOC": []}

        old_ents, doc.ents = doc.ents, []
        for pattern in self.pval_patterns:
            Interpreter._regex_match(pattern, doc, "PVAL")

        for pattern in self.rsid_patterns:
            Interpreter._regex_match(pattern, doc, "RSID")

        for pattern in self.abbrev_pattens:
            Interpreter._regex_match(pattern[1], doc, pattern[0], ignore_case=False)

        # self.__basic_matcher(doc)
        self.__phrase_matcher(doc)

        # Ensure that rule-matched entities override data model entities when needed.
        for ent in old_ents:
            doc.ents += (ent,)

        # Ensure that multi-token entities are merged for extraction and association processing.
        for ent_label in self._Interpreter__entity_labels:
            self._Interpreter__merge_spans(doc, ent_label)

        return doc

    def calculate_sdp(self, phenotype_sents, top_phenotype=None):
        """[Calculates the shortest dependency path for each phenotype/marker/p-value combination,
        returning the shortest for each one.]

        Args:
            phenotype_sents ([list]): [List of SpaCy sent objects containing phenotype entities.]
            top_phenotype (LexiconEntry): Most occurring phenotype within a study.
        Returns:
            [dict]: [Dictionary containing extracted phenotype, marker and p-values
            associated based on SDP calculation.]
        """
        results = []
        # Iterate through each sentence containing a phenotype named entity label
        for sent in phenotype_sents:
            edges = []
            token_indexes = {}
            for token in sent:
                for child in token.children:
                    # Add unique id to token strings.
                    token_text = F"{token}<id{token.idx}>"
                    child_text = F"{child}<id{child.idx}>"
                    if token.idx not in token_indexes.keys():
                        token_indexes[token.idx] = token.i
                    if child.idx not in token_indexes.keys():
                        token_indexes[child.idx] = child.i
                    edges.append(('{0}'.format(token_text),
                                  '{0}'.format(child_text)))

            graph = nx.Graph(edges)
            phenotypes = [x for x in sent.ents if x._.has_trait]
            phenotypes = [x for x in phenotypes if type(x) == Span]
            phenotypes = Interpreter._validate_node_entities(
                phenotypes, graph.nodes) if not top_phenotype else None

            markers = Interpreter._validate_node_entities(
                [x for x in sent.ents if x.label_ == 'RSID' or x.label_ == "GENE"], graph.nodes)
            pvals = Interpreter._validate_node_entities(
                [x for x in sent.ents if x.label_ == 'PVAL'], graph.nodes)

            relations = []

            if phenotypes and markers and pvals:
                for pval in pvals:
                    for s, m, p in self.gc_relations:
                        significance = None
                        marker = None
                        phenotype = None
                        if re.search(s, pval[0].text, flags=re.IGNORECASE):
                            significance = pval[1]
                            marker = next((x[1] for x in markers if x[0].text == m), None)

                            best_pheno_distance = None
                            best_pheno = None
                            for phenotype in phenotypes:
                                is_token = type(phenotype[0]) == Token
                                if is_token:
                                    temp_distance = nx.shortest_path_length(
                                        graph, target=phenotype[1], source=significance)
                                else:
                                    temp_distance = nx.shortest_path_length(
                                        graph, target=F"{phenotype[0][0].text}<id{phenotype[0][0].idx}>",
                                        source=significance)

                                if not best_pheno_distance or temp_distance < best_pheno_distance:
                                    if is_token:
                                        best_pheno = nx.shortest_path(
                                            graph, target=phenotype[1], source=significance)[-1]
                                    else:
                                        best_pheno = nx.shortest_path(
                                            graph, target=F"{phenotype[0][0].text}<id{phenotype[0][0].idx}>",
                                            source=significance)[-1]
                                    best_pheno_distance = temp_distance
                                else:
                                    continue
                            if not pval:
                                continue
                            if best_pheno_distance:
                                phenotype = best_pheno
                        if significance and marker and phenotype:
                            relations.append([significance, marker, phenotype])

                # Validate associations are triples
                pheno_assocs = [x for x in relations if len(x) == 3]
                for (pval, rsid, phenotype) in pheno_assocs:
                    temp_marker = sent.doc[
                        token_indexes[int(rsid[rsid.find("<id") + 3: rsid.find(">", rsid.find("<id") + 3)])]]
                    temp_pval = sent.doc[
                        token_indexes[int(pval[pval.find("<id") + 3: pval.find(">", pval.find("<id") + 3)])]]
                    temp_pheno = sent.doc[token_indexes[
                        int(phenotype[phenotype.find("<id") + 3: phenotype.find(">", phenotype.find("<id") + 3)])]]
                    result_marker = Marker(temp_marker)
                    result_significance = Significance(temp_pval)
                    result_pheno = Phenotype(temp_pheno)
                    results.append(
                        Association(marker=result_marker, significance=result_significance, phenotype=result_pheno))
            elif top_phenotype and markers and pvals:
                for pval in pvals:
                    for s, m, p in self.gc_relations:
                        significance = None
                        marker = None
                        phenotype = top_phenotype
                        if re.search(s, pval[0].text) and top_phenotype["ID"] == p:
                            significance = pval[1]
                            marker = next((x[1] for x in markers if x[0].text == m), None)
                            phenotype = next((x for x in sent.doc.ents if x.label_ == phenotype["ID"]), None)
                        if significance and marker and phenotype:
                            relations.append([significance, marker, phenotype])

                pheno_assocs = [x for x in relations if len(x) == 3]
                for (pval, rsid, phenotype) in pheno_assocs:
                    temp_marker = sent.doc[
                        token_indexes[int(rsid[rsid.find("<id") + 3: rsid.find(">", rsid.find("<id") + 3)])]]
                    temp_pval = sent.doc[
                        token_indexes[int(pval[pval.find("<id") + 3: pval.find(">", pval.find("<id") + 3)])]]
                    temp_pheno = sent.doc[phenotype.start:phenotype.end]
                    result_marker = Marker(temp_marker)
                    result_significance = Significance(temp_pval)
                    result_pheno = Phenotype(temp_pheno)
                    results.append(
                        Association(marker=result_marker, significance=result_significance, phenotype=result_pheno))
        return results


def get_relation(relation, passage, nlp):
    current_datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    if type(relation.phenotype.token) == Token:
        pheno_id = [x.id for x in passage['annotations'] if
                    relation.phenotype.token.idx + passage["offset"] == x.locations[0].offset][0]
    else:
        pheno_id = [x.id for x in passage['annotations'] if
                    relation.phenotype.token.start_char + passage["offset"] == x.locations[0].offset][0]
    marker_id = [x.id for x in passage['annotations'] if
                 relation.marker.token.idx + passage["offset"] == x.locations[0].offset][0]
    significance_id = [x.id for x in passage['annotations'] if
                       relation.significance.token.idx + passage["offset"] == x.locations[0].offset][0]
    phenotype_node = BioC.BioCNode(refid=pheno_id, role="")
    marker_node = BioC.BioCNode(refid=marker_id, role="")
    significance_node = BioC.BioCNode(refid=significance_id, role="")
    bioc_relation = BioC.BioCRelation(id=F"R{nlp.r}",
                                      infons={"type": "GeneticVariant_Trait_Significance", "annotator": "GWASMiner@le.ac.uk",
                                              "updated_at": current_datetime},
                                      nodes=[phenotype_node, marker_node, significance_node])
    return bioc_relation, nlp


def process_study(nlp, study):
    if not study:
        return False
    current_datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    document_relations = []
    results_present = False

    for passage in study['documents'][0]['passages']:
        if passage["infons"]["section_type"].lower() == "results":
            results_present = True
            break
    for passage in study['documents'][0]['passages']:
        # footnotes need to be excluded.
        if results_present and passage["infons"]["section_type"].lower() not in ["abstract", "results",
                                                                                 "discussion", "conclusion"]:
            continue
        passage_text = passage['text']
        doc = nlp.process_corpus(passage_text)
        top_phenotype = nlp.get_ubiquitous_phenotype(passage_text, nlp)
        doc.user_data["top_phenotype"] = top_phenotype
        annotations = nlp.get_entities(doc)
        used_annots = []
        if annotations:
            for annot in annotations:
                loc = BioC.BioCLocation(offset=annot["offset"] + passage["offset"], length=annot["length"])
                if "RSID" not in annot["entity_type"] and "PVAL" not in annot["entity_type"] \
                        and "GENE" not in annot["entity_type"]:
                    genomic_trait = BioC.BioCAnnotation(id=F"T{nlp.t}",
                                                        infons={"type": "trait", "identifier": F"MeSH:{annot['id']}",
                                                                "annotator": "GWASMiner@le.ac.uk",
                                                                "updated_at": current_datetime},
                                                        locations=[loc], text=annot["text"])
                    passage['annotations'].append(genomic_trait)
                    nlp.t += 1
                elif "RSID" in annot["entity_type"]:
                    marker_identifier = BioC.BioCAnnotation(id=F"V{nlp.v}",
                                                            infons={"type": "genetic_variant",
                                                                    "identifier": F"dbSNP:{annot['text']}",
                                                                    "annotator": "GWASMiner@le.ac.uk",
                                                                    "updated_at": current_datetime},
                                                            locations=[loc], text=annot["text"])
                    passage['annotations'].append(marker_identifier)
                    nlp.v += 1
                elif "PVAL" in annot["entity_type"]:
                    p_value = BioC.BioCAnnotation(id=F"S{nlp.s}",
                                                  infons={"type": "significance", "identifier": "PVAL",
                                                          "annotator": "GWASMiner@le.ac.uk",
                                                          "updated_at": current_datetime},
                                                  locations=[loc], text=annot["text"])
                    passage['annotations'].append(p_value)
                    nlp.s += 1
                elif "GENE" in annot["entity_type"]:
                    gene = BioC.BioCAnnotation(id=F"G{nlp.g}",
                                               infons={"type": "gene", "identifier": F"Entrez:{annot['id']}",
                                                       "annotator": "GWASMiner@le.ac.uk",
                                                       "updated_at": current_datetime},
                                               locations=[loc], text=annot["text"])
                    passage['annotations'].append(gene)
                    nlp.g += 1
                used_annots.append(annot["text"])

            relations, uncertain_relations = nlp.extract_phenotypes(doc)
            if relations:
                relations = validate_relations(nlp, relations)
            if uncertain_relations:
                uncertain_relations = validate_relations(nlp, uncertain_relations)
            for relation in relations:
                bioc_relation, nlp = get_relation(relation, passage, nlp)
                document_relations.append(bioc_relation)
                nlp.r += 1
            for relation in uncertain_relations:
                bioc_relation, nlp = get_relation(relation, passage, nlp)
                document_relations.append(bioc_relation)
                nlp.r += 1
    if document_relations:
        study['documents'][0]['relations'] = document_relations
    study, nlp = befree_annotate.get_befree_annotations(study, nlp, current_datetime)
    study = clean_output_annotations(study)
    OutputConverter.output_xml(json.dumps(study, default=BioC.ComplexHandler),
                               F"output/xml/PMC{study['documents'][0]['id']}_result.xml")
    OutputConverter.output_json(study, F"output/json/PMC{study['documents'][0]['id']}_result.json")
    return study, nlp


def generate_pval_regex_strings(input_string: str) -> str:
    if "." in input_string:
        pval_input_int = int(input_string[:input_string.find(".")])
    else:
        pval_input_int = int(input_string[:input_string.lower().find("e")])
    pval_input_range = F"{pval_input_int - 1}{pval_input_int}{pval_input_int + 1}"
    power_digits = input_string[input_string.find('-') + 1:]
    pval = [F"([{pval_input_range}][ \.0-9]*[ xX*×]*10 ?(<sup>)?[-−] ?"
            F"{power_digits[0] + '?' if power_digits[0] == '0' else power_digits[0]}"
            F"{power_digits[1:] + '(</sup>)?)' if len(power_digits) > 1 else '(</sup>)?)'}",
            F"([{pval_input_range}][ \.0-9]*E ?(<sup>)?[-−] ?"
            F"{power_digits[0] + '?' if power_digits[0] == '0' else power_digits[0]}"
            F"{power_digits[1:] + '(</sup>)?)' if len(power_digits) > 1 else '(</sup>)?)'}",
            F"(?:[pP][- \(\)metavalueofcbind]" + "{0,13}"
            + F"[ =<of]*[ ]?)([{pval_input_range}][ \.0-9]*[ xX*×]*10 ?[-−] ?"
              F"{power_digits[0] + '?' if power_digits[0] == '0' else power_digits[0]}"
              F"{power_digits[1:] + ')' if len(power_digits) > 1 else ')'}"]
    return pval


def get_matching_data(input_file: str, bioc_pmcids: list, headers_skipped=False) -> dict:
    gc_data = {}
    with open(input_file, "r", encoding="utf-8") as f_in:
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
    return gc_data


def validate_relations(nlp, relations):
    new_relations = []
    for relation in relations:
        marker = relation.marker
        gene = relation.gene
        pval = relation.significance
        phenotype = relation.phenotype
        # above are mismatched!
        for gc_relation in nlp.gc_relations:
            if (marker and (marker.token.text == gc_relation[1])) or (gene and (gene.token.text == gc_relation[1])):
                phenotype_text = phenotype.token.ent_type_ if type(phenotype.token) != Span else \
                    phenotype.token.ents[0].label_
                if re.match(gc_relation[0], pval.token.text) and phenotype_text == gc_relation[2]:
                    new_relations.append(relation)
                    break

    return new_relations


def main():
    # load bioc pmc ids
    bioc_pmcids = [x.replace(".json", "").replace("_abbreviations", "") for x in listdir("BioC_Studies") if
                   isfile(join("BioC_Studies", x))]

    # retrieve matching data.
    gc_data = get_matching_data("GC_content.tsv", bioc_pmcids)

    lexicon = Ontology.get_master_lexicon()
    nlp = GCInterpreter(lexicon)
    failed_documents = []
    study_processing_times = []
    for pmc_id in gc_data.keys():
        if pmc_id != "PMC5536245":
            continue
        start_time = datetime.now()
        pvals = []
        rsids = []
        mesh_terms = []
        gc_relations = []
        nlp.reset_annotation_identifiers()
        for relation in gc_data[pmc_id]:
            rsid = relation[0]
            mesh_id = relation[2]
            new_pvals = generate_pval_regex_strings(relation[1])
            pvals.extend(new_pvals)
            rsids.append(F"({rsid})") #(?:\[[a-zA-Z0-9]\])?")
            mesh_terms.append(mesh_id)
            for x in new_pvals:
                gc_relations.append([x, rsid, mesh_id])
        nlp.set_ontology_terms([x for x in mesh_terms if x])
        nlp.pval_patterns = pvals
        nlp.rsid_patterns = rsids
        nlp.gc_relations = gc_relations

        study = Experimental.load_bioc_study("BioC_Studies", F"{pmc_id}.json")
        if not study:
            continue

        fulltext = "\n".join([x['text'] for x in study['documents'][0]['passages']])
        altered_text = re.sub(r"(?:\w)(\()", lambda x: x.group().replace("(", " ("), fulltext)

        abbreviations = nlp.get_all_abbreviations(altered_text)
        file_abbrevs = nlp.get_study_abbreviations(F"BioC_Studies/{pmc_id}_abbreviations.json")
        if file_abbrevs:
            abbreviations += file_abbrevs
        nlp.set_abbreviations(abbreviations)  # TODO: Check abbreviation partial entity HPC.

        result, nlp = process_study(nlp, study)
        nlp.clear_saved_study_data()

        study_tables, contains_annotations = GCTableExtractor.parse_tables(F"BioC_Studies/{pmc_id}_tables.json", nlp)

        time_taken = (datetime.now() - start_time).total_seconds()
        study_processing_times.append((pmc_id, time_taken))

        if contains_annotations:
            TableExtractor.output_tables(F"output/{pmc_id}_tables.json", study_tables)
        if not result['documents'][0]['relations'] and not contains_annotations:
            failed_documents.append(pmc_id)

    def avg(times):
        sum = 0
        for i in times:
            sum += i
        return sum / len(times) if times else 0

    print(F"Times taken: {study_processing_times}")
    print(F"Average processing time: {avg([y for x, y in study_processing_times if y < 600])}")


if __name__ == '__main__':
    main()
