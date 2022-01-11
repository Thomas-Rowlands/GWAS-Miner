import json
import re
import sys
from datetime import datetime
from os import listdir
from os.path import isfile, join

import networkx as nx
import spacy
from spacy.matcher import PhraseMatcher, Matcher, DependencyMatcher
from spacy.tokens import Span, Token, Doc

import Ontology
from GWAS_Miner import BioC, OutputConverter, Experimental
from GWAS_Miner.DataStructures import Marker, Significance, Phenotype, Association
from NLP import Interpreter
import config

gc_data = {}
headers_skipped = False


class GCInterpreter(Interpreter):

    def __init__(self, lexicon, ontology_only=False):
        self.lexicon = lexicon
        self.__nlp = spacy.load("en_core_sci_lg", disable=["ner"])
        self.__failed_matches = []
        # self.__nlp.tokenizer.add_special_case(",", [{"ORTH": ","}])
        infixes = self.__nlp.Defaults.infixes + [r'(?!\w)(\()']
        infix_regex = spacy.util.compile_infix_regex(infixes)
        self.__nlp.tokenizer.infix_finditer = infix_regex.finditer
        self.__basic_matcher = None
        self.__phrase_matcher = None
        self.__abbreviation_matcher = PhraseMatcher(self.__nlp.vocab)
        self.__entity_labels = ["MESH", "HPO", "PVAL"]
        self.pval_patterns = []
        self.rsid_patterns = []
        self.abbrev_pattens = []
        self.gc_relations = []
        self.association_patterns = config.pheno_assoc_patterns
        self.__dep_matcher = DependencyMatcher(self.__nlp.vocab, validate=True)
        self.__add_semgrex()
        if not ontology_only:
            self.__add_matchers(lexicon)

    def __add_semgrex(self):
        self.__dep_matcher.add("PHENO_ASSOC", self.association_patterns, on_match=self.__on_dep_match)

    def __on_dep_match(self, matcher, doc, id, matches):
        for match_id, token_ids in matches:
            pattern_name = self.__nlp.vocab[match_id].text
            entities = [doc[x].text for x in token_ids]
            if entities not in doc.user_data["relations"]["PHENO_ASSOC"]:
                print(F"{pattern_name}: {' | '.join(entities)}")
                doc.user_data["relations"]["PHENO_ASSOC"].append(entities)

    def set_abbreviations(self, abbrevs):
        result = []
        for abbrev in abbrevs:
            entry = self.lexicon.get_ordered_lexicons()[0].get_entry_by_term(abbrev[1])
            if entry:
                result.append([entry.identifier, abbrev[0]])

        self.abbrev_pattens = result

    def set_ontology_terms(self, term_ids):
        new_matcher = PhraseMatcher(self.__nlp.vocab, attr="LOWER")
        for lexicon in self.lexicon.get_ordered_lexicons():
            if lexicon.name == "HPO":
                continue
            for term_id in term_ids:
                entry = lexicon.get_entry_by_id(term_id)
                patterns = [entry.name()]
                if len(entry.name()) > 4:
                    patterns.append(Interpreter.get_plural_variation(entry.name()))
                if "," in entry.name():
                    patterns.append(Interpreter.remove_comma_variation(entry.name()))
                for synonym in entry.synonyms():
                    if synonym["name"] not in patterns:
                        patterns.append(synonym["name"])
                        patterns.append(Interpreter.get_plural_variation(synonym["name"]))
                        if "," in synonym["name"]:
                            patterns.append(Interpreter.remove_comma_variation(synonym["name"]))
                patterns = self.__nlp.tokenizer.pipe(patterns)
                new_matcher.add(entry.identifier, patterns, on_match=self.__on_match)
        self.__phrase_matcher = new_matcher

    def __add_matchers(self, lexicon):
        self.__basic_matcher = Matcher(self.__nlp.vocab)
        # self.__basic_matcher.add('marker', [[self.__marker_regex]], on_match=self.__on_match)

        new_matcher = PhraseMatcher(self.__nlp.vocab, attr="LOWER")
        for lexicon in lexicon.get_ordered_lexicons():
            if lexicon.name == "HPO":
                continue
            for entry in lexicon.get_entries():
                patterns = [entry.name()]
                if len(entry.name()) > 4:
                    patterns.append(Interpreter.get_plural_variation(entry.name()))
                if "," in entry.name():
                    patterns.append(Interpreter.remove_comma_variation(entry.name()))
                for synonym in entry.synonyms():
                    if synonym["name"] not in patterns:
                        patterns.append(synonym["name"])
                        if synonym["name"].upper() != synonym["name"]:
                            patterns.append(Interpreter.get_plural_variation(synonym["name"]))
                        if "," in synonym["name"]:
                            patterns.append(Interpreter.remove_comma_variation(synonym["name"]))
                patterns = self.__nlp.tokenizer.pipe(patterns)
                new_matcher.add(entry.identifier, patterns, on_match=self.__on_match)
        self.__phrase_matcher = new_matcher
        # Assign extension getters
        Token.set_extension("matches_ontology", getter=self.ontology_getter)
        Token.set_extension("is_trait", getter=self.is_trait_getter)
        Span.set_extension("has_ontology_term", getter=self.has_ontology_getter)
        Span.set_extension("has_trait", getter=self.has_trait_getter)
        Doc.set_extension("has_ontology_term", getter=self.has_ontology_getter)
        Doc.set_extension("has_trait", getter=self.has_trait_getter)
        # Add matcher patterns for parsing hyphenated and compound words.
        hyphenated_pattern = [{'POS': 'PROPN'}, {
            'IS_PUNCT': True, 'LOWER': '-'}, {'POS': 'VERB'}]
        compound_pattern = [{'POS': 'NOUN', 'DEP': 'compound'}, {
            'POS': 'NOUN', 'DEP': 'compound'}, {'POS': 'NOUN'}]
        self.__basic_matcher.add(
            "JOIN", [hyphenated_pattern, compound_pattern])

    def process_corpus(self, corpus):
        """[Applies tokenization, entity recognition and dependency parsing to the supplied corpus.]

        Args:
            corpus ([string]): [corpus text for information extraction]
            ontology_only (bool, optional): [Only apply ontology term matching to the supplied corpus]. Defaults to False.

        Returns:
            [SpaCy doc object]: [Parsed SpaCy doc object containing the processed input text with entities, tokens and dependencies.]
        """
        # Clean corpus with NLPre parsers
        # parsers = [dedash(), titlecaps(), separate_reference(), unidecoder()]
        # for parser in parsers:
        #     corpus = parser(corpus)

        doc = self.__nlp(corpus)
        doc.user_data["relations"] = {"PHENO_ASSOC": []}

        old_ents, doc.ents = doc.ents, []
        for pattern in self.pval_patterns:
            self.__regex_match(pattern, doc, "PVAL")

        for pattern in self.rsid_patterns:
            self.__regex_match(pattern, doc, "RSID")

        for pattern in self.abbrev_pattens:
            self.__regex_match(pattern[1], doc, pattern[0], ignore_case=False)

        # self.__basic_matcher(doc)
        self.__phrase_matcher(doc)

        # Ensure that rule-matched entities override data model entities when needed.
        for ent in old_ents:
            try:
                doc.ents += (ent,)
            except:  # Default SpaCy entities should never override others.
                continue

        # Ensure that multi-token entities are merged for extraction and association processing.
        for ent_label in self.__entity_labels:
            self.__merge_spans(doc, ent_label)

        # self.__dep_matcher(doc)

        return doc

    def __on_match(self, matcher, doc, i, matches):
        """
        (Event handler) Add matched entity to document entity list if no overlap is caused.
        @param matcher: Matcher object which fired the event
        @param doc: nlp doc object
        @param i: index of the current match
        @param matches: list of matches found by the matcher object
        """
        match_id, start, end = matches[i]
        entity = Span(doc, start, end,
                      label=self.__nlp.vocab.strings[match_id])
        # if entity.label_ in ["RSID", "PVAL"]:
        #     entity.set_extension("is_trait", default=False, force=True)
        #     entity.set_extension("ontology", default=None, force=True)
        # else:
        #     entity.set_extension("ontology", getter=self.get_ent_ontology, force=True)
        #     entity.set_extension("is_trait", default=True, force=True)
        try:
            doc.ents += (entity,)
        except Exception as Ex:
            entities_to_replace = []
            for ent in doc.ents:
                if (start <= ent.start < end) or (start < ent.end <= end):
                    if len(ent) <= len(entity):
                        entities_to_replace.append(ent)
            if entities_to_replace:
                doc.ents = [x for x in doc.ents if x not in entities_to_replace]
                doc.ents += (entity,)
            return  # print(entity.text)  # self.__failed_matches.append(entity.text)

    @staticmethod
    def __regex_match(pattern, doc, label, ignore_case=True):
        chars_to_tokens = {}
        for token in doc:
            for i in range(token.idx, token.idx + len(token.text)):
                chars_to_tokens[i] = token.i
        test = None
        if ignore_case:
            test = re.finditer(pattern, doc.text, flags=re.IGNORECASE)
        else:
            test = re.finditer(pattern, doc.text)
        for match in test:
            start, end = match.span()
            if label == "PVAL":
                span = doc.char_span(start, end, label=label, alignment_mode="expand")
            else:
                span = doc.char_span(start, end, label=label)
            if span is not None:
                try:
                    doc.ents += (span,)
                except:
                    test = None
            else:
                start_token = chars_to_tokens.get(start)
                end_token = chars_to_tokens.get(end)
                if start_token is not None and end_token is not None:
                    span = doc[start_token:end_token + 1]
                    try:
                        doc.ents += (span,)
                    except Exception as e:
                        return  # print(e)

    @staticmethod
    def __merge_spans(doc, entity_label):
        """
        Merge document spans that are part of the same entity.
        @param doc: Processed document
        @param entity_label: Label of entities to merge.
        """
        with doc.retokenize() as retokenizer:
            for ent in [x for x in doc.ents if x.label_ == entity_label]:
                retokenizer.merge(doc[ent.start:ent.end])

    def get_relations(self, doc):
        return self.__dep_matcher(doc)

    def calculate_sdp(self, phenotype_sents, top_phenotype=None):
        """[Calculates the shortest dependency path for each phenotype/marker/p-value combination, returning the shortest for each one.]

        Args:
            phenotype_sents ([list]): [List of SpaCy sent objects containing phenotype entities.]

        Returns:
            [dict]: [Dictionary containing extracted phenotype, marker and p-values associated together based on SDP calculation.]
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

            phenotypes = Interpreter._validate_node_entities(
                [x for x in sent.ents if x._.has_trait], graph.nodes) if not top_phenotype else None
            markers = Interpreter._validate_node_entities(
                [x for x in sent.ents if x.label_ == 'RSID'], graph.nodes)
            pvals = Interpreter._validate_node_entities(
                [x for x in sent.ents if x.label_ == 'PVAL'], graph.nodes)

            phenotype_count = len(phenotypes) if not top_phenotype else None
            marker_count = len(markers)
            pval_count = len(pvals)

            relations = []

            if phenotypes and markers and pvals:
                for pval in pvals:
                    for s, m, p in self.gc_relations:
                        significance = None
                        marker = None
                        phenotype = None
                        if re.search(s, pval[0].text):
                            significance = pval[1]
                            marker = next((x[1] for x in markers if x[0].text == m), None)
                            phenotype = next((x[1] for x in phenotypes if x[0].label_ == p), None)
                        if significance and marker and phenotype:
                            relations.append([significance, marker, phenotype])

                # Validate associations are triples
                pheno_assocs = [x for x in relations if len(x) == 3]
                for (rsid, pval, phenotype) in pheno_assocs:
                    temp_marker = sent.doc[token_indexes[int(rsid[rsid.find("<id") + 3: rsid.find(">", rsid.find("<id") + 3)])]]
                    temp_pval = sent.doc[token_indexes[int(pval[pval.find("<id") + 3: pval.find(">", pval.find("<id") + 3)])]]
                    temp_pheno = sent.doc[token_indexes[int(phenotype[phenotype.find("<id") + 3: phenotype.find(">", phenotype.find("<id") + 3)])]]
                    result_marker = Marker(temp_marker)
                    result_significance = Significance(temp_pval)
                    result_pheno = Phenotype(temp_pheno)
                    results.append(Association(marker=result_marker, significance=result_significance, phenotype=result_pheno))
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
                for (rsid, pval, phenotype) in pheno_assocs:
                    temp_marker = sent.doc[token_indexes[int(rsid[rsid.find("<id") + 3: rsid.find(">", rsid.find("<id") + 3)])]]
                    temp_pval = sent.doc[token_indexes[int(pval[pval.find("<id") + 3: pval.find(">", pval.find("<id") + 3)])]]
                    temp_pheno = sent.doc[phenotype.start:phenotype.end]
                    result_marker = Marker(temp_marker)
                    result_significance = Significance(temp_pval)
                    result_pheno = Phenotype(temp_pheno)
                    results.append(Association(marker=result_marker, significance=result_significance, phenotype=result_pheno))
        return results


def get_ubiquitous_phenotype(fulltext, nlp):
    doc = nlp.process_corpus(fulltext)
    top_phenotypes = nlp.get_phenotype_stats(doc, nlp.lexicon)
    top_phenotype = None
    for pheno in top_phenotypes:
        if not top_phenotype:
            top_phenotype = top_phenotypes[pheno]
            top_phenotype["label"] = pheno
            top_phenotype["ID"] = top_phenotypes[pheno]["ID"]
        elif top_phenotype["Count"] < top_phenotypes[pheno]["Count"]:
            top_phenotype = top_phenotypes[pheno]
            top_phenotype["label"] = pheno
            top_phenotype["ID"] = top_phenotypes[pheno]["ID"]
    return top_phenotype


def process_study(nlp, study):
    if not study:
        return False
    t, m, p, r = 0, 0, 0, 0
    study_fulltext = "\n".join([x['text'] for x in study['documents'][0]['passages']])
    # abbreviations = nlp.get_all_abbreviations(study_fulltext)
    current_datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    document_relations = []
    current_offset = 0
    results_present = False
    top_phenotype = get_ubiquitous_phenotype(study_fulltext, nlp)
    for passage in study['documents'][0]['passages']:
        if passage["infons"]["section_type"].lower() == "results":
            results_present = True
            break
    for passage in study['documents'][0]['passages']:
        if results_present and passage["infons"]["section_type"].lower() not in ["abstract", "results"]:
            continue
        passage_text = passage['text']
        # passage_text = re.sub(r"(?:\w)(\()", lambda x: x.group().replace("(", " ("), passage_text)
        doc = nlp.process_corpus(passage_text)
        doc.user_data["top_phenotype"] = top_phenotype
        annotations = nlp.get_entities(doc, ["MESH", "HPO", "RSID", "PVAL"])
        used_annots = []
        if annotations:
            for annot in annotations:
                loc = BioC.BioCLocation(offset=annot["offset"] + passage["offset"], length=annot["length"])
                if annot["text"] in used_annots:
                    for old_annot in passage['annotations']:
                        if old_annot.text == annot["text"] and loc not in old_annot.locations:
                            old_annot.locations.append(loc)
                if "RSID" not in annot["entity_type"] and "PVAL" not in annot["entity_type"]:
                    genomic_trait = BioC.BioCAnnotation(id=F"T{t}", infons={"type": "trait", "identifier": annot["id"],
                                                                            "annotator": "tr142@le.ac.uk",
                                                                            "updated_at": current_datetime},
                                                        locations=[loc], text=annot["text"])
                    passage['annotations'].append(genomic_trait)
                    t += 1
                elif "RSID" in annot["entity_type"]:
                    marker_identifier = BioC.BioCAnnotation(id=F"M{m}",
                                                            infons={"type": "genomic_marker", "identifier": annot["id"],
                                                                    "annotator": "tr142@le.ac.uk",
                                                                    "updated_at": current_datetime},
                                                            locations=[loc], text=annot["text"])
                    passage['annotations'].append(marker_identifier)
                    m += 1
                elif "PVAL" in annot["entity_type"]:
                    p_value = BioC.BioCAnnotation(id=F"P{p}", infons={"type": "significance", "identifier": annot["id"],
                                                                      "annotator": "tr142@le.ac.uk",
                                                                      "updated_at": current_datetime},
                                                  locations=[loc], text=annot["text"])
                    passage['annotations'].append(p_value)
                    p += 1
                used_annots.append(annot["text"])

            relations, uncertain_relations = nlp.extract_phenotypes(doc)
            pheno_id = None
            for relation in relations:
                # test = doc[]
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
                bioc_relation = BioC.BioCRelation(id=F"R{r}",
                                                  infons={"type": "disease_assoc", "annotator": "tr142@le.ac.uk",
                                                          "updated_at": current_datetime},
                                                  nodes=[phenotype_node, marker_node, significance_node])
                # passage['relations'].append(bioc_relation)
                document_relations.append(bioc_relation)
                r += 1
            for relation in uncertain_relations:
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
                bioc_relation = BioC.BioCRelation(id=F"R{r}",
                                                  infons={"type": "possible_disease_assoc", "annotator": "tr142@le.ac.uk",
                                                          "updated_at": current_datetime},
                                                  nodes=[phenotype_node, marker_node, significance_node])
                # passage['relations'].append(bioc_relation)
                document_relations.append(bioc_relation)
                r += 1
    with open("top_phenotypes.txt", "a+", encoding="utf-8") as fin:
        fin.write(F"{pmc_id}\t{top_phenotype['ID']}\t{top_phenotype['label']}\n")
    if document_relations:
        study['documents'][0]['relations'] = document_relations
    OutputConverter.output_xml(json.dumps(study, default=BioC.ComplexHandler),
                               F"output/PMC{study['documents'][0]['id']}_result.xml")
    return study


# load bioc pmc ids
bioc_pmcids = [x.replace(".json", "") for x in listdir("BioC_Studies") if isfile(join("BioC_Studies", x))]

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

lexicon = Ontology.get_master_lexicon()
nlp = GCInterpreter(lexicon)
failed_documents = []
for pmc_id in gc_data.keys():
    # if pmc_id != "PMC5581217":
    #     continue
    pvals = []
    rsids = []
    mesh_terms = []
    gc_relations = []
    for relation in gc_data[pmc_id]:
        rsid = relation[0]
        mesh_id = relation[2]
        pval_input_range = ""
        pval_input_int = 0
        if "." in relation[1]:
            pval_input_int = int(relation[1][:relation[1].find(".")])
        else:
            pval_input_int = int(relation[1][:relation[1].lower().find("e")])
        pval_input_range = F"{pval_input_int - 1}{pval_input_int}{pval_input_int + 1}"
        power_digits = relation[1][relation[1].find('-') + 1:]
        pval = F"(?:[pP][- \(\)metavalueofcbind]" + "{0,13}" + F"[ =<of]*[ ]?)([{pval_input_range}][ \.0-9]*[ xX*]*10 ?- ?" \
                                                          F"{power_digits[0] + '?' if power_digits[0] == '0' else power_digits[0]}" \
                                                          F"{power_digits[1:] + ')' if len(power_digits) > 1 else ')'}"
        pvals.append(pval)
        rsids.append(F"({rsid})")
        mesh_terms.append(mesh_id)
        gc_relations.append([pval, rsid, mesh_id])
    # nlp.set_ontology_terms(mesh_terms)
    nlp.pval_patterns = pvals
    nlp.rsid_patterns = rsids
    nlp.gc_relations = gc_relations
    study = Experimental.load_bioc_study("BioC_Studies", F"{pmc_id}.json")
    fulltext = "\n".join([x['text'] for x in study['documents'][0]['passages']])
    altered_text = re.sub(r"(?:\w)(\()", lambda x: x.group().replace("(", " ("), fulltext)
    abbreviations = nlp.get_all_abbreviations(altered_text)
    nlp.set_abbreviations(abbreviations)
    result = process_study(nlp, study)
    if not result['documents'][0]['relations']:
        failed_documents.append(pmc_id)
print(failed_documents)
# sys.exit()
