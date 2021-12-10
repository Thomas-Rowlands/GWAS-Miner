import json
import re
from datetime import datetime
from os import listdir
from os.path import isfile, join

import spacy
from spacy.matcher import PhraseMatcher, Matcher
from spacy.tokens import Span

import Ontology
from GWAS_Miner import BioC, OutputConverter, Experimental
from NLP import Interpreter

gc_data = {}
headers_skipped = False


class GCInterpreter(Interpreter):

    def __init__(self, lexicon, ontology_only=False):
        self.lexicon = lexicon
        self.__nlp = spacy.load("en_core_sci_lg", disable=["ner"])
        self.__failed_matches = []
        self.__nlp.tokenizer.add_special_case(",", [{"ORTH": ","}])
        self.__basic_matcher = None
        self.__phrase_matcher = None
        self.__dep_matcher = None
        self.__entity_labels = ["MESH", "HPO"]
        if not ontology_only:
            self.__add_matchers(lexicon)

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
                        patterns.append(Interpreter.get_plural_variation(synonym["name"]))
                        if "," in synonym["name"]:
                            patterns.append(Interpreter.remove_comma_variation(synonym["name"]))
                patterns = self.__nlp.tokenizer.pipe(patterns)
                new_matcher.add(entry.identifier, patterns, on_match=self.__on_match)
        self.__phrase_matcher = new_matcher
        Span.set_extension("ontology", getter=self.get_ent_ontology)
        Span.set_extension("is_trait", getter=self.get_is_trait)
        # Add matcher patterns for parsing hyphenated and compound words.
        hyphenated_pattern = [{'POS': 'PROPN'}, {
            'IS_PUNCT': True, 'LOWER': '-'}, {'POS': 'VERB'}]
        compound_pattern = [{'POS': 'NOUN', 'DEP': 'compound'}, {
            'POS': 'NOUN', 'DEP': 'compound'}, {'POS': 'NOUN'}]
        self.__basic_matcher.add(
            "JOIN", [hyphenated_pattern, compound_pattern])

    def process_corpus(self, corpus, pval_patterns, rsid_patterns):
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

        old_ents, doc.ents = doc.ents, []
        for pattern in pval_patterns:
            self.__regex_match(pattern, doc, "PVAL")

        for pattern in rsid_patterns:
            self.__regex_match(pattern, doc, "RSID")

        self.__basic_matcher(doc)
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
        except:
            print(entity.text)#self.__failed_matches.append(entity.text)

    @staticmethod
    def __regex_match(pattern, doc, label):
        chars_to_tokens = {}
        for token in doc:
            for i in range(token.idx, token.idx + len(token.text)):
                chars_to_tokens[i] = token.i
        for match in re.finditer(pattern, doc.text, flags=re.IGNORECASE):
            start, end = match.span()
            span = doc.char_span(start, end, label=label, alignment_mode="expand")
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
                        print(e)

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


def process_study(nlp, study, pval_patterns, rsid_patterns):
    if not study:
        return False
    t, m, p = 0, 0, 0
    study_fulltext = "\n".join([x['text'] for x in study['documents'][0]['passages']])
    # abbreviations = nlp.get_all_abbreviations(study_fulltext)
    for passage in study['documents'][0]['passages']:
        passage_text = passage['text']
        # if abbreviations:
        #     for abbrev in abbreviations:
        #         passage_text = passage_text.replace(abbrev[0], abbrev[1])

        doc = nlp.process_corpus(passage_text, pval_patterns, rsid_patterns)
        for sent in doc.sents:
            training_sent = [x.label_ for x in sent.ents]
            if training_sent:
                if [x for x in training_sent if x[0] == "D"] and "RSID" in training_sent and "PVAL" in training_sent:
                    training_string = sent.text_with_ws
                    for ent in sent.ents:
                        ent_string = F"<!TRAIT:{ent.text_with_ws}!>" if ent.label_[
                                                                            0] == "D" else F"<!{ent.label_}:{ent.text_with_ws}!>"
                        training_string = training_string.replace(ent.text_with_ws, ent_string)
                    with open("training_input/training_input.txt", "a+", encoding="utf-8") as f_in:
                        f_in.write(training_string + "\n")

        annotations = nlp.get_entities(doc, ["MESH", "HPO", "RSID", "PVAL"])
        used_annots = []
        if annotations:
            for annot in annotations:
                loc = BioC.BioCLocation(offset=annot["offset"] + passage["offset"], length=annot["length"])
                current_datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
                if annot["text"] in used_annots:
                    for old_annot in passage['annotations']:
                        if old_annot.text == annot["text"]:
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

        # relations = nlp.extract_phenotypes(doc)
        relations = None
        if relations:
            print(relations)
    OutputConverter.output_xml(json.dumps(study, default=BioC.ComplexHandler),
                               F"output/PMC{study['documents'][0]['id']}_result.xml")
    return True


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

for pmc_id in gc_data.keys():
    pvals = []
    rsids = []
    mesh_terms = []
    for relation in gc_data[pmc_id]:
        rsid = relation[0]
        mesh_id = relation[2]
        pval = F"(?:p[ =<]*[ ]?)({relation[1][0]}[\.0-9 ]*{relation[1][1:]})"
        pvals.append(pval)
        pvals.append(pval.replace("e", " x 10"))
        pvals.append(pval.replace("e", " [*] 10"))
        pvals.append(pval.replace("e", "x10"))
        pvals.append(pval.replace("e", "[*]10"))
        temp = pvals
        for val in temp:
            pvals.append("".join(val.rsplit("0", 1)))
        rsids.append(F"({rsid})")
        mesh_terms.append(mesh_id)

    nlp.set_ontology_terms(mesh_terms)
    study = Experimental.load_bioc_study("BioC_Studies", F"{pmc_id}.json")
    result = process_study(nlp, study, pvals, rsids)
