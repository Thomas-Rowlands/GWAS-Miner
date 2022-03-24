import itertools
import logging
import re
from datetime import datetime

from spacy.pipeline import merge_entities

import config
import networkx as nx
import spacy
from DataStructures import Marker, Phenotype, Significance, Association, LexiconEntry
from spacy import displacy
from spacy.matcher import Matcher, PhraseMatcher, DependencyMatcher
from spacy.tokens import Span, Token, Doc


class Interpreter:
    def __init__(self, lexicon, ontology_only=False):
        self.lexicon = lexicon
        self.nlp = spacy.load("en_core_sci_lg", disable=["ner"])
        self.nlp.add_pipe("merge_noun_chunks")
        self.__failed_matches = []
        self.nlp.tokenizer.add_special_case(",", [{"ORTH": ","}])
        self.__rsid_regex = {"LOWER": {"REGEX": "((?:[(]?)(rs[0-9]{1,}){1,})"}}  # "(?:rs[0-9]{1,}){1}"}}
        self.__marker_regex = {"TEXT": {"REGEX": r"([ATCG]{1}[a-z]{1,}[0-9]{1,}[ATCG]{1}[a-z]{1,})"}}
        self.__gene_seq_regex = {"TEXT": {"REGEX": "([ ][ACTG]{3,}[ ])"}}
        self.__basic_matcher = None
        self.__phrase_matcher = None
        self.__dep_matcher = None
        self.__logger = logging.getLogger("GWAS Miner")
        self.__entity_labels = ["MESH", "HPO", "PVAL"]
        self.lexicon = lexicon
        self.nlp = spacy.load("en_core_sci_lg", disable=["ner"])
        self.__failed_matches = []
        # self.nlp.tokenizer.add_special_case(",", [{"ORTH": ","}])
        infixes = self.nlp.Defaults.infixes + [r'(?!\w)(\()']
        infix_regex = spacy.util.compile_infix_regex(infixes)
        self.nlp.tokenizer.infix_finditer = infix_regex.finditer
        self.__basic_matcher = None
        self.__phrase_matcher = None
        self.__abbreviation_matcher = PhraseMatcher(self.nlp.vocab, attr="ORTH")
        self.__entity_labels = ["MESH", "HPO", "PVAL"]
        self.pval_patterns = []
        self.rsid_patterns = []
        self.abbrev_pattens = []
        self.t = 0
        self.v = 0
        self.s = 0
        self.g = 0
        self.r = 0
        self.annotations = []
        self.relations = []
        self.association_patterns = config.pheno_assoc_patterns
        if not ontology_only:
            self.__add_matchers(lexicon)

    def __add_matchers(self, lexicon):
        self.__basic_matcher = Matcher(self.nlp.vocab)
        # self.__basic_matcher.add('marker', [[self.__marker_regex]], on_match=self.__on_match)

        new_matcher = PhraseMatcher(self.nlp.vocab, attr="LOWER")
        for lexicon in lexicon.get_ordered_lexicons():
            if lexicon.name == "HPO":
                continue
            for entry in lexicon.get_entries():
                patterns = Interpreter.get_term_variations(entry)
                patterns = self.nlp.tokenizer.pipe(patterns)
                new_matcher.add(entry.identifier, patterns, on_match=self.__on_match)
        self.__phrase_matcher = new_matcher
        # Assign extension getters
        Token.set_extension("matches_ontology", getter=self.ontology_getter)
        Token.set_extension("is_trait", getter=self.is_trait_getter)
        Token.set_extension("has_trait", getter=self.has_trait_getter)
        Span.set_extension("has_ontology_term", getter=self.has_ontology_getter)
        Span.set_extension("has_trait", getter=self.has_trait_getter)
        Span.set_extension("is_trait", getter=self.is_trait_getter)
        Doc.set_extension("has_ontology_term", getter=self.has_ontology_getter)
        Doc.set_extension("has_trait", getter=self.has_trait_getter)
        Doc.set_extension("is_trait", getter=self.is_trait_getter)

    @staticmethod
    def get_term_variations(term: LexiconEntry):
        """
        Calculate every plausible variation of the input term, including synonyms.
        :param term: LexiconEntry object containing the desired term.
        :return: List of string variations.
        """
        patterns = None
        if type(term) == str:
            patterns = [term]
        else:
            patterns = [term.name()]
            if term.synonyms():
                for syn in term.synonyms():
                    patterns.append(syn['name'])
        funcs = [Interpreter.remove_comma_variation, Interpreter.get_hyphenated_variations,
                 Interpreter.get_roman_numeral_variation, Interpreter.get_plural_variation] # TODO: plurals + reverse the roman numerals too!
        new_patterns = []
        combos = [itertools.combinations(funcs, 1), itertools.combinations(funcs, 2),
                  itertools.combinations(funcs, 3), itertools.combinations(funcs, 4)]
        combos = [x for y in combos for x in y]
        for pattern in patterns:
            for combo in combos:
                for func in combo:
                    new_variant = func(pattern)
                    if new_variant and new_variant not in new_patterns:
                        if type(new_variant) == list:
                            new_patterns += new_variant
                        else:
                            new_patterns.append(new_variant)
        new_patterns = list(set(patterns + new_patterns))
        return new_patterns

    @staticmethod
    def get_plural_variation(term: str):
        if term.lower()[-1] == "s":
            return term[:-1]
        else:
            return term + "s"

    @staticmethod
    def get_hyphenated_variations(term: str):
        output = []
        # Remove hyphens.
        if "-" in term:
            output.append(term.replace("-", " "))
        # Add each variation of hyphenations.
        if " " in term:
            location = term.find(" ")
            for i in range(term.count(" ")):
                hyphenated = term[:location] + "-" + term[location + 1:]
                output.append(hyphenated)
                location = term.find(" ", location + 1)
                # no further spaces found.
                if location == -1:
                    break
        return output

    @staticmethod
    def get_roman_numeral_variation(term: str):
        search_result = re.search(r"(\d+)", term)
        result = None
        replacement = ''
        if search_result:
            num = search_result.group(0)
            replacement = num
            # Adapted from solution found here
            # https://www.w3resource.com/python-exercises/class-exercises/python-class-exercise-1.php
            val = [
                1000, 900, 500, 400,
                100, 90, 50, 40,
                10, 9, 5, 4,
                1
            ]
            syb = [
                "M", "CM", "D", "CD",
                "C", "XC", "L", "XL",
                "X", "IX", "V", "IV",
                "I"
            ]
            roman_num = ''
            i = 0
            num = int(num)
            while num > 0:
                for _ in range(num // val[i]):
                    roman_num += syb[i]
                    num -= val[i]
                i += 1
            result = term.replace(replacement, roman_num)
        return result

    @staticmethod
    def remove_comma_variation(term: str):
        if "," in term:
            adjusted_term = term.split(",")
            adjusted_term = adjusted_term[1] + " " + adjusted_term[0]
            adjusted_term = adjusted_term.lstrip()
            return adjusted_term
        else:
            return None

    def add_rule_matcher(self, label, rule):
        self.__basic_matcher.add(label, self.__on_match, rule)

    # def add_study_specific_abbreviations(self, abbrevs):
    #     self.nlp

    def ontology_getter(self, token):
        if token.ent_type_:
            if "HP:" in token.ent_type_:
                return "HPO"
            elif token.ent_type_[0] == "D":
                return "MESH"
            else:
                return False
        else:
            return False

    def has_ontology_getter(self, obj):
        if any(["HP:" in token.ent_type_ or token.ent_type_[0] == "D" for token in obj if token.ent_type_]):
            return True
        else:
            return False

    def is_trait_getter(self, token):
        if type(token) != Token:
            return False
        if token.ent_type_:
            if "HP:" in token.ent_type_ or token.ent_type_[0] == "D":
                return True
            else:
                return False
        else:
            return False

    def has_trait_getter(self, obj):
        if type(obj) != Span:
            return False
        if any(["HP:" in token.ent_type_ or token.ent_type_[0] == "D" for token in obj if token.ent_type_]):
            return True
        else:
            return False

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
                      label=self.nlp.vocab.strings[match_id])
        # if entity.label_ in ["RSID", "PVAL"]:
        #     entity.set_extension("is_trait", default=False, force=True)
        #     entity.set_extension("ontology", default=None, force=True)
        # else:
        #     entity.set_extension("ontology", getter=self.get_ent_ontology, force=True)
        #     entity.set_extension("is_trait", default=True, force=True)
        try:
            doc.ents += (entity,)
        except:
            self.__logger.log(level=0, msg=entity.text)
            # print(entity.text)#self.__failed_matches.append(entity.text)

    @staticmethod
    def __regex_match(pattern, doc, label):
        chars_to_tokens = {}
        for token in doc:
            for i in range(token.idx, token.idx + len(token.text)):
                chars_to_tokens[i] = token.i
        for match in re.finditer(pattern, doc.text, flags=re.IGNORECASE):
            start, end = match.span(1)
            span = None
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
                        print(e)

    def process_corpus(self, corpus, ontology_only=False):
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

        doc = self.nlp(corpus)

        old_ents, doc.ents = doc.ents, []

        #  Additional regex matches unnecessary when limited to ontology entities.
        if not ontology_only:
            for ent_label in config.regex_entity_patterns:
                if isinstance(config.regex_entity_patterns[ent_label], list):
                    for pattern in config.regex_entity_patterns[ent_label]:
                        self.__regex_match(pattern, doc, ent_label)
                        if ent_label not in self.__entity_labels:
                            self.__entity_labels.append(ent_label)
                else:
                    self.__regex_match(config.regex_entity_patterns[ent_label], doc, ent_label)
                    if ent_label not in self.__entity_labels:
                        self.__entity_labels.append(ent_label)

        self.__basic_matcher(doc)
        self.__phrase_matcher(doc)

        # Ensure that rule-matched entities override data model entities when needed.
        if not ontology_only:
            for ent in old_ents:
                try:
                    doc.ents += (ent,)
                except:  # Default SpaCy entities should never override others.
                    continue

        # Ensure that multi-token entities are merged for extraction and association processing.
        for ent_label in self.__entity_labels:
            self.__merge_spans(doc, ent_label)
        doc = merge_entities(doc)

        return doc

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

    @staticmethod
    def __filter_sents_by_entity(sents, entity_list, property_list=[]):
        """
        Remove sentence objects from a list if they do not contain all of the provided entities.
        @param sents: List of sentence objects, with nested lists for OR conditions
        @param entity_list: List of entity label strings to search for
        @return: List of sentence objects containing all of the required entities.
        """
        output = []
        for sent in sents:
            missing_entity = False
            missing_property = False
            ents = sent.ents
            for ent in entity_list:
                if type(ent) == list:
                    found_match = False
                    for or_ent in ent:
                        if or_ent in [x.label_ for x in ents]:
                            found_match = True
                            break
                    if not found_match:
                        missing_entity = True
                        break
                elif ent not in [x.label_ for x in ents]:
                    missing_entity = True
                    break
            for property in property_list:
                if type(property) == list:
                    found_match = False
                    for or_prop in property:
                        if [x.has_extension(or_prop) for x in ents]:
                            found_match = True
                            break
                    if not found_match:
                        missing_property = True
                        break
                elif [x.has_extension(property) for x in ents]:
                    missing_property = True
                    break
            if not missing_entity:
                output.append(sent)
        return output

    @staticmethod
    def _validate_node_entities(ents, nodes):
        output = []
        used_indexes = []
        for item in ents:
            token_count = len([x for x in item if type(x) == Token])
            if token_count > 1:
                split_item = [x.text + F"<id{x.idx}>" for x in item]
                for node in nodes:
                    for token in split_item:
                        if token == str(node):
                            output.append((item, F"{item}<id{item.start_char}>"))
                            break
            else:
                for node in nodes:
                    if F"{item.text}<id{item.start_char}>" == str(node):
                        output.append((item, F"{item}<id{item.start_char}>"))
        return output

    @staticmethod
    def __validate_phenotype_context(token):
        temp = []
        for left in token.lefts:
            temp.append(left)
            if left.n_lefts == 0:
                break
        head = token

    @staticmethod
    def __expand_sentence_dependency_search(token):
        contains_pheno = False
        contains_marker = False
        next_token = None
        old_token = None
        if type(token) == Span:
            next_token = token.root.head
            old_token = token.root
        else:
            next_token = token.head
            old_token = token
        # Expand the list of tokens which are related to the original token (unless containing a subordinating
        # conjunction)
        while old_token != next_token:
            if next_token.pos_ == 'SCONJ':
                return False, False
            old_token = next_token
            next_token = old_token.head
        for child in next_token.subtree:
            if contains_pheno and contains_marker:
                return contains_pheno, contains_marker
            if child._.is_trait:
                contains_pheno = True
            elif child.ent_type_ == "RSID":
                contains_marker = True
        return contains_pheno, contains_marker

    @staticmethod
    def allocate_contiguous_phenotypes(sent):
        """[Identify phenotype entities bordering p-value entities]

        Args:
            sent ([SpaCy sent object]): [NLP processed sentence object]

        Returns:
            [list]: [list containing tuples of phenotypes with contiguous p-value entities]
        """
        result = []
        prev_token = None
        for token in sent:
            if prev_token:
                if "PVAL" in token.ent_type_ and prev_token.ent_type_ == "MESH":
                    result.append(
                        [(prev_token, F"{prev_token.lower_}<id{prev_token.i}>"),
                         (token, F"{token.lower_}<id{token.i}>")])
            prev_token = token
        return result

    @staticmethod
    def trim_brackets(input):
        removed_brackets = 0
        if input[0] == "(":
            input = input[1:]
            removed_brackets += 1
        if input[-1] == ")":
            input = input[:-1]
            removed_brackets += 1
        return input

    @staticmethod
    def get_entities(doc, entity_list):
        result = []
        for ent in doc.ents:
            result.append({
                "entity_type": ent.label_,
                "id": ent.label_[ent.label_.index(":") + 2:] if ":" in ent.label_ else ent.label_,
                "text": ent.text,  # self.trim_brackets(ent.text),
                "offset": ent.start_char,
                "length": ent.end_char - ent.start_char
            })
        return result

    def extract_phenotypes(self, doc):
        """[Extract phenotype to genotype associations from the provided SpaCy doc object]

        Args:
            doc ([SpaCy doc object]): [SpaCy processed document containing entities for extraction]

        Returns:
            [dict]: [Dictionary containing extracted phenotype, marker and p-values associated together based on SDP calculation.]
        """
        phenotype_sents = Interpreter.__filter_sents_by_entity(
            doc.sents, ["PVAL", "RSID"], ["has_ontology"])
        uncertain_sents = Interpreter.__filter_sents_by_entity(doc.sents, ["PVAL", "RSID"])
        results, uncertain_results = [], []
        results = self.calculate_sdp(phenotype_sents)
        uncertain_results = self.calculate_sdp(uncertain_sents, doc.user_data["top_phenotype"])
        filtered_uncertain_results = []
        if results and uncertain_results:
            for association in results:
                match_found = False
                for u_association in uncertain_results:
                    if association.marker.token.text == u_association.marker.token.text and association.significance.token.text == u_association.significance.token.text:
                        match_found = True
                        break
                if not match_found:
                    filtered_uncertain_results += uncertain_results
        else:
            filtered_uncertain_results = uncertain_results
        return results, filtered_uncertain_results

    @staticmethod
    def calculate_sdp(phenotype_sents, top_phenotype=None):
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
            # immediate_relations = Interpreter.allocate_contiguous_phenotypes(sent)

            pheno_assocs = []
            for marker in markers:
                best_pval_distance = None
                best_pval = None
                for pval in pvals:
                    if pval[1] in [y for (x, y) in pheno_assocs]:
                        continue
                    temp_distance = nx.shortest_path_length(
                        graph, source=marker[1], target=pval[1])
                    if not best_pval_distance or temp_distance < best_pval_distance:
                        best_pval = nx.shortest_path(
                            graph, source=marker[1], target=pval[1])[-1]
                        best_pval_distance = temp_distance
                    else:
                        continue
                if not pval:
                    continue
                if best_pval_distance:
                    pheno_assocs.append([marker[1], best_pval])
            if not top_phenotype:
                for pair in pheno_assocs:
                    best_pheno_distance = None
                    best_pheno = None
                    for phenotype in phenotypes:
                        temp_distance = nx.shortest_path_length(
                            graph, target=phenotype[1], source=pair[0])
                        if not best_pheno_distance or temp_distance < best_pheno_distance:
                            best_pheno = nx.shortest_path(
                                graph, target=phenotype[1], source=pair[0])[-1]
                            best_pheno_distance = temp_distance
                        else:
                            continue
                    if not pval:
                        continue
                    if best_pheno_distance:
                        pair.append(best_pheno)

            if not top_phenotype:
                # Validate associations are triples
                pheno_assocs = [x for x in pheno_assocs if len(x) == 3]
                for (rsid, pval, phenotype) in pheno_assocs:
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
            else:
                # Validate associations are doubles
                pheno_assocs = [x for x in pheno_assocs if len(x) == 2]
                for (rsid, pval) in pheno_assocs:
                    temp_marker = sent.doc[
                        token_indexes[int(rsid[rsid.find("<id") + 3: rsid.find(">", rsid.find("<id") + 3)])]]
                    temp_pval = sent.doc[
                        token_indexes[int(pval[pval.find("<id") + 3: pval.find(">", pval.find("<id") + 3)])]]
                    result_marker = Marker(temp_marker)
                    result_significance = Significance(temp_pval)
                    result_pheno = None
                    for ent in sent.doc.ents:
                        if ent.label_ == top_phenotype["ID"]:
                            result_pheno = Phenotype(next(x for x in ent.sent if x.ent_type_ == ent.label_))
                            break

                    if result_pheno:
                        results.append(
                            Association(marker=result_marker, significance=result_significance, phenotype=result_pheno))
        return results

    @staticmethod
    def count_entities(doc, label):
        """[Counts the number of named entities matching the provided label string]

        Args:
            doc ([SpaCy doc object]): SpaCy processed document containing entities to count]
            label ([string]): [Label of entities that will be counted e.g. 'MESH'.]

        Returns:
            [int]: [Total number of entities found matching the supplied label]
        """
        count = 0
        for ent in doc.ents:
            if ent.label_ == label:
                Interpreter.__logger.info((ent.lower_))
                count += 1
        return count

    @staticmethod
    def split_doc_sents(doc):
        """
        Retrieve a list of sentence objects from a Spacy Doc object.
        @param doc: Spacy Doc object.
        @return: List of sentence objects
        """
        return [x for x in list(doc.sents)]

    @staticmethod
    def merge_reliant(doc, word):
        """
        Merge words with their neighbours where appropriate relations are present (e.g. compounds)
        @param doc: The document to which the base word belongs
        @param word: The base word object to search from
        @return: The word accompanied by it's reliant siblings
        """
        result = {"original": str(word), "index": word.i,
                  "l_words": [], "replacement": str(word)}
        while word.n_lefts > 0:
            left_word = doc[word.i - 1]
            if left_word.dep_ == "punct":
                result["replacement"] = F"{left_word}{result['replacement']}"
                result["l_words"].append(left_word.idx)
                left_word = doc[word.i - 2]
            if left_word.dep_ in ["compound", "npadvmod"] and word.n_lefts > 0:
                result[
                    "replacement"] = F"{left_word}{'' if (result['replacement'][0:1] == '-') else ' '}{result['replacement']}"
                result["l_words"].append(left_word.idx)
            word = left_word
        result['l_words'].sort()
        if result['l_words']:
            result['l_words'] = result['l_words'][:1][0]
        return result

    @staticmethod
    def replace_list(doc):
        """
        Replace a list within a sentence with multiple statements e.g. epi-, collagen- and adp-induced platelet...
        @param doc: The sentence string to be changed
        @return: List of new sentences generated from the input statement
        """
        results = []
        temp_word = ""
        temp_word_two = ""
        is_changed = False
        for word in doc:
            if word.dep_ == "amod":  # adjectival modifier
                temp_word = Interpreter.merge_reliant(doc, word)
                temp_word_two = Interpreter.merge_reliant(doc, word.head)
                results.append({"original": str(word), "index": word.idx,
                                "l_words": [temp_word['l_words'], temp_word_two['l_words']],
                                "replacement": F"{temp_word['replacement']} {temp_word_two['replacement']}"})
                is_changed = True
        if not is_changed:
            return None
        new_text = doc.text
        output = {"items": [], "statements": []}
        for result in results:
            output["statements"].append(
                F"{new_text[:results[0]['l_words'][0]]}{result['replacement']}.")
            output['items'].append(result['replacement'])
        return output

    @staticmethod
    def display_structure(sentence_spans, markup_only=False, theme="light"):
        """
        Start running the Displacy visualization of the tokenized sentences identified by the NLP pipeline.
        @param doc: The NLP processed document.
        @param markup_only: If True, returns a HTML string instead of hosting.
        """
        options = None
        if theme.lower() == "light":
            options = {"compact": True}
        else:
            options = {"compact": True, "bg": "#19232D", "color": "white"}
        # __logger.info([doc[match[1]:match[2]] for match in matches])
        if markup_only:
            svgs = []
            for sent in sentence_spans:
                ent_labels = [x.label_ for x in sent.ents]
                ent_traits = [x._.has_trait for x in sent.ents]
                if "RSID" in ent_labels and "PVAL" in ent_labels and ent_traits:
                    svgs.append(displacy.render(sent, style="dep", options=options))
            return svgs
        displacy.serve(sentence_spans, style="dep", options=options)

    @staticmethod
    def display_ents(doc, markup_only=False, theme="light"):
        """
        Start running the Displacy visualization of the named entities recognised by the NLP pipeline.
        @param doc: The NLP processed document.
        @param markup_only: If True, returns a HTML string instead of hosting.
        """
        colors = None
        if theme == "light":
            colors = {"MESH": "rgb(247, 66, 145)", "EFO": "rgb(247, 66, 145)", "HPO": "rgb(147, 66, 245)",
                      "RSID": "rgb(245, 66, 72)",
                      "PVAL": "rgb(102, 255, 51)", "PTYPE": "rgb(51, 102, 255)", "marker": "rgb(0, 255, 204)"}
        else:
            colors = {"MESH": "rgb(156, 25, 82)", "EFO": "rgb(217, 36, 115)", "HPO": "rgb(117, 36, 215)",
                      "RSID": "rgb(140, 22, 26)",
                      "PVAL": "rgb(48, 150, 14)", "PTYPE": "rgb(21, 72, 225)", "marker": "rgb(0, 225, 174)",
                      "TABLE REF": "rgb(61, 38, 6)"}
        options = {"colors": colors}
        if markup_only:
            return displacy.render(doc, style="ent", page=True, options=options, jupyter=False)
        displacy.serve(doc, style="ent", options=options)

    def onto_match(self, doc):
        doc = self.process_corpus(doc, ontology_only=True)
        return [(x.text, x.label_) for x in doc.ents]

    @staticmethod
    def get_phenotype_stats(doc, master_lexicon):
        results = {}
        for ent in doc.ents:
            if ent._.has_trait:  # in ["MESH", "HPO"]:
                entry = master_lexicon.get_lexicon_entry(lexicon_name="MESH", ident=ent.label_)
                # if not entry:
                # print(ent)
                if entry.name() in results:
                    results[entry.name()]["Count"] += 1
                else:
                    results[entry.name()] = {}
                    results[entry.name()]["Count"] = 1
                    results[entry.name()]["Ontology"] = "MESH"
                    results[entry.name()]["ID"] = entry.identifier

        return results

    @staticmethod
    def insert_phrase(abbrevs, token):
        """
        Replaces an abbreviated form of a phrase/name with the complete version, assuming it is provided.
        @param abbrevs: Counter provided by NLPre's identify_parenthetical_phrases function
        @param token: Abbreviation to be replaced
        @return: String containing the complete phrase or None
        """
        for abbrev in abbrevs:
            if abbrev[1] == token:
                result = ""
                for word in abbrev[0]:
                    result += word + " "
                return result[:-1]
        return None

    @staticmethod
    def __check_single_word_abbrev_proto(fulltext, token):
        temp = [F"(\\b[{token[0]}{token[0].lower()}][\w]+)(?:[a-z\s-]*)"]
        for char in token[1:]:
            temp.append(F"(\\b[{char}{char.lower()}][\w]+)(?:[a-z\s-]*)")
        pattern = F"{''.join(temp)}(?:[ ]\({token}[^\w]*\))"
        match = re.search(pattern, fulltext)
        if not match:
            # print(token)
            return None
        result = []
        for group in match.groups():
            result.append(group)
        return " ".join(result)

    @staticmethod
    def __check_single_word_abbrev(fulltext, token):
        """
        Locate the expanded phrase for a single abbreviation
        @param fulltext: Text containing the declaration of the abbreviation
        @param token: Abbreviation to be expanded
        @return: The expanded text for the abbreviation.
        """
        # Identify first letter of the first word for the abbreviation
        result = ""
        first_char = token[0:1]
        # Get the number of occurrences of this letter in abbreviation.
        first_char_count = token.count(first_char)
        # Locate abbreviations in parenthesis
        search_regex = r"([ \-'\n\w0-9]+\n?\({0}[^a-zA-Z]{0,}\))".replace("{0}", re.escape(
            token))  # r"([ \-'\na-zA-Z0-9]+\n?\({0}\))".format(re.escape(token))
        declaration_match = re.search(
            search_regex, fulltext, re.IGNORECASE | re.MULTILINE)
        if declaration_match is None:  # No declaration found for token
            return None
        # First match SHOULD be the declaration of this abbreviation.
        # Split REGEX match to a list of words
        split_sent = declaration_match.group(0).replace(
            " )", ")").replace("( ", "(").replace("\n", "").split(" ")
        found_counter = 0
        found_indexes = []
        i = len(split_sent) - 2  # Indexing + ignore the actual abbreviation
        #  Moving backwards from the abbreviation, count each word in the sentence matching the first character.
        while i >= 0:
            if split_sent[i][0:1] == first_char:
                found_counter += 1
                found_indexes.append(i)
            i -= 1
        #  Add each word following (inclusively) the nth word with a matching first character.
        if first_char_count <= found_counter:
            found_indexes.sort()
            for x in split_sent[found_indexes[(first_char_count - found_counter)]:-1]:
                result += x + " "
        result = result.strip()
        if result:
            return result
        else:
            return False

    @staticmethod
    def replace_all_abbreviations(fulltext, section=None):
        """
        Returns the expanded form of an abbreviation from the fulltext.
        @param fulltext: The document containing the abbreviations and their declarations
        @return: String containing the expanded version of the abbreviation.
        """
        changes = []
        # r"([^(-]\b[a-z]{0,}[A-Z]{2,}[a-z]{0,}\b[^)-])"
        pattern = r"([^ \"',.(-]\b)?([a-z]{0,})([A-Z]{2,})([a-z]{0,})(\b[^;,.'\" )-]?)"
        input_text = None
        if section:
            input_text = section
        else:
            input_text = fulltext
        for match in re.findall(pattern, input_text):
            target = ""
            if type(match) == str:
                target = match
            else:
                target = match[2]
            target = target.strip()
            if target not in [x for [x, y] in changes]:
                changes.append(
                    [target, Interpreter.replace_abbreviations(target, fulltext)])
        for change in changes:
            if change[1]:
                if change[1] != change[0]:
                    input_text = input_text.replace(change[0], change[1])
        # Interpreter.__logger.info(changes) #  Can error due to strange encodings used.
        return Interpreter.__clean_reference_remains(input_text)

    @staticmethod
    def get_all_abbreviations(fulltext):
        """
        Returns a list of all abbreviations and their expanded forms.
        @param fulltext: The document containing the abbreviations and their declarations
        @return: 2D list containing abbreviations and their expanded forms.
        """
        result = []
        # r"([^(-]\b[a-z]{0,}[A-Z]{2,}[a-z]{0,}\b[^)-])"
        # pattern = r"([^ \"',.(-]\b)?([a-z]{0,})([A-Z]{2,})([a-z]{0,})(\b[^;,.'\" )-]?)"
        pattern = r"(?:\()([a-z]{0,3}[A-Z]{2,}[a-z]{0,3})(?:\))"
        input_text = fulltext
        matches = set(re.findall(pattern, input_text))
        for match in matches:
            target = ""
            if type(match) == str:
                target = match
            else:
                continue
            target = target.strip()
            # expanded = Interpreter.__check_single_word_abbrev(fulltext, target)
            expanded = Interpreter.__check_single_word_abbrev_proto(fulltext, target)
            # print(F"Expanded ({target}) - {datetime.now()}")
            if expanded:
                result.append([target, expanded.strip()])
        # Interpreter.__logger.info(changes) #  Can error due to strange encodings used.
        return result

    @staticmethod
    def __clean_reference_remains(text):
        return text.replace("()", "").replace("(, )", "")

    @staticmethod
    def replace_abbreviations(token, fulltext):
        """
        Returns the expanded form of an abbreviation from the fulltext.
        @param token: The abbreviation to be expanded.
        @param fulltext: The document containing the declaration of the abbreviation.
        @return: String containing the expanded version of the abbreviation.
        """
        # Remove all preceding and trailing white space.
        doc = token.strip()
        if len(doc) < 5:
            result = Interpreter.__check_single_word_abbrev(
                fulltext, doc.upper())
            if result:
                return result
            else:
                return doc
        else:
            return doc
