import logging
import re

import config
import networkx as nx
import spacy
from DataStructures import Marker
# from nlpre import dedash, titlecaps, separate_reference, unidecoder, identify_parenthetical_phrases, replace_acronyms
from spacy import displacy
from spacy.matcher import Matcher
from spacy.matcher import PhraseMatcher
from spacy.tokens import Span


class Interpreter:
    def __init__(self, lexicon, ontology_only=False):
        self.__nlp = spacy.load("en_core_sci_md", disable=["ner"])
        self.__failed_matches = []
        self.__nlp.tokenizer.add_special_case(",", [{"ORTH": ","}])
        self.__rsid_regex = [{"TEXT": {"REGEX": "(?:rs[0-9]{1,}){1}"}}]
        self.__marker_regex = [
            {"TEXT": {"REGEX": r"([ATCG]{1}[a-z]{1,}[0-9]{1,}[ATCG]{1}[a-z]{1,})"}}]
        self.__gene_seq_regex = [{"TEXT": {"REGEX": "([ ][ACTG]{3,}[ ])"}}]
        self.__basic_matcher = None
        self.__phrase_matchers = []
        self.__logger = logging.getLogger("GWAS Miner")
        self.__entity_labels = ["MeSH", "HPO"]
        if not ontology_only:
            self.__add_matchers(lexicon)

    def __add_matchers(self, lexicon):
        self.__basic_matcher = Matcher(self.__nlp.vocab)
        self.__basic_matcher.add('RSID', self.__on_match, self.__rsid_regex)
        self.__basic_matcher.add('marker', self.__on_match, self.__marker_regex)
        for entry in lexicon.keys():
            new_matcher = PhraseMatcher(self.__nlp.vocab, attr="LOWER")
            sub_list = sorted(lexicon[entry].keys(), reverse=True)
            for sub_entry in sub_list:
                patterns = list(self.__nlp.tokenizer.pipe(
                    lexicon[entry][sub_entry]))
                new_matcher.add(entry, self.__on_match, *patterns)
            self.__phrase_matchers.append(new_matcher)

    def add_rule_matcher(self, label, rule):
        self.__basic_matcher.add(label, self.__on_match, rule)

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
        try:
            doc.ents += (entity,)
        except:
            self.__failed_matches.append(entity.text)

    @staticmethod
    def __regex_match(pattern, doc, label):
        for match in re.finditer(pattern, doc.text, flags=re.IGNORECASE):
            start, end = match.span()
            span = doc.char_span(start, end, label=label)
            if span is not None:
                try:
                    doc.ents += (span,)
                except:
                    test = None

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

        doc = self.__nlp(corpus)

        old_ents, doc.ents = doc.ents, []

        #  Additional regex matches unnecessary when limited to ontology entities.
        if not ontology_only:
            for ent_label in config.regex_entity_patterns:
                if isinstance(config.regex_entity_patterns[ent_label], list):
                    for pattern in config.regex_entity_patterns[ent_label]:
                        self.__regex_match(pattern, doc, ent_label)
                        self.__entity_labels.append(ent_label)
                else:
                    self.__regex_match(config.regex_entity_patterns[ent_label], doc, ent_label)
                    self.__entity_labels.append(ent_label)

        # Add matcher patterns for parsing hyphenated and compound words.
        hyphenated_pattern = [{'POS': 'PROPN'}, {
            'IS_PUNCT': True, 'LOWER': '-'}, {'POS': 'VERB'}]
        compound_pattern = [{'POS': 'NOUN', 'DEP': 'compound'}, {
            'POS': 'NOUN', 'DEP': 'compound'}, {'POS': 'NOUN'}]
        self.__basic_matcher.add(
            "JOIN", None, hyphenated_pattern, compound_pattern)
        self.__basic_matcher(doc)

        for matcher in self.__phrase_matchers:
            matcher(doc)

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
    def __filter_sents_by_entity(sents, entity_list):
        """
        Remove sentence objects from a list if they do not contain all of the provided entities.
        @param sents: List of sentence objects, with nested lists for OR conditions
        @param entity_list: List of entity label strings to search for
        @return: List of sentence objects containing all of the required entities.
        """
        output = []
        for sent in sents:
            missing_entity = False
            ents = [x.label_ for x in sent.ents]
            for ent in entity_list:
                if type(ent) == list:
                    found_match = False
                    for or_ent in ent:
                        if or_ent in ents:
                            found_match = True
                            break
                    if not found_match:
                        missing_entity = True
                        break
                elif ent not in ents:
                    missing_entity = True
                    break
            if not missing_entity:
                output.append(sent)
        return output

    @staticmethod
    def __validate_node_entities(ents, nodes):
        output = []
        used_indexes = []
        for item in ents:
            for node in nodes:
                if F"{item.lower_}<id{item.start}>" == str(node):
                    output.append((item, F"{item}<id{item.start}>"))
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
        contains_pval = False
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
        test = []
        for child in next_token.subtree:
            test.append(child)
            if child.ent_type_ == "PVAL":
                contains_pval = True
            elif child.ent_type_ == "RSID":
                contains_marker = True
        return contains_pval, contains_marker

    def allocate_contiguous_phenotypes(self, sent):
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
                if "PVAL" in token.ent_type_ and prev_token.ent_type_ == "MeSH":
                    result.append(
                        [(prev_token, F"{prev_token.lower_}<id{prev_token.i}>"),
                         (token, F"{token.lower_}<id{token.i}>")])
            prev_token = token
        return result

    def extract_phenotypes(self, doc):
        """[Extract phenotype to genotype associations from the provided SpaCy doc object]

        Args:
            doc ([SpaCy doc object]): [SpaCy processed document containing entities for extraction]

        Returns:
            [dict]: [Dictionary containing extracted phenotype, marker and p-values associated together based on SDP calculation.]
        """
        phenotype_sents = Interpreter.__filter_sents_by_entity(
            doc.sents, [["MeSH", "HPO"], ["PVAL", "PVAL-G"], "RSID"])
        results = self.calculate_sdp(phenotype_sents)
        return results

    def calculate_sdp(self, phenotype_sents):
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
            for token in sent:
                # TODO: Re-label duplicates with some easily but uniquely identifiable id.
                for child in token.children:
                    token_text = F"{token.lower_}<id{token.i}>"
                    child_text = F"{child.lower_}<id{child.i}>"
                    edges.append(('{0}'.format(token_text),
                                  '{0}'.format(child_text)))

            graph = nx.Graph(edges)

            phenotypes = Interpreter.__validate_node_entities(
                [x for x in sent.ents if x.label_ == 'MeSH'], graph.nodes)
            markers = Interpreter.__validate_node_entities(
                [x for x in sent.ents if x.label_ == 'RSID'], graph.nodes)
            pvals = Interpreter.__validate_node_entities(
                [x for x in sent.ents if x.label_ == 'PVAL'], graph.nodes)

            phenotype_count = len(phenotypes)
            marker_count = len(markers)
            pval_count = len(pvals)
            immediate_relations = self.allocate_contiguous_phenotypes(sent)

            for phenotype in phenotypes:
                # if phenotype.lower_ == 'olfactory receptors': #DEBUGGING ONLY
                #     Interpreter.display_structure(sent)
                contains_marker = False
                contains_pval = False
                subtree_count = 0
                for child in phenotype[0][0].subtree:
                    if child.ent_type_ == "PVAL":
                        contains_pval = True
                    elif child.ent_type_ == "RSID":
                        contains_marker = True
                    subtree_count += 1

                if not contains_marker and not contains_pval:
                    contains_pval, contains_marker = Interpreter.__expand_sentence_dependency_search(
                        phenotype[0][0])
                if not contains_pval and not contains_marker:
                    continue

                # if phenotype[0].lower_ == "hematocrit":
                #     Interpreter.display_structure(sent)

                # Maximum length of dependency path for association.
                marker_distance = 4
                pval_distance = 4
                rsid = None
                pval = None
                best_marker_distance = marker_distance + 1
                best_pval_distance = pval_distance + 1
                for marker in markers:  # Check shortest dependency path between each marker and each phenotype & RSID.
                    temp_distance = nx.shortest_path_length(
                        graph, source=phenotype[1].lower(), target=marker[1].lower())
                    if temp_distance <= marker_distance and temp_distance < best_marker_distance:
                        marker_distance = temp_distance
                        rsid = nx.shortest_path(
                            graph, source=phenotype[1].lower(), target=marker[1].lower())[-1]
                        best_marker_distance = temp_distance
                    else:
                        continue
                if not rsid:  # RSID must be present for an association to be made.
                    continue
                for pvalue in pvals:
                    temp_distance = nx.shortest_path_length(
                        graph, source=rsid, target=pvalue[1].lower())
                    if temp_distance <= pval_distance and temp_distance < best_pval_distance:
                        pval_distance = temp_distance
                        pval = nx.shortest_path(
                            graph, source=rsid, target=pvalue[1].lower())[-1]
                        best_pval_distance = temp_distance
                    else:
                        continue
                if not pval:
                    continue
                result_marker = Marker(rsid[:rsid.find("<id")])
                result_marker.phenotype = phenotype[0].lower_
                result_marker.misc_p_val = pval[:pval.find("<id")]
                results.append(result_marker)
                # if phenotype[0].lower_ not in results:
                #     results[phenotype[0].lower_] = []
                # results[phenotype[0].lower_].append([rsid[:rsid.find("<id")]])
                # results[phenotype[0].lower_].append([pval[:pval.find("<id")]])
        return results

    @staticmethod
    def count_entities(doc, label):
        """[Counts the number of named entities matching the provided label string]

        Args:
            doc ([SpaCy doc object]): SpaCy processed document containing entities to count]
            label ([string]): [Label of entities that will be counted e.g. 'MeSH'.]

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
            colors = {"MESH": "rgb(217, 36, 115)", "EFO": "rgb(217, 36, 115)", "HPO": "rgb(117, 36, 215)",
                      "RSID": "rgb(215, 36, 42)",
                      "PVAL": "rgb(72, 225, 21)", "PTYPE": "rgb(21, 72, 225)", "marker": "rgb(0, 225, 174)"}
        options = {"colors": colors}
        if markup_only:
            return displacy.render(doc, style="ent", page=True, options=options, jupyter=False)
        displacy.serve(doc, style="ent", options=options)

    def onto_match(self, doc):
        doc = self.process_corpus(doc, ontology_only=True)
        return [(x.text, x.label_) for x in doc.ents]

    @staticmethod
    def get_phenotype_stats(doc, lexicon):

        results = {}
        for ent in doc.ents:
            if ent.label_ in ["MeSH", "HPO"]:
                if ent.lower_ in results:
                    results[ent.lower_]["Count"] += 1
                else:
                    results[ent.lower_] = {}
                    results[ent.lower_]["Count"] = 1
                    results[ent.lower_]["Ontology"] = ent.label_
                    #results[ent.lower_]["ID"]

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
        split_sent = declaration_match.group(0).lower().replace(
            " )", ")").replace("( ", "(").split(" ")
        found_counter = 0
        found_indexes = []
        i = len(split_sent) - 2  # Indexing + ignore the actual abbreviation
        #  Moving backwards from the abbreviation, count each word in the sentence matching the first character.
        while i >= 0:
            if split_sent[i][0:1].lower() == first_char.lower():
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
                    input_text = input_text.replace(change[0], F" {change[1]} ")
        # Interpreter.__logger.info(changes) #  Can error due to strange encodings used.
        return input_text

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
        # if len(doc) < 5:  # Abbreviations are more likely to be less than 5 characters long, to avoid noise.
        #     # Prepare the document for processing and try a first pass with the NLPre library.
        #     new_doc = doc.upper()
        #     fulltext = fulltext.upper()
        #     fulltext = fulltext.replace(doc, new_doc)
        #     abbrevs = identify_parenthetical_phrases()(fulltext)
        #     result = Interpreter.insert_phrase(abbrevs, new_doc)
        #     if result is None:
        #         result = Interpreter.__check_single_word_abbrev(fulltext, new_doc)
        #     return result
        # else:
        #     return doc
