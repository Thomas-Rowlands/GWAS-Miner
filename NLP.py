import itertools
import re
from pprint import pprint
import spacy
from nlpre import dedash, titlecaps, separate_reference, unidecoder, identify_parenthetical_phrases, replace_acronyms
from spacy import displacy
from spacy.matcher.matcher import Matcher
from spacy.matcher.phrasematcher import PhraseMatcher
from spacy.tokens.span import Span
from Utility_Functions import Utility
import networkx as nx
import string


class Interpreter:
    __failed_matches = []
    __nlp = spacy.load("en_core_web_md", disable=["ner"])
    __nlp.tokenizer.add_special_case(",", [{"ORTH": ","}])
    __rsid_regex = [{"TEXT": {"REGEX": "(?:rs[0-9]{1,}){1}"}}]
    __p_value_regex = r"((\(?\bp[ -=<]{1,}(val{1,}[ue]{0,})?[ <≥=×xX-]{0,}[ \(]?\d+[\.]?[\d]{0,}[-^*() \d×xX]{0,}))"
    __p_value_regex_inline = r"(\d?\.?\d[ ]?[*×xX]{1}[ ]?\d{1,}[ ]?-\d{1,})"
    __p_type_regex = r"(\(?GEE\)?)|(\(?FBAT\)?)"
    __SNP_regex = [{"TEXT": {"REGEX": r"([ATCG]{1}[a-z]{1,}[0-9]{1,}[ATCG]{1}[a-z]{1,})"}}]
    __gene_seq_regex = [{"TEXT": {"REGEX": "([ ][ACTG]{3,}[ ])"}}]
    __table_ref_regex = r"(table[- ]{0,}\d{1,})"
    __basic_matcher = None
    __phrase_matchers = []

    def __init__(self, lexicon, ontology_only=False):
        if not ontology_only:
            self.__add_matchers(lexicon)

    def __add_matchers(self, lexicon):
        self.__basic_matcher = Matcher(self.__nlp.vocab, validate=True)
        self.__basic_matcher.add('RSID', self.__on_match, self.__rsid_regex)
        self.__basic_matcher.add('SNP', self.__on_match, self.__SNP_regex)
        for entry in lexicon.keys():
            for sub_entry in sorted(lexicon[entry].keys(), reverse=True):
                new_matcher = PhraseMatcher(self.__nlp.vocab, attr="LOWER")
                patterns = list(self.__nlp.tokenizer.pipe(lexicon[entry][sub_entry]))
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
        entity = Span(doc, start, end, label=self.__nlp.vocab.strings[match_id])
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
        # Clean corpus with NLPre parsers
        parsers = [dedash(), titlecaps(), separate_reference(), unidecoder()]
        for parser in parsers:
            corpus = parser(corpus)

        doc = self.__nlp(corpus)

        old_ents, doc.ents = doc.ents, []
        #  Additional regex matches unnecessary when limited to ontology entities.
        if not ontology_only:
            self.__regex_match(self.__table_ref_regex, doc, "TABLE")
            self.__regex_match(self.__p_value_regex_inline, doc, "PVAL")
            self.__regex_match(self.__p_value_regex, doc, "PVAL-G")
            self.__regex_match(self.__p_type_regex, doc, "PTYPE")

        hyphenated_pattern = [{'POS': 'PROPN'}, {'IS_PUNCT': True, 'LOWER': '-'}, {'POS': 'VERB'}]
        compound_pattern = [{'POS': 'NOUN', 'DEP': 'compound'}, {'POS': 'NOUN', 'DEP': 'compound'}, {'POS': 'NOUN'}]
        self.__basic_matcher.add("JOIN", None, hyphenated_pattern, compound_pattern)
        self.__basic_matcher(doc)

        new_doc = ""
        # for i in range(len(sentence_spans)):
        #     change = Interpreter.replace_list(sentence_spans[i])
        #     if change:
        #         for statement in change['statements']:
        #             new_doc = F"{new_doc} {statement}"
        #     else:
        #         new_doc = F"{new_doc} {sentence_spans[i]}"
        for matcher in self.__phrase_matchers:
            matcher(doc)
        if not ontology_only:
            for ent in old_ents:
                try:
                    doc.ents += (ent,)
                except:  # Default SpaCy entities should never override others.
                    continue
        self.__merge_spans(doc, "PVAL")
        self.__merge_spans(doc, "MeSH")
        return doc

    @staticmethod
    def __merge_spans(doc, entity_label):
        with doc.retokenize() as retokenizer:
            for ent in [x for x in doc.ents if x.label_ == entity_label]:
                retokenizer.merge(doc[ent.start:ent.end])

    @staticmethod
    def __filter_sents_by_entity(sents, entity_list):
        """
        Remove sentence objects from a list if they do not contain all of the provided entities.
        @param sents: List of sentence objects
        @param entity_list: List of entity label strings to search for
        @return: List of sentence objects containing all of the required entities.
        """
        output = []
        for sent in sents:
            missing_entity = False
            for ent in entity_list:
                ents = [x.label_ for x in sent.ents]
                if ent not in ents:
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
                if item.lower_ == str(node) or item.lower_ == str(node):
                    duplicate_node_indexes = Utility.retrieve_value_indexes(str(node), list(nodes))
                    for idx in duplicate_node_indexes:
                        current_idx = None
                        try:
                            current_idx = output.index((item, idx))
                        except(ValueError):
                            current_idx = None
                        if current_idx:
                            continue
                        else:
                            output.append((item, list(nodes).index(node)))
                            break
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
        contains_snp = False
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
                contains_snp = True
        return contains_pval, contains_snp

    def extract_phenotypes(self, doc):
        output = []
        phenotype_sents = Interpreter.__filter_sents_by_entity(doc.sents, ["MeSH", "PVAL", "RSID"])
        new_phenotype_sents = []
        # Split any sentences containing duplicate phenotypes up.
        # for sent in phenotype_sents:
        #     duplicate_found = False
        #     for token in sent:
        #         duplicate_count = len([x for x in sent if x.lower_ == token.lower_ and token.lower_ in [i.lower_ for i in sent.ents]])
        #         if duplicate_count > 1:
        #             duplicate_found = True
        #             length = len(token.lower_)
        #             for dupe_token in [x for x in sent if x.lower_ == token.lower_]:
        #                 start = dupe_token.sent_start
        #                 new_sent = F"{sent.text[:start]} <dupe> {sent.text[start + length:]}"
        #                 new_sent = self.process_corpus(new_sent)
        #                 new_phenotype_sents.append([dupe_token, new_sent])
        #     if not duplicate_found:
        #         new_phenotype_sents.append(sent)


        results = {}
        # Iterate through each sentence containing a phenotype named entity label
        for sent in phenotype_sents:
            edges = []
            for token in sent:
                for child in token.children:# TODO: Re-label duplicates with some easily but uniquely identifiable id.
                    edges.append(('{0}'.format(token.lower_), '{0}'.format(child.lower_)))
            graph = nx.Graph(edges)

            phenotypes = Interpreter.__validate_node_entities([x for x in sent.ents if x.label_ == 'MeSH'], graph.nodes)
            snps = Interpreter.__validate_node_entities([x for x in sent.ents if x.label_ == 'RSID'], graph.nodes)
            pvals = Interpreter.__validate_node_entities([x for x in sent.ents if x.label_ == 'PVAL'], graph.nodes)

            combinations = [phenotypes, snps, pvals]

            phenotype_count = len(phenotypes)
            snp_count = len(snps)
            pval_count = len(pvals)

            #combinations = list(itertools.product(*combinations))

            for phenotype in phenotypes:
                # if phenotype.lower_ == 'olfactory receptors':
                #     Interpreter.display_structure(sent)
                contains_snp = False
                contains_pval = False
                subtree_count = 0
                for child in phenotype[0].subtree:
                    if child.ent_type_ == "PVAL":
                        contains_pval = True
                    elif child.ent_type_ == "RSID":
                        contains_snp = True
                    subtree_count += 1
                if not contains_snp and not contains_pval:
                    contains_pval, contains_snp = Interpreter.__expand_sentence_dependency_search(phenotype[0])
                if not contains_pval and not contains_snp:
                    continue
                #if phenotype.lower_ == "olfactory receptor":
                    # test = Interpreter.__validate_phenotype_context(phenotype)
                    # Interpreter.display_structure(sent)
                # if not Interpreter.__validate_phenotype_context(phenotype):
                #     continue
                snp_distance = 4
                rsid = None
                pval_distance = 4
                pval = None
                for snp in snps:
                    temp_distance = nx.shortest_path_length(graph, source=phenotype[0].lower_, target=snp[0].lower_)
                    if temp_distance <= snp_distance:
                        snp_distance = temp_distance
                        rsid = nx.shortest_path(graph, source=phenotype[0].lower_, target=snp[0].lower_)[-1]
                    else:
                        continue
                if not rsid:
                    continue
                for pvalue in pvals:
                    temp_distance = nx.shortest_path_length(graph, source=rsid, target=pvalue[0].lower_)
                    if temp_distance <= pval_distance:
                        pval_distance = temp_distance
                        pval = nx.shortest_path(graph, source=rsid, target=pvalue[0].lower_)[-1]
                    else:
                        continue
                if not pval:
                    continue
                results[phenotype[0]] = [rsid]
                results[phenotype[0]].append(pval)
            # for (pheno, snp, pval) in combinations:
            # results.append({"Phenotype": pheno, "SNP": snp, "PVAL": pval,
            #                 "pheno>snp": {"Distance": nx.shortest_path_length(graph, source=pheno, target=snp),
            #                               "Path": nx.shortest_path(graph, source=pheno, target=snp)},
            #                 "snp>pval": {"Distance": nx.shortest_path_length(graph, source=snp, target=pval),
            #                              "Path": nx.shortest_path(graph, source=snp, target=pval)}
            #                 }
            #                )

        return results


    @staticmethod
    def count_entities(doc, label):
        count = 0
        for ent in doc.ents:
            if ent.label_ == label:
                print(ent.lower_)
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
        result = {"original": str(word), "index": word.i, "l_words": [], "replacement": str(word)}
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
            output["statements"].append(F"{new_text[:results[0]['l_words'][0]]}{result['replacement']}.")
            output['items'].append(result['replacement'])
        return output

    @staticmethod
    def display_structure(sentence_spans):
        """
        Start running the Displacy visualization of the tokenized sentences identified by the NLP pipeline.
        @param doc: The NLP processed document.
        """
        options = {"compact": True}
        # pprint([doc[match[1]:match[2]] for match in matches])
        displacy.serve(sentence_spans, style="dep", options=options)

    @staticmethod
    def display_ents(doc):
        """
        Start running the Displacy visualization of the named entities recognised by the NLP pipeline.
        @param doc: The NLP processed document.
        """
        colors = {"MESH": "rgb(247, 66, 145)", "HP": "rgb(147, 66, 245)", "RSID": "rgb(245, 66, 72)",
                  "PVAL": "rgb(102, 255, 51)", "PTYPE": "rgb(51, 102, 255)", "SNP": "rgb(0, 255, 204)"}
        options = {"colors": colors}
        displacy.serve(doc, style="ent", options=options)

    def onto_match(self, doc):
        doc = self.process_corpus(doc, ontology_only=True)
        return [(x.text, x.label_) for x in doc.ents]

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
        first_char_count = token.count(first_char)  # Get the number of occurrences of this letter in abbreviation.
        # Locate abbreviations in parenthesis
        search_regex = r"([ \-'\n\w0-9]+\n?\({0}[^a-zA-Z]{0,}\))".replace("{0}", re.escape(token))  # r"([ \-'\na-zA-Z0-9]+\n?\({0}\))".format(re.escape(token))
        declaration_match = re.search(search_regex, fulltext, re.IGNORECASE | re.MULTILINE)
        if declaration_match is None:  # No declaration found for token
            return None
        # First match SHOULD be the declaration of this abbreviation.
        # Split REGEX match to a list of words
        split_sent = declaration_match.group(0).lower().replace(" )", ")").replace("( ", "(").split(" ")
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
    def replace_all_abbreviations(fulltext):
        """
        Returns the expanded form of an abbreviation from the fulltext.
        @param fulltext: The document containing the abbreviations and their declarations
        @return: String containing the expanded version of the abbreviation.
        """
        changes = []
        pattern = r"([^ \"',.(-]\b)?([a-z]{0,})([A-Z]{2,})([a-z]{0,})(\b[^;,.'\" )-]?)"#r"([^(-]\b[a-z]{0,}[A-Z]{2,}[a-z]{0,}\b[^)-])"
        for match in re.findall(pattern, fulltext):
            target = ""
            if type(match) == str:
                target = match
            else:
                target = match[2]
            target = target.strip()
            if target not in [x for [x, y] in changes]:
                changes.append([target, Interpreter.replace_abbreviations(target, fulltext)])
        for change in changes:
            if change[1]:
                fulltext = fulltext.replace(change[0], F" {change[1]} ")
        pprint(changes)
        return fulltext

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
        if len(doc) < 5:  # Abbreviations are more likely to be less than 5 characters long, to avoid noise.
            # Prepare the document for processing and try a first pass with the NLPre library.
            new_doc = doc.upper()
            fulltext = fulltext.upper()
            fulltext = fulltext.replace(doc, new_doc)
            abbrevs = identify_parenthetical_phrases()(fulltext)
            result = Interpreter.insert_phrase(abbrevs, new_doc)
            if result is None:
                result = Interpreter.__check_single_word_abbrev(fulltext, new_doc)
            return result
        else:
            return doc
