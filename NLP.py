import re
from pprint import pprint
import spacy
from nlpre import dedash, titlecaps, separate_reference, unidecoder, identify_parenthetical_phrases, replace_acronyms
from spacy import displacy
from spacy.matcher.matcher import Matcher
from spacy.matcher.phrasematcher import PhraseMatcher
from spacy.tokens.span import Span
import string


class Interpreter:
    __failed_matches = []
    __nlp = spacy.load("en_core_web_sm")
    __nlp.tokenizer.add_special_case(",", [{"ORTH": ","}])
    __rsid_regex = [{"TEXT": {"REGEX": "(?:rs[0-9]{1,}){1}"}}]
    __p_value_regex = r"((\(?\bp[ -=<]{1,}(val{1,}[ue]{0,})?[ <≥=×xX-]{0,}[ \(]?\d+[\.]?[\d]{0,}[-^*() \d×xX]{0,}))"
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
                doc.ents += (span,)

    def __process_corpus(self, corpus, ontology_only=False):
        # Clean corpus with NLPre parsers
        parsers = [dedash(), titlecaps(), separate_reference(), unidecoder()]
        for parser in parsers:
            corpus = parser(corpus)

        doc = self.__nlp(corpus)
        old_ents, doc.ents = doc.ents, []
        #  Additional regex matches unnecessary when limited to ontology entities.
        if not ontology_only:
            self.__regex_match(self.__table_ref_regex, doc, "TABLE")
            self.__regex_match(self.__p_value_regex, doc, "PVAL")
            self.__regex_match(self.__p_type_regex, doc, "PTYPE")
            self.__basic_matcher(doc)
        for matcher in self.__phrase_matchers:
            matcher(doc)
        if not ontology_only:
            for ent in old_ents:
                try:
                    doc.ents += (ent,)
                except:  # Default SpaCy entities should never override others.
                    continue

        return doc

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
    def check_lists(doc):
        results = []
        temp_word = ""
        temp_word_two = ""
        for word in doc:
            if word.dep_ == "amod":
                temp_word = Interpreter.merge_reliant(doc, word)
                temp_word_two = Interpreter.merge_reliant(doc, word.head)
                results.append({"original": str(word), "index": word.idx, "l_words": [temp_word['l_words'], temp_word_two['l_words']], "replacement": F"{temp_word['replacement']} {temp_word_two['replacement']}"})
        new_text = doc.text
        test = []
        for result in results:
            test.append(F"{new_text[:result['l_words'][0]]} {result['replacement']}.")
            #TODO: Remove the additional terms which stack up on every subsequent list addition.
        #new_text = F"{new_text[:temp_word_two['index']]} {new_text[temp_word_two['index'] + len(temp_word_two['original']):]}"
        #print(new_text)
        pprint(test)

    def display_structure(self, doc):
        """
        Start running the Displacy visualization of the tokenized sentences identified by the NLP pipeline.
        @param doc: The document to process for tokenized sentence structures.
        """
        doc = doc.replace(",", "").replace(" and", "")

        hyphenated_pattern = [{'POS': 'PROPN'}, {'IS_PUNCT': True}, {'POS': 'VERB'}]
        compound_pattern = [{'POS': 'NOUN', 'DEP': 'compound'}, {'POS': 'NOUN', 'DEP': 'compound'}, {'POS': 'NOUN'}]

        matcher = Matcher(self.__nlp.vocab)
        matcher.add("JOIN", None, hyphenated_pattern, compound_pattern)
        doc = self.__nlp(doc)
        matches = matcher(doc)
        spans = []
        for match_id, start, end in matches:
            string_id = self.__nlp.vocab.strings[match_id]  # Get string representation
            span = doc[start:end]  # The matched span
            spans.append(span)

        Interpreter.check_lists(doc)

        sentence_spans = [x for x in list(doc.sents)]
        options = {"compact": True}
        # pprint([doc[match[1]:match[2]] for match in matches])
        displacy.serve(sentence_spans, style="dep", options=options)

    def display_ents(self, doc):
        """
        Start running the Displacy visualization of the named entities recognised by the NLP pipeline.
        @param doc: The document to process for named entities.
        """
        doc = self.__process_corpus(doc)
        colors = {"MESH": "rgb(247, 66, 145)", "HP": "rgb(147, 66, 245)", "RSID": "rgb(245, 66, 72)",
                  "PVAL": "rgb(102, 255, 51)", "PTYPE": "rgb(51, 102, 255)", "SNP": "rgb(0, 255, 204)"}
        options = {"colors": colors}
        pprint(self.__failed_matches)
        displacy.serve(doc, style="ent", options=options)

    def onto_match(self, doc):
        doc = self.__process_corpus(doc, ontology_only=True)
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
        search_regex = r"([ \-'\na-zA-Z0-9]+\n?\({0}\))".format(re.escape(token))
        declaration_match = re.search(search_regex, fulltext, re.IGNORECASE | re.MULTILINE)
        if declaration_match is None:  # No declaration found for token
            return None
        # First match SHOULD be the declaration of this abbreviation.
        # Split REGEX match to a list of words
        split_sent = declaration_match.group(0).lower().split(" ")
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
        if result[0:1] == " ":
            result = result[1:]
        return result[:-1]

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
