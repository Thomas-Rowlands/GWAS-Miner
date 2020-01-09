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
    __rsid_regex = [{"TEXT": {"REGEX": "(?:rs[0-9]{1,}){1}"}}]
    __p_value_regex = r"((\(?\bp[ -=<]{1,}(val{1,}[ue]{0,})?[ <≥=×xX-]{0,}[ \(]?\d+[\.]?[\d]{0,}[-^*() \d×xX]{0,}))"
    __p_type_regex = r"(\(?GEE\)?)|(\(?FBAT\)?)"
    __SNP_regex = [{"TEXT": {"REGEX": r"([ATCG]{1}[a-z]{1,}[0-9]{1,}[ATCG]{1}[a-z]{1,})"}}]
    __gene_seq_regex = [{"TEXT": {"REGEX": "([ ][ACTG]{3,}[ ])"}}]
    __table_ref_regex = r"(table[- ]{0,}\d{1,})"
    __basic_matcher = None
    __phrase_matchers = []

    def __init__(self, lexicon):
        self.__add_matchers(lexicon)

    def __add_matchers(self, lexicon):
        self.__basic_matcher = Matcher(self.__nlp.vocab)
        self.__basic_matcher.add('RSID', self.__on_match, self.__rsid_regex)
        self.__basic_matcher.add('SNP', self.__on_match, self.__SNP_regex)
        for entry in lexicon.keys():
            for sub_entry in sorted(lexicon[entry].keys(), reverse=True):
                new_matcher = PhraseMatcher(self.__nlp.vocab, attr="LOWER")
                patterns = list(self.__nlp.tokenizer.pipe(lexicon[entry][sub_entry]))
                new_matcher.add(entry, self.__on_match, *patterns)
                self.__phrase_matchers.append(new_matcher)

    def __on_match(self, matcher, doc, i, matches):
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

    def display_structure(self, doc):
        doc = self.__process_corpus(doc)
        sentence_spans = [x for x in list(doc.sents) if "TABLE" in [i.label_ for i in x.ents]]
        displacy.serve(sentence_spans, style="dep")

    def display_ents(self, doc):
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
    def check_single_word_abbrev(fulltext, token):
        result = ""
        first_char = token[0:1]
        first_char_count = token.count(first_char)
        search_regex = r"([ \-'\na-zA-Z0-9]+\n?\({0}\))".format(re.escape(token))
        declaration_match = re.search(search_regex, fulltext, re.IGNORECASE | re.MULTILINE)
        if declaration_match is None:
            return None
        split_sent = declaration_match.group(0).lower().split(" ")
        found_counter = 0
        found_indexes = []
        i = len(split_sent) - 2  # Indexing + ignore the actual abbreviation
        while i >= 0:
            if split_sent[i][0:1].lower() == first_char.lower():
                found_counter += 1
                found_indexes.append(i)
            i -= 1
        if first_char_count <= found_counter:
            found_indexes.sort()
            for x in split_sent[found_indexes[(first_char_count - found_counter)]:-1]:
                result += x + " "
        if result[0:1] == " ":
            result = result[1:]
        return result[:-1]

    @staticmethod
    def replace_abbreviations(doc, fulltext):
        doc = doc.lstrip()
        doc = doc.rstrip()
        if len(doc) < 5:
            new_doc = doc.upper()
            fulltext = fulltext.upper()
            fulltext = fulltext.replace(doc, new_doc)
            abbrevs = identify_parenthetical_phrases()(fulltext)
            result = Interpreter.insert_phrase(abbrevs, new_doc)
            if result is None:
                result = Interpreter.check_single_word_abbrev(fulltext, new_doc)
            return result
        else:
            return doc
