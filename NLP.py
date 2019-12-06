import re
from pprint import pprint
import spacy
from nlpre import dedash, titlecaps, separate_reference, unidecoder
from spacy import displacy
from spacy.matcher.matcher import Matcher
from spacy.matcher.phrasematcher import PhraseMatcher
from spacy.tokens.span import Span


class PhenoNLP:
    failed_matches = []
    nlp = spacy.load("en_core_web_sm")
    rsid_regex = [{"TEXT": {"REGEX": "(?:rs[0-9]{1,}){1}"}}]
    p_value_regex = r"((\(?\bp[ -=<]{1,}(val{1,}[ue]{0,})?[ <≥=×xX-]{0,}[ \(]?\d+[\.]?[\d]{0,}[-^*() \d×xX]{0,}))"
    p_type_regex = r"(\(?GEE\)?)|(\(?FBAT\)?)"
    SNP_regex = [{"TEXT": {"REGEX": r"([ATCG]{1}[a-z]{1,}[0-9]{1,}[ATCG]{1}[a-z]{1,})"}}]
    gene_seq_regex = [{"TEXT": {"REGEX": "([ ][ACTG]{3,}[ ])"}}]
    basic_matcher = None
    phrase_matcher = None

    def __init__(self, tagging_data):
        self.__add_matchers(tagging_data)

    def __add_matchers(self, tagging_data):
        self.basic_matcher = Matcher(self.nlp.vocab)
        self.basic_matcher.add('RSID', self.__on_match, self.rsid_regex)
        self.basic_matcher.add('SNP', self.__on_match, self.SNP_regex)
        self.phrase_matcher = PhraseMatcher(self.nlp.vocab, attr="LOWER")
        for entry in tagging_data.keys():
            patterns = list(self.nlp.tokenizer.pipe(tagging_data[entry]))
            self.phrase_matcher.add(entry, self.__on_match, *patterns)

    def __on_match(self, matcher, doc, i, matches):
        match_id, start, end = matches[i]
        entity = Span(doc, start, end, label=self.nlp.vocab.strings[match_id])
        try:
            doc.ents += (entity,)
        except:
            self.failed_matches.append(entity.text)

    def __regex_match(self, pattern, doc, label):
        for match in re.finditer(pattern, doc.text, flags=re.IGNORECASE):
            start, end = match.span()
            span = doc.char_span(start, end, label=label)
            if span is not None:
                doc.ents += (span,)

    def process_corpus(self, corpus):
        # Clean corpus with NLPre parsers
        parsers = [dedash(), titlecaps(), separate_reference(), unidecoder()]
        for parser in parsers:
            corpus = parser(corpus)

        doc = self.nlp(corpus)
        old_ents, doc.ents = doc.ents, []
        self.__regex_match(self.p_value_regex, doc, "PVAL")
        self.__regex_match(self.p_type_regex, doc, "PTYPE")
        self.basic_matcher(doc)
        self.phrase_matcher(doc)

        for ent in old_ents:
            try:
                doc.ents += (ent,)
            except:
                continue

        colors = {"MESH": "rgb(247, 66, 145)", "HP": "rgb(147, 66, 245)", "RSID": "rgb(245, 66, 72)",
                  "PVAL": "rgb(102, 255, 51)", "PTYPE": "rgb(51, 102, 255)", "SNP": "rgb(0, 255, 204)"}
        options = {"colors": colors}
        # sentence_spans = list(doc.sents)
        # displacy.serve(sentence_spans, style="dep")
        pprint(self.failed_matches)
        displacy.serve(doc, style="ent", options=options)
