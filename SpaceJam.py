import csv
import re
import sys

import spacy
from spacy.matcher.matcher import Matcher
from spacy.tokens.span import Span
from spacy.util import filter_spans
from spacy.vocab import Vocab
from spacy.matcher import PhraseMatcher
from spacy import displacy
from nlpre import titlecaps, dedash, identify_parenthetical_phrases
from nlpre import replace_acronyms, replace_from_dictionary
# Load English tokenizer, tagger, parser, NER and word vectors
import pprint

nlp = spacy.load("en_core_web_sm")


def on_match(matcher, doc, i, matches):
    match_id, start, end = matches[i]
    entity = Span(doc, start, end, label=nlp.vocab.strings[match_id])
    doc.ents += (entity,)

def _regex_match(pattern, doc, label):
    for match in re.finditer(pattern, doc.text):
        start, end = match.span()
        span = doc.char_span(start, end, label=label)
        if span is not None:
            doc.ents += (span,)

def process_text(data, tagging_data):
    # spacy.cli.download("en_core_web_sm")
    corpus = ""
    for text in data:
        corpus = corpus + " " + text
    ABBR = identify_parenthetical_phrases()(corpus)
    parsers = [dedash(), titlecaps(), replace_acronyms(ABBR), replace_from_dictionary(prefix="MeSH_")]
    for parser in parsers:
        corpus = parser(text)

    HP_patterns = [nlp(text) for text in tagging_data]
    rsid_regex = [{"TEXT": {"REGEX": "(?:rs[0-9]{1,}){1}"}}]
    p_value_regex = r"((\(?\bp[ -=<]{1,}(val{1,}[ue]{0,})?[ <≥=×xX-]{0,}[ \(]?\d+[\.]?[\d]{0,}[-^*() \d×xX]{0,}))"
    p_type_regex = r"(\(?GEE\)?)|(\(?FBAT\)?)"
    SNP_regex = [{"TEXT": {"REGEX": r"([ATCG]{1}[a-z]{1,}[0-9]{1,}[ATCG]{1}[a-z]{1,})"}}]
    mesh_regex = [{"TEXT": {"REGEX": r"(\bMeSH_[^ ]{1,})"}}]
    print(HP_patterns)

    phrase_matcher = PhraseMatcher(nlp.vocab)
    basic_matcher = Matcher(nlp.vocab)

    phrase_matcher.add('HP', on_match, *HP_patterns)
    basic_matcher.add('MeSH', on_match, mesh_regex)
    basic_matcher.add('RSID', on_match, rsid_regex)
    basic_matcher.add('SNP', on_match, SNP_regex)

    doc = nlp(corpus)
    old_ents = doc.ents
    doc.ents = []

    _regex_match(p_value_regex, doc, "PVAL")
    _regex_match(p_type_regex, doc, "PTYPE")
    matches = basic_matcher(doc)
    matches = phrase_matcher(doc)
    # for ent in old_ents:
    #     try:
    #         doc.ents += (ent,)
    #         print("Added new ent " + ent)
    #     except:
    #         print("Error adding ent ")
    #         continue
    colors = {"MESH": "rgb(247, 66, 145)", "HP": "rgb(147, 66, 245)", "RSID": "rgb(245, 66, 72)",
              "PVAL": "rgb(102, 255, 51)", "PTYPE": "rgb(51, 102, 255)", "SNP": "rgb(0, 255, 204)"}
    options = {"colors": colors}
    #sentence_spans = list(doc.sents)
    #displacy.serve(sentence_spans, style="dep")
    displacy.serve(doc, style="ent", options=options)
