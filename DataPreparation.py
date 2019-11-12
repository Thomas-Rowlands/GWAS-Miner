import re
import sys

from lxml import etree

from Study import Study
from nltk.stem.porter import PorterStemmer
from nltk.stem.wordnet import WordNetLemmatizer
from nltk.corpus import stopwords
import nltk.tag
import nltk.data
import nltk.chunk
nltk.download('stopwords')
nltk.download('wordnet')
nltk.download('punkt')


class BackoffTagger:
    def __init__(self):
        self._taggers = [nltk.data.load('taggers/maxent_treebank_pos_tagger/PY3/english.pickle')]


class PreProcessing:
    import pandas
    from string import punctuation
    import re

    def __init__(self, study_data, tagging_data):
        """
        Preprocesses free text data ready for analyses.
        :param study_data: Unprocessed corpus text (numpy array).
        :param tagging_data: Ontology terms to use for POS tagging
        """
        self.study_data = self.pandas.DataFrame(data=study_data, columns=["content"])
        self.study_data["word_count"] = self.study_data["content"].apply(lambda x: len(str(x).split(" ")))
        self.common_words = self.pandas.Series(''.join(self.study_data["content"]).split()).value_counts()[:20]
        self.uncommon_words = self.pandas.Series(''.join(self.study_data["content"]).split()).value_counts()[-20:]
        self.stop_words = nltk.corpus.stopwords
        self.stem = nltk.PorterStemmer()  # Not currently used but may be soon
        self.lem = nltk.WordNetLemmatizer()
        self.stop_words = set(nltk.corpus.stopwords.words('english') + list(self.punctuation))
        self.corpus = self.prepare_corpus()
        self.default_tagger = BackoffTagger()
        self.tagging_data = tagging_data
        self.phenotype_tagger = self.build_semantic_taggers()
        self.tagged_text = self.tag_text()

    def add_stop_words(self, new_words):
        """
        Add a new set of stop words
        :param new_words: List of stop words to append to the current set.
        :return:
        """
        self.stop_words = self.stop_words.union(new_words)

    def prepare_corpus(self):
        """
        Processes the current corpus text, cleaning & preparing the text.
        :return: Pre-processed corpus text.
        """
        corpus = []
        for i in range(len(self.study_data.index)):
            # Remove punctuations
            text = self.re.sub('[^a-zA-Z]', ' ', self.study_data['content'][i])
            # Convert to lowercase
            text = text.lower()
            # remove tags
            text = self.re.sub("&lt;/?.*?&gt;", " &lt;&gt; ", text)
            # remove special characters and digits <---------Separate processing of numeric data?
            text = self.re.sub("(\\d|\\W)+", " ", text)
            # Convert to list from string
            text = text.split()
            text = [self.lem.lemmatize(word) for word in text] #if word not in self.stop_words]
            text = " ".join(text)
            corpus.append(text)
        return corpus

    def build_semantic_taggers(self):
        """
        Uses the MeSH & HPO terms to create a new tagger that backs off to the NLTK tagger.
        :return:
        """
        uni_tagger = nltk.UnigramTagger(model=self.tagging_data, backoff=self.default_tagger)
        bi_tagger = nltk.BigramTagger(model=self.tagging_data, backoff=uni_tagger)
        tri_tagger = nltk.TrigramTagger(model=self.tagging_data, backoff=bi_tagger)
        n4_tagger = nltk.NgramTagger(model=self.tagging_data, backoff=tri_tagger, n=4)
        n5_tagger = nltk.NgramTagger(model=self.tagging_data, backoff=n4_tagger, n=5)
        n6_tagger = nltk.NgramTagger(model=self.tagging_data, backoff=n5_tagger, n=6)
        n7_tagger = nltk.NgramTagger(model=self.tagging_data, backoff=n6_tagger, n=7)
        n8_tagger = nltk.NgramTagger(model=self.tagging_data, backoff=n7_tagger, n=8)
        n9_tagger = nltk.NgramTagger(model=self.tagging_data, backoff=n8_tagger, n=9)
        return n9_tagger

    def tag_text(self):
        """
        Attaches POS tags to each word in the corpus
        :return: Tagged version of the corpus
        """
        text = ""
        for i in self.corpus:
            text = text + i
        text = nltk.word_tokenize(text)
        return self.phenotype_tagger.tag(text)

    def chunk_text(self):
        pattern = "<DT>?<JJ>*<NN.*>"
       # regexp_pattern = tag_pattern2re_pattern(pattern)



    # Retrieves rs identifiers from input string
    @staticmethod
    def get_rs_identifiers(text):
        result = re.findall("(?:rs[0-9]{1,}){1}", text)
        if not result:
            return None
        else:
            return result

    def strip_xml(self, pmid, xml):
        study = Study()
        tree = etree.fromstring(xml, etree.get_default_parser())
        study.pmid = pmid

        #  Abstract

        tmp = tree.xpath("//abstract//*/text()")
        abstract = ""
        for i in tmp:
            abstract = abstract + " " + i
        study.abstract = abstract

        #  Article Title

        study.title = tree.xpath("//article-title[1]/text()")[0]

        #  Authors

        study.authors = []
        for i in range(int(tree.xpath("count(//contrib[@contrib-type='author'])"))):
            author = [tree.xpath("//contrib[@contrib-type='author'][" + str(i + 1) + "]//text()")[0],
                      tree.xpath("//contrib[@contrib-type='author'][" + str(i + 1) + "]//text()")[1],
                      tree.xpath("//contrib[@contrib-type='author'][" + str(i + 1) + "]//text()")[-1:][0]]
            study.authors.append(author)

        #  Sections

        sections = []
        for i in range(int(tree.xpath("count(//body//sec)"))):  # Iterate through each section title
            content = ""
            section_text = tree.xpath("//body/sec[" + str(i + 1) + "]//*[not(self::title)]/text()")
            section_title = tree.xpath("//body/sec[" + str(i + 1) + "]/title/text()")
            for tmp in section_title:
                content = content + tmp
            section_title = content
            content = ""
            for tmp in section_text:
                content = content + tmp
            section_text = content
            section = [section_title,
                       section_text]
            sections.append(section)
        study.sections = sections

        #  SNP Extraction

        snips = []
        titleSNP = PreProcessing.get_rs_identifiers(study.title)
        if titleSNP:
            for i in range(len(titleSNP)):
                snips.append(titleSNP[i])
        abstractSNP = PreProcessing.get_rs_identifiers(study.abstract)
        if abstractSNP:
            for i in range(len(abstractSNP)):
                snips.append(abstractSNP[i])
        for (title, body) in study.sections:
            results = PreProcessing.get_rs_identifiers(body)
            if results:
                for i in range(len(results)):
                    snips.append(results[i])
        study.snps = snips

        #  Back sections

        acknowledgements = ""
        content = tree.xpath("//back/ack//*[not(self::title)]/text()")
        for tmp in content:
            acknowledgements = acknowledgements + tmp
        study.acknowledgements = acknowledgements

        #  Citations pending

        return study
