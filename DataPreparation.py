import pprint
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
from pprint import pprint
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
        self.lem = nltk.WordNetLemmatizer()
        self.stop_words = set(nltk.corpus.stopwords.words('english') + list(self.punctuation))
        self.corpus = self.prepare_corpus()
        self.default_tagger = BackoffTagger()
        self.tagging_data = tagging_data
        self.phenotype_tagger = self.build_semantic_taggers()
        self.tagged_text = self.tag_text()
        self.chunked_text = self.chunk_text()

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
            # text = self.re.sub('[^a-zA-Z0-9]', ' ', self.study_data['content'][i])
            # Convert to lowercase
            text = self.study_data['content'][i].lower()
            # remove tags
            # text = self.re.sub("&lt;/?.*?&gt;", " &lt;&gt; ", text)
            # remove special characters
            text = self.re.sub("(\\W)+", " ", text)
            # Convert to list from string
            text = text.split()
            text = [self.lem.lemmatize(word) for word in text]  # if word not in self.stop_words]
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
        rsid_tagger = nltk.RegexpTagger(regexps=[(r'(?:rs[0-9]{1,}){1}', 'RSID')], backoff=tri_tagger)

        # facet joint arthrosis
        return rsid_tagger

    def tag_text(self):
        """
        Attaches POS tags to each word in the corpus
        :return: Tagged version of the corpus
        """
        text = ""
        for i in self.corpus:
            text = text + " " + i
        text = nltk.word_tokenize(text)
        return self.phenotype_tagger.tag(text)

    def chunk_text(self):
        pattern = r"""NP: {<DT|PP\$>?<JJ>*<NN>}
                      {<NNP>+}"""
        chonk = nltk.RegexpParser(pattern)
        result = chonk.parse(self.tagged_text)
        return result

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
        xml = re.sub("</", " </", xml)  # Ensure white space is present between nested tags
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

        #  SNP Extraction   <------REDUNDANT with new NLP Regex tags

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

        #  Tables
        tables = tree.xpath("//table")
        snps_with_values = {}
        for elem in tables:#Iterate through tables in study
            header = elem.xpath(".//thead/tr//td")
            i = 1
            fbat_p_val_col = None
            gee_p_val_col = None
            snp_col = None
            phenotype_col = None
            for cell in header:#Identify purpose of each column
                text = cell.xpath(".//text()")[0]
                p_type = re.search(r"(p{1}([- ]?val(ue)?))", text, flags=re.IGNORECASE)
                if p_type:
                    if "GEE" in p_type.string:
                        gee_p_val_col = i
                        i += 1
                        continue
                    elif "FBAT" in p_type.string:
                        fbat_p_val_col = i
                        i += 1
                        continue
                if "snp" in text.lower():
                    snp_col = i
                    i += 1
                    continue
                if "phenotype" == text.lower().replace("* ", ""):
                    phenotype_col = i
                    i += 1
                    continue
                i += 1
            if fbat_p_val_col or gee_p_val_col:#Check if table contains associated p-vals
                body = elem.xpath(".//tbody")
                rows = body[0].xpath(".//tr")
                for row in rows:
                    fbat_p_value = "".join(row.xpath(".//td["+str(fbat_p_val_col)+"]//text()"))
                    gee_p_value = "".join(row.xpath(".//td["+str(gee_p_val_col)+"]//text()"))
                    phenotype = "".join(row.xpath(".//td["+str(phenotype_col)+"]//text()"))
                    snp = "".join(row.xpath(".//td["+str(snp_col)+"]//text()"))
                    if phenotype not in snps_with_values.keys():
                        snps_with_values[phenotype] = []
                    snps_with_values[phenotype].append({"FBAT":fbat_p_value, "GEE":gee_p_value, "SNP":snp})
            else:
                continue
        pprint(snps_with_values)

        #  Citations pending

        return study
