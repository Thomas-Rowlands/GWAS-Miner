class PreProcessing:
    import pandas
    from nltk.stem.porter import PorterStemmer
    from nltk.stem.wordnet import WordNetLemmatizer
    from nltk.corpus import stopwords
    from string import punctuation
    import re

    def __init__(self, study_data):
        self.study_data = self.pandas.DataFrame(data=study_data, columns=["id", "abstract"])
        self.study_data["word_count"] = self.study_data["abstract"].apply(lambda x: len(str(x).split(" ")))
        self.common_words = self.pandas.Series(''.join(self.study_data["abstract"]).split()).value_counts()[:20]
        self.uncommon_words = self.pandas.Series(''.join(self.study_data["abstract"]).split()).value_counts()[-20:]
        self.stop_words = self.stopwords
        self.stem = self.PorterStemmer()#Not currently used but may be soon
        self.lem = self.WordNetLemmatizer()
        self.stop_words = set(self.stopwords.words('english') + list(self.punctuation))
        self.corpus = self.process_corpus()

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
            text = self.re.sub('[^a-zA-Z]', ' ', self.study_data['abstract'][i])
            # Convert to lowercase
            text = text.lower()
            # remove tags
            text = self.re.sub("&lt;/?.*?&gt;", " &lt;&gt; ", text)
            # remove special characters and digits <---------Separate processing of numeric data?
            text = self.re.sub("(\\d|\\W)+", " ", text)
            # Convert to list from string
            text = text.split()
            text = [self.lem.lemmatize(word) for word in text if not word in self.stop_words]
            text = " ".join(text)
            corpus.append(text)
        return corpus
