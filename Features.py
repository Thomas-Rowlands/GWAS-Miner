class Extraction:
    from sklearn.feature_extraction.text import CountVectorizer
    from sklearn.feature_extraction.text import TfidfTransformer
    from scipy.sparse import coo_matrix

    def __init__(self, stop_words, corpus):
        self.cv = self.CountVectorizer(max_df=0.8,stop_words=stop_words, max_features=10000, ngram_range=(1,3))
        self.X = self.cv.fit_transform(corpus)
        self.tfidf_transformer = self.TfidfTransformer(smooth_idf=True, use_idf=True)
        self.tfidf_transformer.fit(self.X)
        self.feature_names = self.cv.get_feature_names()

    # Most frequently occurring words
    def get_top_n_words(self, corpus, n=None):
        vec = self.CountVectorizer().fit(corpus)
        bag_of_words = vec.transform(corpus)
        sum_words = bag_of_words.sum(axis=0)
        words_freq = [(word, sum_words[0, idx]) for word, idx in vec.vocabulary_.items()]
        words_freq = sorted(words_freq, key=lambda x: x[1], reverse=True)
        return words_freq[:n]

    # Most frequently occurring Bi-grams
    def get_top_n2_words(self, corpus, n=None):
        vec1 = self.CountVectorizer(ngram_range=(2, 2), max_features=2000).fit(corpus)
        bag_of_words = vec1.transform(corpus)
        sum_words = bag_of_words.sum(axis=0)
        words_freq = [(word, sum_words[0, idx]) for word, idx in vec1.vocabulary_.items()]
        words_freq = sorted(words_freq, key=lambda x: x[1], reverse=True)
        return words_freq[:n]

    # Most frequently occurring Tri-grams
    def get_top_n3_words(self, corpus, n=None):
        vec1 = self.CountVectorizer(ngram_range=(3, 3), max_features=2000).fit(corpus)
        bag_of_words = vec1.transform(corpus)
        sum_words = bag_of_words.sum(axis=0)
        words_freq = [(word, sum_words[0, idx]) for word, idx in vec1.vocabulary_.items()]
        words_freq = sorted(words_freq, key=lambda x: x[1], reverse=True)
        return words_freq[:n]

    def sort_coo(self, coo_matrix):
        tuples = zip(coo_matrix.col, coo_matrix.data)
        return sorted(tuples, key=lambda x: (x[1], x[0]), reverse=True)

    def extract_topn_from_vector(self, feature_names, sorted_items, topn=10):
        # feature names and tf-idf score of top n items

        # use only topn items from vector
        sorted_items = sorted_items[:topn]

        score_vals = []
        feature_vals = []

        # word index and corresponding tf-idf score
        for idx, score in sorted_items:
            # keep track of feature name and its corresponding score
            score_vals.append(round(score, 3))
            feature_vals.append(feature_names[idx])

        # create a tuple of feature,score
        # results = zip(feature_vals, score_vals)
        results = {}
        for idx in range(len(feature_vals)):
            results[feature_vals[idx]] = score_vals[idx]

        return results

    def process_corpus(self):
        results = []
        for i in range(len(self.corpus)):
            doc = self.corpus[i]

            # generate tf-idf for the given document
            tf_idf_vector = self.tfidf_transformer.transform(self.cv.transform([doc]))

            sorted_items = self.sort_coo(tf_idf_vector.tocoo())

            # extract only the top n; n here is 10
            keywords = self.extract_topn_from_vector(self.feature_names, sorted_items, 10)

            # now store the results
            results.append([doc, [keywords]])
        return results

