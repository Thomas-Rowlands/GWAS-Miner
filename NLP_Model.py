# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 08:47:22 2019

@author: Thomas Rowlands

Based heavily on the tutorial by Sowmya Vivek
Found here: https://medium.com/analytics-vidhya/automated-keyword-extraction-from-articles-using-nlp-bfd864f41b34

"""
passw="PASSWORD_HERE"
from nltk.corpus import wordnet
from nltk.book import *
import nltk
#from itertools import chain
#from rdflib import Graph

#g = Graph()
#g.parse("C:\\mesh2019.nt", format="nt")
#
#query = """ PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
#PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
#PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
#PREFIX owl: <http://www.w3.org/2002/07/owl#>
#PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
#PREFIX mesh: <http://id.nlm.nih.gov/mesh/>
#PREFIX mesh2019: <http://id.nlm.nih.gov/mesh/2019/>
#PREFIX mesh2018: <http://id.nlm.nih.gov/mesh/2018/>
#PREFIX mesh2017: <http://id.nlm.nih.gov/mesh/2017/>
#
#SELECT DISTINCT ?x ?y
#FROM <http://id.nlm.nih.gov/mesh>
#WHERE { ?x a meshv:SCR_Disease . ?x rdfs:label ?y .}
#ORDER BY ?y
#"""
#result = g.query(query)




import mysql.connector
import pandas
import numpy as np
import sys
from nltk.stem.porter import PorterStemmer
from nltk.stem.wordnet import WordNetLemmatizer
import re
import nltk
from nltk.corpus import stopwords
from nltk.tokenize import RegexpTokenizer
from string import punctuation


host = "localhost"
user = "root"
db = "ontology_partial"
con = mysql.connector.connect(host=host, user=user, passwd=passw, database=db)
#
#Retrieve MeSH terms
query = """SELECT meshID, termID, termName
	FROM mesh_term
    WHERE meshID in 
		(SELECT meshID
			FROM mesh_termtree
			WHERE LEFT(treeID,1) = 'C');
                        """
cursor = con.cursor()
cursor.execute(query)

testMesh = {"meshID":[], "termID":[], "termName":[]}
meshTerms = []
for (meshID, termID, termName) in cursor:
    testMesh["meshID"].append(meshID)
    testMesh["termName"].append(termName)
    testMesh["termID"].append(termID)
    meshTerms.append(termName)
    
#Retrieve HPO terms
query = """SELECT hpoID, name
        	FROM hpo_term"""

cursor.execute(query)

testHPO = {"hpoID":[], "name":[]}
hpoTerms = []
for (hpoID, name) in cursor:
    testHPO["hpoID"].append(hpoID)
    testHPO["name"].append(name)
    hpoTerms.append(name)


#Retrieve HPO synonyms
query = """SELECT hpoID, synonymText
            FROM hpo_synonym
                        """
cursor.execute(query)

testHPOSyn = {"hpoID":[], "synText":[]}
hpoSynonyms = []
for (hpoID, synText) in cursor:
    testHPOSyn["hpoID"].append(hpoID)
    testHPOSyn["synText"].append(synText)
    hpoSynonyms.append(synText)

#Retrieve HPO to Mesh mappings
query = """SELECT hpoID, meshID
        	FROM hpo2mesh
                        """
cursor.execute(query)

hpo2mesh = {"hpoID":[], "meshID":[]}
for (hpoID, meshID) in cursor:
    hpo2mesh["hpoID"].append(hpoID) 
    hpo2mesh["meshID"].append(meshID)

##Get Abstracts from dataset
query = "SELECT identifier, StudyAbstract FROM study_partial.study"
cursor.execute(query)
#add header row
tmp = []

for (identifier, studyAbstract) in cursor:
    tmp.append([identifier, studyAbstract])
    
cursor.close()
con.close()
dataset = np.array(tmp)
dataset = pandas.DataFrame(data=dataset,columns=["id", "abstract"])

dataset['word_count'] = dataset['abstract'].apply(lambda x:len(str(x).split(" ")))

print(dataset.word_count.describe())

freq = pandas.Series(''.join(dataset["abstract"]).split()).value_counts()[:20]

print("Common words: \n")
print(freq)

inFreq = pandas.Series(''.join(dataset["abstract"]).split()).value_counts()[-20:]

print("Uncommon words: \n")
print(inFreq)

lem = WordNetLemmatizer()
stem = PorterStemmer()
#List of stop words
stop_words = set(stopwords.words('english') + list(punctuation))
#Custom stop words
#new_words = set(test)
#terms = new_words
#new_words = ["using", "show", "result", "large", "also",
#             "iv", "one", "two", "new", "previously", "shown"]
#stop_words = stop_words.union(new_words)


corpus = []
for i in range(len(dataset.index)):
    #Remove punctuations
    text = re.sub('[^a-zA-Z]', ' ', dataset['abstract'][i])
    #Convert to lowercase
    text = text.lower()
    #remove tags
    text = re.sub("&lt;/?.*?&gt;"," &lt;&gt; ", text)
    #remove special characters and digits
    text = re.sub("(\\d|\\W)+", " ", text)
    #Convert to list from string
    text = text.split()
    #Stemming
    ps=PorterStemmer()
    #Lemmatisation
    lem = WordNetLemmatizer()
    text = [lem.lemmatize(word) for word in text if not word in stop_words]
    text = " ".join(text)
    corpus.append(text)

#word cloud
from os import path
from PIL import Image
from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator

import matplotlib.pyplot as plt


#wordcloud = WordCloud(background_color='white',
#                      stopwords=stop_words,
#                      max_words=100,
#                      max_font_size=50,
#                      random_state=42).generate(str(corpus))
#
#
#fig = plt.figure(1)
#plt.imshow(wordcloud)
#plt.axis('off')
##plt.show()
#fig.savefig("test.png", dpi=900)

#word count vector creation
from sklearn.feature_extraction.text import CountVectorizer

cv = CountVectorizer(max_df=0.8, stop_words=stop_words, max_features=10000, ngram_range=(1,3))
X = cv.fit_transform(corpus)

#Most frequently occuring words
def get_top_n_words(corpus, n=None):
    vec = CountVectorizer().fit(corpus)
    bag_of_words = vec.transform(corpus)
    sum_words = bag_of_words.sum(axis=0)
    words_freq = [(word, sum_words[0, idx]) for word, idx in vec.vocabulary_.items()]
    words_freq = sorted(words_freq, key=lambda x: x[1], reverse=True)
    return words_freq[:n]

#Convert most freq words to dataframe for plotting bar plot
top_words = get_top_n_words(corpus, n=20)
top_df = pandas.DataFrame(top_words, columns=["Word","Freq"])
    
#Barplot of most freq words
import seaborn as sns

g = sns.barplot(x="Word", y="Freq", data=top_df)
g.set_xticklabels(g.get_xticklabels(), rotation=25)
plt.savefig("freq_words.png", dpi=900)

#Most frequently occuring Bi-grams
def get_top_n2_words(corpus, n=None):
    vec1 = CountVectorizer(ngram_range=(2,2), max_features=2000).fit(corpus)
    bag_of_words = vec1.transform(corpus)
    sum_words = bag_of_words.sum(axis=0)
    words_freq = [(word, sum_words[0, idx]) for word, idx in vec1.vocabulary_.items()]
    words_freq = sorted(words_freq, key = lambda x: x[1], reverse=True)
    return words_freq[:n]

top2_words = get_top_n2_words(corpus, n=20)
top2_df = pandas.DataFrame(top2_words)
top2_df.columns = ["Bi-gram", "Freq"]
print(top2_df)

#Barplot of most freq Bi-grams
h = sns.barplot(x="Bi-gram", y="Freq", data=top2_df)
h.set_xticklabels(h.get_xticklabels(), rotation=45)
plt.savefig("bi-grams.png", dpi=900)


#Most frequently occuring Tri-grams
def get_top_n3_words(corpus, n=None):
    vec1 = CountVectorizer(ngram_range=(3,3), max_features=2000).fit(corpus)
    bag_of_words = vec1.transform(corpus)
    sum_words = bag_of_words.sum(axis=0)
    words_freq = [(word, sum_words[0, idx]) for word, idx in vec1.vocabulary_.items()]
    words_freq = sorted(words_freq, key = lambda x: x[1], reverse=True)
    return words_freq[:n]

top3_words = get_top_n3_words(corpus, n=20)
top3_df = pandas.DataFrame(top3_words, columns=["Tri-gram", "Freq"])
print(top3_df)

#Barplot of most freq Tri-grams
j = sns.barplot(x="Tri-gram", y="Freq", data=top3_df)
j.set_xticklabels(j.get_xticklabels(), rotation=45)
plt.savefig("tri-grams.png", dpi=900)

from sklearn.feature_extraction.text import TfidfTransformer

tfidf_transformer = TfidfTransformer(smooth_idf=True, use_idf=True)
tfidf_transformer.fit(X)

#get feature names
feature_names = cv.get_feature_names()



#Function for sorting tf_idf in descending order
from scipy.sparse import coo_matrix

def sort_coo(coo_matrix):
    tuples = zip(coo_matrix.col, coo_matrix.data)
    return sorted(tuples, key=lambda x: (x[1], x[0]), reverse=True)

def extract_topn_from_vector(feature_names, sorted_items, topn=10):
    #feature names and tf-idf score of top n items
    
    #use only topn items from vector
    sorted_items = sorted_items[:topn]
    
    score_vals = []
    feature_vals = []
    
    #word index and corresponding tf-idf score
    for idx, score in sorted_items:
        #keep track of feature name and its corresponding score
        score_vals.append(round(score, 3))
        feature_vals.append(feature_names[idx])
    
    #create a tuple of feature,score
    #results = zip(feature_vals, score_vals)
    results = {}
    for idx in range(len(feature_vals)):
        results[feature_vals[idx]] = score_vals[idx]
    
    return results

#fetch document for which keywords needs to be extracted
db = "study_keywords"
con = mysql.connector.connect(host=host, user=user, passwd=passw, database=db)
cursor = con.cursor()    
for i in range(len(corpus)):
    doc=corpus[i]

    #generate tf-idf for the given document
    tf_idf_vector = tfidf_transformer.transform(cv.transform([doc]))

    sorted_items = sort_coo(tf_idf_vector.tocoo())

    #extract only the top n; n here is 10
    keywords = extract_topn_from_vector(feature_names, sorted_items, 10)
    
    #now print the results
    print("\nAbstract:")
    print(doc)
    print("\nKeywords:")
    for k in keywords:
        print(k, keywords[k])
    sys.exit()
    keys = ""
    #for k in keywords:
    #    for t in terms:
    #        if k == t["meshName"] or k == t["synonym"]:
    #            keys = keys + k + ","
    #if keys == "":
    #    keys = "NONE."
    query = "INSERT INTO study_keywords_ft (keywords, study_identifier) VALUES ('"+keys[:-1]+"', '"+dataset["id"][i]+"')"
    cursor.execute(query)
    con.commit()


    
cursor.close()
con.close()


