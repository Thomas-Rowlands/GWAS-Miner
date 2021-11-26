import kindred
import os
import GWASMiner

if __name__ == '__main__':
    print("Setting up output directory")
    annotatedDir = os.path.join("training_output", 'annotated_relations')
    if not os.path.isdir(annotatedDir):
        os.makedirs(annotatedDir)
    unannotatedDir = os.path.join("training_output", 'missing_relations')
    if not os.path.isdir(unannotatedDir):
        os.makedirs(unannotatedDir)

    print("Loading and parsing corpus:")

    with open(args.corpus) as f:
        corpus = kindred.Corpus(f.read())

    parser = kindred.Parser()
    parser.parse(corpus)

    print(
        "Splitting the corpus into sentences so that we can save any annotated sentences and don't need to annotate it all")
    sentenceCorpus = corpus.splitIntoSentences()

    print("Loading wordlists:")
    lexicon = GWASMiner.__prepare_ontology_data().get_lexicon_by_name("MESH")
    wordlistDict = {}
    for wordlist in args.wordlists.split(','):
        assert os.path.isfile(wordlist), 'Unable to access file: %s' % wordlist
        entityType = os.path.splitext(os.path.basename(wordlist))[0]
        wordlistDict[entityType] = wordlist
        print("  %s - %s" % (entityType, wordlist))

    assert len(wordlistDict) == 2, "This annotation tool currently only handles two entity relations of different types"

    wordlistLookup = kindred.EntityRecognizer.loadWordlists(wordlistDict, idColumn=0, termsColumn=0)

    print("Annotating entities in corpus with wordlists")
    entityRecognizer = kindred.EntityRecognizer(wordlistLookup)
    entityRecognizer.annotate(sentenceCorpus)

    print("Finding all candidate relations")
    acceptedEntityTypes = wordlistDict
    candidateBuilder = kindred.CandidateBuilder(entityCount=len(wordlistDict),
                                                acceptedEntityTypes=[tuple(sorted(wordlistDict.keys()))])
    candidateRelations = candidateBuilder.build(sentenceCorpus)

    print("Time to through some of the candidate relations and annotate some...")
    annotatedCorpus, unannotatedCorpus = kindred.manuallyAnnotate(sentenceCorpus, candidateRelations)

    print("\nSaving annotated corpus of %d sentences (with relations that you have just annotated)" % len(
        annotatedCorpus.documents))
    kindred.save(annotatedCorpus, 'standoff', annotatedDir)

    print("Saving unannotated corpus of %d sentences (which you did not review)" % len(unannotatedCorpus.documents))
    kindred.save(unannotatedCorpus, 'standoff', unannotatedDir)