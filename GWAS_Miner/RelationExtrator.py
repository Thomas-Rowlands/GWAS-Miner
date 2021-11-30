import re

import kindred
import os
import GWASMiner


# Overrides kindred's sentence generator method for single sentence documents.
class CustomParser(kindred.Parser):
    def _sentencesGenerator(self, text):
        parsed = self.nlp(text)
        sentence = None
        for token in parsed:
            if sentence is None:
                if not sentence is None:
                    yield sentence
                sentence = []
            sentence.append(token)

        if not sentence is None and len(sentence) > 0:
            yield sentence


if __name__ == '__main__':
    print("Setting up output directory")
    annotatedDir = os.path.join("training_output", 'annotated_relations')
    if not os.path.isdir(annotatedDir):
        os.makedirs(annotatedDir)
    unannotatedDir = os.path.join("training_output", 'missing_relations')
    if not os.path.isdir(unannotatedDir):
        os.makedirs(unannotatedDir)

    print("Loading and parsing corpus:")

    corpus = kindred.Corpus()
    with open("training_input/training_input.txt") as f:
        temp = f.read().split("\n")
        for line in temp:
            ents = []
            for match in re.finditer(pattern=r"(<!TRAIT:[0-9a-zA-Z ]{1,}!>)", string=line):
                e = kindred.Entity("TRAIT", str(line[match.start():match.end()]), [match.span()], sourceEntityID="TRAIT")
                # doc.addEntity(e)
                ents.append(e)
            for match in re.finditer(pattern=r"(<!RSID:[0-9a-zA-Z ]{1,}!>)", string=line):
                e = kindred.Entity("RSID", str(line[match.start():match.end()]), [match.span()], sourceEntityID="RSID")
                # doc.addEntity(e)
                ents.append(e)
            for match in re.finditer(pattern=r"(<!SIGNIFICANCE:[0-9a-zA-Z ]{1,}!>)", string=line):
                e = kindred.Entity("SIGNIFICANCE", str(line[match.start():match.end()]), [match.span()], sourceEntityID="SIGNIFICANCE")
                # doc.addEntity(e)
                ents.append(e)
            corpus.addDocument(kindred.Document(text=line, entities=ents))

    parser = CustomParser(model="en_core_sci_lg")
    # parser = kindred.Parser(model="en_core_sci_lg")
    parser.parse(corpus)

    print("Finding all candidate relations")
    candidateBuilder = kindred.CandidateBuilder(entityCount=2,
                                                acceptedEntityTypes=[("RSID", "TRAIT")])
    candidateRelations = candidateBuilder.build(corpus)

    print("Time to through some of the candidate relations and annotate some...")
    annotatedCorpus, unannotatedCorpus = kindred.manuallyAnnotate(corpus, candidateRelations)

    print("\nSaving annotated corpus of %d sentences (with relations that you have just annotated)" % len(
        annotatedCorpus.documents))
    kindred.save(annotatedCorpus, 'standoff', annotatedDir)

    print("Saving unannotated corpus of %d sentences (which you did not review)" % len(unannotatedCorpus.documents))
    kindred.save(unannotatedCorpus, 'standoff', unannotatedDir)
