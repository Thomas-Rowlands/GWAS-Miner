class Study:
    def __init__(self, title=None, abstract=None, authors=None, snps=None, concepts=None, p_values=None, results=None,
                 pmid=None, gwas_id=None, sections=None, acknowledgements=None, citations=None, tables=None):
        self.title = title
        self.abstract = abstract
        self.authors = authors
        self.snps = snps
        self.concepts = concepts
        self.p_values = p_values
        self.results = results
        self.pmid = pmid
        self.gwas_id = gwas_id
        self.sections = sections
        self.acknowledgements = acknowledgements
        self.citations = citations
        self.tables = tables