class BioCNode:
    refid = None
    role = None

    def __init__(self, refid=None, role=None):
        self.refid = refid
        self.role = role

    def jsonable(self):
        return self.__dict__


class BioCRelation:
    id = None
    infons = None
    nodes = None

    def __init__(self, id=None, infons=None, nodes=None):
        self.id = id
        self.infons = infons
        self.nodes = nodes

    def jsonable(self):
        return self.__dict__


class BioCLocation:
    offset = None
    length = None

    def __init__(self, offset=None, length=None):
        self.offset = offset
        self.length = length

    def jsonable(self):
        return self.__dict__


class BioCAnnotation:

    def __init__(self, id=None, infons=None, locations=None, text=None, length=None):
        self.id = id
        self.infons = infons
        self.text = text
        self.locations = locations

    def jsonable(self):
        return self.__dict__


class BioCSentence:
    infons = None
    offset = None
    text = None
    annotations = None
    relations = None

    def __init__(self, infons=None, offset=None, text=None, annotations=None, relations=None):
        self.infons = infons
        self.offset = offset
        self.text = text
        self.annotations = annotations
        self.relations = relations

    def jsonable(self):
        return self.__dict__


class BioCPassage:
    infons = None
    offset = None
    text = None
    sentences = None
    annotations = None
    relations = None

    def __init__(self, infons=None, offset=None, text=None, sentences=None, annotations=None, relations=None):
        self.infons = infons
        self.offset = offset
        self.sentences = sentences
        self.annotations = annotations
        self.relations = relations

    def jsonable(self):
        return self.__dict__


class BioCDocument:
    id = None
    infons = None
    passages = None

    def __init__(self, id=None, infons=None, passages=None):
        self.id = id
        self.infons = infons
        self.passages = passages

    def jsonable(self):
        return self.__dict__


class BioCCollection:
    source = None
    date = None
    key = None
    infons = None
    documents = None

    def __init__(self, source=None, date=None, key=None, infons=None, documents=None):
        self.source = source
        self.date = date
        self.key = key
        self.infons = infons
        self.documents = documents

    def jsonable(self):
        return self.__dict__


def ComplexHandler(Obj):
    if hasattr(Obj, 'jsonable'):
        return Obj.jsonable()
    else:
        raise TypeError('Object of type %s with value of %s is not JSON serializable' % (type(Obj), repr(Obj)))
