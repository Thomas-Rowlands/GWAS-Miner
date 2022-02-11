from datetime import datetime


def get_bioc_annotations(doc, used_annots, offset, t, m, p, r, table_elem_id=None, table_cell_id=None):
    current_datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    annotations = []
    output_annotations = []
    for ent in doc.ents:
        annotations.append({
            "entity_type": ent.label_,
            "id": ent.label_[ent.label_.index(":") + 2:] if ":" in ent.label_ else ent.label_,
            "text": ent.text,  # self.trim_brackets(ent.text),
            "offset": ent.start_char,
            "length": ent.end_char - ent.start_char,
            "table_element_id": table_elem_id,
            "table_cell_id": table_cell_id
        })
    if annotations:
        for annot in annotations:
            loc = BioCLocation(offset=annot["offset"] + offset, length=annot["length"], table_element=annot['table_element_id'], table_cell_id=annot['table_cell_id'])
            if annot["text"] in used_annots:
                for old_annot in output_annotations:
                    if old_annot.text == annot["text"] and loc not in old_annot.locations:
                        old_annot.locations.append(loc)
            if "RSID" not in annot["entity_type"] and "PVAL" not in annot["entity_type"]:
                genomic_trait = BioCAnnotation(id=F"T{t}", infons={"type": "trait", "identifier": annot["id"],
                                                                   "annotator": "tr142@le.ac.uk",
                                                                   "updated_at": current_datetime},
                                               locations=[loc], text=annot["text"])
                output_annotations.append(genomic_trait)
                t += 1
            elif "RSID" in annot["entity_type"]:
                marker_identifier = BioCAnnotation(id=F"M{m}",
                                                   infons={"type": "genomic_marker", "identifier": annot["id"],
                                                           "annotator": "tr142@le.ac.uk",
                                                           "updated_at": current_datetime},
                                                   locations=[loc], text=annot["text"])
                output_annotations.append(marker_identifier)
                m += 1
            elif "PVAL" in annot["entity_type"]:
                p_value = BioCAnnotation(id=F"P{p}", infons={"type": "significance", "identifier": annot["id"],
                                                             "annotator": "tr142@le.ac.uk",
                                                             "updated_at": current_datetime},
                                         locations=[loc], text=annot["text"])
                output_annotations.append(p_value)
                p += 1
            used_annots.append(annot['text'])
    return output_annotations, used_annots, t, m, p, r


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
    table_element = None
    table_cell_id = None

    def __init__(self, offset=None, length=None, table_element=None, table_cell_id=None):
        self.offset = offset
        self.length = length
        self.table_element = table_element
        self.table_cell_id = table_cell_id

    def jsonable(self):
        return self.__dict__


class BioCAnnotation:

    def __init__(self, id=None, infons=None, locations=None, text=None):
        self.text = text
        self.infons = infons
        self.id = id
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
