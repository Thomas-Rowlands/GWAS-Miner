import json
import re
from datetime import datetime

from Exceptions import TableTypeError
import BioC

table_significance_pattern = r""


def output_tables(destination, tables):
    try:
        with open(destination, "w", encoding="utf-8") as fout:
            json.dump(tables, fout, default=ComplexHandler, ensure_ascii=False, indent=4)
    except IOError as ie:
        print(F"IO Error occurred: {ie}")
    except Exception as e:
        print(F"An unknown error occurred: {e}")


def process_tables(nlp, tables):
    t, m, p, r = 0, 0, 0, 0
    for table in tables:
        nlp = table.add_spacy_docs(nlp)
        table.annotations = nlp.annotations
        table.relations = nlp.relations
        nlp.annotations = []
        nlp.relations = []
        if table.table_type:
            t, m, p, r = table.assign_annotations(t, m, p, r)
    return tables, nlp


def get_cell_entity_annotation(nlp, ent, table_element, cell_id, element_offset=0):
    current_datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    offset = ent.start_char + element_offset
    length = ent.end_char - ent.start_char
    loc = BioC.BioCLocation(offset=offset, length=length, table_cell_id=cell_id, table_element=table_element)
    entity_type = None
    entity_identifier = None
    if ent.label_ == "RSID":
        entity_type = "genetic_variant"
        entity_identifier = F"dbSNP:{ent.text}"
    elif ent.label_ == "PVAL":
        entity_type = "significance"
        entity_identifier = "PVAL"
    else:
        entity_type = "trait"
        entity_identifier = F"MeSH:{ent.label_}"
    entity_id = None
    if entity_type == "trait":
        entity_id = F"T{nlp.t}"
    elif entity_type == "genetic_variant":
        entity_id = F"V{nlp.v}"
    else:
        entity_id = F"S{nlp.s}"
    annotation = BioC.BioCAnnotation(id=entity_id, infons={"type": entity_type, "identifier": entity_identifier,
                                                           "annotator": "GWASMiner@le.ac.uk",
                                                           "updated_at": current_datetime},
                                     locations=[loc], text=ent.text)
    return annotation, nlp


def parse_tables(file_input, nlp):
    tables = []
    tables_data = None
    contains_annotations = False
    try:
        with open(file_input, "r", encoding="utf-8") as fin:
            tables_data = json.load(fin)
        for doc in tables_data["documents"]:  # each doc is a table
            table_id = doc['id']
            title = None
            title_offset = None
            caption = None
            caption_offset = None
            content_offset = None
            columns = []
            table_sections = []
            footer = None
            footer_offset = None
            doc['annotations'] = []
            for passage in doc['passages']:
                if "section_title_2" in passage['infons'].keys():
                    print(F"Second table title found in: {file_input}")
                if passage['infons']['section_title_1'] == 'table_title':
                    title = passage['text']
                    title_offset = passage['offset']
                elif passage['infons']['section_title_1'] == 'table_caption':
                    caption = passage['text']
                    caption_offset = passage['offset']
                elif passage['infons']['section_title_1'] == 'table_footer':
                    footer = passage['text']
                    footer_offset = passage['offset']
                elif passage['infons']['section_title_1'] == 'table_content':
                    content_offset = passage['offset']
                    col_row = TableRow()
                    for col in passage['column_headings']:
                        column_cell = TableCell(col['cell_id'], str(col['cell_text']))
                        col_row.cells.append(column_cell)
                    columns.append(col_row)
                    for section in passage['data_section']:
                        table_section = TableSection(title=section['table_section_title_1'])
                        for row in section['data_rows']:
                            data_row = TableRow()
                            for cell in row:
                                data_cell = TableCell(cell['cell_id'], str(cell['cell_text']))
                                data_row.cells.append(data_cell)
                            if abs(len(data_row.cells) - len(columns[0].cells)) > 1:
                                print(F"{file_input} contains potentially problematic cell counts!")
                            table_section.rows.append(data_row)
                        table_sections.append(table_section)

            table = Table(title, table_id, title_offset, content_offset, columns,
                          table_sections, caption, caption_offset, footer, footer_offset)
            tables.append(table)
    except IOError as ie:
        print(F"No tables file found for {file_input.replace('_tables.json', '')}")
    if tables_data:
        annotated_tables, nlp = process_tables(nlp, tables)
        for table in annotated_tables:
            if table.annotations:
                contains_annotations = True
                for bioc_table in tables_data["documents"]:
                    for passage in bioc_table["passages"]:
                        passage["relations"] = []
                    if bioc_table["id"] == table.table_id:
                        bioc_table["annotations"] = table.annotations
                        bioc_table["relations"] = table.relations
                        break

    # if contains_annotations:
    #     print(F"Table(s) have annotation(s) in: {file_input}")
    # else:
    #     print(F"No table annotations found in: {file_input}")
    # if not contains_annotations:
    #     print(F"No table annotations found in: {file_input}")
    return tables_data, contains_annotations


class Table:
    # Table contains a list of variants and significances for a single trait
    TYPE_MARKER_LIST = 1
    # Table contains a list of traits and significances for a single variant
    TYPE_TRAIT_LIST = 2
    # Table contains a list of traits as well as markers and significances
    TYPE_TRAIT_AND_MARKER_LIST = 3
    # Table columns are on the left rather than the top.
    TYPE_REVERSED_COLUMNS = 4

    COLUMN_TRAIT = 1
    COLUMN_SIGNIFICANCE = 2
    COLUMN_MARKER = 3

    trait_strings = [r"(?:[ ]|^)(trait)(?:[^\w]|[\s]|\b)", r"(?:[ ]|^)(disease)(?:[^\w]|[\s]|\b)",
                     r"(?:[ ]|^)(phenotype)(?:[^\w]|[\s]|\b)"]
    significance_strings = [r"(?:[ ]|^)(significance)(?:[^\w]|[\s]|\b)",
                            r"(?:[ ]|^)(p-?val[ue]{0,2}s?)(?:[^\w]|[\s]|\b)", r"(?:[ ]|^)(p)(?:[^\w]|[\s]|\b)"]
    marker_strings = [r"(?:[ ]|^)(rsids?)(?:[^\w]|[\s]|\b)", r"(?:[ ]|^)(markers?)(?:[^\w]|[\s]|\b)",
                      r"(?:[ ]|^)(variants?)(?:[^\w]|[\s]|\b)", r"(?:[ ]|^)(snps?)(?:[^\w]|[\s]|\b)"]

    def __init__(self, title, table_id, title_offset, content_offset, column_cells, data_sections,
                 caption_text, caption_offset, footer_text, footer_offset):
        """
        Table object from parsed JSON input.
        """
        self.title = title
        self.table_id = table_id
        self.title_offset = title_offset
        self.title_doc = None
        self.content_offset = content_offset
        self.column_rows = column_cells
        self.caption_text = caption_text
        self.caption_offset = caption_offset
        self.caption_doc = None
        self.footer_text = footer_text
        self.footer_offset = footer_offset
        self.footer_doc = None
        self.annotations = []
        self.relations = []
        self.table_type = None
        self.has_super_rows = False
        self.column_types = []
        self.data_sections = data_sections
        self.caption_ents = []
        self.section_title_ents = []
        self.footer_ents = []
        self.title_ents = []
        self.contains_trait = False
        self.contains_marker = False
        self.contains_significance = False
        self.section_ents = []
        self.passages = []

    def __set_table_type(self):
        # Check column types
        self.contains_marker, self.contains_trait, self.contains_pval = False, False, False
        if self.title_doc:
            if self.title_doc.ents and any(x for x in self.title_doc.ents if x._.is_trait or x._.has_trait):
                self.contains_trait = True
                self.title_ents.append(self.COLUMN_TRAIT)

        if self.caption_doc:
            if self.caption_doc.ents:
                if any(x for x in self.caption_doc.ents if x._.is_trait or x._.has_trait):
                    self.contains_trait = True
                    self.caption_ents.append(self.COLUMN_TRAIT)
                if any(x for x in self.caption_doc.ents if x.label_ == "RSID"):
                    self.contains_marker = True
                    self.caption_ents.append(self.COLUMN_MARKER)

        if self.footer_doc and not self.contains_trait:
            if self.footer_doc.ents and any(x for x in self.footer_doc.ents if x._.is_trait or x._.has_trait):
                self.contains_trait = True
                self.footer_ents.append(self.COLUMN_TRAIT)

        for i in range(len(self.data_sections)):
            section = self.data_sections[i]
            self.section_ents.append([])
            if section.doc:
                self.has_super_rows = True
                if not self.contains_trait or not self.contains_marker:
                    if section.doc.ents and any(x for x in section.doc.ents if x._.is_trait or x._.has_trait):
                        self.contains_trait = True
                        self.section_ents[i].append(self.COLUMN_TRAIT)
                    if section.doc.ents and any(x for x in section.doc.ents if x.label_ == "RSID"):
                        self.contains_marker = True
                        self.section_ents[i].append(self.COLUMN_MARKER)

        for row in self.column_rows:
            for i in range(len(row.cells)):
                cell = row.cells[i]
                if cell.doc:
                    if cell.doc.ents:
                        entity = cell.doc.ents[0]
                        if entity.label_ == "RSID":
                            self.column_types.insert(i, Table.COLUMN_MARKER)
                            self.contains_marker = True
                        elif entity.label_ == "PVAL":
                            self.column_types.insert(i, Table.COLUMN_SIGNIFICANCE)
                            self.contains_significance = True
                        elif entity._.is_trait or entity._.has_trait:
                            self.column_types.insert(i, Table.COLUMN_TRAIT)
                            self.contains_trait = True
                    if len(self.column_types) != i + 1:
                        self.column_types.insert(i, 0)  # no entity found
                    cell_text = [cell.text.lower()] if not cell.text.find("|") else cell.text.lower().split(
                        "|")
                    for text in [x for x in cell_text if x]:
                        if [x for x in Table.trait_strings if re.search(x, text)]:
                            if isinstance(self.column_types[i], list):
                                self.column_types[i].append(Table.COLUMN_TRAIT)
                            else:
                                self.column_types[i] = [self.column_types[i], Table.COLUMN_TRAIT]
                            self.contains_trait = True
                        elif [x for x in Table.marker_strings if re.search(x, text)]:
                            if isinstance(self.column_types[i], list):
                                self.column_types[i].append(Table.COLUMN_MARKER)
                            else:
                                self.column_types[i] = [self.column_types[i], Table.COLUMN_MARKER]
                            self.contains_marker = True
                        elif [x for x in Table.significance_strings if re.search(x, text)] \
                                or text.lower().strip() == "p":
                            if isinstance(self.column_types[i], list):
                                self.column_types[i].append(Table.COLUMN_SIGNIFICANCE)
                            else:
                                self.column_types[i] = [self.column_types[i], Table.COLUMN_SIGNIFICANCE]
                            self.contains_significance = True
                    if len(self.column_types) != i + 1:
                        self.column_types.append(0)  # default added in case of no entities.
                else:
                    self.column_types.append(0)  # blank cell
                if isinstance(self.column_types[i], list):  # Clean list values
                    desirables = [x for x in self.column_types[i] if x != 0]
                    desirables_len = len(desirables)
                    if desirables_len == 0:
                        self.column_types[i] = 0
                    else:
                        self.column_types[i] = [x for x in self.column_types[i] if x != 0]
        if self.contains_significance and self.contains_trait and self.contains_marker:
            self.table_type = Table.TYPE_TRAIT_LIST
        elif self.contains_significance and self.contains_marker:
            self.table_type = Table.TYPE_MARKER_LIST

    def add_spacy_docs(self, nlp):
        if self.caption_text:
            self.caption_doc = nlp.process_corpus(self.caption_text, )
        if self.footer_text:
            self.footer_doc = nlp.process_corpus(self.footer_text, )
        if self.title:
            if isinstance(self.title, list):
                self.title = self.title[0]
            self.title_doc = nlp.process_corpus(self.title, )
        if self.column_rows:
            for row in self.column_rows:
                for cell in row.cells:
                    cell.add_spacy_docs(nlp)
        if self.data_sections:
            for section in self.data_sections:
                section.add_spacy_docs(nlp)
                for row in section.rows:
                    for cell in row.cells:
                        cell.add_spacy_docs(nlp)
        self.__set_table_type()
        return nlp

    def __get_cell_annotation(self, cell, col_type, t, m, p):
        annotation = None
        if self.COLUMN_TRAIT == col_type:
            if cell.doc and cell.doc.ents:
                annotation = BioC.convert_cell_to_annotation(cell.doc, t,
                                                             "table_content",
                                                             cell.id)
                t += 1

        elif self.COLUMN_MARKER == col_type:
            if cell.doc and cell.doc.ents:
                annotation = BioC.convert_cell_to_annotation(cell.doc, m,
                                                             "table_content",
                                                             cell.id)
                m += 1
        elif self.COLUMN_SIGNIFICANCE == col_type:
            if cell.doc:
                new_entity = cell.doc[:]
                new_entity.label_ = "PVAL"
                try:
                    cell.doc.ents += (new_entity,)
                    annotation = BioC.convert_cell_to_annotation(cell.doc, p,
                                                                 "table_content",
                                                                 cell.id)
                    p += 1
                except Exception as e:
                    print(e)
                    return annotation, t, m, p
        return annotation, t, m, p

    def __get_row_annotations(self, row, col_types, t, m, p, r):
        row_annotations = []
        for i in range(len(row.cells)):
            cell = row.cells[i]
            column_types = [col_types[i]] if isinstance(col_types[i], int) \
                else col_types[i]
            for col_type in set(column_types):
                annotation, t, m, p = self.__get_cell_annotation(cell, col_type, t, m, p)
                if annotation:
                    row_annotations.append(annotation)
        row_relations, r = self.__get_row_relations(row_annotations, r)
        return row_annotations, row_relations, t, m, p, r

    def __get_row_relations(self, row_annotations, r):
        row_relations = []
        row_annotations_dict = {
            "traits": [x for x in row_annotations if x.infons['type'] == "trait"],
            "pvals": [x for x in row_annotations if x.infons['type'] == "significance"],
            "variants": [x for x in row_annotations if x.infons['type'] == "genetic_variant"]
        }

        if self.table_type == self.TYPE_TRAIT_LIST:
            row_relations = self.__get_trait_list_row_relations(row_annotations_dict, r)
        elif self.table_type == self.TYPE_MARKER_LIST:
            row_relations = self.__get_marker_list_row_relations(row_annotations_dict, r)
        elif self.table_type == self.TYPE_TRAIT_AND_MARKER_LIST:
            row_relations = self.__get_trait_marker_list_row_relations(row_annotations_dict, r)

        return row_relations, r

    def __get_trait_list_row_relations(self, row_annotations_dict, r):
        row_relations = []
        trait_count, pval_count, variant_count = self.__get_row_annotation_counts(row_annotations_dict)

        return row_relations, r

    def __get_marker_list_row_relations(self, row_annotations_dict, r):
        row_relations = []
        trait_count, pval_count, variant_count = self.__get_row_annotation_counts(row_annotations_dict)
        if trait_count > 1 or pval_count > 1 or variant_count > 1:
            raise TableTypeError(self)
        trait_annot = row_annotations_dict['traits'][0]
        variant_annot = row_annotations_dict['variants'][0]
        pval_annot = row_annotations_dict['pvals'][0]
        row_relations.append(
            BioC.BioCRelation()
        )
        return row_relations, r

    def __get_trait_marker_list_row_relations(self, row_annotations_dict, r):
        row_relations = []
        trait_count, pval_count, variant_count = self.__get_row_annotation_counts(row_annotations_dict)
        return row_relations, r

    @staticmethod
    def __get_row_annotation_counts(row_annotations_dict):
        trait_count, pval_count, variant_count = len(row_annotations_dict['traits']), len(
            row_annotations_dict['pvals']), len(row_annotations_dict['variants'])
        return trait_count, pval_count, variant_count

    def assign_gc_annotations(self, t, m, p, r, used_annots, gc_relations):
        docs_iter = [(self.title_doc, self.title_offset, "table_title"),
                     (self.caption_doc, self.caption_offset, "table_caption"),
                     (self.footer_doc, self.footer_offset, "table_footer")]
        used_annots = []
        for (doc, offset, elem) in docs_iter:
            if doc and doc.ents:
                used_annots, t, m, p, r = BioC.get_bioc_annotations(doc, used_annots, offset, t, m, p, r, elem)

        for row in self.column_rows:
            row_annotations, row_relations, t, m, p, r = self.__get_gc_row_annotations(row, self.column_types, t, m, p,
                                                                                       r)
            used_annots += row_annotations

        for section in self.data_sections:
            for row in section.rows:
                row_annotations, row_relations, t, m, p, r = self.__get_gc_row_annotations(row, self.column_types, t, m,
                                                                                           p, r)
                used_annots += row_annotations
        self.annotations = used_annots
        return t, m, p, r, used_annots

    def assign_annotations(self, t, m, p, r):
        docs_iter = [(self.title_doc, self.title_offset, "table_title"),
                     (self.caption_doc, self.caption_offset, "table_caption"),
                     (self.footer_doc, self.footer_offset, "table_footer")]
        used_annots = []
        for (doc, offset, elem) in docs_iter:
            if doc and doc.ents:
                used_annots, t, m, p, r = BioC.get_bioc_annotations(doc, used_annots, offset, t, m, p, r, elem)

        for row in self.column_rows:
            row_annotations, row_relations, t, m, p, r = self.__get_row_annotations(row, self.column_types, t, m, p, r)
            used_annots += row_annotations

        for section in self.data_sections:
            for row in section.rows:
                row_annotations, row_relations, t, m, p, r = self.__get_row_annotations(row, self.column_types, t, m,
                                                                                        p, r)
                used_annots += row_annotations
        self.annotations = used_annots
        return t, m, p, r

    def jsonable(self):
        output_dict = self.__dict__
        del output_dict['title_doc']
        del output_dict['footer_doc']
        del output_dict['caption_doc']
        return output_dict


class TablePassage:
    def __init__(self, title):
        self.title = title
        self.sections = []
        self.annotations = []
        self.relations = []

    def jsonable(self):
        output_dict = self.__dict__
        del output_dict['title']
        return output_dict


class TableSection:
    def __init__(self, rows=None, title=None, title_offset=None):
        if rows is None:
            rows = []
        self.rows = rows
        self.title = str(title)
        self.title_offset = title_offset
        self.doc = None

    def add_spacy_docs(self, nlp):
        if self.title:
            self.doc = nlp.process_corpus(self.title, )

    def jsonable(self):
        output_dict = self.__dict__
        del output_dict['doc']
        return output_dict


class TableRow:
    def __init__(self, cells=None):
        if cells is None:
            cells = []
        self.cells = cells

    def jsonable(self):
        return self.__dict__


class TableCell:
    def __init__(self, id, text):
        self.id = id
        self.text = text
        self.doc = None

    def add_spacy_docs(self, nlp):
        if self.text:
            self.doc = nlp.process_corpus(self.text, )

    def jsonable(self):
        output_dict = self.__dict__
        del output_dict['doc']
        return output_dict


def ComplexHandler(Obj):
    if hasattr(Obj, 'jsonable'):
        return Obj.jsonable()
    else:
        raise TypeError('Object of type %s with value of %s is not JSON serializable' % (type(Obj), repr(Obj)))
