import json
import re

import BioC

table_significance_pattern = r""


def output_tables(destination, tables):
    try:
        with open(destination, "w", encoding="utf-8") as fout:
            json.dump(tables, fout, default=ComplexHandler, ensure_ascii=False)
    except IOError as ie:
        print(F"IO Error occurred: {ie}")
    except Exception as e:
        print(F"An unknown error occurred: {e}")


def process_tables(nlp, tables):
    t, m, p, r = 0, 0, 0, 0
    used_annots = []
    for table in tables:
        table.add_spacy_docs(nlp)
        if table.table_type:
            t, m, p, r, used_annots = table.assign_annotations(t, m, p, r, used_annots)
    return tables


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
                        table_section = TableSection()
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
        annotated_tables = process_tables(nlp, tables)
        for table in annotated_tables:
            if table.annotations:
                contains_annotations = True
                break
        for table in annotated_tables:
            if table.annotations:
                for bioc_table in tables_data["documents"]:
                    if bioc_table["id"] == table.table_id:
                        for annot in table.annotations:
                            bioc_table['annotations'].append(annot)
                        break
    # if contains_annotations:
    # print(F"Table(s) have annotation(s) in: {file_input}")
    # else:
    #     print(F"No table annotations found in: {file_input}")
    if not contains_annotations:
        print(F"No table annotations found in: {file_input}")
    return tables_data, contains_annotations


class Table:
    TYPE_MARKER_LIST = 1
    TYPE_TRAIT_LIST = 2
    TYPE_TRAIT_AND_MARKER_LIST = 3
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
        self.annotations = None
        self.table_type = None
        self.column_types = []
        self.data_sections = data_sections

    def __set_table_type(self):
        # Check column types
        contains_marker, contains_trait, contains_pval = False, False, False
        if self.title_doc:
            if self.title_doc.ents and any(x for x in self.title_doc.ents if x._.is_trait or x._.has_trait):
                contains_trait = True
        if self.caption_doc:
            if self.caption_doc.ents and any(x for x in self.caption_doc.ents if x._.is_trait or x._.has_trait):
                contains_trait = True
        if self.footer_doc:
            if self.footer_doc.ents and any(x for x in self.footer_doc.ents if x._.is_trait or x._.has_trait):
                contains_trait = True

        for row in self.column_rows:
            for i in range(len(row.cells)):
                cell = row.cells[i]
                if cell.doc:
                    if cell.doc.ents:
                        entity = cell.doc.ents[0]
                        if entity.label_ == "RSID":
                            self.column_types.insert(i, Table.COLUMN_MARKER)
                            contains_marker = True
                        elif entity.label_ == "PVAL":
                            self.column_types.insert(i, Table.COLUMN_SIGNIFICANCE)
                            contains_pval = True
                        elif entity._.is_trait or entity._.has_trait:
                            self.column_types.insert(i, Table.COLUMN_TRAIT)
                            contains_trait = True
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
                            contains_trait = True
                        elif [x for x in Table.marker_strings if re.search(x, text)]:
                            if isinstance(self.column_types[i], list):
                                self.column_types[i].append(Table.COLUMN_MARKER)
                            else:
                                self.column_types[i] = [self.column_types[i], Table.COLUMN_MARKER]
                            contains_marker = True
                        elif [x for x in Table.significance_strings if re.search(x, text)] \
                                or text.lower().strip() == "p":
                            if isinstance(self.column_types[i], list):
                                self.column_types[i].append(Table.COLUMN_SIGNIFICANCE)
                            else:
                                self.column_types[i] = [self.column_types[i], Table.COLUMN_SIGNIFICANCE]
                            contains_pval = True
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
        if contains_pval and contains_trait and contains_marker:
            self.table_type = Table.TYPE_TRAIT_LIST
        elif contains_pval and contains_marker:
            self.table_type = Table.TYPE_MARKER_LIST

    def add_spacy_docs(self, nlp):
        if self.caption_text:
            self.caption_doc = nlp.process_corpus(self.caption_text)
        if self.footer_text:
            self.footer_doc = nlp.process_corpus(self.footer_text)
        if self.title:
            if isinstance(self.title, list):
                self.title = self.title[0]
            self.title_doc = nlp.process_corpus(self.title)
        if self.column_rows:
            for row in self.column_rows:
                for cell in row.cells:
                    cell.add_spacy_docs(nlp)
        if self.data_sections:
            for section in self.data_sections:
                for row in section.rows:
                    for cell in row.cells:
                        cell.add_spacy_docs(nlp)
        self.__set_table_type()

    def assign_annotations(self, t, m, p, r, used_annots):
        docs_iter = [(self.title_doc, self.title_offset, "table_title"),
                     (self.caption_doc, self.caption_offset, "table_caption"),
                     (self.footer_doc, self.footer_offset, "table_footer")]
        data_cells_iter = [(self.data_sections, self.content_offset, "table_content")]
        column_cells_iter = [(self.column_rows, self.content_offset, "table_content")]
        used_annots = []
        for (doc, offset, elem) in docs_iter:
            if doc and doc.ents:
                used_annots, t, m, p, r = BioC.get_bioc_annotations(doc, used_annots, offset, t, m, p, r, elem)
        for row in self.column_rows:
            for cell in row.cells:
                if cell.doc and cell.doc.ents:
                    used_annots, t, m, p, r = BioC.convert_cell_to_annotation(cell.doc, used_annots,
                                                                              self.content_offset, t, m,
                                                                              p, r, "table_content", cell.id)
        for section in self.data_sections:
            for row in section.rows:
                for i in range(len(row.cells)):
                    cell = row.cells[i]
                    column_types = [self.column_types[i]] if isinstance(self.column_types[i], int) \
                        else self.column_types[i]
                    for type in set(column_types):
                        if Table.COLUMN_TRAIT == type or Table.COLUMN_MARKER == type:
                            if cell.doc and cell.doc.ents:
                                used_annots, t, m, p, r = BioC.convert_cell_to_annotation(cell.doc, used_annots,
                                                                                          self.content_offset,
                                                                                          t, m, p, r, "table_content",
                                                                                          cell.id)
                        elif Table.COLUMN_SIGNIFICANCE == type:
                            if cell.doc:
                                new_entity = cell.doc[:]
                                new_entity.label_ = "PVAL"
                                try:
                                    cell.doc.ents += (new_entity,)
                                    used_annots, t, m, p, r = BioC.convert_cell_to_annotation(cell.doc, used_annots,
                                                                                              self.content_offset,
                                                                                              t, m, p, r,
                                                                                              "table_content",
                                                                                              cell.id)
                                except Exception as e:
                                    print(e)
                                    return t, m, p, r, used_annots

        self.annotations = used_annots
        return t, m, p, r, used_annots

    def jsonable(self):
        output_dict = self.__dict__
        del output_dict['title_doc']
        del output_dict['footer_doc']
        del output_dict['caption_doc']
        return output_dict


class TableSection:
    def __init__(self, rows=None):
        if rows is None:
            rows = []
        self.rows = rows

    def jsonable(self):
        return self.__dict__


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
            self.doc = nlp.process_corpus(self.text)

    def jsonable(self):
        output_dict = self.__dict__
        del output_dict['doc']
        return output_dict


def ComplexHandler(Obj):
    if hasattr(Obj, 'jsonable'):
        return Obj.jsonable()
    else:
        raise TypeError('Object of type %s with value of %s is not JSON serializable' % (type(Obj), repr(Obj)))
