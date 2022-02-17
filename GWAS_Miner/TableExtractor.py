import json

import BioC


def output_tables(destination, tables):
    try:
        with open(destination, "w", encoding="utf-8") as fout:
            json.dump(tables, fout, default=ComplexHandler)
    except IOError as ie:
        print(F"IO Error occurred: {ie}")
    except Exception as e:
        print(F"An unknown error occurred: {e}")


def process_tables(nlp, tables):
    t, m, p, r = 0, 0, 0, 0
    used_annots = []
    for table in tables:
        table.add_spacy_docs(nlp)
        t, m, p, r, used_annots = table.assign_annotations(t, m, p, r, used_annots)
    return tables


def parse_tables(file_input, nlp):
    tables = []
    tables_data = None
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
            rows = []
            footer = None
            footer_offset = None
            doc['annotations'] = []
            for passage in doc['passages']:
                if passage['infons']['section_title_1'] == 'table_title': # TODO: Could be more titles, specify in annotations which one.
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
                    for col in passage['column_headings']:
                        column_cell = TableCell(col['cell_id'], str(col['cell_text']))
                        columns.append(column_cell)
                    for row in passage['data_section'][0]['data_rows']:
                        data_row = []
                        for cell in row:
                            data_cell = TableCell(cell['cell_id'], str(cell['cell_text']))
                            data_row.append(data_cell)
                        if len(data_row) != len(columns):
                            print(F"{file_input} contains potentially problematic cell counts!")
                        rows.append(data_row)

            table = Table(title, table_id, title_offset, content_offset, columns,
                          rows, caption, caption_offset, footer, footer_offset)
            tables.append(table)
    except IOError as ie:
        print(F"No tables file found for {file_input.replace('_tables.json', '')}")
    if tables_data:
        annotated_tables = process_tables(nlp, tables)
        for table in annotated_tables:
            if table.annotations:
                # print(file_input)
                for bioc_table in tables_data["documents"]:
                    if bioc_table["id"] == table.table_id:
                        for annot in table.annotations:
                            bioc_table['annotations'].append(annot)
                        break
    return tables_data


class Table:
    TYPE_MARKER_LIST = 1
    TYPE_TRAIT_LIST = 2
    COLUMN_TRAIT = 1
    COLUMN_SIGNIFICANCE = 2
    COLUMN_MARKER = 3
    trait_strings = ["trait", "disease", "phenotype"]
    significance_strings = ["p-value", "significance", "p value", "p"]
    marker_strings = ["rsid", "marker", "variant", "snp"]

    def __init__(self, title, table_id, title_offset, content_offset, column_cells, data_cells,
                 caption_text, caption_offset, footer_text, footer_offset):
        """
        Table object from parsed JSON input.
        """
        self.title = title
        self.table_id = table_id
        self.title_offset = title_offset
        self.title_doc = None
        self.content_offset = content_offset
        self.column_cells = column_cells
        self.data_cells = data_cells
        self.caption_text = caption_text
        self.caption_offset = caption_offset
        self.caption_doc = None
        self.footer_text = footer_text
        self.footer_offset = footer_offset
        self.footer_doc = None
        self.annotations = None
        self.table_type = None
        self.column_types = []

    def __set_table_type(self):
        # Check column types
        contains_marker, contains_trait, contains_pval = False, False, False
        if self.title_doc:
            if self.title_doc.ents and any(x for x in self.title_doc.ents if x._.is_trait):
                contains_trait = True
        if self.caption_doc:
            if self.caption_doc.ents and any(x for x in self.caption_doc.ents if x._.is_trait):
                contains_trait = True

        if self.footer_doc:
            if self.footer_doc.ents and any(x for x in self.footer_doc.ents if x._.is_trait):
                contains_trait = True

        for cell in self.column_cells:
            if cell.doc:
                if cell.doc.ents:
                    entity = cell.doc.ents[0]
                    if entity._.is_trait:
                        self.column_types.append(Table.COLUMN_TRAIT)
                        contains_trait = True
                    elif entity.label_ == "RSID":
                        self.column_types.append(Table.COLUMN_MARKER)
                        contains_marker = True
                    elif entity.label_ == "PVAL":
                        self.column_types.append(Table.COLUMN_SIGNIFICANCE)
                        contains_pval = True
                    else:
                        self.column_types.append(0)  # no entity found
                else:
                    cell_text = [cell.text.lower()] if not cell.text.find("||") else cell.text.lower().split("||")
                    for text in cell_text:
                        if [x for x in Table.trait_strings if x == cell_text]:
                            self.column_types.append(Table.COLUMN_TRAIT)
                            contains_trait = True
                        elif [x for x in Table.marker_strings if x == cell_text]:
                            self.column_types.append(Table.COLUMN_MARKER)
                            contains_marker = True
                        elif [x for x in Table.significance_strings if x == cell_text]:
                            self.column_types.append(Table.COLUMN_SIGNIFICANCE)
                            contains_pval = True
                        else:
                            self.column_types.append(0)  # no entity found
            else:
                self.column_types.append(0)  # blank cell
        if contains_pval and contains_trait and contains_marker:
            self.table_type = Table.TYPE_TRAIT_LIST
        elif contains_pval and contains_marker:
            self.table_type = Table.TYPE_MARKER_LIST
        if self.table_type and self.table_type > 0:
            print(self.table_type)

    def add_spacy_docs(self, nlp):
        if self.caption_text:
            self.caption_doc = nlp.process_corpus(self.caption_text)
        if self.footer_text:
            self.footer_doc = nlp.process_corpus(self.footer_text)
        if self.title:
            if type(self.title) == list:
                self.title = self.title[0]
            self.title_doc = nlp.process_corpus(self.title)
        if self.column_cells:
            for cell in self.column_cells:
                cell.add_spacy_docs(nlp)
        if self.data_cells:
            for row in self.data_cells:
                for cell in row:
                    cell.add_spacy_docs(nlp)
        self.__set_table_type()

    def assign_annotations(self, t, m, p, r, used_annots):
        docs_iter = [(self.title_doc, self.title_offset, "table_title"),
                     (self.caption_doc, self.caption_offset, "table_caption"),
                     (self.footer_doc, self.footer_offset, "table_footer")]
        cells_iter = [(self.column_cells, self.content_offset, "table_content"),
                      (self.data_cells, self.content_offset, "table_content")]
        annots, used_annots = [], []
        for (doc, offset, elem) in docs_iter:
            if doc and doc.ents:
                annots, used_annots, t, m, p, r = BioC.get_bioc_annotations(doc, used_annots, offset, t, m, p, r, elem)
        for (rows, offset, elem) in cells_iter:
            for row in rows:
                if type(row) == TableCell:
                    if row.doc and row.doc.ents:
                        annots, used_annots, t, m, p, r = BioC.get_bioc_annotations(row.doc, used_annots, offset, t, m,
                                                                                    p, r, elem, row.id)
                else:
                    for i in range(len(row)):
                        cell = row[i]
                        if self.column_types[i] == Table.COLUMN_TRAIT or self.column_types[i] == Table.COLUMN_MARKER:
                            if cell.doc and cell.doc.ents:
                                annots, used_annots, t, m, p, r = BioC.get_bioc_annotations(cell.doc, used_annots, offset,
                                                                                            t, m, p, r, elem, cell.id)
                        elif self.column_types[i] == Table.COLUMN_SIGNIFICANCE:
                            if cell.doc:
                                cell.doc.ents += cell.doc[:]
                                annots, used_annots, t, m, p, r = BioC.get_bioc_annotations(cell.doc, used_annots, offset,
                                                                                            t, m, p, r, elem, cell.id)

        self.annotations = annots
        return t, m, p, r, used_annots

    def jsonable(self):
        output_dict = self.__dict__
        del output_dict['title_doc']
        del output_dict['footer_doc']
        del output_dict['caption_doc']
        return output_dict


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
