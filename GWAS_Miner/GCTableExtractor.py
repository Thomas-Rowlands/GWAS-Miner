import json
from datetime import datetime

import BioC
from TableExtractor import Table, get_cell_entity_annotation, TableRow, TableSection, TableCell, TablePassage

table_significance_pattern = r""


def process_tables(nlp, tables):
    for table in tables:
        nlp = table.add_spacy_docs(nlp)
        if table.table_type:
            nlp = table.get_gc_annotations(nlp)
        table.annotations = nlp.annotations
        table.relations = nlp.relations
        nlp.annotations = []
        nlp.relations = []
    return tables, nlp


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
                table_passage = TablePassage(passage['infons']['section_title_1'])
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
                        if 'table_section_title_1' in section.keys():
                            return None, False
                        section_title = section['text'] if 'text' in section.keys() else None
                        title_offset = section['offset'] if section_title else None
                        table_section = TableSection(title=section_title, title_offset=title_offset)
                        for row in section['data_rows']:
                            data_row = TableRow()
                            for cell in row:
                                data_cell = TableCell(cell['cell_id'], str(cell['cell_text']))
                                data_row.cells.append(data_cell)
                            if abs(len(data_row.cells) - len(columns[0].cells)) > 1:
                                print(F"{file_input} contains potentially problematic cell counts!")
                            table_section.rows.append(data_row)
                        table_sections.append(table_section)
                table_passage.sections = table_sections
            table = GCTable(title, table_id, title_offset, content_offset, columns,
                            table_sections, caption, caption_offset, footer, footer_offset)
            tables.append(table)
    except IOError:
        print(F"No tables file found for {file_input.replace('_tables.json', '')}")
    if tables_data:
        annotated_tables, nlp = process_tables(nlp, tables)
        for table in annotated_tables:
            if table.annotations:
                contains_annotations = True
                for bioc_table in tables_data["documents"]:
                    if bioc_table["id"] == table.table_id:
                        for passage in bioc_table["passages"]:
                            passage["annotations"] = []
                            used_ids = []
                            passage["relations"] = []
                            for annotation in table.annotations:
                                if annotation.locations[0].table_element == passage['infons']['section_title_1']:
                                    passage["annotations"].append(annotation)
                                    used_ids.append(annotation.id)
                            for relation in table.relations:
                                if [x for x in relation.nodes if x.refid in used_ids]:
                                    passage["relations"].append(relation)
                    bioc_table["relations"] = []
    return tables_data, contains_annotations


class GCTable(Table):

    def __init__(self, title, table_id, title_offset, content_offset, column_cells, data_sections,
                 caption_text, caption_offset, footer_text, footer_offset):
        super().__init__(title, table_id, title_offset, content_offset, column_cells, data_sections,
                         caption_text, caption_offset, footer_text, footer_offset)

    def get_gc_annotations(self, nlp):
        current_datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
        annotations = []
        relations = []
        section_num = 1
        for section in self.data_sections:
            contains_trait, contains_variant, contains_significance = False, False, False
            for row in section.rows:
                row_contains_variant = False
                row_contains_significance = False
                row_annotations = []
                for cell in row.cells:
                    if cell.doc and cell.doc.ents:
                        entity = cell.doc.ents[0]
                        if entity.label_ == "PVAL":
                            contains_significance = True
                            row_contains_significance = True
                        elif entity.label_ == "RSID":
                            contains_variant = True
                            row_contains_variant = True
                        else:
                            contains_trait = True
                        annotation, nlp = get_cell_entity_annotation(nlp, entity, "table_content", cell.id)
                        row_annotations.append(annotation)
                if not contains_trait:
                    if self.COLUMN_TRAIT in self.section_ents[section_num - 1]:
                        entity = [x for x in section.doc.ents if x._.is_trait or x._.has_trait][0]
                        annotation, nlp = get_cell_entity_annotation(nlp, entity, "section_title",
                                                                     None, section.title_offset)
                        row_annotations.append(annotation)
                        contains_trait = True
                    elif self.COLUMN_TRAIT in self.caption_ents:
                        entity = [x for x in self.caption_doc.ents if x._.is_trait or x._.has_trait][0]
                        annotation, nlp = get_cell_entity_annotation(nlp, entity, "table_caption", None)
                        row_annotations.append(annotation)
                        contains_trait = True
                    elif self.COLUMN_TRAIT in self.footer_ents:
                        entity = [x for x in self.footer_doc.ents if x._.is_trait or x._.has_trait][0]
                        annotation, nlp = get_cell_entity_annotation(nlp, entity, "table_footer", None)
                        row_annotations.append(annotation)
                        contains_trait = True
                if row_contains_variant and row_contains_significance: #len(row_annotations) == 3:
                    nlp.v += 1
                    if contains_trait:
                        nlp.t += 1
                    nlp.s += 1
                    annotations.extend(row_annotations)
                    if contains_trait:
                        relations.append(BioC.BioCRelation(id=F"R{nlp.r}",
                                                           infons={"type": "disease_assoc",
                                                                   "annotator": "GWASMiner@le.ac.uk",
                                                                   "updated_at": current_datetime},
                                                           nodes=[BioC.BioCNode(refid=x.id, role="") for x in
                                                                  row_annotations]))
                        nlp.r += 1
            section_num += 1
        nlp.annotations += annotations
        nlp.relations += relations
        return nlp

    def jsonable(self):
        output_dict = self.__dict__
        del output_dict['title_doc']
        del output_dict['footer_doc']
        del output_dict['caption_doc']
        return output_dict
