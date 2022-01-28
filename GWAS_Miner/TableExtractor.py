import json

import BioC


def parse_tables(file_input):
    tables = []
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
            for passage in doc['passages']:
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
    return tables


class Table:
    def __init__(self, title, table_id, title_offset, content_offset, column_cells, data_cells,
                 caption_text, caption_offset, footer_text, footer_offset):
        """
        Table object from parsed JSON input.
        """
        self.title = title
        self.table_id = table_id
        self.title_offset = title_offset
        self.content_offset = content_offset
        self.column_cells = column_cells
        self.data_cells = data_cells
        self.caption_text = caption_text
        self.caption_offset = caption_offset
        self.footer_text = footer_text
        self.footer_offset = footer_offset


class TableCell:
    def __init__(self, id, text):
        self.id = id
        self.text = text
