from collections import OrderedDict
from lxml import etree


def convert_to_list(num):
    result = []
    for i in range(num):
        result.append(1)
    return result


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


class Table:
    def __init__(self, xml, targets=None, table_num=None):
        """
        Table data structure for processing/storing PubMed table XML/HTML.
        :param xml: The xml element containing data for the PubMed table
        :param targets: (Optional) List of heading strings to identify column indexes for
        :param table_num: (Optional) Integer identifier for the table
        """
        self.xml = xml
        self.p_values = None
        self.snps = None
        self.targets = targets
        self.table_num = table_num
        self.map, self.data = Table.__map_table(xml, table_num=self.table_num)
        self.rows = [x for x in self.data["body"]]
        self.headings = [y for y in self.data["header"]]
        self.target_indexes = None
        if self.targets:
            self.target_indexes = Table.__get_target_headings(self.targets)

    def set_targets(self, targets):
        """
        Set target headers for this table instance
        :param targets: List of string headers to search for
        :return:
        """
        self.targets = targets
        self.target_indexes = Table.__get_target_headings(table=self.data, target_headings=self.targets)

    @staticmethod
    def __validate_table(elem):
        """
        (Internal) Extract table data to a dictionary structure whilst duplicating spanned headers for ease of access.
        :param elem: The XML element containing the table data
        :return: Dictionary containing the table structure & data.
        """
        #  Make table array from XML table
        header = elem.xpath(".//thead")[0]
        header_rows = header.xpath(".//tr")
        table = {"header": []}
        for row in header_rows:
            row_cells = row.xpath(".//td")
            row_data = []
            for cell in row_cells:
                content = cell.xpath(".//text()")[0]
                span = cell.xpath(".//@colspan")
                if span:  # Duplicate cells where they span to align columns
                    for i in range(int(span[0])):
                        row_data.append(content)
                else:
                    row_data.append(content)
            table["header"].append(row_data)
        body = elem.xpath(".//tbody")[0]
        body_rows = body.xpath(".//tr")
        table["body"] = []
        for row in body_rows:
            row_cells = row.xpath(".//td")
            row_data = []
            for cell in row_cells:
                content = cell.xpath(".//text()")[0]
                row_data.append(content)
            table["body"].append(row_data)
        return table

    @staticmethod
    def __map_table(elem, table_num=None):
        """
        (Internal) Maps the table's original structure only, using 1s to represent cells and nested arrays
        representing spans.
        :param elem: The XML element containing the table data
        :param table_num: (Optional) The table's integer identifier
        :return: 2 Dictionaries, the table map and the validated table data
        """
        header = elem.xpath(".//thead")[0]
        body = elem.xpath(".//tbody")[0]
        header_row_count = int(elem.xpath("count(.//thead//tr)"))
        body_row_count = int(elem.xpath("count(.//tbody//tr)"))
        table_map = {"h": [], "b": []}
        #  Step 1 -> Basic mapping of table structure with column span values
        for i in range(int(header_row_count)):
            current_header_row = header.xpath(".//tr[" + str(i + 1) + "]")[0]
            table_map["h"].append([])
            for o in range(int(current_header_row.xpath("count(.//td)"))):
                current_cell_span = current_header_row.xpath(".//td[" + str(o + 1) + "]//@colspan")
                if not current_cell_span:
                    current_cell_span = 1
                else:
                    current_cell_span = int(current_cell_span[0])
                table_map["h"][i].append(current_cell_span)
        for i in range(int(body_row_count)):
            current_body_row = body.xpath(".//tr[" + str(i + 1) + "]")[0]
            table_map["b"].append([])
            for o in range(int(current_body_row.xpath("count(.//td)"))):
                current_cell_span = current_body_row.xpath(".//td[" + str(o + 1) + "]//@colspan")
                if not current_cell_span:
                    current_cell_span = 1
                else:
                    current_cell_span = int(current_cell_span[0])
                table_map["b"][i].append(current_cell_span)
        # Step 2 -> Process table map to split spanned headings
        for row in table_map["h"]:
            for i in range(len(row)):
                if row[i] > 1:
                    row[i] = convert_to_list(row[i])
        for row in table_map["b"]:
            for i in range(len(row)):
                if row[i] > 1:
                    row[i] = convert_to_list(row[i])
        table = Table.__validate_table(elem)
        return table_map, table

    @staticmethod
    def __get_target_headings(table, target_headings):
        """
        (Internal) Identifies indexes of columns which match the target_headings strings
        :param table: Dictionary of table data
        :param target_headings: List of headings to look for.
        :return: List of indexes for each identified headings
        """
        target_indexes = {}
        for target in target_headings:
            tmp = []
            # Process headers
            for i in range(len(table["header"])):  # Row
                for o in range(len(table["header"][i])):  # Col
                    cell_value = table["header"][i][o].lower().replace("* ", "")
                    if target == cell_value:
                        tmp.append(o)
                        continue
            if not tmp:
                target_indexes[target] = None
            else:
                target_indexes[target] = tmp
        return target_indexes


class SNP:
    def __init__(self):
        self.rs_identifier = None
        self.gee_p_val = None
        self.fbat_p_val = None
        self.misc_p_val = None
        self.phenotype = None
        self.internal_marker = None
