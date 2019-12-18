import pprint
import re
from collections import OrderedDict
from lxml import etree


def convert_to_list(num):
    result = []
    for i in range(num):
        result.append(1)
    return result


class Study:
    def __init__(self):
        self.title = None
        self.abstract = None
        self.authors = None
        self.snps = None
        self.concepts = None
        self.p_values = None
        self.results = None
        self.pmid = None
        self.gwas_id = None
        self.sections = None
        self.acknowledgements = None
        self.citations = None
        self.tables = None


class Table:
    def __init__(self, xml, targets=None, table_num=None):
        """
        Table data structure for processing/storing PubMed table XML/HTML.
        :param xml: The xml element containing data for the PubMed table
        :param targets: (Optional) List of heading strings to identify column indexes for
        :param table_num: (Optional) Integer identifier for the table
        """
        self.xml = xml
        self.caption = self.__get_caption(xml)
        self.p_values = None
        self.snps = []
        self.targets = targets
        self.table_num = table_num
        self.map, self.data = Table.__map_table(xml, table_num=self.table_num)
        self.rows = [x for x in self.data["body"]]
        self.headings = [y for y in self.data["header"]]
        self.target_indexes = None
        if self.targets:
            self.target_indexes = Table.__get_target_headings(self.targets)
        self.__get_snps()

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
                content_list = cell.xpath(".//text()")
                content = ""
                for text in content_list:
                    content += text
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
                content_list = cell.xpath(".//text()")
                content = ""
                if content_list:
                    for text in content_list:
                        content += text
                else:
                    content = content_list
                row_data.append(content)
            table["body"].append(row_data)
        return table

    def __get_table_column_types(self):
        valuable_fields = {"Phenotypes": [], "GEE": [], "FBAT": [], "MISC_PVAL": [], "SNP": []}
        acceptable_threshold = 85
        for i in range(len(self.headings)):
            for o in range(len(self.headings[i])):
                #  Test each body cell value in column
                rsid_count = 0  # Contains RSID, should be a SNP column.
                phenotype_count = 0  # High chance of being a phenotype descriptor of some sort.
                integer_count = 0  # Likely to be a quantity or other non-desirable value
                p_val_count = 0  # Likely to contain a p-value
                value_test_count = len(self.rows)
                for l in range(value_test_count):
                    cell_value = None
                    try:  # Malformed tables can contain horizontal rules which will cause an Index Error.
                        cell_value = str(self.rows[l][o])
                    except IndexError as e:
                        continue
                    if re.search(r"(?:rs[0-9]{1,}){1}", cell_value, re.IGNORECASE):
                        rsid_count += 1
                    elif re.fullmatch(r"(^[0-9]{1,}[ ]?$)", cell_value, re.IGNORECASE):
                        integer_count += 1
                    elif re.search(r"(\d+\.?\d?[ ]?[×xX][ ]?\d+-\d[\(]?\d?[\)]?)|(\d\.\d+ ?)$",
                                   cell_value, re.IGNORECASE):
                        p_val_count += 1
                    elif cell_value == " ":
                        continue
                    else:
                        phenotype_count += 1
                is_rsid = (rsid_count * (100 / value_test_count)) > acceptable_threshold
                is_phenotype = (phenotype_count * (100 / value_test_count)) > acceptable_threshold
                is_integer = (integer_count * (100 / value_test_count)) > acceptable_threshold
                is_p_val = (p_val_count * (100 / value_test_count)) > acceptable_threshold
                heading = self.headings[i][o].lower()
                if "snp" in heading:
                    if is_rsid:
                        valuable_fields["SNP"].append([i, o])
                elif "gee" in heading:
                    if is_p_val:
                        valuable_fields["GEE"].append([i, o])
                elif "fbat" in heading:
                    if is_p_val:
                        valuable_fields["FBAT"].append([i, o])
                elif "p-val" in heading:
                    if is_p_val:
                        valuable_fields["MISC_PVAL"].append([i, o])
                elif "phenotype" in heading:
                    if "p-val" not in heading:
                        if is_phenotype:
                            valuable_fields["Phenotypes"].append([i, o])
                elif heading and is_p_val:
                    back_counter = i - 1
                    while back_counter >= 0:
                        if "p-val" in self.headings[back_counter][o].lower():
                            valuable_fields["GEE"].append([i, o, heading])
                            break
                        else:
                            back_counter -= 1
        print(self.table_num)
        pprint.pprint(valuable_fields)
        return valuable_fields

    @staticmethod
    def __strip_pval(text):
        match = re.search(r"(\d+\.?\d?[ ]?[×xX][ ]?\d+-\d[\(]?\d?[\)]?)|(\d\.\d+ ?)$", text,
                          re.IGNORECASE)
        if match:
            return match.group()
        else:
            return None

    @staticmethod
    def __strip_rsid(text):
        match = re.search(r"(?:rs[0-9]{1,}){1}", text, re.IGNORECASE)
        if match:
            return match.group()
        else:
            return None

    def __get_snps(self):
        """
        Assigns genetic markers with their RS identifier and associated P-values & phenotype.
        @return:
        """
        table_targets = self.__get_table_column_types()
        for i in range(len(self.rows)):
            new_snp = SNP()
            is_snp_added = False
            if table_targets["Phenotypes"]:
                new_snp.phenotype = self.rows[i][table_targets["Phenotypes"][0][1]]
            if table_targets["GEE"]:
                if len(table_targets["GEE"]) > 1:
                    for entry in table_targets["GEE"]:
                        if len(entry) == 3:
                            if entry[2].replace(" ", ""):
                                new_snp = SNP()
                                new_snp.phenotype = entry[2]
                                new_snp.gee_p_val = Table.__strip_pval(self.rows[i][entry[1]])
                                if table_targets["SNP"]:
                                    new_snp.rs_identifier = Table.__strip_rsid(self.rows[i][table_targets["SNP"][0][1]])
                                self.snps.append(new_snp)
                                is_snp_added = True
                else:
                    new_snp.gee_p_val = Table.__strip_pval(self.rows[i][table_targets["GEE"][0][1]])

            if table_targets["FBAT"]:
                new_snp.fbat_p_val = Table.__strip_pval(self.rows[i][table_targets["FBAT"][0][1]])
            if table_targets["SNP"]:
                new_snp.rs_identifier = Table.__strip_rsid(self.rows[i][table_targets["SNP"][0][1]])
            if not is_snp_added:
                self.snps.append(new_snp)

    @staticmethod
    def __get_caption(elem):
        """
        Retrieves the accompanying table caption text.
        @param elem: XML root from which to retrieve the caption.
        @return: String containing caption text.
        """
        caption = elem.xpath(".//preceding-sibling::caption//p//text()")[0]
        return caption

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
    def __init__(self, rs_identifier=None):
        self.rs_identifier = rs_identifier
        self.gee_p_val = None
        self.fbat_p_val = None
        self.misc_p_val = None
        self.phenotype = None
        self.internal_marker = None


class MasterLexicon:
    def __init__(self, vocab=None):
        self.vocab = None
        self.master = None
        if vocab:
            self.parse(vocab)

    def parse(self, vocab):
        self.vocab = vocab
        self.master = self.__separate_lexicon()
        return self.master

    def __separate_lexicon(self):
        """
        (Internal) Divide the input vocabulary by the number of tokens in each term, for each ontology dict entry.
        @return: Dictionary containing the master lexicon generated.
        """
        master = {}
        for ontology in self.vocab.keys():
            master[ontology] = {}
            for phrase in self.vocab[ontology]:
                space_count = phrase.count(" ")
                if space_count in master[ontology]:
                    master[ontology][space_count].append(phrase)
                else:
                    master[ontology][space_count] = [phrase]
        return master

    def vocab_count(self, key):
        """
        Calculate the number of vocabularies in the master lexicon for a specified key.
        @return: Integer count of master lexicon vocabularies.
        """
        return len(self.master[key].keys())
