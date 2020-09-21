import logging
import re

logger = logging.getLogger("GWAS Miner")


def convert_to_list(num):
    result = []
    for i in range(num):
        result.append(1)
    return result


class Study:
    def __init__(self, study_json, study_tables_json):
        self.title = None
        self.abstract = ""
        self.introduction = ""
        self.discussion = ""
        self.methods = ""
        self.results = ""
        self.acknowledgements = ""
        self.citations = ""
        self.authors = None
        self.__snps = None
        self.concepts = None
        self.p_values = None
        self.pmid = None
        self.gwas_id = None
        self.sections = []
        self.__tables = []
        self.table_text = ""
        self.original = None
        self.__populate_study(study_json)
        self.__load_tables(study_tables_json)

    def __populate_study(self, study_json):
        """
        Parses the study JSON data to populate the study object.
        """
        self.title = list(study_json.keys())[0]
        for section in study_json[self.title]:
            if "abstract" in section[0].lower():
                self.abstract = self.abstract + section[2]
            elif "acknowledgements" in section[0].lower():
                self.acknowledgements = self.acknowledgements + section[2]
            elif "discussion" in section[0].lower():
                self.discussion = self.discussion + section[2]
            elif "introduction" in section[0].lower():
                self.introduction = self.introduction + section[2]
            elif "results" in section[0].lower():
                self.results = self.results + section[2]
            elif "methods" in section[0].lower():
                self.methods = self.methods + section[2]
            self.sections.append([section[0], section[2]])

    def __load_tables(self, study_tables_json):
        """
        Parses the results tables from the accompanying JSON file.
        """
        for table in study_tables_json["tables"]:
            new_table = Table(table)
            self.__tables.append(new_table)
            text = new_table.convert_to_text()
            if text:
                self.table_text += text
            else:
                self.__tables = None
                return None

    def get_fulltext(self):
        """
        Generates a free text string for the study excluding acknowledgements, citations and tables.
        @return: String containing the full text study
        """
        if self.sections:
            result = F"{self.title}.\n{self.abstract}\n"
            for section in self.sections:
                result += F" {section[1]}\n"
            result += F"\n{self.table_text}"
            return result
        else:
            return None

    def append_snp(self, new_snp):
        """
        Append a new unique SNP object to the study
        @param new_snp: SNP object to append to the study
        @return: True on a successful append, False if a duplicate is found
        """
        duplicate_found = False
        for snp in self.__snps:
            if snp.__dict__ == new_snp.__dict__:
                duplicate_found = True
                break
        if duplicate_found:
            return False
        else:
            self.__snps.append(new_snp)
            return True

    def set_snps(self, new_snps):
        """
        Setter for the list of SNP objects of this study.
        @param new_snps: List of snps to assign to the study.
        @return: True on success
        """
        self.__snps = new_snps
        return True

    def get_snps(self):
        """
        Getter for all SNP objects for this study
        @return: List of SNP objects
        """
        return self.__snps


class TableSection:
    def __init__(self, table_section):
        self.name = table_section["section_name"]
        self.rows = [x for x in table_section["results"]]


class Table:
    def __init__(self, data, targets=None):
        """
        Table data structure
        :param caption: (Optional) String representing the table's caption
        :param targets: (Optional) List of heading strings to identify column indexes for
        :param table_num: (Optional) Integer identifier for the table
        """
        self.data = data
        self.caption = self.data["title"]
        self.p_values = None
        self.snps = []
        self.targets = targets
        self.table_num = self.data["identifier"]
        self.sections = [TableSection(x) for x in self.data["section"]]
        self.rows = self.__get_rows()
        self.columns = [x for x in self.data["columns"]]
        self.target_indexes = None
        # if self.targets: ###removed while testing sentence conversion.
        #     self.target_indexes = Table.__get_target_headings(self.targets)
        # self.__get_snps()

    def __get_rows(self):
        rows = []
        for section in self.sections:
            for row in section.rows:
                rows.append(row)
        return rows

    def set_targets(self, targets):
        """
        Set target headers for this table instance
        :param targets: List of string headers to search for
        :return:
        """
        self.targets = targets
        self.target_indexes = Table.__get_target_headings(table=self.data, target_headings=self.targets)

    def __get_table_column_types(self):
        valuable_fields = {"Phenotypes": [], "GEE": [], "FBAT": [], "MISC_PVAL": [], "SNP": []}
        acceptable_threshold = 80
        value_test_count = len(self.rows)
        data = [x for x in [i for i in self.columns] if x]
        if not data:
            return valuable_fields
        for i in range(len(self.columns)):
            if self.columns[i] == '' or self.columns[i].count(self.columns[i][0]) == len(self.columns[i]):
                continue
            #  Test each body cell value in column
            rsid_count = 0  # Contains RSID, should be a marker column.
            phenotype_count = 0  # High chance of being a phenotype descriptor of some sort.
            integer_count = 0  # Likely to be a quantity or other non-desirable value
            p_val_count = 0  # Likely to contain a p-value
            blank_cells = 0
            for o in range(value_test_count):
                cell_value = str(self.rows[o][i])
                if re.search(r"(?:rs[0-9]{1,}){1}", cell_value, re.IGNORECASE):
                    rsid_count += 1
                elif re.fullmatch(r"(^[0-9]{1,}[ ]?$)", cell_value, re.IGNORECASE):
                    integer_count += 1
                elif re.search(r"(\d+\.?\d?[ ]?[×xX*][ ]?\d+-\d[\(]?\d?[\)]?)|(\d\.\d+ ?)$",
                               cell_value, re.IGNORECASE):
                    p_val_count += 1
                elif not cell_value.replace(" ", ""):
                    blank_cells += 1
                else:
                    phenotype_count += 1
            tested_cells = value_test_count - blank_cells
            is_rsid = False
            is_phenotype = False
            is_p_val = False
            if tested_cells > 0:
                is_rsid = (rsid_count * (100 / value_test_count)) > acceptable_threshold
                is_phenotype = (phenotype_count * (100 / value_test_count)) > acceptable_threshold
                is_p_val = (p_val_count * (100 / value_test_count)) > acceptable_threshold
            heading = self.columns[i].lower()
            if "snp" in heading:
                if is_rsid:
                    valuable_fields["SNP"].append(i)
            elif "gee" in heading:
                if is_p_val:
                    valuable_fields["GEE"].append(i)
            elif "fbat" in heading:
                if is_p_val:
                    valuable_fields["FBAT"].append(i)
            elif "p-val" in heading:
                if is_p_val:
                    valuable_fields["MISC_PVAL"].append(i)
            elif "phenotype" in heading or "trait" in heading:
                if "p-val" not in heading:
                    if is_phenotype:
                        valuable_fields["Phenotypes"].append([i, heading])  # Do not remove heading!
        if not valuable_fields["Phenotypes"] or not valuable_fields["SNP"]:
            return None
        elif not valuable_fields["GEE"] and not valuable_fields["FBAT"] and not valuable_fields["MISC_PVAL"]:
            return None
        # Debugging Only
        # logger.info(self.caption)
        # logger.info(self.table_num)
        # logger.info(valuable_fields)
        return valuable_fields

    @staticmethod
    def __strip_pval(text):
        match = re.search(r"(\d+\.?\d?[ ]?[×xX*][ ]?\d+-\d[\(]?\d?[\)]?)|(\d\.\d+ ?)$", text,
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
        if not table_targets:
            return
        for i in range(len(self.rows)):
            new_snp = SNP()
            is_snp_added = False
            if table_targets["Phenotypes"]:
                pheno_val = self.rows[i][table_targets["Phenotypes"][0][0]]
                if not pheno_val.replace(" ", ""):
                    back_counter = i
                    while back_counter >= 0:  # Phenotypes can be stated once but used for multiple rows.
                        prev_val = self.rows[back_counter][table_targets["Phenotypes"][0][0]].replace(" ", "")
                        if not prev_val:
                            back_counter -= 1
                            continue
                        else:
                            pheno_val = self.rows[back_counter][table_targets["Phenotypes"][0][0]]
                            break
                new_snp.phenotype = pheno_val
            if table_targets["GEE"]:
                if len(table_targets["GEE"]) > 1:
                    for entry in table_targets["GEE"]:
                        if len(entry) == 3:
                            if entry[2].replace(" ", ""):
                                new_snp = SNP()
                                new_snp.phenotype = entry[2]
                                new_snp.gee_p_val = Table.__strip_pval(self.rows[i][entry[1]])
                                if table_targets["SNP"]:
                                    new_snp.rs_identifier = Table.__strip_rsid(self.rows[i][table_targets["SNP"][0]])
                                self.snps.append(new_snp)
                                is_snp_added = True
                else:
                    new_snp.gee_p_val = Table.__strip_pval(self.rows[i][table_targets["GEE"][0]])

            if table_targets["FBAT"]:
                new_snp.fbat_p_val = Table.__strip_pval(self.rows[i][table_targets["FBAT"][0]])
            if table_targets["SNP"]:
                new_snp.rs_identifier = Table.__strip_rsid(self.rows[i][table_targets["SNP"][0]])
            if not is_snp_added:
                self.snps.append(new_snp)

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

    def convert_to_text(self):
        output = ""
        try:
            for i in range(len(self.rows)):
                output += F"{self.caption.replace('.', ',')}, "
                for o in range(len(self.columns)):
                    column_value = self.columns[o]
                    cell_value = self.rows[i][o]
                    if not column_value:
                        column_value = "<!blank!>"
                    if not cell_value:
                        cell_value = "<!blank!>"
                    output += F"{column_value} is {cell_value} and "
                output = output[:-5]
                output += ".\n"
        except IndexError as ie:
            print(F"Index error found in table:{self.table_num}.\nPlease investigate this publication table JSON for "
                  F"discrepancies.")
            return None
        return output.replace("<!>", "")


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
