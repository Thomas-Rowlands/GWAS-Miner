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
        self.authors = None
        self.__markers = []
        self.concepts = None
        self.p_values = None
        self.pmid = None
        self.gwas_id = None
        self.sections = []
        self.__tables = []
        self.__table_text = ""
        self.original = None
        self.__add_core_sections()
        self.__populate_study(study_json)
        self.__load_tables(study_tables_json)

    def __populate_study(self, study_json):
        """
        Parses the study JSON data to populate the study object.
        """
        misc = []
        self.title = list(study_json.keys())[0]
        for input_section in study_json[self.title]:
            for section in self.sections:
                match_found = False
                if section.get_name() in input_section[0].lower():
                    match_found = True
                    section.add_text(input_section[2])
                    break
                if section.get_name() == "Misc":
                    section.add_text(input_section[2])

    def __add_core_sections(self):
        self.sections.append(StudySection("Abstract", weighting=10))
        self.sections.append(StudySection("Introduction", weighting=3))
        self.sections.append(StudySection("Discussion", weighting=8))
        self.sections.append(StudySection("Results", weighting=10))
        self.sections.append(StudySection("Methods", weighting=5))
        self.sections.append(StudySection("Acknowledgements", weighting=1))
        self.sections.append(StudySection("Misc", weighting=5))

    def __load_tables(self, study_tables_json):
        """
        Parses the results tables from the accompanying JSON file.
        """
        for table in study_tables_json["tables"]:
            new_table = Table(table)
            self.__tables.append(new_table)
            text = new_table.get_text()
            if text:
                self.__table_text += text
            else:
                self.__tables = None
                return None

    def get_tables(self):
        return self.__tables

    def get_fulltext(self):
        """
        Generates a formatted free text string for visualisation.
        @return: String containing the full text study
        """
        result = F"{self.title} \n <new_line><new_line> \n"
        for section in self.sections:
            result += F"{section.get_text()} \n <new_line><new_line> \n"
        result += self.get_table_text()
        return result

    def get_formatted_fulltext(self):
        """
        Generates a free text string for the study excluding acknowledgements, citations and tables.
        @return: String containing the full text study
        """
        result = F"{self.title} \n\n"
        for section in self.sections:
            result += F"{section.get_text()} \n\n"
        result += self.get_table_text()
        return result

    def get_table_text(self):
        result = ""
        for table in self.__tables:
            result += F"<new_line><new_line> Table {table.table_num} \n <new_line><new_line> \n"
            result += F"{table.caption} \n <new_line><new_line> \n {table.get_text()}"
        return result  # self.__table_text

    def get_formatted_table_text(self):
        return self.__table_text

    def append_marker(self, new_marker):
        """
        Append a new unique marker object to the study
        @param new_marker: marker object to append to the study
        @return: True on a successful append, False if a duplicate is found
        """
        if not isinstance(new_marker, list):
            new_marker = [new_marker]
        for input_marker in new_marker:
            duplicate_found = False
            for marker in self.__markers:
                if marker.__dict__ == input_marker.__dict__:
                    duplicate_found = True
                    break
            if duplicate_found:
                continue
            else:
                self.__markers.append(input_marker)

    def set_markers(self, new_markers):
        """
        Setter for the list of marker objects of this study.
        @param new_markers: List of markers to assign to the study.
        @return: True on success
        """
        self.__markers = new_markers
        return True

    def get_markers(self):
        """
        Getter for all marker objects for this study
        @return: List of marker objects
        """
        return self.__markers


class StudySection:
    def __init__(self, name, text="", weighting=0):
        self.__name = name.lower()
        self.__text = text.lower()
        self.__weighting = weighting

    def set_name(self, name):
        self.__name = name

    def get_name(self):
        return self.__name

    def add_text(self, text):
        self.__text += text

    def set_text(self, text):
        self.__text = text

    def get_text(self):
        return self.__text

    def set_weighting(self, value):
        self.__weighting = value

    def get_weighting(self):
        return self.__weighting


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
        self.markers = []
        self.targets = targets
        self.table_num = self.data["identifier"]
        self.sections = [TableSection(x) for x in self.data["section"]]
        self.rows = self.__get_rows()
        self.columns = [x for x in self.data["columns"]]
        self.target_indexes = None
        self.__text = self.__convert_to_text()

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
        valuable_fields = {"Phenotypes": [], "GEE": [], "FBAT": [], "MISC_PVAL": [], "marker": []}
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
            if "marker" in heading:
                if is_rsid:
                    valuable_fields["marker"].append(i)
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
        if not valuable_fields["Phenotypes"] or not valuable_fields["marker"]:
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

    def __get_markers(self):
        """
        Assigns genetic markers with their RS identifier and associated P-values & phenotype.
        @return:
        """
        table_targets = self.__get_table_column_types()
        if not table_targets:
            return
        for i in range(len(self.rows)):
            new_marker = Marker()
            is_marker_added = False
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
                new_marker.phenotype = pheno_val
            if table_targets["GEE"]:
                if len(table_targets["GEE"]) > 1:
                    for entry in table_targets["GEE"]:
                        if len(entry) == 3:
                            if entry[2].replace(" ", ""):
                                new_marker = Marker()
                                new_marker.phenotype = entry[2]
                                new_marker.gee_p_val = Table.__strip_pval(self.rows[i][entry[1]])
                                if table_targets["marker"]:
                                    new_marker.rs_identifier = Table.__strip_rsid(
                                        self.rows[i][table_targets["marker"][0]])
                                self.markers.append(new_marker)
                                is_marker_added = True
                else:
                    new_marker.gee_p_val = Table.__strip_pval(self.rows[i][table_targets["GEE"][0]])

            if table_targets["FBAT"]:
                new_marker.fbat_p_val = Table.__strip_pval(self.rows[i][table_targets["FBAT"][0]])
            if table_targets["marker"]:
                new_marker.rs_identifier = Table.__strip_rsid(self.rows[i][table_targets["marker"][0]])
            if not is_marker_added:
                self.markers.append(new_marker)

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

    def __convert_to_text(self):
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

    def get_text(self):
        return self.__text


class Marker:
    def __init__(self, rs_identifier=None, gee_p_val=None, fbat_p_val=None, misc_p_val=None, phenotype=None,
                 internal_marker=None, weight=0):
        self.rs_identifier = rs_identifier
        self.gee_p_val = gee_p_val
        self.fbat_p_val = fbat_p_val
        self.misc_p_val = misc_p_val
        self.phenotype = phenotype
        self.internal_marker = internal_marker
        self.weight = weight


class MeshDescriptor:
    def __init__(self):
        self.name = ""
        self.ui = ""
        self.concepts = []
        self.tree_nums = []


class MeshConcept:
    def __init__(self):
        self.name = ""
        self.ui = ""
        self.is_preferred = ''
        self.terms = []


class MeshTerm:
    def __init__(self):
        self.name = ""
        self.ui = ""
        self.is_preferred = ''


class LexiconEntry:
    def __init__(self, identifer, name, mesh_descriptor=None):
        self.identifier = identifer
        self.__mesh_descriptor = mesh_descriptor
        self.__name = name
        self.__synonyms = []
        self.__token_size = name.count(" ")

    def get_token_size(self):
        return self.__token_size

    def set_name(self, name):
        self.__name = name
        self.__token_size = name.count(" ")

    def name(self):
        return self.__name

    def add_synonym(self, id, name):
        if id not in [x["id"] for x in self.__synonyms]:
            self.__synonyms.append({"id": id, "name": name})

    def remove_synonym(self, synonym):
        self.__synonyms.remove(synonym)

    def synonyms(self):
        return self.__synonyms


class Lexicon:
    def __init__(self, name):
        self.name = name
        self.__entries = []
        self.__identifiers = []
        self.__longest_term = 0

    def add_entry(self, entry):
        if isinstance(entry, LexiconEntry):
            self.__entries.append(entry)
            if entry.get_token_size() > self.__longest_term:
                self.__longest_term = entry.get_token_size()
            self.__identifiers.append(entry.identifier)
        else:
            raise TypeError("entry input must be of type LexiconEntry")

    def remove_entry(self, entry):
        if isinstance(entry, LexiconEntry):
            self.__entries.remove(entry)
        else:
            raise TypeError("entry input must be of type LexiconEntry")

    def get_entries(self):
        return self.__entries

    def add_descriptor(self, descriptor):
        if isinstance(descriptor, MeshDescriptor):
            self.__descriptors.append(descriptor)
        else:
            raise TypeError("descriptor input must be of type MeshDescriptor")

    def remove_descriptor(self, descriptor):
        if isinstance(descriptor, MeshDescriptor):
            self.__descriptors.remove(descriptor)
        else:
            raise TypeError("descriptor input must be of type MeshDescriptor")

    def add_descriptors(self, descriptors):
        if isinstance(descriptors, list):
            for descriptor in descriptors:
                self.add_descriptor(descriptor)
        else:
            raise TypeError("descriptors input must be of type List containing MeshDescriptor items")

    def identifier_used(self, identifier):
        return identifier in self.__identifiers

    def assign_synonym(self, identifier, name):
        entry = next((x for x in self.__entries if x.identifier == identifier), None)
        if entry:
            entry.add_synonym(name)

    def get_entry_by_term(self, term):
        for entry in self.__entries:
            if entry.name().lower() == term:
                return entry
            if term in [x.lower() for x in entry.synonyms()]:
                return entry
        return None


class MasterLexicon:
    def __init__(self):
        self.__lexicons = []
        self.__priority_order = {}

    def add_lexicon(self, new_lexicon):
        if isinstance(new_lexicon, Lexicon):
            if new_lexicon not in self.__lexicons:
                self.__lexicons.append(new_lexicon)
            else:
                logger.info("Lexicon already present in MasterLexicon object.")

    def remove_lexicon(self, target_lexicon):
        if not isinstance(target_lexicon, Lexicon):
            raise TypeError("Input must be of type Lexicon")
        for lexicon in self.__lexicons:
            if lexicon is target_lexicon:
                self.__lexicons.remove(lexicon)
                return True
        return False

    def get_lexicon_by_name(self, name):
        for lexicon in self.__lexicons:
            if lexicon.name == name:
                return lexicon

    def set_priority_order(self, new_priority):
        """

        :param new_priority: dict containing lexicon name (case insensitive) -> priority integer pairs
            (lower numbers have higher priority)
        :return: True if successful else False.
        """
        if not isinstance(new_priority, dict):
            raise TypeError("Input must be a dictionary")
        if len(new_priority.keys()) != len(self.__lexicons):
            raise IndexError("Dictionary key length does not match the number of current lexicons.")
        self.__priority_order = {k: v for k, v in sorted(new_priority.items(), key=lambda item: item[1])}
        return True

    def get_priority_order(self):
        return self.__priority_order

    def get_ordered_lexicons(self):
        ordered_lexicons = []
        for key in self.__priority_order:
            ordered_lexicons.append(self.get_lexicon_by_name(key))
        return ordered_lexicons

    def get_lexicon_entry(self, term, lexicon_name):
        lexicon = self.get_lexicon_by_name(lexicon_name)
        return lexicon.get_entry_by_term(term)
