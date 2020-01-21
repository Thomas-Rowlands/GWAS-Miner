import pprint
import re
import sys

from lxml import etree
from lxml.etree import tostring

from DataStructures import Study, SNP, Table
from nltk.stem.porter import PorterStemmer
from nltk.stem.wordnet import WordNetLemmatizer
from nltk.corpus import stopwords
import nltk.tag
import nltk.data
import nltk.chunk
from pprint import pprint

nltk.download('stopwords')
nltk.download('wordnet')
nltk.download('punkt')


class BackoffTagger:
    def __init__(self):
        self._taggers = [nltk.data.load('taggers/maxent_treebank_pos_tagger/PY3/english.pickle')]


class PreProcessing:
    import pandas
    from string import punctuation
    import re

    # regexp_pattern = tag_pattern2re_pattern(pattern)

    # Retrieves rs identifiers from input string
    @staticmethod
    def get_rs_identifiers(text):
        result = re.findall("(?:rs[0-9]{1,}){1}", text)
        if not result:
            return None
        else:
            return result

    @staticmethod
    def __get_authors(tree):
        results = []
        for i in range(int(tree.xpath("count(//contrib[@contrib-type='author'])"))):
            author = [tree.xpath("//contrib[@contrib-type='author'][" + str(i + 1) + "]//text()")[0],
                      tree.xpath("//contrib[@contrib-type='author'][" + str(i + 1) + "]//text()")[1],
                      tree.xpath("//contrib[@contrib-type='author'][" + str(i + 1) + "]//text()")[-1:][0]]
            results.append(author)
        return results

    @staticmethod
    def __get_abstract(tree):
        tmp = tree.xpath("//abstract//*/text()")
        abstract = ""
        for i in tmp:
            abstract = abstract + " " + i
        return abstract

    @staticmethod
    def __get_title(tree):
        return tree.xpath("//article-title[1]/text()")[0]

    @staticmethod
    def __get_sections(tree):
        sections = []
        tables = tree.xpath("//table")
        for t in tables:
            children = t.getchildren()
            for child in children:
                t.remove(child)
        for i in range(int(tree.xpath("count(//body//sec)"))):  # Iterate through each section title
            content = ""
            section_text = tree.xpath("//body/sec[" + str(i + 1) + "]//*[not(self::title)]/text()")
            section_title = tree.xpath("//body/sec[" + str(i + 1) + "]/title/text()")
            for tmp in section_title:
                content = content + tmp
            section_title = content
            content = ""
            for tmp in section_text:
                content = content + tmp
            section_text = content
            section = [section_title,
                       section_text]
            sections.append(section)
        return sections

    @staticmethod
    def __divide_table(xml_table, row_nums):
        """
        Splits an XML table into multiple XML tables at specified break points
        @param xml_table: Lxml Element containing table body elements
        @param row_nums: List of Integers representing the row numbers to split the table at.
        @return: List of XML tables created by splitting the input table.
        """
        for i in row_nums:
            if type(i) != list:
                print(xml_table.xpath(".//tbody//tr[position() = " + str(i) + "]//text()"))
            else:
                print("Between: \n")
                print(xml_table.xpath(".//tbody//tr[position() = " + str(i[0]) + "]//text()"))
                print("\nAnd:")
                print(xml_table.xpath(".//tbody//tr[position() = " + str(i[1]) + "]//text()"))

        result = []
        ignore_row_index = 0
        first_table_xml = "<table><thead>"
        for x in xml_table.xpath(".//thead//tr"):
            test = tostring(x, encoding="unicode", method="xml")
            first_table_xml += tostring(x, encoding="unicode", method="xml")
        if type(row_nums[0]) != list:
            if row_nums[0] <= 3:
                ignore_row_index = row_nums[0]
                for elem in xml_table.xpath(".//tbody//tr[position() <" + str(row_nums[0]) + "]"):
                    first_table_xml += tostring(elem, encoding="unicode", method="xml")
                # first_table_xml += "</thead>"
        else:
            if row_nums[0][1] <= 3:
                for elem in xml_table.xpath(".//tbody//tr[position() <" + str(row_nums[0][1]) + "]"):
                    first_table_xml += tostring(elem, encoding="unicode", method="xml")
        first_table_xml += "</thead><tbody>"
        for x in xml_table.xpath(".//tbody//tr[position() < " + str(row_nums[1][0]) + "]"):
            first_table_xml += tostring(x, encoding="unicode", method="xml")
        first_table_xml += "</tbody></table>"
        result.append(first_table_xml)
        row_nums = row_nums[1:]

        for i in range(len(row_nums)):
            new_table_xml = "<table><thead>"
            for elem in xml_table.xpath(".//tbody//tr[position() >=" + str(row_nums[i][0]) + " and position() <= " +
                                        str(row_nums[i][1]) + "]"):
                new_table_xml += tostring(elem, encoding="unicode", method="xml")

            new_table_xml += "</thead><tbody>"
            if i == (len(row_nums) - 1):
                for elem in xml_table.xpath(".//tbody//tr[position() > " + str(row_nums[i][1]) + "]"):
                    new_table_xml += tostring(elem, encoding="unicode", method="xml")
            else:
                for elem in xml_table.xpath(".//tbody//tr[position() > " + str(row_nums[i][1]) +
                                            "and position() < " + str(row_nums[i + 1][0]) + "]"):
                    new_table_xml += tostring(elem, encoding="unicode", method="xml")

            new_table_xml += "</tbody></table>"
            result.append(new_table_xml)
        return result

    @staticmethod
    def __sort_ranges(values):
        """
        Group list of ascending integers by neighbouring values.
        @param values: List of integers in ascending order
        @return: List containing lone integers and lists of ranges for neighbouring values (e.g. [2, [27-29]]
        """
        result = []
        counter = 0
        for i in range(len(values)):
            if i == (len(values) - 1):
                if counter == 0:
                    result.append(values[i])
                else:
                    result.append([values[i - counter], values[i]])
                break
            first = values[i]
            second = values[i + 1]
            if (first + 2) >= second >= (first - 2):
                counter += 1
                continue
            elif counter == 0:
                result.append(values[i])
            else:
                result.append([values[i - counter], values[i]])
                counter = 0
        return result

    @staticmethod
    def __get_caption(tree):
        return str(tree.xpath("../caption//p/text()")[0])

    @staticmethod
    def __validate_caption(caption):
        if "candidate gene" in caption.lower():
            return False
        else:
            return True

    @staticmethod
    def __get_tables(tree):
        tables = tree.xpath("//table")
        new_tables = []
        results = []
        captions = []
        snps = []
        table_count = 1
        for table in tables:
            caption = PreProcessing.__get_caption(table)
            captions.append(caption)
            if not PreProcessing.__validate_caption(caption):
                continue
            header_cell_count = table.xpath("count(.//thead//td)")
            separator_indexes = []
            body_rows = table.xpath(".//tbody//tr")
            row_index = 1
            is_prev_row_rule = False
            is_prev_row_header = True
            for row in body_rows:
                test_content = row.xpath(".//text()")
                span_count = 0
                if row.xpath(".//hr"):
                    if is_prev_row_header:
                        separator_indexes.append(row_index)
                    is_prev_row_rule = True
                    row_index += 1
                    continue
                span_count = row.xpath(".//td[1]//@colspan")
                if span_count:
                    span_count = int(span_count[0])
                else:
                    span_count = 0
                bold_count = row.xpath("count(.//td//bold)")
                if (is_prev_row_rule or is_prev_row_header) and (bold_count > (header_cell_count - (span_count - 1))):
                    separator_indexes.append(row_index)
                    is_prev_row_header = True
                    row_index += 1
                    continue
                row_index += 1
                is_prev_row_rule = False
                is_prev_row_header = False
                table_count += 1
            if separator_indexes:
                ranges = PreProcessing.__sort_ranges(separator_indexes)
                is_two_dimensional = True
                for i in ranges[1:]:
                    if type(i) == list:
                        is_two_dimensional = False
                        break
                if is_two_dimensional:
                    new_tables.append(table)
                    continue
                new_tables.append([x for x in PreProcessing.__divide_table(table, ranges)])
            else:
                new_tables.append(table)
        #  end of experiment
        table_num = 1
        caption_counter = 0
        for table in new_tables:
            parsed_table = None
            if type(table) == list:
                for i in table:
                    parsed_table = Table(i, table_num=table_num, caption=captions[caption_counter])
                    results.append(parsed_table)
                    if parsed_table.snps:
                        for snp in parsed_table.snps:
                            snps.append(snp)
            else:
                parsed_table = Table(table, table_num=table_num, caption=captions[caption_counter])
                results.append(parsed_table)
                if parsed_table.snps:
                    for snp in parsed_table.snps:
                        snps.append(snp)
            table_num += 1
            caption_counter += 1
        return results, snps

    @staticmethod
    def strip_xml(pmid, xml):
        study = Study()
        xml = re.sub("</", " </", xml)  # Ensure white space is present between nested tags
        study.original = xml
        parser = etree.XMLParser(encoding='utf-8')
        tree = etree.fromstring(xml, parser)
        study.pmid = pmid

        study.abstract = PreProcessing.__get_abstract(tree)
        study.title = PreProcessing.__get_title(tree)
        study.authors = PreProcessing.__get_authors(tree)
        table_data = PreProcessing.__get_tables(tree)
        study.tables = table_data[0]
        study.set_snps(table_data[1])
        study.sections = PreProcessing.__get_sections(tree)

        #  Back sections
        acknowledgements = ""
        content = tree.xpath("//back/ack//*[not(self::title)]/text()")
        for tmp in content:
            acknowledgements = acknowledgements + tmp
        study.acknowledgements = acknowledgements

        #  Tables

        #  Citations pending

        return study
