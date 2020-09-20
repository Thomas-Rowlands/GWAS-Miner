import re
from lxml import etree
from lxml.etree import tostring
from DataStructures import Study, Table
from Utility_Functions import Utility
import CharacterDealer
import logging


class PreProcessing:
    __logger = logging.getLogger("Phenotype Finder")

    # Retrieves rs identifiers from input string
    @staticmethod
    def get_rs_identifiers(text):
        """
        Checks to see if an rs identifier is present within the input string
        @param text: string to check for rs identifiers
        @return: rs identifier if identified, otherwise None
        """
        result = re.findall("(?:rs[0-9]{1,}){1}", text)
        if not result:
            return None
        else:
            return result

    @staticmethod
    def __get_authors(tree):
        """
        Retrieve list of authors identified in the authors section of a publication.
        @param tree: etree parsed xml root element
        @return: list of author lists containing their name and associated details.
        """
        results = []
        for i in range(int(tree.xpath("count(//contrib[@contrib-type='author'])"))):
            author = [tree.xpath("//contrib[@contrib-type='author'][" + str(i + 1) + "]//text()")[0],
                      tree.xpath(
                          "//contrib[@contrib-type='author'][" + str(i + 1) + "]//text()")[1],
                      tree.xpath("//contrib[@contrib-type='author'][" + str(i + 1) + "]//text()")[-1:][0]]
            results.append(author)
        return results

    @staticmethod
    def __get_abstract(tree):
        """
        Retrieve the abstract section from a pubmed publication
        @param tree: etree parsed xml root element
        @return: String contaning the abstract section of the publication.
        """
        tmp = tree.xpath("//abstract//*/text()")
        abstract = ""
        # xpath can sometimes return a list with empty strings, typically if other elements are within it.
        for i in tmp:
            abstract = abstract + " " + i
        return abstract

    @staticmethod
    def __get_title(tree):
        """
        Retrieve the title text of the pubmed publication.
        @param tree: etree parsed xml root element
        @return: String containing the title of the publication.
        """
        return Utility.expand_xpath_output(tree.xpath("//article-title[1]/text()"))

    @staticmethod
    def __get_sections(tree):
        """
        Divides the publication text by section using XPath
        @param tree: the etree parsed XML root element
        @return: list containing lists of section titles along with their bodies of text
        """
        sections = []
        tables = tree.xpath("//table")
        for t in tables:
            children = t.getchildren()
            for child in children:
                t.remove(child)
        # Iterate through each section title
        for i in range(int(tree.xpath("count(//body//sec)"))):
            content = ""
            section_text = tree.xpath(
                "//body/sec[" + str(i + 1) + "]//*[not(self::title)]/text()")
            section_title = tree.xpath(
                "//body/sec[" + str(i + 1) + "]/title/text()")
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
        # Debugging only
        # for i in row_nums:
        #     if type(i) != list:
        #         PreProcessing.__logger.info(xml_table.xpath(
        #             ".//tbody//tr[position() = " + str(i) + "]//text()"))
        #     else:
        #         PreProcessing.__logger.info("Between: \n")
        #         PreProcessing.__logger.info(xml_table.xpath(
        #             ".//tbody//tr[position() = " + str(i[0]) + "]//text()"))
        #         PreProcessing.__logger.info("\nAnd:")
        #         PreProcessing.__logger.info(xml_table.xpath(
        #             ".//tbody//tr[position() = " + str(i[1]) + "]//text()"))
        result = []
        first_table_xml = "<table><thead>"
        for x in xml_table.xpath(".//thead//tr"):
            test = tostring(x, encoding="unicode", method="xml")
            first_table_xml += tostring(x, encoding="unicode", method="xml")
        if type(row_nums[0]) != list:
            if row_nums[0] <= 3:
                for elem in xml_table.xpath(".//tbody//tr[position() <" + str(row_nums[0]) + "]"):
                    first_table_xml += tostring(elem,
                                                encoding="unicode", method="xml")
        else:
            if row_nums[0][1] <= 3:
                for elem in xml_table.xpath(".//tbody//tr[position() <" + str(row_nums[0][1]) + "]"):
                    first_table_xml += tostring(elem,
                                                encoding="unicode", method="xml")
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
                new_table_xml += tostring(elem,
                                          encoding="unicode", method="xml")
            new_table_xml += "</thead><tbody>"
            if i == (len(row_nums) - 1):
                for elem in xml_table.xpath(".//tbody//tr[position() > " + str(row_nums[i][1]) + "]"):
                    new_table_xml += tostring(elem,
                                              encoding="unicode", method="xml")
            else:
                for elem in xml_table.xpath(".//tbody//tr[position() > " + str(row_nums[i][1]) +
                                            "and position() < " + str(row_nums[i + 1][0]) + "]"):
                    new_table_xml += tostring(elem,
                                              encoding="unicode", method="xml")

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
        """
        Retrieve the caption under the provided xml element
        @param tree: etree parsed xml element
        @return: string containing the caption text
        """
        caption = Utility.expand_xpath_output(tree.xpath("../caption//p/text()"))
        if caption == None:
            caption = Utility.expand_xpath_output(tree.xpath("//table/../../caption//text()"))
        return caption

    @staticmethod
    def __validate_caption(caption):
        """
        Check the caption text for signs of containing data we definitely do not want (e.g. candidate genes).
        @param caption:
        @return:
        """
        invalidation_list = ["candidate gene", "in/near", "genes"]
        for item in invalidation_list:
            if item in caption.lower():
                return False
        return True

    @staticmethod
    def __get_tables(tree):
        """
        Parses the results tables and retrieves SNPs with their associated p-value and phenotype.
        @param tree: etree parsed XML document containing results tables.
        @return: List of results tables and a list of SNP objects.
        """
        tables = tree.xpath("//table")
        new_tables = []
        results = []
        captions = []
        snps = []
        table_count = 1
        for table in tables:
            # Retrieve caption for each table with validation.
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
            # PubMed XML docs' typical format separates data & headings through horizontal rules.
            # Separator_indexes used to log each point in a table where a separate table could be made.
            # Table columns will sometimes change mid-table multiple times, changing the datatypes in each column
            # so separating the tables is a MUST!
            for row in body_rows:
                span_count = 0
                if row.xpath(".//hr"):
                    if is_prev_row_header:
                        separator_indexes.append(row_index)
                    is_prev_row_rule = True
                    row_index += 1
                    continue
                # Span value is used to track column ordering compared to the headers, as they can split headers into
                # multiple columns of data.
                span_count = row.xpath(".//td[1]//@colspan")
                if span_count:
                    span_count = int(span_count[0])
                else:
                    span_count = 0
                # Bold count can also be used to indicate the presence of a row of new headers.
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
            # Divide up the table into a new table for each new detected header change.
            if separator_indexes:
                ranges = PreProcessing.__sort_ranges(separator_indexes)
                is_two_dimensional = True
                # Check if the new header row is actually made up of 2 or more rows together.
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
        table_num = 1
        caption_counter = 0
        # Parse the newly defined tables and store their data as Table objects.
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
    def __char_encoding_filter(input_string):
        """
        Removes all problematic characters from the input string.
        @param input_string: The string to replace
        @return: The input string with problematic characters replaced where appropriate, or removed.
        """
        spaces = [x.lower_ for x in CharacterDealer.spaces]
        for space in spaces:
            input_string = input_string.replace(space, " ")
        return input_string

    @staticmethod
    def strip_xml(pmid, xml):
        """
        Divides the XML document into study sections and tables
        @param pmid: pubmed identifier for the document
        @param xml: the xml document to parse
        @return: Study object containing the split sections of the publication.
        """
        study = Study()
        # xml = __char_encoding_filter(xml)
        xml = re.sub("</", " </", xml)  # Ensure white space is present between nested tags
        study.original = xml
        parser = etree.XMLParser(encoding='utf-8')
        tree = etree.fromstring(xml, parser)
        study.pmid = pmid
        study.abstract = PreProcessing.__get_abstract(tree)
        study.title = PreProcessing.__get_title(tree)
        study.authors = PreProcessing.__get_authors(tree)

        #  Back sections
        acknowledgements = ""
        content = tree.xpath("//back/ack//*[not(self::title)]/text()")
        for tmp in content:
            acknowledgements = acknowledgements + tmp
        study.acknowledgements = acknowledgements

        #  Tables
        table_data = PreProcessing.__get_tables(tree)
        study.__tables = table_data[0]
        study.set_snps(table_data[1])

        # Sections
        study.sections = PreProcessing.__get_sections(tree)

        #  Citations pending

        return study
