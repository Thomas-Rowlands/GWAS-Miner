import re
import sys
from collections import OrderedDict
from pprint import pprint

from lxml import etree

from DataStructures import SNP, Table

def parse_table_data(elem, table_num=None):
    table = Table(elem, table_num=table_num)
    print("\n TABLE - " + str(table_num) + "\n")
    # end of experimental block

    # for row_num in range(header_row_count):
    #     current_row = header.xpath(".//tr[" + str(row_num + 1) + "]//td")
    #     cell_count = int(header.xpath("count(.//tr[" + str(row_num + 1) + "]//td)"))
    #     i = 1
    #     for cell in current_row:  # Identify purpose of each column
    #         col_span = cell.xpath(".//@colspan")
    #         text = cell.xpath(".//text()")[0]
    #         p_type = re.search(r"(p{1}([- ]?val(ue)?))", text, flags=re.IGNORECASE)
    #         if p_type:
    #             if "GEE" in p_type.string and (gee_p_val_col is None):
    #                 gee_p_val_col = i
    #                 i += 1
    #                 continue
    #             elif "FBAT" in p_type.string:
    #                 fbat_p_val_col = i
    #                 i += 1
    #                 continue
    #             else:
    #                 misc_p_val_col = i
    #                 i += 1
    #                 continue
    #         if "snp" in text.lower():
    #             snp_col = i
    #             i += 1
    #             continue
    #         if "phenotype" == text.lower().replace("* ", ""):
    #             phenotype_col = i
    #             i += 1
    #             continue
    #         i += 1
    # if (fbat_p_val_col or gee_p_val_col) and snp_col:  # Check if table contains associated p-vals
    #     body = elem.xpath(".//tbody")
    #     rows = body[0].xpath(".//tr")
    #     for o in range(len(rows)):
    #         fbat_p_value = "".join(rows[o].xpath(".//td[" + str(fbat_p_val_col) + "]//text()")).replace(" ", "")
    #         gee_p_value = "".join(rows[o].xpath(".//td[" + str(gee_p_val_col) + "]//text()")).replace(" ", "")
    #         phenotype = "".join(rows[o].xpath(".//td[" + str(phenotype_col) + "]//text()")).replace(" ", "")
    #         snp = "".join(rows[o].xpath(".//td[" + str(snp_col) + "]//text()")).replace(" ", "")
    #         misc_p_value = "".join(rows[o].xpath(".//td[" + str(misc_p_val_col) + "]//text()")).replace(" ", "")
    #         if not re.search(r"(?:rs[0-9]{1,}){1}", snp):  # SNP validation
    #             print("INVALID SNP: " + snp + " table " + str(table_num))
    #             if table_num == 5:
    #                 sys.exit("Column: " + str(snp_col))
    #             continue
    #         elif (len(gee_p_value) < 4) and (len(fbat_p_value) < 4):  # P-validation
    #             continue
    #         elif not re.search(p_pattern, gee_p_value) and (len(fbat_p_value) < 4):
    #             print("IGNORED: " + gee_p_value + " ___ " + fbat_p_value)
    #             continue
    #         elif not re.search(p_pattern, fbat_p_value) and (len(gee_p_value) < 4):
    #             print("IGNORED: " + gee_p_value + " ___ " + fbat_p_value)
    #             continue
    #         newSNP = SNP()
    #         newSNP.rs_identifier = snp
    #         newSNP.gee_p_val = gee_p_value
    #         newSNP.fbat_p_val = fbat_p_value
    #         newSNP.misc_p_val = misc_p_value
    #         newSNP.phenotype = phenotype
    #         snps.append(newSNP)
    return table
