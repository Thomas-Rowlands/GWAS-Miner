from lxml import etree

class Table:
    def __init__(self):
        self.columns = None
        self.column_count = None
        self.rows = None
        self.row_count = None
        self.p_values = None
        self.snps = None
        self.header_row = None

    def parse_table(self, xml):
        tree = etree.fromstring(xml, etree.get_default_parser())
        self.row_count = tree.xpath("count(//tbody//tr)")
        self.column_count = tree.xpath("count(//thead/tr//td)")
        self.columns = [x for x in tree.xpath("//thead/tr//td/text()")]
