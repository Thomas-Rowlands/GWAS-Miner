class TableTypeError(Exception):
    def __init__(self, table, message="Table type does not match row annotation layout"):
        self.table = table
        self.message = message
        super().__init__(self.message)
