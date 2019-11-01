import mysql.connector


class Connection:
    user = ""
    password = ""
    db = ""
    server = ""
    con = None
    cursor = None

    def __init__(self, user, password, db, server):
        self.user = user
        self.password = password
        self.db = db
        self.server = server

    def connect(self):
        """
        Connect using the assigned connection string.
        :return:
        """
        try:
            self.con = mysql.connector.connect(host=self.server, user=self.user, passwd=self.password, database=self.db)
            return True
        except:
            return False

    def disconnect(self):
        """
        Disconnect from the server.
        :return:
        """
        try:
            self.con.close()
            return True
        except:
            return False

    def query(self, querystring):
        """
        Execute a query using the connection object.
        :param querystring: String containing the SQL query.
        :return: Cursor object containing results.
        """
        self.cursor = self.con.cursor()
        try:
            self.cursor.execute(querystring)
        except:
            return None
        return self.cursor

    def insert(self, table, cols, vals):
        """
        Insert a new record in the current DB
        :param table: Target data table string
        :param cols: Target columns, array
        :param vals: Target values, array
        :return: Boolean representing success or failure
        """
        columns = ""
        values = ""
        for col in cols:
            if type(col) == str:
                columns = columns + "'" + str(col) + "',"
            else:
                columns = columns + col + ", "
        columns = columns[:-1]
        for val in vals:
            if type(val) == str:
                values = values + "'" + str(val) + "',"
            else:
                values = values + val + ", "
        values = values[:-1]
        query_string = "INSERT INTO " + table + " (" + columns + ") VALUES (" + values + ")"
        self.cursor = self.con.cursor()
        self.cursor.execute(query_string)
        self.cursor.commit()

    def set_db(self, database):
        """
        Change the database and reconnect
        :param database: Name of the DB
        :return:
        """
        self.disconnect()
        self.db = database
        self.connect()
