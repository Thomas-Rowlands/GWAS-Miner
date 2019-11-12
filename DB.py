import mysql.connector


class Connection:

    def __init__(self, user, password, db, server, port):
        self.user = user
        self.password = password
        self.db = db
        self.server = server
        self.cursor = None
        self.con = None
        self.port = port

    def connect(self):
        """
        Connect using the assigned connection string.
        :return:
        """
        try:
            self.con = mysql.connector.connect(host=self.server, user=self.user, passwd=self.password, database=self.db,
                                               port=self.port)
            self.cursor = self.con.cursor()
            return True
        except mysql.connector.Error as err:
            print(err.msg)
            return False
        except:
            return False

    def disconnect(self):
        """
        Disconnect from the server.
        :return:
        """
        try:
            self.con.close()
            self.cursor = None
            return True
        except:
            return False

    def query(self, querystring):
        """
        Execute a query using the connection object.
        :param querystring: String containing the SQL query.
        :return: Cursor object containing results.
        """
        if self.connect():
            try:
                self.cursor.execute(querystring)
                result_cursor = self.cursor
                self.disconnect()
                return result_cursor
            except mysql.connector.Error as err:
                print(err.msg)
                return None
        else:
            print("Failed to connect to DB")
            return None

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
            values = values + "%s,"
        values = values[:-1]
        query_string = "INSERT INTO " + table + " (" + columns + ") VALUES (" + values + ")"
        cursor = self.con.cursor(prepared=True)
        cursor.execute(query_string, vals)
        self.con.commit()
        self.disconnect()

    def set_db(self, database):
        """
        Change the database and reconnect
        :param database: Name of the DB
        :return:
        """
        self.disconnect()
        self.db = database
        self.connect()
