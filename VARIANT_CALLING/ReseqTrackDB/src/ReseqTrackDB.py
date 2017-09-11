'''
Created on 25 Oct 2016

@author: ernesto
'''

import os
import subprocess
from datetime import datetime
import pymysql

class ReseqTrackDB(object):
    '''
    Class representing a Reseqtrack MySQL DB
    '''

    def __init__(self, host=None, user=None, port=0, pwd=None, db=None):
        '''
        Constructor

        Class variables
        ---------------
        host : str
            DB hostname
        user : str
            DB username
        port : str
            DB port
        pwd : str
            DB password
        db : str
            database to connect

        '''

        self.host = host
        self.user = user
        self.port = port
        self.pwd = pwd
        self.db = db

        db = pymysql.connect(host=self.host, user=self.user, passwd=self.pwd,
                             port=self.port, db=self.db)

        self.db = db

    def fetch_file_by_filename(self, name):
        '''
        Fetch a file by its basename

        Parameters
        ----------
        name : str
            Basename (without the path) in the DB

        Returns
        -------
        A File object

        '''

        query = "SELECT * FROM file WHERE name like %s"

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query, ['%' + name])
            result_set = cursor.fetchall()
            if not result_set:
                raise Exception("More than one file retrieved"
                                " from DB using filename: %s" % name)
            if not result_set:
                raise Exception("No file retrieved"
                                " from DB using filename: %s" % name)
            for row in result_set:
                return File(dbID=row['file_id'], path=row['name'], type=row['type'])
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the file: {}".format(err))

    def fetch_file_by_url(self, url):
        '''
        Fetch a file by its url (absolute path)

        Parameters
        ----------
        url : str
            Absolute url (absolute path) in the DB

        Returns
        -------
        A File object

        '''

        query = "SELECT * FROM file WHERE name = '{0}'".format(url)
        print(query)

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            result_set = cursor.fetchall()
            if len(result_set) > 1:
                raise Exception("More than one file"
                                " retrieved from DB using url: %s" % url)
            if not result_set:
                raise Exception("No file retrieved"
                                " from DB using url: %s" % url)
            for row in result_set:
                return File(dbID=row['file_id'], path=row['name'], type=row['type'])
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the file: {}".format(err))

    def fetch_file_by_id(self, id):
        '''
        Fetch a file by its internal dbID

        Parameters
        ----------
        id : int
            internal dbID for a file

        Returns
        -------
        A File object

        '''

        query = "SELECT * FROM file WHERE file_id= %d" % id

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            if cursor.rowcount == 0:
                raise Exception("No file retrieved "
                                "from DB using dbID: %d" % id)
            if cursor.rowcount > 1:
                raise Exception("More than one file "
                                "retrieved from DB using dbID: %d" % id)
            result_set = cursor.fetchall()
            for row in result_set:
                return File(dbID=row['file_id'], path=row['name'], type=row['type'])
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the file: {}".format(err))

    def fetch_files_by_type(self, type):
        '''
        Fetch files by their type

        Parameters
        ----------
        type : str
            File type

        Returns
        -------
        A list of File objects

        '''

        query = "SELECT * FROM file WHERE type= '%s'" % type

        list_of_files = []

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            if cursor.rowcount == 0:
                raise Exception("No file retrieved "
                                "from DB using type: %s" % type)
            result_set = cursor.fetchall()
            for row in result_set:
                list_of_files.append(File(dbID=row['file_id'], path=row['name'], type=row['type']))
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the file: {}".format(err))

        return list_of_files

    def fetch_collection_by_name(self, name):
        '''
        Fetch a collection by its name

        Parameters
        ----------
        name : str
            Collection name

        Returns
        -------
        A Collection object

        '''

        query = "SELECT * FROM collection WHERE name= '%s'" % name

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            result_set = cursor.fetchall()
            if len(result_set) > 1:
                raise Exception("More than one collection "
                                "retrieved from DB using name: %s" % name)
            if not result_set:
                raise Exception("No collection retrieved "
                                "from DB using name: %s" % name)
            for row in result_set:
                return Collection(dbID=row['collection_id'], name=row['name'],
                                  type=row['type'], table_name=row['table_name'])
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the collection: {}".format(err))

    def fetch_attribute_by_params(self, table_name, other_id, attribute_name):
        '''
        Fetch an Attribute objects using the following params:

        Parameters
        ----------

        table_name : str
        other_id : int
        attribute_name : str

        This function is useful if one wants to check the presence of a particular
        attribute in the database in order to prevent duplication errors

        Returns
        -------
        An Attribute object

        '''

        query = "SELECT * FROM attribute WHERE table_name= '%s' and other_id = "\
        "%d and attribute_name = '%s'" % (table_name, other_id, attribute_name)

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            result_set = cursor.fetchall()
            if len(result_set) > 1:
                raise Exception("More than one attribute retrieved "
                                "from DB using table_name-other_id-"
                                "attribute_name: %s-%d-%s" % (table_name,
                                                              other_id,
                                                              attribute_name))
            for row in result_set:
                return Attribute(table_name=row['table_name'], name=row['attribute_name'],
                                 value=row['attribute_value'], other_id=row['other_id'])
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the attributes: {}".format(err))

    def fetch_attributes_by_other_id(self, other_id):
        '''
        Fetch a list of Attribute objects using their 'other_id'

        Parameters
        ----------
        other_id : int

        Returns
        -------
        A list of Attributes objects

        '''

        attributes_list = []

        query = "SELECT * FROM attribute WHERE other_id= '%d'" % other_id

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            result_set = cursor.fetchall()
            if not result_set:
                raise Exception("No attribute "
                                "records retrieved "
                                "from DB using other_id: %d"\
                                % other_id)
            for row in result_set:
                attributes_list.append(Attribute(table_name=row['table_name'],
                                                 name=row['attribute_name'],
                                                 value=row['attribute_value'],
                                                 other_id=row['other_id']))
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the attributes: {}".format(err))

        return attributes_list

    def fetch_attributes_by_name(self, name):
        '''
        Fetch a list of Attribute objects using their 'attribute_name'

        Parameters
        ----------
        name : attribute_name in the DB

        Returns
        -------
        A list of Attributes objects

        '''

        attributes_list = []

        query = "SELECT * FROM attribute WHERE attribute_name = '%s'" % name

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            result_set = cursor.fetchall()
            if not result_set:
                raise Exception("No attribute records retrieved "
                                "from DB using attribute_name: %d" % name)
            for row in result_set:
                attributes_list.append(Attribute(table_name=row['table_name'],
                                                 name=row['attribute_name'],
                                                 value=row['attribute_value'],
                                                 other_id=row['other_id']))
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the attributes: {}".format(err))

        return attributes_list

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()

class File(object):
    '''
    Class to represent a file in the ReseqTrack DB
    '''

    def __init__(self, path, type, dbID=None, md5=None,
                 size=None, host_id=1, withdrawn=0, created=None):
        '''
        Create a File object

        Class members
        -------------
        path : str, Required
               Path to the file
        type : str, Required
               Type of the file
        dbID : int, Optional
               Internal db identifier
        name : str, Optional
               File name (without the path)
        md5 : str, Optional
              md5 of the file
        host_id : int, Optional
                  id of the host of this file. Default=1
        withdrawn : int (0 or 1), Optional
                    is the file withdrawn?. Default=0
        created : str representing the date, Optional
                  Date of creation of this file

        '''
        self.path = path
        self.type = type
        self.dbID = dbID

        #set the name class member
        self.name = os.path.split(path)[1]

        self.md5 = md5
        self.host_id = host_id
        self.withdrawn = withdrawn
        self.created = created
        self.size = size

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()

    def change_file_type(self, reseqdb, new_type):
        '''
        Update file type of this object in the DB.

        Parameters
        ----------
        reseqdb : MySQL db connection object  pointing to a reseqtrack db
        new_type : Str, required
                   New file type
        '''

        update_sql = "update file set type=\"%s\" where file_id=%d" % (new_type, self.dbID)

        try:
            cursor = reseqdb.db.cursor()
            # Execute the SQL command
            cursor.execute(update_sql)
            # Commit your changes in the database
            reseqdb.db.commit()
        except pymysql.Error as exc:
            print(exc[0], exc[1])
            # Rollback in case there is any error
            reseqdb.rollback()

    def calc_md5(self):
        '''
        Calculate the md5 sum of a file

        Returns
        -------
        The md5sum string
        '''
        command = "md5sum %s"% self.path

        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()

        md5sum = stdout.decode("utf-8").split('  ')[0]
        stderr = stderr.decode("utf-8")

        if stderr:
            print(stderr)

        return md5sum

    def store(self, reseqdb, do_md5=False):
        '''
        Store this File in the DB.

        Parameters
        ----------
        reseqdb : MySQL db connection object  pointing to a reseqtrack db
        do_md5 : Optional, Default=False
                 Calculate the md5 on the file

        '''

        if do_md5 is True:
            md5sum = self.calc_md5()
            self.md5 = md5sum

        if self.size is None:
            size = os.path.getsize(self.path)
            self.size = size

        if self.created is None:
            self.created = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        # insert new File in the 'file' table
        sql_insert_attr = "INSERT INTO file (file_id,name,md5,type,size,host_id,"\
                          "withdrawn,created) VALUES (NULL,'%s','%s','%s','%s',"\
                          "'%s','%s','%s')" % (self.path, self.md5, self.type,
                                               self.size, self.host_id,
                                               self.withdrawn, self.created)

        try:
            cursor = reseqdb.db.cursor()
            # Execute the SQL command
            cursor.execute(sql_insert_attr)
            # Commit your changes in the database
            reseqdb.db.commit()
        except pymysql.Error as e:
            print(e[0], e[1])
            # Rollback in case there is any error
            reseqdb.rollback()

    def move(self, reseqdb, newpath, do_md5=False):
        '''
        Move this File. It will move this file to a new location

        Parameters
        ----------
        reseqdb : MySQL db connection object  pointing to a reseqtrack db, Required
        newpath : str, Required
                  new location for file
        do_md5 : Optional, Default=False
                 Calculate the md5 on the file
        '''

        oldpath = self.path

        if do_md5 is True:
            md5sum = self.calc_md5()
            self.md5 = md5sum

        if self.size is None:
            size = os.path.getsize(self.path)
            self.size = size

        if self.created is None:
            self.created = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        self.path = newpath

        # insert new File in the 'file' table
        sql_insert_attr = "INSERT INTO file (file_id,name,md5,type,size,"\
        "host_id,withdrawn,created) VALUES (NULL,'%s','%s','%s','%s','%s','%s','%s')"\
         % (self.path, self.md5, self.type, self.size, self.host_id, self.withdrawn, self.created)

        try:
            cursor = reseqdb.db.cursor()
            # Execute the SQL command
            cursor.execute(sql_insert_attr)
            # Commit your changes in the database
            reseqdb.db.commit()
        except pymysql.Error as exc:
            print(exc[0], exc[1])
            # Rollback in case there is any error
            reseqdb.rollback()

        os.rename(oldpath, newpath)

class Collection(object):
    '''
    Class to represent a Collection of entities in the ReseqTrack DB
    '''

    def __init__(self, name, type, table_name, others_dbIDs=None, dbID=None):
        '''
        Create a Collection object

        Class variables
        ---------------
        name : str, Required
               Name of the collection
        type : str, Required
               Type of the collection
        others : list
            List of objects that will be part of the collection
        others_dbIDs : list
            List of dbIDs corresponding to the objects that will be part of the collection
        table_name : str
            This param controls to what table the objects in 'others'
            point to (i.e. Collection, File, Run)
        dbID : int, optional
            Internal dbID for this collection
        '''
        self.dbID = dbID
        self.name = name
        self.type = type
        self.table_name = table_name

        if others_dbIDs is None:
            others_dbIDs = []
        # convert to list
        if not isinstance(others_dbIDs, list):
            others_dbIDs = [others_dbIDs]
        self.others_dbIDs = others_dbIDs

        if table_name not in ['file', 'sample']:\
        raise Exception(("table_name %s is not valid."
                         "Valid table_names are file,sample")) % table_name

    def fetch_others_dbIDs(self, reseqdb):
        '''
        Fetch the dbIDs of the objects referenced by this collection.

        Parameters
        ----------
        reseqdb : MySQL db connection object  pointing to a reseqtrack db
        '''

        select_sql_query = "select other_id from collection_group where "\
        "collection_id = %d" % (self.dbID)

        try:
            cursor = reseqdb.db.cursor()
            cursor.execute(select_sql_query)
            if cursor.rowcount == 0:
                raise Exception("No entry exists in "
                                "table %s with id: %d"\
                                 % (self.table_name,
                                    self.dbID))
            other_ids = [row[0] for row in cursor.fetchall()]
            self.others_dbIDs = other_ids
            reseqdb.db.commit()
            cursor.close() 
            return other_ids
        except reseqdb.db.Error as err:
            print("Something went wrong in the select query: {}".format(err))

    def store(self, reseqdb):
        '''
        Store a collection in the DB. It will also update the 'collection_group' table accordingly

        Parameters
        ----------
        reseqdb : MySQL db connection object  pointing to a reseqtrack db

        '''

        # insert new Collection in the 'collection' table
        sql_insert_collection = "INSERT INTO collection (collection_id, name, type, table_name) "\
        "VALUES (NULL, '%s', '%s', '%s')" % (self.name, self.type, self.table_name)

        try:
            cursor = reseqdb.db.cursor()
            # Execute the SQL command
            cursor.execute(sql_insert_collection)
            # get id of the inserted collection
            self.dbID = cursor.lastrowid
            # Commit your changes in the database
            reseqdb.db.commit()
        except MySQLdb.Error as exc:
            print(exc[0], exc[1])
            # Rollback in case there is any error
            reseqdb.rollback()

        # Now, update also the 'collection_group' table
        for id in self.others_dbIDs:
            sql_insert_collection_group = "INSERT INTO collection_group (collection_id,other_id) "\
            "VALUES ('%d','%d')" % (self.dbID, id)
            try:
                cursor = reseqdb.db.cursor()
                # Execute the SQL command
                cursor.execute(sql_insert_collection_group)
                # Commit your changes in the database
                reseqdb.db.commit()
            except MySQLdb.Error as exc:
                print(exc[0], exc[1])
                # Rollback in case there is any error
                reseqdb.rollback()

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()

class Attribute(object):
    '''
    Class to represent an Attribute in the ReseqTrack DB
    '''

    def __init__(self, table_name, other_id, name, value, units=None):
        '''
        Create an Attribute object

        Class members
        ---------------
        table_name : str
            Attribute object points to this table_name
        other_id : int
            id number in table_name
        name : str
            Attribute name
        value : str or float
            Value assigned to this attribute name
        units : str, optional
        '''
        self.table_name = table_name
        self.other_id = other_id
        self.name = name
        if self.__numeric_type(value) == "int":
            self.value = int(value)
        elif self.__numeric_type(value) == "float":
            self.value = float(value)
        else:
            self.value = value
        if units is None:
            self.units = 'NULL'

        if table_name not in ['experiment', 'sample', 'run', 'file', 'collection']:
            raise Exception(("table_name %s is not valid. Valid "
                             "table_names are file,collection") % table_name)

    def __numeric_type(self, x):
        '''
        Check if the given value is either a str or number. If it is a number
        then check if it is float or int
        '''

        type = ""
        try:
            a = float(x)
            if a.is_integer() is True:
                type = "int"
            else:
                type = "float"
        except ValueError:
            type = "str"
        return type

    def store(self, reseqdb):
        '''
        Store this object in the DB.

        Parameters
        ----------
        reseqdb : MySQL db connection object  pointing to a reseqtrack db

        '''

        #check if object referenced by self.other_id exists in self.table_name
        select_sql_query = "select * from %s where %s_id = %d" % (self.table_name,
                                                                  self.table_name,
                                                                  self.other_id)

        try:
            cursor = reseqdb.db.cursor()
            cursor.execute(select_sql_query)
            result_set = cursor.fetchall()
            if not result_set:
                raise Exception("No entry exists in table "
                                "%s with id: %d" % (self.table_name,
                                                    self.other_id))
            cursor.close()
            reseqdb.db.commit()
        except reseqdb.db.Error as err:
            print("Something went wrong in the select query: {}".format(err))

        # insert new Attribute in the 'attribute' table
        sql_insert_attr = "INSERT INTO attribute (attribute_id,table_name,"\
        "other_id,attribute_name,attribute_value,attribute_units) VALUES (NULL,"\
        " '%s', '%d', '%s','%s','%s')" % (self.table_name, self.other_id,
                                          self.name, self.value, self.units)

        try:
            cursor = reseqdb.db.cursor()
            # Execute the SQL command
            cursor.execute(sql_insert_attr)
            # Commit your changes in the database
            reseqdb.db.commit()
        except pymysql.Error as e:
            print(e[0], e[1])
            # Rollback in case there is any error
            reseqdb.rollback()

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()
