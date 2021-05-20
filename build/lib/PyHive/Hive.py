'''
Created on 01 Aug 2019

@author: ernesto lowy
'''

import pymysql
import pdb

class HiveDB(object):
    '''
    Class representing a Hive MySQL DB
    '''

    def __init__(self, host=None, user=None, port=0, pwd=None, db=None):
        '''
        Constructor

        Parameters
        ----------
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

    def fetch_jobs_by_analysis_id(self,id,limit=None):
        '''
        Fetch a list of jobs by its analysis_id

        Parameters
        ----------
        id: int
            Analysis_data_id used in the job table
        limit: int
               Limit the number of jobs returned
               Default=None

        Returns
        -------
        List of Job objects        
        '''

        query = "SELECT * FROM job WHERE analysis_id= %d " % id

        if limit is not None:
            query+= "LIMIT {0}".format(limit)

        jobs=[]

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            if cursor.rowcount == 0:
                raise Exception("No Job retrieved "
                                "from DB using analysis_id: %d" % id)
            result_set = cursor.fetchall()
            for row in result_set:
                jobs.append(Job(job_id=row['job_id'], analysis_id=row['analysis_id'],
                                input_id=row['input_id'], prev_job_id=row['prev_job_id'],
                                status=row['status']))
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the job: {}".format(err))
        
        return jobs

    def fetch_job_byid(self,id):
        '''
        Fetch a certain job by ids id

        Parameters
        ----------
        id: int
            Job id

        Returns
        -------
        Job object
        '''
        query = "SELECT * FROM job WHERE job_id= %d" % id

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            if cursor.rowcount == 0:
                raise Exception("No Job retrieved "
                                "from DB using job_id: %d" % id)
            if cursor.rowcount > 1:
                raise Exception("More than 1 Job retrieved "
                                "from DB using job_id: %d" % id)
            result_set = cursor.fetchall()
            for row in result_set:
                return Job(job_id=row['job_id'], analysis_id=row['analysis_id'],
                           input_id=row['input_id'], prev_job_id=row['prev_job_id'],
                           status=row['status'])
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the job: {}".format(err))

    def fetch_analysisdata_byid(self,id):
        '''
        Fetch a certain analysis_data by its id

        Parameters
        ----------
        id: int
            Analysis data id

        Returns
        -------
        AnalysisData object
        '''

        query = "SELECT * FROM analysis_data WHERE analysis_data_id= %d" % id

        try:
            cursor = self.db.cursor(pymysql.cursors.DictCursor)
            cursor.execute(query)
            if cursor.rowcount == 0:
                raise Exception("No analysis_data retrieved "
                                "from DB using analysis_data_id: %d" % id)
            result_set = cursor.fetchall()
            for row in result_set:
                return AnalysisData(analysis_data_id=row['analysis_data_id'], md5sum=row['md5sum'],
                        data=row['data'])
            cursor.close()
            self.db.commit()
        except self.db.Error as err:
            print("Something went wrong when getting the job: {}".format(err))

        return jobs
    
class AnalysisData(object):
    '''
    Class representing an AnalysisData in the Hive DB
    '''
    def __init__(self, analysis_data_id, md5sum, data):
        '''
        Constructor

        Class variables
        ---------------
        analysis_data_id: int
                          Analysis_data_id
        md5sum: str
                MD5sum 
        data: str
        '''
        self.analysis_data_id = analysis_data_id
        self.md5sum = md5sum
        self.data = data

class Job(object):
    '''
    Class representing a Job in the Hive DB
    '''
    def __init__(self, job_id, analysis_id, input_id, prev_job_id, status):
        '''
        Constructor

        Class variables
        ---------------
        job_id: int
                Job id
        analysis_id: int
                     Analysis id
        input_id: str
                  String representing the input id
        prev_job_id: int
                     Prev job id
        status: str
                Status of job
        '''

        self.job_id = job_id
        self.analysis_id = analysis_id
        self.input_id = input_id
        self.prev_job_id = prev_job_id
        self.status = status
