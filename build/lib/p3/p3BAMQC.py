'''
Created on 27 Jan 2017

@author: ernesto
'''
import pandas as pd


class p3BAMQC(object):
    '''
    Class representing a spreadsheet located at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx
    containing information on the BAM QC done for the p3
    '''
    
    def __init__(self,filepath):
        '''
        Constructor

        Parameters
        ----------
        
        filepath: str
             Path to the spreadsheet
        book: ExcelFile object
        '''
        
        self.filepath=filepath
        # Import the excel file and call it xls_file
        xls_file = pd.ExcelFile(filepath)
        self.book=xls_file
        
    def get_final_qc_results(self,group):
        '''
        Method to get the sheet corresponding to 'Final QC Results'
        
        Parameters
        ----------
        group: str, 'low coverage' or 'exome'
                    Get the data for the low coverage or exome worksheet
        
        Returns
        -------
        A data.frame object 
        '''
        
        sheet = self.book.parse('Final QC Results',skiprows=1,index_col=[0,1])
        
        new_column_names=['VerifyBam_Omni_Free','VerifyBam_Affy_Free','VerifyBam_Omni_Chip','VerifyBam_Affy_Chip','Indel_Ratio','Passed_QC']
        
        df=""
        
        if group=="low coverage":
            df=sheet.iloc[:,6:12]
        elif group=="exome":
            df=sheet.iloc[:,0:6]
        
        df.columns=new_column_names
        return df           