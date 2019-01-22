'''
Created on 27 Feb 2017

@author: ernesto
'''
import pandas as pd
import numpy as np
import os
import pdb
import pickle
import gc
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import Imputer

class MLclassifier(object):
    '''
    Class to filter a VCF using a supervised machine learning binary classifier. This class
    relies on a truth set (for example the GIAB consensus call set) that can be used to train
    the model
    '''
    def __init__(self, fitted_model=None, score=None, bcftools_folder=None):
        '''
        Constructor

        Parameters
        ----------
        fitted_model : filename, optional
                       Path to file containing the serialized fitted model
        score : float, optional
                Score of the fitted model
        bcftools_folder : str, optional
                          Path to folder containing the bcftools binary
        '''
        self.fitted_model = fitted_model
        self.bcftools_folder = bcftools_folder

    def __chunk_preprocessing(self, chunk, is_valid=None):
        '''
        Private method to preprocess each Data frame chunk

        Parameters
        ----------
        chunk : dataframe
        is_valid : int, optional
                    Value to add in the 'is_valid' column

        Returns
        -------
        A preprocessed dataframe for a particular chunk
        
        '''

        # remove NA values
        chunk.dropna(inplace=True)

        # normalization of the different features
        feature_names=chunk.columns
        std_scale = preprocessing.StandardScaler().fit(chunk[feature_names])
        std_array = std_scale.transform(chunk[feature_names])

        # 'preprocessing' returns a NumPy array, so we need to transform to a Pandas data frame:
        aDF_std=pd.DataFrame(data=std_array,columns=feature_names)
        
        # Now, let's add the column with the status of the variant( is_valid=0 or is_valid=1) to the normalized data frame
        aDF_std.insert(loc=0, column='is_valid', value=is_valid)

        return aDF_std

    def train(self, tp_annotations, fp_annotations, outprefix, test_size=0.25):
        '''
        Function to train the binary classifier using a gold standart call set

        Parameters
        ----------
        tp_annotations : filename
                         Path to file with the variant annotations derived from the call set with the True positives
        fp_annotations : filename
                         Path to file with the variant annotations derived from the call set with the False positives
        outprefix : str
                    String used as the prefix for the fitted model
        test_size : float
                    Fraction of the initial call set that will be used for assessing the model
                    Default: 0.25
        
        Returns
        -------
        filename
                 Path to serialized fitted model
        '''
        # check if tp_annotations and fp_annotations have the same columns and get columns names
        DF_TP_columns=pd.read_csv(tp_annotations, sep="\t", na_values=['.'], nrows=1).columns
        DF_FP_columns=pd.read_csv(fp_annotations, sep="\t", na_values=['.'], nrows=1).columns

        if DF_TP_columns.equals((DF_FP_columns)) is False: raise("Indices in the passed dataframes are not equal")

        # create 2 dataframes from tsv files skipping the 2 first columns, as it is assumed that the 1st is 'chr' and 2nd is 'pos'
        DF_TP=pd.read_csv(tp_annotations, sep="\t", na_values=['.'], usecols=[i for i in range(2, len(DF_TP_columns))])
        DF_FP=pd.read_csv(fp_annotations, sep="\t", na_values=['.'], usecols=[i for i in range(2, len(DF_FP_columns))])

        #assign outcome=1 if TP and 0 if FP
        DF_TP=DF_TP.assign(is_valid=1)
        DF_FP=DF_FP.assign(is_valid=0)

        #now, we combine the 2 dataframes
        frames = [DF_TP,DF_FP]
        DF = pd.concat(frames)

        #we have several columns with NA values, we will impute the missing values with the median
        imputer = Imputer(strategy="median")
        imputer.fit(DF)
        X = imputer.transform(DF)

        #transforming back to a DF
        DF_tr = pd.DataFrame(X, columns=DF.columns)
        
        #Normalization of the different features
        feature_names=DF_tr.columns.drop(['is_valid'])
        std_scale = preprocessing.StandardScaler().fit(DF_tr[feature_names])
        std_array = std_scale.transform(DF_tr[feature_names])

        aDF_std=pd.DataFrame(data=std_array,columns=feature_names)
        aDF_std.insert(loc=0, column='is_valid', value=DF_tr['is_valid'].values)

        #now, we fit a ML model
        predictors=aDF_std[feature_names]
        outcome=aDF_std[['is_valid']]

        x_train, x_test, y_train, y_test = train_test_split(predictors, outcome, test_size=0.25)
        logisticRegr = LogisticRegression(verbose=1)
        logisticRegr.fit(x_train, y_train)

        #and now we assess the fitted model performance
        score = logisticRegr.score(x_test, y_test)
        print("Score for the logistic regression fitted model is: {0}".format(score))
        self.score=score

        # Model persistence
        outfile = outprefix+".sav"
        pickle.dump(logisticRegr, open(outfile, 'wb'))

        self.fitted_model=outfile

        return outfile

    def predict(self, outprefix, annotation_f, filter_label='MLFILT', cutoff=0.8):
        '''
        Function to apply a serialized logistic regression model on a file containing the annotations for each site
        and to predict if the variant is real

        Note. Sites with missing annotations will be imputed with the median of that annotation

        Parameters
        ----------
        outprefix: str
                   String used as the prefix for the fitted model
        annotation_f: filename
                      Path to file with the sites and annotations that will be classified
        filter_label: str, Optional
                      String with the label used for FILTER that can be used with 'bcftools annotate'
                      in order to annotate the VCF file. Default='MLFILT'
        cutoff : float, Optional
                 Cutoff value used for deciding if a variant is a TP. Sites with a prob_1<cutoff will
                 be considered as FP and the FILTER column will be the string in the 'filter_label' option.
                 If prob_1>=cutoff then the FILTER column will be PASS. Default=0.8

        Returns
        -------
        filename
                 Path to table file with predictions
        '''
        imputer = Imputer(strategy="median")

        # load the serialized model
        loaded_model = pickle.load(open(self.fitted_model,'rb'))

        outfile='{0}.tsv'.format(outprefix)

        DF_all=pd.read_csv(annotation_f,sep="\t",na_values=['.'],index_col=False)
        DF_all_num = DF_all.drop("# [1]CHROM", axis=1)

        # we will impute the missing values by the median
        imputer.fit(DF_all_num)
        X = imputer.transform(DF_all_num)
        DF_all_imp = pd.DataFrame(X, columns=DF_all.columns.drop(['# [1]CHROM']))

        # normalization
        feature_names=DF_all_imp.columns.drop(['[2]POS'])

        std_scale = preprocessing.StandardScaler().fit(DF_all_imp[feature_names])
        std_array = std_scale.transform(DF_all_imp[feature_names])

        # Now, we can calculate the probabilities for each of the categories of the dependent variable:
        predictions_probs = loaded_model.predict_proba(std_array)

        #decide if it is a TP or a FP
        filter_outcome=[]
        for i in predictions_probs[:,1]:
            if i >= cutoff: 
                filter_outcome.append('PASS')
            else:
                filter_outcome.append(filter_label)

        final_df = pd.DataFrame({
            '#CHR': DF_all['# [1]CHROM'],
            'POS': DF_all['[2]POS'].astype(int),
            'FILTER': filter_outcome,
            'prob_TP': [ round(elem, 2) for elem in predictions_probs[:,1]]})
        # change order of columns
        final_df=final_df[['#CHR','POS','FILTER','prob_TP']]

        final_df.to_csv(outfile, sep='\t', header=True, index= False)

        return outfile
