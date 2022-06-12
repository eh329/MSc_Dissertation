import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from random import randint


def help():
    print('''

To use this module, you must save a copy of the following files in a directory on their own:
'mutational_fingerprints_for sample.csv'
'signatures_probabilities.txt'
'CosmicGenomeScreensMutantExport.tsv'
'sigs_pres.csv'
 
Once that is done go through the signatures.py script and change directory = '/home/s/sk/skw24/Fingerprints/' to point
to where you have saved the files.

Key functions: 
    signature(samples):
        average signature as percentage
    other_samples(cancer,samples):
        Given a list of samples in a particular cancer, will return the list of all
        samples that are in the cancer but not in the list
    main_cancer(samples):
        the cancer corresponding to the samples
    breakdown(sig = '',samples = '', show = False, others = False, reduced = False,cancer='',sigs_to_check=[]):   
        Input: 
        Provide preferably a list of samples using the keyword argument samples
        The function loads all the signatures and picks out those for the samples in the list or uses the signature given.
        Finds the mean and then decomposes them using non-negative matrix factorisation.
        Using the keyword show - you can display this decomposition.
        If you use the keyword others = True 
        then it automatically finds the samples in your cancer
        that are not in your list and calculates the breakdown for them
        reduced uses just those signatures that are known to exist in your cancer - and also 
        includes sig12 because it may be connected with TOP2A.
        Alternatively you can add your cancer manually by using the keyword cancer, 
        and you can also specify which signatures to check
        by using the keyword sigs_to_check.
    
    Data.info also gives you information on the different signatures
''')

class Data:
    directory = 'D:/Project/cosmic2/'
    
    def info(n):
        '''Information in the form of a dictionary on the nth signature'''
        if '_info' not in Data.__dict__.keys():
            with open(Data.directory+'signature_descriptions.json','r') as f:
                Data._info = json.load(f)   
        return Data._info[str(n)]
    
    _samples = {}
    
    def sigs():
        if '_sigs' not in Data.__dict__.keys():
            Data._sigs= pd.DataFrame.from_csv(Data.directory + 'mutational_fingerprints_for sample.csv')
        return Data._sigs

    def all_mutations():
        """Cosmic"""
        if '_all_mutations' not in Data.__dict__.keys():
            print('loading all mutations... this may take a while', end = '')
            Data._all_mutations = pd.DataFrame.from_csv(Data.directory+'CosmicGenomeScreensMutantExport.tsv',sep = '\t')
            print('...loaded')
        return Data._all_mutations
    
    def primaries():
        if '_primaries' not in Data.__dict__.keys():
            Data._primaries = set(Data.all_mutations()['Primary site'])
        return Data._primaries
    
    def others():
       if '_others' not in Data.__dict__.keys():
           Data._others = set(Data.all_mutations()['Site subtype 1'])
       return Data._others

    def samples(cancer):
        prim = Data.all_mutations().loc[
            Data.all_mutations()['Primary site']==cancer]
        sec = Data.all_mutations().loc[
            Data.all_mutations()['Site subtype 1']==cancer]
        return list(pd.concat([prim,sec])['Sample name'])

    
class MutSig:
    '''Given a pandas series or np array of 96 mutational frequencies, decomposes them into the 30 mutation signatures
    INPUT:
    np.array of shape (96)
    optional max_iter  - default maximum iteration is 50

    FUNCTIONS:

        show()
        bar graph of mutational signature decomposition

    ATTRIBUTES:

        decompose - np.array of shape 30 = decomposition into mutational signatures
        reconstruction - reconstruction of original vector using decomposition
        error  - Frobenius error


    '''    
    
     def sigs():
        if '_sigs' not in MutSig.__dict__.keys(): 
            sigs = pd.DataFrame.from_csv(Data.directory+'signatures_probabilities.txt', sep = '\t')
            cols = sigs.columns
            MutSig._sigs = sigs[cols[2:32]].T.as_matrix()
            
        return MutSig._sigs
