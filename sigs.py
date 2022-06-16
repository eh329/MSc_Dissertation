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

    def sigs_pres():
        if '_sigs_pres' not in MutSig.__dict__.keys():
            #MutSig._sigs_pres = pd.DataFrame.from_csv(Data.directory+'signatures_present_primary_cancers.csv')
            MutSig._sigs_pres = pd.DataFrame.from_csv(Data.directory+'sigs_pres.csv')
        return MutSig._sigs_pres
    
    cancer_sigs = lambda cancer: list(MutSig.sigs_pres()[cancer][MutSig.sigs_pres()[cancer]==1].index)
    
    def bump_up(short_array, index):

        a0 = pd.Series(short_array)
        a0.index = index

        extra = pd.Series([0 for i in range(30)])
        extra.index = range(1,31)
        a1 = a0+extra

        return a1.fillna(0).as_matrix()
    
    convert_sigs_to_pres = lambda mylist: [int(i in mylist) for i in range(30)]
    
    def __init__(self,s, max_iter = 50, cancer = '', sigs_to_check =[], ticks_to_show=[]):
        self.ticks_to_show = ticks_to_show
        self.sigs_to_check = sigs_to_check
        self.cancer = cancer
        if sigs_to_check!=[]:
            pres = MutSig.convert_sigs_to_pres(sigs_to_check)
            self.sigs = np.array([MutSig.sigs()[i] for i in range(30) if pres[i]==1])

        elif cancer=='':
            self.sigs = MutSig.sigs()
        else:
            pres = list(MutSig.sigs_pres()[cancer])
            self.sigs = np.array([MutSig.sigs()[i] for i in range(30) if pres[i]==1])

        Ms = np.matmul(self.sigs,s)
        MTM  = np.matmul(self.sigs,self.sigs.T)
        w = np.array([np.random.ranf() for i in range(self.sigs.shape[0])])
        ws = [w]
        self.changes =[1] 
        i=0
        while i<max_iter and self.changes[-1]>0.00001:
            i+=1
            w = w*Ms/np.matmul(MTM,w)
            ws.append(w)
            self.changes.append(sum((ws[-1]/ws[-1].sum()-ws[-2]/ws[-2].sum())**2))
        self.decompose = ws[-1]
        self.reconstruction = np.matmul(self.sigs.T,self.decompose)
        self.error =0.5*((self.reconstruction/self.reconstruction.sum()-s/s.sum())**2).sum()
        self.ws = ws
        
    def show(self):
        if self.cancer !='':
            self.long = MutSig.bump_up(self.decompose,MutSig.cancer_sigs(self.cancer))
        else:
            self.long = self.decompose
        if self.sigs_to_check!=[]:
            the_range=self.sigs_to_check
        else:
            the_range=range(1,31)
        if self.ticks_to_show!=[]:
            ticks_to_show = self.ticks_to_show
        else:
            ticks_to_show = the_range
        plt.bar(the_range,self.long)
        x=plt.xticks(ticks_to_show,ticks_to_show)


        plt.show()

# we can reduce Data.sigs to any set of samples and take the mean
def signature(samples):
    '''average signature as percentage'''
    s = Data.sigs().loc[list(samples)].mean()
    return s/s.sum()*100

def other_samples(cancer,samples):
    '''Given a list of samples in a particular cancer, will return the list of all
    samples that are in the cancer but not in the list'''
    all_s = Data.samples(cancer)
    return list(set(all_s)-set(samples))

main_cancer  = lambda s0:list(Data.all_mutations().loc[
    Data.all_mutations()['Sample name']==s0[0]]['Primary site'])[0]

def breakdown(sig = '',samples = '', show = False, others = False, reduced = False, cancer='',sigs_to_check=[], ticks_to_show=[]): 
    ''' Input: either provide a signature or a list of samples using the keyword arguments sig,samples
        The function loads all the signatures and picks out those for the samples in the list or uses the signature given.
        Finds the mean and then decomposes them using non-negative matrix factorisation.
        Using the keyword show - you can display this decomposition.
        If you use the keyword others = True then it automatically finds the samples in your cancer
        that are not in your list and calculates the breakdown for them.
        
        You can add your cancer manually by using the keyword cancer, and you can also specify which signatures to check
        by using the keyword sigs_to_check. You can also specify which ticks to show using keyword ticks_to_show
        '''

    if sig !='':
        sg=sig
    
    elif others == False:
        if cancer=='':
            cancer= main_cancer(samples)
        sg = signature(samples)

    else:
        if cancer=='':
            cancer = main_cancer(samples)
        s = other_samples(cancer,samples)
        sg = signature(list(s))
    
    if reduced:
        m = MutSig(sg,cancer = cancer,sigs_to_check=sigs_to_check,ticks_to_show=ticks_to_show)
    else:
        m = MutSig(sg,sigs_to_check=sigs_to_check,ticks_to_show=ticks_to_show)
    
    if show:
        m.show()
    return m.decompose

def change(samples, reduced = False,sigs_to_check = [],ticks_to_show=[]):
    '''Shows the breakdown for the samples and for the list of samples in your cancer that are
    not in your list. Returns a dictionary with keys as follows:
    breakdown,others_breakdown,difference. You can also specify a list of ticks to show using the keyword ticks_to_show'''
    
    a = breakdown(samples = samples,others = False,show = True,
                  reduced = reduced,sigs_to_check = sigs_to_check,ticks_to_show=ticks_to_show )
    b = breakdown(samples = samples,others = True,show = True,
                  reduced = reduced,sigs_to_check =sigs_to_check,ticks_to_show=ticks_to_show)
    c = pd.Series((a-b))
    c.index = MutSig.cancer_sigs(main_cancer(samples))
    d = c.sort_values(ascending = False)
    return {'breakdown':a,
            'others_breakdown':b,
            'difference':d}

def find_main_sigs(sample):
    '''Given a specific sample, the function breaks it down 10 times and 
    each time identifies those signatures that make up more than 20% of the overall signal. It then
    returns as a list just those signatures that turn up at least half the time.'''
    main=[]
    for i in range(10):
        b = pd.Series(breakdown(samples = [sample]))
        b.index = range(1,31)
        main+=list(b[b>20].index)
    mcv = pd.Series(main).value_counts()
    sigs = list(mcv[mcv>4].index)
    return sigs

cancer_samples = lambda p: Data.all_mutations()[Data.all_mutations()['Primary site']==p]['Sample name']

def bad_samples(gene,cancer):
    a = pd.concat([Data.all_mutations().loc[Data.all_mutations()['Primary site']==cancer],
                   Data.all_mutations().loc[Data.all_mutations()['Site subtype 1']==cancer]])
    b=a.loc[gene].drop_duplicates()
    c = b.loc[b['FATHMM prediction']=='PATHOGENIC']
    return list(c['Sample name'].drop_duplicates())
