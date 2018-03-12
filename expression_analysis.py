
# coding: utf-8

# In[2]:

import numpy as np
import yaml
import os
from os.path import join 
import pybedtools as pb
import gffutils
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import matplotlib.pyplot as plt
import metaseq
import pandas as pd
from metaseq.results_table import ResultsTable, DESeq2Results


# In[3]:

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12


# In[4]:

tbl_f = 'TLX3vsRAG_DESeq2-results.txt'

tbl = ResultsTable(tbl_f, import_kwargs=dict(index_col=0))
tbl.data['TLX3-mean'] = (tbl.data['TLX3.1_1']+tbl.data['TLX3.1_5']+tbl.data['TLX3.1_P'])/3.
tbl.data['RAG-mean'] = (tbl.data['R2.RAG1W.RAG1']+tbl.data['RAGS.RAGZ']+tbl.data['RAGZ'])/3.


# In[84]:

list ={
    'Samsn1-003':'ENSMUST00000114239',
    'Samsn1-001':'ENSMUST00000114240',
    'Samsn1-004':'ENSMUST00000139691',
    'Samsn1-007':'ENSMUST00000150248',
    'Samsn1-006':'ENSMUST00000147447',
    '1700041M19Rik-003':'ENSMUST00000208105',
    '1700041M19Rik-001':'ENSMUST00000190908',
    '1700041M19Rik-004':'ENSMUST00000207519',
    '1700041M19Rik-002':'ENSMUST00000190317',
    '1700041M19Rik-005':'ENSMUST00000207924',
    'Nrip1-001':'ENSMUST00000121927',
    'Nrip1-003':'ENSMUST00000054178',
    'Nrip1-006':'ENSMUST00000140483',
    'Nrip1-005':'ENSMUST00000141182',
    'Nrip1-004':'ENSMUST00000128871',
    'Nrip1-002':'ENSMUST00000145649',
    'D16Ertd472e-001':'ENSMUST00000114220',
    'D16Ertd472e-002':'ENSMUST00000114219',
    'D16Ertd472e-003':'ENSMUST00000114218',
    'Cxadr-001':'ENSMUST00000023572',
    'Cxadr-002':'ENSMUST00000114229',
    'mir99a':'ENSMUST00000083596',
    'mir125b':'ENSMUST00000083538',
    'Usp25':'ENSMUST00000023580',
    'Mir99ahg-001':'ENSMUST00000193546',
    'Mir99ahg-002':'ENSMUST00000179255',
    'Btg3-001':'ENSMUST00000023570',
    'Btg3-002':'ENSMUST00000148124'
    }


d = {'Gene name':list.keys(), 
     'Gene':list.values()}


# In[85]:

expr_df = pd.DataFrame(data=d, index=d['Gene'])
expr_df.index.name = 'Gene'

#expr_df


# In[86]:

sel_df = tbl.data.ix[expr_df['Gene']][['log2FoldChange', 'TLX3-mean', 'RAG-mean']]


# In[90]:

tot_df = pd.concat([sel_df,expr_df.drop(['Gene'], axis=1)], axis=1)
tot_df.sort_values(by='Gene name')


# In[ ]:



