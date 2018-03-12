#!/usr/bin/env python3

import numpy as np
import yaml
import os
from os.path import join 
import pybedtools as pb
import gffutils
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import matplotlib.pyplot as plt

import pandas as pd
import gseapy as gp


tbl = pd.read_table(join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt'), index_col=0)


tbl = tbl[(tbl.padj < 0.05)].dropna()

# === Load gene names 
names = pd.read_table("tracks/UCSC_mm9_transcripID_to_geneSymbol.sort.txt", 
                           index_col=0,
                           names=["Geneid", "NAME"])


names = names.loc[tbl.index]
assert names.shape[0] == tbl.shape[0]

tbl_n=names.join(tbl, how ='right')



# === Run GSEA
tbl_n['NAME']=tbl_n['NAME'].str.upper()
tbl_c = tbl_n.copy() 

tbl_c.index=tbl_n['NAME']


tbl_c = tbl_c.groupby(tbl_c.index).agg({'NAME': 'first',
                                'R2.RAG1W.RAG1':sum,
                                'RAGS.RAGZ':sum,
                                'RAGZ':sum,
                                'TLX3.1_1':sum,
                                'TLX3.1_5':sum,
                                'TLX3.1_P':sum,
                                'TAP':sum,
                                'TAP1B':sum,
                                'TAP2B':sum})

tbl_c = tbl_c[['NAME',
            'R2.RAG1W.RAG1',
            'RAGS.RAGZ',
            'RAGZ',
            'TLX3.1_1',
            'TLX3.1_5',
            'TLX3.1_P']]
            #~ 'TAP',
            #~ 'TAP1B',
            #~ 'TAP2B']]

classi = ['RAG','RAG','RAG','TLX3','TLX3','TLX3'] #,'TLX3','TLX3','TLX3']

gs_res = gp.gsea.call(data=tbl_c, 
                        gene_sets= 'tracks/GSEA_gene_sets/ENCODE-Histone_Mod2015_T_ALL.gmt', #gene_sets='KEGG_2016', 
                        cls=classi,
                        min_size=10,
                        max_size = 3000,
                        top_enrich = 250,
                        permutation_type='gene_set',#~ permutation_type='phenotype',
                        outdir='GSEA/ENCODE-Histone_Mod2015_T_ALL')


