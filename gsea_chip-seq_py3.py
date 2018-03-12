#!/usr/bin/env python2

import numpy as np
import yaml
import os
from os.path import join 
import pybedtools as pb
import gffutils
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import matplotlib.pyplot as plt
#~ import metaseq
import pandas as pd
import gseapy as gp

#~ from metaseq.results_table import ResultsTable, DESeq2Results


# Paths
#tracks/annot_tracks/references
#~ ref_dir = 'tracks/annot_tracks/referencesmm9'
#~ data_dir = 'tracks'
#~ rel_path = "/home/sergio/media"


# Files list

#~ file_t = open(join(data_dir,"TLX3_list.yaml"), "r")
#~ file_r = open(join(data_dir,"RAG_list.yaml"), "r")

#~ tlx_lst = yaml.load(file_t)
#~ rag_lst = yaml.load(file_r)


# Martix of signals
#~ features3K, arrays3K = metaseq.persistence.load_features_and_arrays(prefix='tlx-list3K')
#~ tlx3_mean=arrays3K['tlx-TLX3'].mean(axis=1)


tbl = pd.read_table(join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt'), index_col=0)

#~ tbl = tbl.reindex_to(features3K, attribute='transcript_id')
#~ tbl.data['tlx3_ChiP-seq_3Kb'] = tlx3_mean

tbl = tbl[(tbl.padj < 0.05)].dropna()

# === Load gene names 
names = pd.read_table("tracks/UCSC_mm9_transcripID_to_geneSymbol.sort.txt", 
                           index_col=0,
                           names=["Geneid", "NAME"])


names = names.loc[tbl.index]
assert names.shape[0] == tbl.shape[0]

tbl_n=names.join(tbl, how ='right')

#~ tbl_s = tbl_n.sort_values('tlx3_ChiP-seq_3Kb', axis=0, ascending=False)




#~ gm = ['ChiP_seq_TLX3_TLX3top500', 'ChiP_seq_TLX3_TLX3top500']  + list(tbl_s['NAME'].head(500).dropna().str.upper())

#~ with open('tracks/ChiP_seq_TLX3_TLX3top500.gmt', 'w') as fp:
    #~ fp.write("\t".join(gm))


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
                        gene_sets= 'tracks/ChiP_seq_TLX3_TLX3top500.gmt', #gene_sets='KEGG_2016', 
                        cls=classi,
                        max_size = 2000,
                        top_enrich = 50,
                        permutation_type='gene_set',#~ permutation_type='phenotype',
                        outdir='GSEA/gsea_ChiP_seq_TLX3_TLX3top500')


