## Minimal test of GSEAPY
## Dev version of GSEAPY, run in Python3 

import gseapy as gp
import numpy as np
from os.path import join 
import pandas as pd


tbl = pd.read_table(join('tracks', 'TLX3vsRAG-results_genes.txt'), index_col=0)

tbl = tbl[(tbl.padj < 0.05)].dropna()
names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           index_col=0,
                           header=0,
                           names=['GeneID', 'TransID', 'Gene_name'])

names = names.drop('TransID', axis=1).drop_duplicates()
names = names.loc[tbl.index]

tbn = pd.concat([names,tbl], axis=1)

tbn = tbn.drop(['baseMean', 'log2FoldChange', 'lfcSE', 'stat',  'pvalue',  'padj'], axis=1)

tbn['Gene_name'] = tbn['Gene_name'].str.upper()
tbn = tbn.reset_index(drop=True)

## ==== GSEAPY
classes = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']
g_set = 'Reactome_2013'
gs_res = gp.gsea(data=tbn, gene_sets = g_set, outdir='_test', cls = classes)





