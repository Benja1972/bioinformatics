#!/usr/bin/env python3

import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt


# === Load table 
tbl = pd.read_table(join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt'), index_col=0)

tbl = tbl[(tbl.padj < 0.05)].dropna()

# === Load gene names 
names = pd.read_table("tracks/UCSC_mm9_transcripID_to_geneSymbol.sort.txt", 
                           index_col=0,
                           names=["Geneid", "NAME"])


names = names.loc[tbl.index]
assert names.shape[0] == tbl.shape[0]


tbl_raw = tbl[['R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ',
               'TLX3.1_1','TLX3.1_5','TLX3.1_P',
               'TAP','TAP1B','TAP2B']]


tbl_n=names.join(tbl_raw, how ='right')

tbl_n['NAME']=tbl_n['NAME'].str.upper()

# === Run GSEA
#~ classi = ['WT','WT','WT','MUT','MUT','MUT','MUT','MUT','MUT']
classi = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']

#tbl_g = tbl_n.copy()#head(23737)
#tbl_g = tbl_n.copy()


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

g_set='KEGG_2015'
out_dir = 'GSEA/gsea_'+g_set
gs_res = gp.gsea.call(data=tbl_c, 
                      gene_sets=g_set,#gene_sets= 'tracks/GSEA_gene_sets/c2.cp.kegg.v6.0.symbols.gmt',#gene_sets='KEGG_2016', 
                      cls=classi,
                      graph_num = 100,
                      #~ permutation_type='phenotype',
                      outdir=out_dir)

# === Pictures

gsea_results = gs_res.reset_index().sort_values('fdr',axis=0,ascending=True)

with plt.style.context('ggplot'):
    gsea_results.head(40).plot.bar(y='fdr',x='Term', figsize=(12, 6),fontsize=12)

plt.savefig(out_dir+'/'+g_set+'.pdf', format="pdf")
plt.show() 



#~ with plt.style.context('ggplot'):
    #~ gsea_results = gs_res.reset_index()  
    #~ gsea_results.head(5).plot.barh(y='fdr',x='Term',fontsize=12)

   


