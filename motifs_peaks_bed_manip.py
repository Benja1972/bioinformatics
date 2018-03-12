#!/usr/bin/env python3

import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt


tlx_peak = pb.BedTool('tracks/TLX3_TLX3_peaks_100.bed')
rnx = pb.BedTool('tracks/RUNX_motif.bed')
ets = pb.BedTool('tracks/ETS_motif.bed')
#enh_fantm = pb.BedTool('tracks/annot_tracks/enhancers_FANTOM.bed')



# in all in TLX-peaks
test1 = (tlx_peak + rnx + ets).sort().merge()
test2 = (ets + tlx_peak + rnx).sort().merge() 
test3 = (rnx + ets + tlx_peak).sort().merge()
test4 = (tlx_peak + ets + rnx).sort().merge()


#~ print(test1.count())
#~ print(test2.count())
#~ print(test3.count())
#~ print(test4.count())

test1.saveas('tracks/TLX3_peaks_RUNX_ETS.bed')

genes = test1.closest('tracks/genes.bed')
genes.saveas('tracks/TLX3_peaks_RUNX_ETS_genes.bed')



gs = genes.to_dataframe()
nm = list(gs['thickStart'].str.upper())
nm = list(set(nm))
gmt = {'TLX-RUNX-ETS':nm}



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
            'TLX3.1_P',
            'TAP',
            'TAP1B',
            'TAP2B']]

classi = ['WT','WT','WT','MUT','MUT','MUT','MUT','MUT','MUT']

gs_res = gp.gsea.call(data=tbl_c, 
                       gene_sets= 'tracks/TLX3_peaks_RUNX_ETS_genes.gmt',#gene_sets='KEGG_2016', 
                      cls=classi,
                      permutation_type='gene_set',#~ permutation_type='phenotype',
                      outdir='gsea_TLX-RUNX-ETS')

# === Pictures

gsea_results = gs_res.reset_index().sort_values('fdr',axis=0,ascending=True)

with plt.style.context('ggplot'):
    gsea_results.head(40).plot.bar(y='fdr',x='Term', figsize=(12, 6),fontsize=12)

plt.show() 








#~ print(genes.head())

#~ print(tlx_peak.count())test

#~ st89_peaks.saveas('st89_peaks.bed')
#~ st3_peaks.saveas('st3_peaks.bed')






#~ st6_tlx =  pb.BedTool('TLX3_14_E6_sorted.bed')
#~ st6_rag =  pb.BedTool('RAG_14_E6_sorted.bed')

#~ # Regions in TLX3 but not in RAG
#~ st6_tlx_not_rag = st6_tlx - st6_rag
#~ st6_rag_not_tlx = st6_rag - st6_tlx

#~ #st6_tlx_not_rag.saveas('state6_tlx_not_rag.bed', trackline="track name='State 6 TLX3 -- not RAG' color=128,0,0")
#~ st6_tlx_not_rag.saveas('state6_tlx_not_rag.bed')

#~ print(st6_tlx_not_rag.count())
#~ print(st6_rag_not_tlx.count())
