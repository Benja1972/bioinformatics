#!/usr/bin/env python3

import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

#~ enh_fantm = pb.BedTool('tracks/annot_tracks/enhancers_FANTOM.bed')
#~ enh_chrhmm = pb.BedTool('tracks/TLX3_6_FE_E4_sorted.bed')
#~ tss_3kb = pb.BedTool('tracks/annot_tracks/references/mm9/mm9_tsses-3kb.gtf')

prc2 = pb.BedTool('tracks/TLX3_15_PRC2_E13E14.bed').merge(d=100)
tss2Kb = pb.BedTool('tracks/annot_tracks/RefSeqTSS2kb.mm9.bed')
tlx_peak = pb.BedTool('tracks/TLX3_TLX3_peaks_100.bed')


tss2kb_prc2 = (tss2Kb+prc2).sort().merge(d=50)
tss2kb_prc2_tlx = (tss2Kb+prc2+tlx_peak).sort()





tss2kb_prc2.saveas('tracks/TSS2Kb_PRC2.bed')
tss2kb_prc2_tlx.saveas('tracks/TSS2Kb_PRC2_tlx.bed')



# in all in TLX-peaks
#~ test1 = (tlx_peak + enh_fantm).sort()#.merge()
#~ test2 = (enh_fantm + tlx_peak).sort()#.merge() 
#~ test3 = (tlx_peak + enh_fantm + enh_chrhmm).sort()#.merge()
#~ tlx_tss_3kb = (tlx_peak+tss_3kb).sort()
#~ tlx_chrhmm = (tlx_peak+enh_chrhmm - tlx_tss_3kb).sort()





#~ tlx_tss_3kb_gt = pd.read_table('results01/TLX3_TSS3Kb-gene.txt', header=1, names=['Genes','Regions'])
#~ tlx_chrhmm_gt = pd.read_table('results01/TLX3_chrhmm-gene.txt', header=1, names=['Genes','Regions'])

#~ tlx_tss_3kb_gn = list(tlx_tss_3kb_gt['Genes'].str.upper())
#~ tlx_chrhmm_gn = list(tlx_chrhmm_gt['Genes'].str.upper())

#~ def intersect(a, b):
    #~ return list(set(a) & set(b))

#~ tss3k_chrhmm_gn = intersect(tlx_tss_3kb_gn, tlx_chrhmm_gn)
#~ tss3k_chrhmm_gm = ['TLX_TSS3Kb_ChHMM', 'TLX_TSS3Kb_ChHMM']  + tss3k_chrhmm_gn


#~ with open('tracks/TLX_TSS3Kb_ChHMM.gmt', 'w') as fp:
    #~ fp.write("\t".join(tss3k_chrhmm_gm))

#~ gen_tb = pd.read_table('results01/TLX3peaks_FANTOM_enhChrHMM-gene.txt', header=1, names=['Genes','Regions'])
#~ gm = ['TLX_FANTOM_ChHMM', 'TLX_FANTOM_ChHMM']  + list(gen_tb['Genes'].str.upper())

#~ with open('tracks/TLX_FANTOM_ChHMM.gmt', 'w') as fp:
    #~ fp.write("\t".join(gm))




# === Load table 
#~ tbl = pd.read_table(join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt'), index_col=0)

#~ tbl = tbl[(tbl.padj < 0.05)].dropna()

# === Load gene names 
#~ names = pd.read_table("tracks/UCSC_mm9_transcripID_to_geneSymbol.sort.txt", 
                           #~ index_col=0,
                           #~ names=["Geneid", "NAME"])


#~ names = names.loc[tbl.index]
#~ assert names.shape[0] == tbl.shape[0]


#~ tbl_raw = tbl[['R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ',
               #~ 'TLX3.1_1','TLX3.1_5','TLX3.1_P',
               #~ 'TAP','TAP1B','TAP2B']]


#~ tbl_n=names.join(tbl_raw, how ='right')

#~ tbl_n['NAME']=tbl_n['NAME'].str.upper()



# === Run GSEA
#~ tbl_c = tbl_n.copy() 

#~ tbl_c.index=tbl_n['NAME']

#gnc = list(gen_tb['Genes'].str.upper())
#tbl_c = tbl_c.loc[gen_tb['Genes'].str.upper()].dropna()
#tbl_cc=tbl_cc.dropna()

#~ tbl_c = tbl_c.groupby(tbl_c.index).agg({'NAME': 'first',
                                #~ 'R2.RAG1W.RAG1':sum,
                                #~ 'RAGS.RAGZ':sum,
                                #~ 'RAGZ':sum,
                                #~ 'TLX3.1_1':sum,
                                #~ 'TLX3.1_5':sum,
                                #~ 'TLX3.1_P':sum,
                                #~ 'TAP':sum,
                                #~ 'TAP1B':sum,
                                #~ 'TAP2B':sum})

#~ tbl_c = tbl_c[['NAME',
            #~ 'R2.RAG1W.RAG1',
            #~ 'RAGS.RAGZ',
            #~ 'RAGZ',
            #~ 'TLX3.1_1',
            #~ 'TLX3.1_5',
            #~ 'TLX3.1_P']]
            #~ 'TAP',
            #~ 'TAP1B',
            #~ 'TAP2B']]

#~ classi = ['RAG','RAG','RAG','TLX3','TLX3','TLX3'] #,'TLX3','TLX3','TLX3']

#~ gs_res = gp.gsea.call(data=tbl_c, 
                        #~ gene_sets= 'tracks/TLX_TSS3Kb_ChHMM.gmt', #gene_sets='KEGG_2016', 
                        #~ cls=classi,
                        #~ max_size = 2000,
                        #~ permutation_type='gene_set',#~ permutation_type='phenotype',
                        #~ outdir='gsea_TLX_TSS3Kb_ChHMM')

# === Pictures

#~ gsea_results = gs_res.reset_index().sort_values('fdr',axis=0,ascending=True)

#~ with plt.style.context('ggplot'):
    #~ gsea_results.head(40).plot.bar(y='fdr',x='Term', figsize=(12, 6),fontsize=12)

#~ plt.show() 




#------------------------------------------------------------------------

#~ genes = test1.closest('tracks/genes.bed')
#~ genes.saveas('tracks/TLX3_peaks_RUNX_ETS_genes.bed')



#~ gs = genes.to_dataframe()
#~ nm = list(gs['thickStart'].str.upper())
#~ nm = list(set(nm))
#~ gmt = {'TLX-RUNX-ETS':nm}


# -----------------------------------------------------------------
#~ b1 = [1,2,3,4,5,9,11,15]
#~ b2 = [4,5,6,7,8]
#~ b3 = [val for val in b1 if val in b2]
#~ or

#~ def intersect(a, b):
    #~ return list(set(a) & set(b))

#~ print intersect(b1, b2)

