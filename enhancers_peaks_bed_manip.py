#!/usr/bin/env python3

import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt


## === Bed manipulation
#~ tlx_peak = pb.BedTool('tracks/TLX3_TLX3_peaks_100.bed')
#~ enh_fantm = pb.BedTool('tracks/annot_tracks/enhancers_FANTOM.bed')
#~ enh_chrhmm = pb.BedTool('tracks/TLX3_6_FE_E4_sorted.bed')
#~ tss_3kb = pb.BedTool('tracks/annot_tracks/references/mm9/mm9_tsses-3kb.gtf')

#~ tlx_tss_3kb = (tlx_peak+tss_3kb).sort()
#~ tlx_chrhmm = (tlx_peak+enh_chrhmm - tlx_tss_3kb).sort()

#~ tlx_tss_3kb.saveas('tracks/TLX3_TSS3Kb5.bed')
#~ tlx_chrhmm.saveas('tracks/TLX3_chrhmm5.bed')

## == Proper read .narrowPeak files
 #~ pk = pd.read_table('TLX3_TLX3_peaks.narrowPeak',sep='\t', names=['Chr', 
    #~ 'start', 'end', 'name', 'score', 'strand','FoldChange','-log10pvalue',
    #~ '-log10qvalue', 'smmt_rel_to_start'])

#~ top = pk[pk['score']>1000].sort_values('score', axis=0, ascending=False)


## === Gene lists manipulation (tables comes from GREAT)
tlx_tss_3kb_gt = pd.read_table('gene_lists/tlx_enh_tss3k/TLX3_TSS3Kb-gene.txt', header=1, names=['Genes','Regions'])
tlx_chrhmm_gt = pd.read_table('gene_lists/tlx_enh_tss3k/TLX3_enh_hmm6-gene.txt', header=1, names=['Genes','Regions'])

tlx_tss_3kb_gn = list(tlx_tss_3kb_gt['Genes'].str.upper())
tlx_chrhmm_gn = list(tlx_chrhmm_gt['Genes'].str.upper())

#~ def intersect(a, b):
    #~ return list(set(a) & set(b))

#~ tss3k_chrhmm_gn = intersect(tlx_tss_3kb_gn, tlx_chrhmm_gn)
#~ tss3k_chrhmm_gm = ['TLX_TSS3Kb_Enh_HMM6', 'TLX_TSS3Kb_Enh_HMM6'] + tss3k_chrhmm_gn

#~ with open('tracks/TLX_TSS3Kb_ChHMM.gmt', 'w') as fp:
    #~ fp.write("\t".join(tss3k_chrhmm_gm))

tss = set(tlx_tss_3kb_gn)

enh = set(tlx_chrhmm_gn)

from matplotlib_venn import venn2
venn2([tss,enh], set_labels = ('Tlx3 in Enhancer', 'Tlx3 in Promoter'))

plt.show()







#~ # === Load table 
#~ tbl = pd.read_table(join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt'), index_col=0)

#~ tbl = tbl[(tbl.padj < 0.05)].dropna()

#~ # === Load gene names 
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

