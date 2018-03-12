#!/usr/bin/env python3

import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt


def log2p1(x):
    return np.log2(x + 1)


def write_gmt(st, name, path=''):
    gmt = [name, name] + list(st)
    with open(join(path,name+'.gmt'), 'w') as fp:
        fp.write("\t".join(gmt))


## === Gene lists manipulation (tables comes from GREAT)
## =====================================================
tlx_tss_3kb_gt = pd.read_table('gene_lists/tlx_enh_tss3k/TLX3_TSS3Kb-gene.txt', 
                                header=1, 
                                names=['Genes','Regions'])
tlx_chrhmm_gt = pd.read_table('gene_lists/tlx_enh_tss3k/TLX3_enh_hmm6-gene.txt', 
                                header=1, 
                                names=['Genes','Regions'])

tlx_tss_3kb_gn = list(tlx_tss_3kb_gt['Genes'].str.upper())
tlx_chrhmm_gn = list(tlx_chrhmm_gt['Genes'].str.upper())


tss = set(tlx_tss_3kb_gn)

enh = set(tlx_chrhmm_gn)

#~ from matplotlib_venn import venn2
#~ venn2([tss,enh], set_labels = ('Tlx3 in Enhancer', 'Tlx3 in Promoter'))


tss_and_enh = tss & enh
tss_or_enh = tss | enh
tss_notin_enh = tss - enh
enh_notin_tss = enh - tss
diff_enh_tss = enh ^ tss


# === Run GSEApy Enrichr
gp.enrichr(gene_list=list(tss_and_enh), 
            description='tss_and_enh', 
            gene_sets='GO_Biological_Process_2017b', 
            outdir='gene_lists/tlx_enh_tss3k/'+'tss_and_enh_GO_BP_2017b')

gp.enrichr(gene_list=list(tss_or_enh), 
            description='tss_or_enh', 
            gene_sets='GO_Biological_Process_2017b', 
            outdir='gene_lists/tlx_enh_tss3k/'+'tss_or_enh_GO_BP_2017b')

gp.enrichr(gene_list=list(tss_notin_enh), 
            description='tss_notin_enh', 
            gene_sets='GO_Biological_Process_2017b', 
            outdir='gene_lists/tlx_enh_tss3k/'+'tss_notin_enh_GO_BP_2017b')

gp.enrichr(gene_list=list(enh_notin_tss), 
            description='enh_notin_tss', 
            gene_sets='GO_Biological_Process_2017b', 
            outdir='gene_lists/tlx_enh_tss3k/'+'enh_notin_tss_GO_BP_2017b')

gp.enrichr(gene_list=list(diff_enh_tss), 
            description='diff_enh_tss', 
            gene_sets='GO_Biological_Process_2017b', 
            outdir='gene_lists/tlx_enh_tss3k/'+'diff_enh_tss_GO_BP_2017b')


gp.enrichr(gene_list=list(tss_and_enh), 
            description='tss_and_enh', 
            gene_sets='KEGG_2016', 
            outdir='gene_lists/tlx_enh_tss3k/'+'tss_and_enh_KEGG_2016')

gp.enrichr(gene_list=list(tss_or_enh), 
            description='tss_or_enh', 
            gene_sets='KEGG_2016', 
            outdir='gene_lists/tlx_enh_tss3k/'+'tss_or_enh_KEGG_2016')

gp.enrichr(gene_list=list(tss_notin_enh), 
            description='tss_notin_enh', 
            gene_sets='KEGG_2016', 
            outdir='gene_lists/tlx_enh_tss3k/'+'tss_notin_enh_KEGG_2016')

gp.enrichr(gene_list=list(enh_notin_tss), 
            description='enh_notin_tss', 
            gene_sets='KEGG_2016', 
            outdir='gene_lists/tlx_enh_tss3k/'+'enh_notin_tss_KEGG_2016')

gp.enrichr(gene_list=list(diff_enh_tss), 
            description='diff_enh_tss', 
            gene_sets='KEGG_2016', 
            outdir='gene_lists/tlx_enh_tss3k/'+'diff_enh_tss_KEGG_2016')



#plt.show()






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



