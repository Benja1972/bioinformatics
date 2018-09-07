#!/usr/bin/env python2

import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt
from metaseq.results_table import ResultsTable, DESeq2Results

def log2p1(x):
    return np.log2(x + 1)


def write_gmt(st, name, path=''):
    gmt = [name, name] + list(st)
    with open(join(path,name+'.gmt'), 'w') as fp:
        fp.write("\t".join(gmt))

# Save or not
SV = False

## === Bed manipulation ======================
## ===========================================
#~ tlx_peak = pb.BedTool('tracks/TLX3_TLX3_peaks_100.bed')
#~ enh_fantm = pb.BedTool('tracks/annot_tracks/enhancers_FANTOM.bed')
#~ enh_chrhmm = pb.BedTool('tracks/TLX3_6_FE_E4_sorted.bed')
#~ tss_3kb = pb.BedTool('tracks/annot_tracks/references/mm9/mm9_tsses-3kb.gtf')

#~ tlx_tss_3kb = (tlx_peak+tss_3kb).sort()
#~ tlx_chrhmm = (tlx_peak+enh_chrhmm - tlx_tss_3kb).sort()

#~ tlx_tss_3kb.saveas('tracks/TLX3_TSS3Kb5.bed')
#~ tlx_chrhmm.saveas('tracks/TLX3_chrhmm5.bed')




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

from matplotlib_venn import venn2
venn2([tss,enh], set_labels = ('Tlx3 in Enhancer', 'Tlx3 in Promoter'))


tss_and_enh = tss & enh
tss_or_enh = tss | enh
tss_notin_enh = tss - enh
enh_notin_tss = enh - tss
diff_enh_tss = enh ^ tss

if SV:
    write_gmt(tss_and_enh,   'tss_and_enh',  'gene_lists/tlx_enh_tss3k')
    write_gmt(tss_or_enh,    'tss_or_enh',   'gene_lists/tlx_enh_tss3k')
    write_gmt(tss_notin_enh, 'tss_notin_enh','gene_lists/tlx_enh_tss3k')
    write_gmt(enh_notin_tss, 'enh_notin_tss','gene_lists/tlx_enh_tss3k')
    write_gmt(diff_enh_tss,  'diff_enh_tss', 'gene_lists/tlx_enh_tss3k')

# === Load expression table 
tbl = pd.read_table(join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt'), index_col=0)


# Filter genes (Note: this filter remove microRNA expression)
tbl = tbl[(tbl.padj < 0.05)].dropna()


# === Load gene names 
names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           index_col=1,
                           header=0,
                           names=['GeneID', 'TransID', 'Gene_name'])

names = names.loc[tbl.index]
assert names.shape[0] == tbl.shape[0]

tbl=names.join(tbl, how ='right')

## === Expresion analysis
import RNA_expression_processing as rn

tbn = tbl[['Gene_name', 'R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P', 'padj']]

classes = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']

groups ={
    'TSS and Enh':tss_and_enh,
    'Only TSS':tss_notin_enh,
    'ALL:  TSS or Enh':tss_or_enh,
    'Only Enh':enh_notin_tss,
    'In Enh or TSS but not common':diff_enh_tss
}

if SV:
    rn.write_dic2gmt(groups, 
                     name='TLX3pk_Enh_TSS2Kb', 
                     path='gene_lists/tlx_enh_tss3k')


for nm, gns in groups.iteritems():
    gl = list(gns)
    rn.scatter(tbn,'TLX3', 'RAG', classes=classes, 
                n_top=5, geneList=gl, ttl=nm)
    rn.express(tbn,'TLX3', 'RAG', classes=classes, 
                n_top=30,  geneList=gl, ttl=nm) 


plt.show()






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



