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


# == Write gene set in gmt format
def write_gmt(st, name, path=''):
    gmt = [name, name] + list(st)
    with open(join(path,name+'.gmt'), 'w') as fp:
        fp.write("\t".join(gmt))

# == Load gene set in gmt format
def read_gmt(name, path=''):
    with open(join(path,name)) as f:
         gene_list = f.read().split()[2:]
    return list(gene_list)





## === Gene lists manipulation (tables comes from GREAT)
## =====================================================
#~ tlx_tss_3kb_gt = pd.read_table('gene_lists/tlx_enh_tss3k/TLX3_TSS3Kb-gene.txt', 
                                #~ header=1, 
                                #~ names=['Genes','Regions'])
#~ tlx_chrhmm_gt = pd.read_table('gene_lists/tlx_enh_tss3k/TLX3_enh_hmm6-gene.txt', 
                                #~ header=1, 
                                #~ names=['Genes','Regions'])

#~ tlx_tss_3kb_gn = list(tlx_tss_3kb_gt['Genes'].str.upper())
#~ tlx_chrhmm_gn = list(tlx_chrhmm_gt['Genes'].str.upper())

supp_l = read_gmt('T-ALL_suppressor.gmt', 'tracks/GSEA_gene_sets')
supp_l = [x.upper() for x in supp_l]

enh_l = read_gmt('TLX_FANTOM_ChHMM.gmt', 'tracks')
enh_l = [x.upper() for x in enh_l]


supp = set(supp_l)

enh = set(enh_l)

from matplotlib_venn import venn2
venn2([supp,enh], set_labels = ('T-ALL suppressors ', 'Enhancers and FANTOM'))


supp_and_enh = supp & enh
print supp_and_enh
#~ tss_or_enh = tss | enh
#~ tss_notin_enh = tss - enh
#~ enh_notin_tss = enh - tss
#~ diff_enh_tss = enh ^ tss

#~ write_gmt(tss_and_enh,   'tss_and_enh',  'gene_lists/tlx_enh_tss3k')
#~ write_gmt(tss_or_enh,    'tss_or_enh',   'gene_lists/tlx_enh_tss3k')
#~ write_gmt(tss_notin_enh, 'tss_notin_enh','gene_lists/tlx_enh_tss3k')
#~ write_gmt(enh_notin_tss, 'enh_notin_tss','gene_lists/tlx_enh_tss3k')
#~ write_gmt(diff_enh_tss,  'diff_enh_tss', 'gene_lists/tlx_enh_tss3k')

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

rn.rna_scatter(tbl, list(supp), 8, 'T-ALL suppressors')
rn.rna_expression(tbl, list(supp), 27, 'T-ALL suppressors')

#~ rn.rna_scatter(tbl, list(tss_notin_enh), 5, 'Only TSS')
#~ rn.rna_expression(tbl, list(tss_notin_enh), 30, 'Only TSS')

#~ rn.rna_scatter(tbl, list(tss_or_enh), 5, 'ALL:  TSS or Enh')
#~ rn.rna_expression(tbl, list(tss_or_enh), 30, 'ALL:  TSS or Enh')

#~ rn.rna_scatter(tbl, list(enh_notin_tss), 5,  'Only Enh')
#~ rn.rna_expression(tbl, list(enh_notin_tss), 30,  'Only Enh')


#~ rn.rna_scatter(tbl, list(diff_enh_tss), 5,  'In Enh ot TSS but not common')
#~ rn.rna_expression(tbl, list(diff_enh_tss), 30,  'In Enh ot TSS but not common')

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



