#!/usr/bin/env python3

import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt
#~ from metaseq.results_table import ResultsTable, DESeq2Results

def log2p1(x):
    return np.log2(x + 1)

## === Bed manipulation -- peaks common in RAG and TLX3
tlx_peak = pb.BedTool('tracks/TLX3_TLX3_peaks_100.bed')
rag_peak = pb.BedTool('tracks/RAG_TLX3_repl1_peaks.narrowPeak')


rag_only = (rag_peak - tlx_peak).sort()

rag_tlx = (tlx_peak+rag_peak).sort().merge()

rag_tlx.saveas('tracks/RAG_tlx_TLX3_tlx_common.bed')

genes = pb.BedTool('tracks/annot_tracks/genes.bed')
rg = rag_tlx.closest(genes)
rg_df = rg.to_dataframe()
nm = list(rg_df['thickStart'].str.upper())


## == Proper read .narrowPeak files
pk = pd.read_table('tracks/RAG_TLX3_repl1_peaks.narrowPeak',sep='\t', names=['Chr', 
    'start', 'end', 'name', 'score', 'strand','FoldChange','-log10pvalue',
    '-log10qvalue', 'smmt_rel_to_start'])

shr = 1000
top = pk[pk['score']>shr].sort_values('score', axis=0, ascending=False)

# Write selected back to simple .bed file
top[['Chr','start','end']].to_csv('tracks/RAG_Tlx3_peaks_sc_gr1000.bed', 
                                    sep='\t', 
                                    index=False,
                                    header=False)


# ===== Figures
# Tweak some font settings so the results look nicer
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14

font = {'family': 'arial',
        'color':  'darkblue',
        'weight': 'normal',
        'size': 16,
        }

with plt.style.context('default'):
    fig, ax = plt.subplots()
    val, bins, patches = ax.hist(pk['score'], 
                                bins=10,
                                facecolor='#1F77B4', 
                                edgecolor='#1F77B4')
    ax.set_xlabel('MACS2 score')
    ax.set_ylabel('Number of peaks')
    ax.set_title(r'TLX3 peaks in RAG distribution by score')
    ax.text(50, max(val)-5, 'Total = ' + str(len(pk)), fontdict=font)
    ax.text(shr+10, max(val)-5, 'N(score>'+str(shr)+') = ' + str(len(top)), fontdict=font)
    ax.axvline(x=shr, color='darkred', linewidth=2)
    fig.tight_layout()
    #sns.distplot(pk['score'], bins=150, hist=True, norm_hist=True)



# ==== Genes expression
# ==============================
# === Load expression table 
tbl = pd.read_table(join('tracks', 'TLX3vsRAG-results_genes.txt'), index_col=0)


# Filter genes (Note: this filter remove microRNA expression)
tbl = tbl[(tbl.padj < 0.05)].dropna()


# === Load gene names 
names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           index_col=0,
                           header=0,
                           names=['GeneID', 'TransID', 'Gene_name'])


names = names.drop('TransID', axis=1).drop_duplicates()
names = names.loc[tbl.index]
assert names.shape[0] == tbl.shape[0]

tbl=names.join(tbl, how ='right')

## === Expresion analysis
import RNA_expression_processing as rn

rn.rna_expression(tbl, len(nm), nm, ttl='Genes common for \n RAG-tlx and TLX3-tlx peaks')

rn.rna_scatter(tbl, 8, nm, ttl='Genes common for \n RAG-tlx and TLX3-tlx peaks')
rn.rna_volcano(tbl, 8, nm, ttl='Genes common for \n RAG-tlx and TLX3-tlx peaks')

plt.show()






## === Gene lists manipulation (tables comes from GREAT)
#~ tlx_tss_3kb_gt = pd.read_table('gene_lists/tlx_enh_tss3k/TLX3_TSS3Kb-gene.txt', header=1, names=['Genes','Regions'])
#~ tlx_chrhmm_gt = pd.read_table('gene_lists/tlx_enh_tss3k/TLX3_enh_hmm6-gene.txt', header=1, names=['Genes','Regions'])

#~ tlx_tss_3kb_gn = list(tlx_tss_3kb_gt['Genes'].str.upper())
#~ tlx_chrhmm_gn = list(tlx_chrhmm_gt['Genes'].str.upper())

#~ def intersect(a, b):
    #~ return list(set(a) & set(b))

#~ tss3k_chrhmm_gn = intersect(tlx_tss_3kb_gn, tlx_chrhmm_gn)
#~ tss3k_chrhmm_gm = ['TLX_TSS3Kb_Enh_HMM6', 'TLX_TSS3Kb_Enh_HMM6'] + tss3k_chrhmm_gn

#~ with open('tracks/TLX_TSS3Kb_ChHMM.gmt', 'w') as fp:
    #~ fp.write("\t".join(tss3k_chrhmm_gm))

#~ tss = set(tlx_tss_3kb_gn)

#~ enh = set(tlx_chrhmm_gn)

#~ from matplotlib_venn import venn2
#~ venn2([tss,enh], set_labels = ('Tlx3 in Enhancer', 'Tlx3 in Promoter'))









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



