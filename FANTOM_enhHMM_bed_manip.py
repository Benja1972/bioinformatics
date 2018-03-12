#!/usr/bin/env python3

import subprocess
import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt


chrhmm = pb.BedTool('tracks/TLX3_6_FE_E4_sorted.bed')
fantm = pb.BedTool('tracks/annot_tracks/enhancers_FANTOM.bed')
tlx_peak = pb.BedTool('tracks/TLX3_TLX3_peaks_100.bed')

all_enh = chrhmm.cat(fantm) # Python way of merging
#all_enh = pb.BedTool('tracks/annot_tracks/enhHMM_FANTOM.bed')


chrhmm_tlx = (tlx_peak + chrhmm)#.merge()
all_enh_tlx = (tlx_peak+all_enh)
fantm_tlx = (tlx_peak + fantm)

#~ all_chr_fnt =(chrhmm+fantm).merge()

print(chrhmm_tlx.count(), fantm_tlx.count())


chrhmm_fantm = (fantm_tlx + chrhmm_tlx)#.merge()
chrhmm_only = (chrhmm_tlx - fantm_tlx)#.merge()
fantm_only = (fantm_tlx - chrhmm_tlx)#.merge()

#chrhmm_only.saveas('tracks/enh_tlx_notFantom.bed')


# == Filtering
sc = 40
ch_onl = chrhmm_only.to_dataframe()
fnt_onl = fantm_only.to_dataframe()
ch_fnt = chrhmm_fantm.to_dataframe()
alll =   all_enh_tlx.to_dataframe()
ch_all = chrhmm_tlx.to_dataframe()

top_ch =  ch_onl[ch_onl['score']>sc].drop(['score'], axis = 1)
top_fnt = fnt_onl[fnt_onl['score']>sc].drop(['score'], axis = 1)
top_ch_fnt = ch_fnt[ch_fnt['score']>sc].drop(['score'], axis = 1)
top_alll = alll[alll['score']>sc].drop(['score'], axis = 1)
top_ch_all = ch_all[ch_all['score']>sc].drop(['score'], axis = 1)

# === Back to bed, pythonic way
hmm_onl = pb.BedTool.from_dataframe(top_ch, 
                                    outfile='proj_enh/HMMenh_tlx_notFantom_topScore40.bed')

hmm_all = pb.BedTool.from_dataframe(top_ch_all, 
                                    outfile='proj_enh/HMMenh_tlx_all_topScore40.bed')

fnt_onl = pb.BedTool.from_dataframe(top_fnt, 
                                    outfile='proj_enh/Fantom_tlx_notHMM_topScore40.bed')

hmm_fnt = pb.BedTool.from_dataframe(top_ch_fnt, 
                                    outfile='proj_enh/Fantom_andHMM_tlx_topScore40.bed')

all_all = pb.BedTool.from_dataframe(top_alll, 
                                    outfile='proj_enh/Enh_all_tlx_topScore40.bed')

#~ top_ch_fnt[['chrom','start','end', 'name']].to_csv('tracks/Fantom_andHMM_tlx_topScore40.bed', 
                                        #~ sep='\t', 
                                        #~ index=False,
                                        #~ header=False)

## === get DAN seqs from bed coordinates
#~ mm9 = "/home/sergio/media/NAS4/PFlab/TLX3_project/ChiP-Seq/references/mm9/chromFa/mm9.fa"

#~ # --- 
#~ bd = "tracks/HMMenh_tlx_notFantom_topScore40.bed"
#~ fa = "tracks/HMMenh_tlx_notFantom_topScore40.fa"
#~ bd2fa = "bedtools getfasta -fi {} -bed {} -fo {}".format(mm9,bd,fa)

#~ print('Running process ........ \n')
#~ print(bd2fa)

#~ subprocess.call(['bash','-c', bd2fa])

#~ # --- 
#~ bd = "tracks/Fantom_tlx_notHMM_topScore40.bed"
#~ fa = "tracks/Fantom_tlx_notHMM_topScore40.fa"
#~ bd2fa = "bedtools getfasta -fi {} -bed {} -fo {}".format(mm9,bd,fa)

#~ print('Running process ........ \n')
#~ print(bd2fa)

#~ subprocess.call(['bash','-c', bd2fa])

#~ # --- 
#~ bd = "tracks/Fantom_andHMM_tlx_topScore40.bed"
#~ fa = "tracks/Fantom_andHMM_tlx_topScore40.fa"
#~ bd2fa = "bedtools getfasta -fi {} -bed {} -fo {}".format(mm9,bd,fa)

#~ print('Running process ........ \n')
#~ print(bd2fa)

#~ subprocess.call(['bash','-c', bd2fa])

## ====== Graph
from matplotlib_venn import venn2
venn2(subsets=(len(top_ch), len(top_fnt), len(top_ch_fnt)),  
        set_labels = ('Enh_HMM with TLX3(sc>40)', 'FANTOM_enh with TLX3(sc>40)'))






## === Gene lists manipulation (tables comes from GREAT)
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


tbn = tbl[['Gene_name', 'R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P', 'padj']]


## === Expresion analysis
classes = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']

import RNA_expression_processing as rn


pref = 'proj_enh/'
cas = ['HMMenh_tlx_notFantom_topScore40',
        'HMMenh_tlx_all_topScore40',
        'Fantom_tlx_notHMM_topScore40',
        'Fantom_andHMM_tlx_topScore40',
        'Enh_all_tlx_topScore40'
        ]

for sml in cas:
    gt = pd.read_table(pref+sml+'-gene.txt', 
                        header=1, 
                        names=['Genes','Regions'])

    gt['Genes'] = gt['Genes'].str.upper()
    gn = list(gt['Genes'].str.upper())

    topN,_,_,f1  = rn.express(tbn, 'TLX3', 'RAG', 
                        classes=classes, 
                        n_top=40, 
                        geneList=gn,  
                        ttl='Genes for '+sml)
    f1.savefig('proj_enh/' +sml+'_expr.pdf', bbox_inches='tight')
    up, dn, ax2 = rn.scatter(tbn, 'TLX3', 'RAG', 
            classes=classes, 
            n_top=5, 
            geneList=gn,  
            ttl='Genes for '+sml)
    f2 = plt.gcf()
    f2.savefig('proj_enh/' +sml+'_scatt.pdf', bbox_inches='tight')

plt.show()


#~ enh_notFANTOM_gt = pd.read_table('gene_lists/enh_tlx_notFantom_topScore40-gene.txt', 
                                #~ header=1, 
                                #~ names=['Genes','Regions'])

#~ enh_notFANTOM_rt = pd.read_table('gene_lists/enh_tlx_notFantom_topScore40-region.txt', 
                                #~ header=1, 
                                #~ names=['Regions','Genes'])

#~ enh_notFANTOM_gt['Genes'] = enh_notFANTOM_gt['Genes'].str.upper()

#~ enh_notFANTOM_gn = list(enh_notFANTOM_gt['Genes'].str.upper())

# === Load expression table 
#~ tbl = pd.read_table(join('tracks', 'TLX3vsRAG-results_genes.txt'), index_col=0)


# Filter genes (Note: this filter remove microRNA expression)
#~ tbl = tbl[(tbl.padj < 0.05)].dropna()


# === Load gene names 
#~ names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           #~ index_col=0,
                           #~ header=0,
                           #~ names=['GeneID', 'TransID', 'Gene_name'])


#~ names = names.drop('TransID', axis=1).drop_duplicates()
#~ names = names.loc[tbl.index]
#~ assert names.shape[0] == tbl.shape[0]

#~ tbl=names.join(tbl, how ='right')


#~ tbn = tbl[['Gene_name', 'R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P', 'padj']]


## === Expresion analysis
#~ classes = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']

#~ import RNA_expression_processing as rn

#~ topN  = rn.express(tbn, 'TLX3', 'RAG', 
                    #~ classes=classes, 
                    #~ n_top=40, 
                    #~ geneList=enh_notFANTOM_gn,  
                    #~ ttl='Genes for enhancers (not FANTOM)\n with top score TLX3(sc>40)')

#~ up, dn = rn.updn(tbn, 'TLX3', 'RAG', classes=classes, geneList=enh_notFANTOM_gn)


#~ up, dn = rn.scatter(tbn, 'TLX3', 'RAG', 
            #~ classes=classes, 
            #~ n_top=5, 
            #~ geneList=enh_notFANTOM_gn,  
            #~ ttl='Genes  enhancers')





#~ plt.show()

# === Save to file
#~ up_l = list(up.index)
#~ ds_u = 'UP_enh_notFANTOM_score40'
#~ with open(join('Enrichr',ds_u+'.txt'), 'w') as fp:
    #~ fp.write("\n".join(up_l))


#~ dn_l = list(dn.index)
#~ ds_d = 'DN_enh_notFANTOM_score40'
#~ with open(join('Enrichr',ds_d+'.txt'), 'w') as fp:
    #~ fp.write("\n".join(dn_l))

# === Run GSEApy Enrichr

#~ gs = 'GO_Biological_Process_2017b'
#~ gp.enrichr(gene_list=list(up.index), 
            #~ description=ds_u, 
            #~ gene_sets= gs, 
            #~ outdir='Enrichr/'+ds_u+'_'+gs)


#~ gs = 'KEGG_2016'
#~ gp.enrichr(gene_list=list(up.index), 
            #~ description=ds_u, 
            #~ gene_sets= gs, 
            #~ outdir='Enrichr/'+ds_u+'_'+gs)


#~ gs = 'GO_Biological_Process_2017b'
#~ gp.enrichr(gene_list=list(dn.index), 
            #~ description=ds_d, 
            #~ gene_sets= gs, 
            #~ outdir='Enrichr/'+ds_d+'_'+gs)


#~ gs = 'KEGG_2016'
#~ gp.enrichr(gene_list=list(dn.index), 
            #~ description=ds_d, 
            #~ gene_sets= gs, 
            #~ outdir='Enrichr/'+ds_d+'_'+gs)






