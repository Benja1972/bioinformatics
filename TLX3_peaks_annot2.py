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

## == Annotation file , 'V5' = score/10 
df = pd.read_csv('tracks/TLX3_TLX3_peaks_100_annot.csv',index_col=0)

shr = 500/10.
top = df[df['V5']>shr] #.sort_values('score', axis=0, ascending=False)

intr  = top[top['annotation'].str.contains("Intron")]
promot= top[top['annotation'].str.contains("Promoter")]
distal= top[top['annotation'].str.contains("Distal Intergenic")]
exon  = top[top['annotation'].str.contains("Exon")] 


rst = len(top) - len(intr) - len(promot) - len(distal) - len(exon)
print len(intr), len(promot), len(distal), len(exon), rst




# ===== Figures
# Tweak some font settings so the results look nicer
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14

lab = ['Intron', 'Promoter', 'Distal Intergenic', 'Exons', 'Others']

colors = ['gold', 'lightcoral', 'yellowgreen', 'orange', 'lightskyblue']
explode = (0, 0.1, 0, 0, 0) 

fig, ax = plt.subplots()
patches, texts, autotexts = ax.pie([len(intr), len(promot), len(distal), len(exon), rst],  
                                    labels=lab,
                                    colors=colors, 
                                    explode=explode,
                                    shadow=True,
                                    autopct='%1.1f%%')
plt.axis('equal')


# Make the labels  easier to read.
for t in texts:
    t.set_size('large')

font = {'family': 'arial',
        'color':  'darkblue',
        'weight': 'normal',
        'size': 16,
        }

ax.set_title('Peaks with score >'+ str(int(shr*10)) + ', Total =' +str(len(top)), fontdict=font)


# ==== Genes expression
# ==============================


# === Load expression table ---------------------------
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

tbn = tbl[['Gene_name', 'R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P', 'padj']]

classes = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']

# ---------------------------------------------------

#~ intrs = [promot, intr, distal]


intrs = {'Promoter':promot,'Intron': intr, 'Distal': distal}


for nm, df in intrs.iteritems():
    gn_l = list(df['SYMBOL'].str.upper().unique())

    topN, up, dn, ttl =rn.express(tbn,'TLX3', 'RAG', 
                                classes=classes, 
                                n_top=100, 
                                geneList=gn_l, 
                                ttl='Genes of  peaks (score >' + str(int(shr*10)) + ')  in '+nm)



    # Figures

    ttl_l = len(ttl)
    up_l = len(up)
    dn_l = len(dn)

    unch = len(ttl) - len(up) - len(dn)
    nexp = len(gn_l)-len(ttl)

    lab = ['Up-regulated', 'Down-regulated', 'Unchanged', 'Unknown expr']

    colors = ['red', 'blue', 'grey','lightgrey' ]
    explode = (0.00, 0.00, 0, 0) 

    #print 'Promoters = ', up, dn, unch, nexp

    with open('Enrichr/UP_TLX_'+nm+ '_sc500.txt', 'w') as fp:
            fp.write("\n".join(list(up.index)))


    with open('Enrichr/DN_TLX_'+nm+ '_sc500.txt', 'w') as fp:
            fp.write("\n".join(list(dn.index)))


    fig2, ax2 = plt.subplots()
    patches, texts, autotexts = ax2.pie([up_l, dn_l, unch, nexp],  
                                        labels=lab,
                                        colors=colors, 
                                        explode=explode,
                                        shadow=True,
                                        autopct='%1.1f%%')
    plt.axis('equal')


    # Make the labels  easier to read.
    for t in texts:
        t.set_size('large')


    ax2.set_title('Gene of peaks (score >' + str(int(shr*10)) + ') in '+nm+ ', Total =' +str(len(gn_l)), fontdict=font)




plt.show()


# GSEAPY plots

#~ df = pd.read_csv('Enrichr/KEGG_2016_table-UP_TLXp_sc500.txt', sep='\t')

#~ fg  = gp.plot.dotplot(df,scale=150,figsize=(6,6))

#~ plt.show()
