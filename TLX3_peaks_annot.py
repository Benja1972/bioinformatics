#!/usr/bin/env python3

#~ import os 
import subprocess
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

# ===== Figures
# Tweak some font settings so the results look nicer
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14


# Read annotated peaks
df = pd.read_csv('tracks/TLX3_TLX3_peaks_100_annot.csv',index_col=0)

intr = df[df['annotation'].str.contains("Intron")]
promot=df[df['annotation'].str.contains("Promoter")]
distal=df[df['annotation'].str.contains("Distal Intergenic")]
exon  =df[df['annotation'].str.contains("Exon")] 

rst = len(df) - len(intr) - len(promot) - len(distal) - len(exon)
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

ax.set_title('Peaks statistics, Total =' +str(len(df)), fontdict=font)


# ==== Genes expression
# ==============================

gn_l = list(promot['SYMBOL'].str.upper().unique())



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

tbn = tbl[['Gene_name', 'R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P', 'padj']]

classes = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']


topN, up, dn, ttl =rn.express(tbn,'TLX3', 'RAG', 
                            classes=classes, 
                            n_top=100, 
                            geneList=gn_l, 
                            ttl='Genes of TLX3-tlx peaks in promoters')



# Figures

unch = ttl - up - dn
nexp = len(gn_l)-ttl

lab = ['Up-regulated', 'Down-regulated', 'Unchanged', 'Not-expressed']

colors = ['gold', 'lightcoral', 'lightskyblue','grey' ]
explode = (0.07, 0.04, 0, 0) 

print 'Promoters = ', up, dn, unch, nexp

fig2, ax2 = plt.subplots()
patches, texts, autotexts = ax2.pie([up, dn, unch, nexp],  
                                    labels=lab,
                                    colors=colors, 
                                    explode=explode,
                                    shadow=True,
                                    autopct='%1.1f%%')
plt.axis('equal')


# Make the labels  easier to read.
for t in texts:
    t.set_size('large')


ax2.set_title('Peaks in promoter, Total =' +str(len(gn_l)), fontdict=font)




plt.show()

## === Save to .bed file
#~ promot[['seqnames','start','end']].to_csv('tracks/Tlx3_peaks_100_promoter.bed', 
                                    #~ sep='\t', 
                                    #~ index=False,
                                    #~ header=False)

#~ distal[['seqnames','start','end']].to_csv('tracks/Tlx3_peaks_100_distal.bed', 
                                    #~ sep='\t', 
                                    #~ index=False,
                                    #~ header=False)

#~ intr[['seqnames','start','end']].to_csv('tracks/Tlx3_peaks_100_intron.bed', 
                                    #~ sep='\t', 
                                    #~ index=False,
                                    #~ header=False)


#~ exon[['seqnames','start','end']].to_csv('tracks/Tlx3_peaks_100_exon.bed', 
                                    #~ sep='\t', 
                                    #~ index=False,
                                    #~ header=False)

## === get DAN seqs from bed coordinates
#~ mm9 = "/home/sergio/media/NAS4/PFlab/TLX3_project/ChiP-Seq/references/mm9/chromFa/mm9.fa"

#~ # --- distal
#~ bd = "tracks/Tlx3_peaks_100_distal.bed"
#~ fa = "tracks/Tlx3_peaks_100_distal.fa"
#~ bd2fa = "bedtools getfasta -fi {} -bed {} -fo {}".format(mm9,bd,fa)

#~ print('Running process ........ \n')
#~ print(bd2fa)

#~ subprocess.call(['bash','-c', bd2fa])


#~ # ---- promoter
#~ bd = "tracks/Tlx3_peaks_100_promoter.bed"
#~ fa = "tracks/Tlx3_peaks_100_promoter.fa"
#~ bd2fa = "bedtools getfasta -fi {} -bed {} -fo {}".format(mm9,bd,fa)

#~ print('Running process ........ \n')
#~ print(bd2fa)

#~ subprocess.call(['bash','-c', bd2fa])

#~ # ---- intron
#~ bd = "tracks/Tlx3_peaks_100_intron.bed"
#~ fa = "tracks/Tlx3_peaks_100_intron.fa"
#~ bd2fa = "bedtools getfasta -fi {} -bed {} -fo {}".format(mm9,bd,fa)

#~ print('Running process ........ \n')
#~ print(bd2fa)

#~ subprocess.call(['bash','-c', bd2fa])

