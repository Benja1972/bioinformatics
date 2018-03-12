#!/usr/bin/env python3

import numpy as np
import yaml
from os.path import join 
import pandas as pd
#import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
from metaseq.results_table import ResultsTable, DESeq2Results
from mpl_toolkits.axes_grid1 import ImageGrid



def log2p1(x):
    return np.log2(x + 1)


# === Load expression table 
tbl = pd.read_table(join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt'), index_col=0)

# === Load gene set to analyse
#~ gene_set = pd.read_table("tracks/miRNA_T-ALL.txt", names=['Gene_name'])

# == Load gene set in gmt format
with open('tracks/TLX_TSS3Kb_ChHMM.gmt') as f:
     gene_list = f.read().split()[2:]


# Filter genes (Note: this filter remove microRNA expression)
tbl = tbl[(tbl.padj < 0.05)].dropna()

tbl['TLX3-mean'] = log2p1((tbl['TLX3.1_1']+tbl['TLX3.1_5']+tbl['TLX3.1_P'])/3.)
tbl['TAP-mean'] = log2p1((tbl['TAP']+tbl['TAP1B']+tbl['TAP2B'])/3.)
tbl['RAG-mean'] = log2p1((tbl['R2.RAG1W.RAG1']+tbl['RAGS.RAGZ']+tbl['RAGZ'])/3.)

# === Load gene names 
names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           index_col=1,
                           header=0,
                           names=['GeneID', 'TransID', 'Gene_name'])

names = names.loc[tbl.index]
assert names.shape[0] == tbl.shape[0]

tbl=names.join(tbl, how ='right')

# === Strip expression table to gene list
# df.loc[df['column_name'] == some_value]
# gs_tbl = tbl.loc[tbl['Gene_name'].isin(gene_set['Gene_name'])]

# for gmt gene names should UPPERCASE
tbl['Gene_name']=tbl['Gene_name'].str.upper()
gs_tbl = tbl.loc[tbl['Gene_name'].isin(gene_list)]



# ==== Gene set expression =======================
# ================================================


gs_tbl = ResultsTable(gs_tbl)


up = (gs_tbl.log2FoldChange > 1.).values
dn = (gs_tbl.log2FoldChange < -1.).values
un = ((gs_tbl.log2FoldChange >= -1) & (gs_tbl.log2FoldChange <= 1)).values

upr = gs_tbl[(gs_tbl.log2FoldChange > 1.)].sort_values('log2FoldChange', axis=0, ascending=False)
dnr = gs_tbl[(gs_tbl.log2FoldChange < -1.)].sort_values('log2FoldChange', axis=0, ascending=False)

print 'Total counts = ', len(gs_tbl)
print 'Over-expressed counts = ', len(upr)
print 'Under-expressed counts = ', len(dnr)
print 'Unchanged counts = ', len(gs_tbl) - len(dnr) - len(upr)

# Show names of genes in graph
num_g = 10
up10 = upr[:num_g+1]

dn10 = dnr[-num_g:]

# = Top Differently expressed
gs_tbl.data['FCmod'] = abs(gs_tbl.log2FoldChange)
top_upd = gs_tbl#.sort_values('FCmod', axis=0, ascending=False)

top_upd = top_upd.sort_values('log2FoldChange', axis=0, ascending=False)
top_upd = top_upd[['Gene_name', 'log2FoldChange','R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P']]

top_upd.set_index('Gene_name', inplace='True')


top_upd = top_upd.groupby(top_upd.index).mean().sort_values('log2FoldChange', axis=0, ascending=False)

#~ top_upd.drop('log2FoldChange', axis=1, inplace=True)
#print list(top_upd.index)
top_upd = top_upd[['R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P']]




fig31 = plt.figure(figsize=(6, 12))
ax21 = fig31.add_subplot(111)
sns.heatmap(log2p1(top_upd[:50]), cmap='RdBu_r',  linewidths=0.004) #cmap=plt.cm.RdBu_r/'seismic', linewidths=0.5, square=True)
for txt in ax21.get_yticklabels():
        txt.set_rotation(0)
        #txt.set_fontsize(7)
for txt in ax21.get_xticklabels():
        txt.set_rotation(90)

fig51 = plt.figure(figsize=(6, 12))
ax51 = fig51.add_subplot(111)
sns.heatmap(log2p1(top_upd[-50:]), cmap='RdBu_r',  linewidths=0.004) #cmap=plt.cm.RdBu_r/'seismic', linewidths=0.5, square=True)
for txt in ax51.get_yticklabels():
        txt.set_rotation(0)
        #txt.set_fontsize(7)
for txt in ax51.get_xticklabels():
        txt.set_rotation(90)


# ============= Scatterplot 
ax = gs_tbl.scatter(
    x='RAG-mean',
    y='TLX3-mean',
    #----------------
    #~ xfunc=log2p1,
    #~ yfunc=log2p1,
    xlab='RAG, log2(FPKM + 1)',
    ylab='TLX3, log2(FPKM + 1)',
    one_to_one=dict(color='k', linestyle=':'),
    marginal_histograms=False,
    # add the "unchanged" label
    general_kwargs=dict(marker='.', color='0.5', alpha=0.7, s=10, label='unchanged'),

    genes_to_highlight=[
    (
     up,
     dict(
          color='#da3b3a', alpha=0.8,
          marginal_histograms=True,
          xhist_kwargs=dict(bins=50, linewidth=0),
          yhist_kwargs=dict(bins=50, linewidth=0),

          # add label
          label='upregulated',
          )
     ),
    (
     dn,
     dict(
          color='#00748e', alpha=0.8,
          marginal_histograms=True,
          xhist_kwargs=dict(bins=50, linewidth=0),
          yhist_kwargs=dict(bins=50, linewidth=0),

          # add label
          label='downregulated'
          )
     )
    ],
);




# Plot text of top genes
for i in range(num_g):
    ax.text(up10['RAG-mean'][i],up10['TLX3-mean'][i],up10['Gene_name'][i].upper(), color='#da3b3a')
    ax.text(dn10['RAG-mean'][i],dn10['TLX3-mean'][i],dn10['Gene_name'][i].upper(), color='#00748e')

# Get handles and labels, and then reverse their order
handles, legend_labels = ax.get_legend_handles_labels()
handles = handles[::-1]
legend_labels = legend_labels[::-1]

# Draw a legend using the flipped handles and labels.
leg = ax.legend(handles,
          legend_labels,

          # These values may take some tweaking.
          # By default they are in axes coordinates, so this means
          # the legend is slightly outside the axes.
          loc=(1.01, 1.05),

          # Various style fixes to default legend.
          fontsize=9,
          scatterpoints=1,
          borderpad=0.1,
          handletextpad=0.05,
          frameon=False,
          title=' Transcripts',
          );

# Adjust the legend title after it's created
leg.get_title().set_weight('bold')

top_axes = gs_tbl.marginal.top_hists[-1]
top_axes.set_title('Differential expression, TLX3 vs RAG')

for ax in gs_tbl.marginal.top_hists:
    ax.set_ylabel('No.\ntranscripts', rotation=0, ha='right', va='center', size=8)

for ax in gs_tbl.marginal.right_hists:
    ax.set_xlabel('No.\ntranscripts', rotation=-90, ha='left', va='top', size=8)

# ==================================================================













#~ gs_raw = gs_tbl[['Gene_name', 'R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ',
               #~ 'TLX3.1_1','TLX3.1_5','TLX3.1_P',
               #~ 'TAP','TAP1B','TAP2B']]

#~ gs_raw.set_index('Gene_name', inplace='True')


#~ fig41 = plt.figure(figsize=(5, 5))
#~ ax41 = fig41.add_subplot(111)
#~ sns.heatmap(log2p1(gs_raw[:20]), cmap='RdBu_r', linewidths=0.005)
#~ for txt in ax41.get_yticklabels():
        #~ txt.set_rotation(0)
        #~ #txt.set_fontsize(7)
#~ for txt in ax41.get_xticklabels():
        #~ txt.set_rotation(90)

plt.show()





#~ tbl_n=names.join(tbl_raw, how ='right')

#~ tbl_n['NAME']=tbl_n['NAME'].str.upper()

#~ # === Run GSEA
#~ classi = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']

#~ #tbl_g = tbl_n.copy()#head(23737)
#~ #tbl_g = tbl_n.copy()


#~ tbl_c = tbl_n.copy() 

#~ tbl_c.index=tbl_n['NAME']



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

#~ g_set='KEGG_2015'
#~ out_dir = 'GSEA/gsea_'+g_set
#~ gs_res = gp.gsea.call(data=tbl_c, 
                      #~ gene_sets=g_set,#gene_sets= 'tracks/GSEA_gene_sets/c2.cp.kegg.v6.0.symbols.gmt',#gene_sets='KEGG_2016', 
                      #~ cls=classi,
                      #~ graph_num = 100,
                      #~ #permutation_type='phenotype',
                      #~ outdir=out_dir)

#~ # === Pictures

#~ gsea_results = gs_res.reset_index().sort_values('fdr',axis=0,ascending=True)

#~ with plt.style.context('ggplot'):
    #~ gsea_results.head(40).plot.bar(y='fdr',x='Term', figsize=(12, 6),fontsize=12)

#~ plt.savefig(out_dir+'/'+g_set+'.pdf', format="pdf")
#~ plt.show() 



#~ with plt.style.context('ggplot'):
    #~ gsea_results = gs_res.reset_index()  
    #~ gsea_results.head(5).plot.barh(y='fdr',x='Term',fontsize=12)

   


