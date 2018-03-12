
import numpy as np
import yaml
import os
from os.path import join 
import pybedtools as pb
import gffutils
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import matplotlib.pyplot as plt
import seaborn as sns
import metaseq
import pandas as pd
from metaseq.results_table import ResultsTable, DESeq2Results
from mpl_toolkits.axes_grid1 import ImageGrid

def log2p1(x):
    return np.log2(x + 1)


# Tweak some font settings so the results look nicer
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12


# RNA-seq expression tables

tbl_f = join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt')

tbl = ResultsTable(tbl_f, import_kwargs=dict(index_col=0))
tbl.data['TLX3-mean'] = log2p1((tbl.data['TLX3.1_1']+tbl.data['TLX3.1_5']+tbl.data['TLX3.1_P'])/3.)
tbl.data['TAP-mean'] = log2p1((tbl.data['TAP']+tbl.data['TAP1B']+tbl.data['TAP2B'])/3.)
tbl.data['RAG-mean'] = log2p1((tbl.data['R2.RAG1W.RAG1']+tbl.data['RAGS.RAGZ']+tbl.data['RAGZ'])/3.)

tbl = ResultsTable(tbl[(tbl.padj < 0.05)].dropna())

gene_names = pd.read_table("tracks/UCSC_mm9_transcripID_to_geneSymbol.sort.txt", 
                            index_col=0,
                            names=["Geneid", "Gene_name"])



# Specific gene name, ex. NOTCH1
tbl = tbl.join(gene_names)


notch = tbl.loc[tbl['Gene_name']=='Notch1'][['R2.RAG1W.RAG1',
                                            'RAGS.RAGZ',
                                            'RAGZ',
                                            'TLX3.1_1',
                                            'TLX3.1_5',
                                            'TLX3.1_P',
                                            'TAP',
                                            'TAP1B',
                                            'TAP2B']]


fig41 = plt.figure(figsize=(5, 7))
ax41 = fig41.add_subplot(111)
sns.heatmap(log2p1(notch), cmap='RdBu_r', linewidths=0.005)
for txt in ax41.get_yticklabels():
        txt.set_rotation(0)
        #~ txt.set_fontsize(7)
for txt in ax41.get_xticklabels():
        txt.set_rotation(90)

# -------------------------



up = (tbl.log2FoldChange > 1.).values
dn = (tbl.log2FoldChange < -1.).values
un = ((tbl.log2FoldChange >= -1) & (tbl.log2FoldChange <= 1)).values



upr = tbl[(tbl.log2FoldChange > 1.)].sort_values('log2FoldChange', axis=0, ascending=False)
dnr = tbl[(tbl.log2FoldChange < -1.)].sort_values('log2FoldChange', axis=0, ascending=False)
#~ unr = tbl[(tbl.log2FoldChange >= -2) & (tbl.log2FoldChange <= 2)]

tbl.data['FCmod'] = abs(tbl.log2FoldChange)
top_upd = tbl.sort_values('FCmod', axis=0, ascending=False)
top_upd = top_upd[:200].sort_values('log2FoldChange', axis=0, ascending=False)
top_upd = log2p1(top_upd[['R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P']])

top_names = gene_names.loc[top_upd.index]
top_upd=top_upd.join(top_names)
top_updR = top_upd.set_index('Gene_name')

print list(top_updR.index)

fig31 = plt.figure(figsize=(6, 18))
ax21 = fig31.add_subplot(111)
sns.heatmap(top_updR, cmap='plasma',  linewidths=0.004) #cmap=plt.cm.RdBu_r/'seismic', linewidths=0.5, square=True)
for txt in ax21.get_yticklabels():
        txt.set_rotation(0)
        txt.set_fontsize(7)
for txt in ax21.get_xticklabels():
        txt.set_rotation(90)


up_xpr = upr[:40][['RAG-mean','TLX3-mean','TAP-mean']]

up_names = gene_names.loc[up_xpr.index]


with open('tracks/TLX3vsRAG_up1000.txt', 'wb') as fp:
    up600 = list(gene_names.loc[upr[:1020].index]['Gene_name'])
    up600 = [x for x in up600 if str(x) != 'nan']
    fp.write("\n".join(up600))

with open('tracks/TLX3vsRAG_dn1000.txt', 'wb') as fp:
    dn600 = list(gene_names.loc[dnr[-1020:].index]['Gene_name'])
    dn600 = [x for x in dn600 if str(x) != 'nan']
    fp.write("\n".join(dn600))


up_xpr=up_xpr.join(up_names)
up_xprR = up_xpr.set_index('Gene_name')


# Show names of genes in graph
num_g = 10
up10 = upr[:num_g+1][['log2FoldChange','padj','RAG-mean','TLX3-mean','TAP-mean']]
up10 = up10.join(gene_names.loc[up10.index]).dropna()

dn10 = dnr[-num_g:][['log2FoldChange','padj','RAG-mean','TLX3-mean','TAP-mean']]
dn10 = dn10.join(gene_names.loc[dn10.index]).dropna()

## ---- Pictures --------


fig = plt.figure(figsize=(10, 6))


## TLX3 gene
#~ ENSMUST00000037746 : TLX3
#~ print tbl.data.loc['ENSMUST00000037746']
with plt.style.context('seaborn-talk'):
    #fig, ax1 = plt.subplots()
    ax1 = fig.add_subplot(122)
    expr_tlx = tbl.data.loc['ENSMUST00000037746'][['RAG-mean','TLX3-mean','TAP-mean']]
    cells = ('RAG','TLX3','TAP')
    x_pos = np.arange(len(cells))
    ax1.bar(x_pos, expr_tlx,align='center', color=['green','red','orange'])
    ax1.set_title('TLX3 gene expression')
    ax1.set_ylabel('log2(FPKM + 1)')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(cells)


ax2 = fig.add_subplot(121)
sns.heatmap(up_xprR, ax=ax2, cmap='rainbow', linewidths=0.5)

for txt in ax2.get_yticklabels():
        txt.set_rotation(0)

    
#~ with plt.style.context('seaborn-talk'):    
    #~ tbl.scatter(
        #~ x='RAG-mean',
        #~ y='TLX3-mean')
        #---------------- 
        #~ #xfunc=log2p1,
        #~ #yfunc=log2p1)

# Scatterplot of control vs knockdown FPKM
ax = tbl.scatter(
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
    general_kwargs=dict(marker='.', color='0.5', alpha=0.2, s=5, label='unchanged'),

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
#~ for i in range(num_g):
    #~ ax.text(up10['RAG-mean'][i],up10['TLX3-mean'][i],up10['Gene_name'][i].upper(), color='#da3b3a')
    #~ ax.text(dn10['RAG-mean'][i],dn10['TLX3-mean'][i],dn10['Gene_name'][i].upper(), color='#00748e')

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

top_axes = tbl.marginal.top_hists[-1]
top_axes.set_title('Differential expression, TLX3 vs RAG')

for ax in tbl.marginal.top_hists:
    ax.set_ylabel('No.\ntranscripts', rotation=0, ha='right', va='center', size=8)

for ax in tbl.marginal.right_hists:
    ax.set_xlabel('No.\ntranscripts', rotation=-90, ha='left', va='top', size=8)



## === Volcano plot

def nlog10(x):
    return -np.log10(x)

up2 = (tbl.log2FoldChange > 2.).values
dn2 = (tbl.log2FoldChange < -2.).values
un2 = ((tbl.log2FoldChange >= -2) & (tbl.log2FoldChange <= 2)).values

ax7 = tbl.scatter(
    x='log2FoldChange', 
    y = 'padj', 
    yfunc=nlog10,
    xlab='log FC',
    ylab='-log10($p_{adj}$)',
    general_kwargs=dict(marker='.', color='0.5', alpha=0.2, s=5, label='unchanged'),
    genes_to_highlight=[
    (up2,dict(color='#da3b3a', alpha=0.8,label='upregulated')),
    (dn2,dict(color='#00748e', alpha=0.8,label='downregulated'))
    ])

ax7.set_title('Volcano plot, TLX3 vs RAG')
#~ ax7.axhline(2000, color='#da3b3a', linewidth=1)
ax7.axvline(-2, color='#da3b3a', linewidth=1)
ax7.axvline(2, color='#da3b3a', linewidth=1)

# Plot text of top genes
for i in range(num_g):
    ax7.text(up10['log2FoldChange'][i],nlog10(up10['padj'][i]+1e-300),up10['Gene_name'][i].upper(), color='#da3b3a')
    ax7.text(dn10['log2FoldChange'][i],nlog10(dn10['padj'][i]+1e-300),dn10['Gene_name'][i].upper(), color='#00748e')


# Get handles and labels, and then reverse their order
handles, legend_labels = ax7.get_legend_handles_labels()
handles = handles[::-1]
legend_labels = legend_labels[::-1]

# Draw a legend using the flipped handles and labels.
leg = ax7.legend(handles,
          legend_labels,
          loc=(.95, .95),
          fontsize=9,
          scatterpoints=1,
          borderpad=0.1,
          handletextpad=0.05, 
          frameon=False,
          title='Transcripts',
          );

leg.get_title().set_weight('bold')

#plt.savefig("results01/VolcannoRNA_seq_Expess.pdf", format="pdf")

## ===== TCRalpha expression
#~ ENSMUST00000103740  2.257974  11.859227  11.729335           Trac
#~ ENSMUST00000103567  1.191116   6.974967   0.843533          Trav1
#~ ENSMUST00000103638  0.344003   5.884605   3.590531        Trav6-3


with plt.style.context('seaborn-talk'):
    fig2 = plt.figure(figsize=(16, 6))
    ax5 = fig2.add_subplot(131)
    ac =up_xprR.loc['Trac'][['RAG-mean','TLX3-mean','TAP-mean']]
    ax5.bar(x_pos, ac,align='center', color=['green','red','orange'])
    ax5.set_title('Trac gene expression')
    ax5.set_ylabel('log2(FPKM + 1)')
    ax5.set_xticks(x_pos)
    ax5.set_xticklabels(cells)
    
    ax3 = fig2.add_subplot(132)
    aa =up_xprR.loc['Trav6-3'][['RAG-mean','TLX3-mean','TAP-mean']]
    ax3.bar(x_pos, aa,align='center', color=['green','red','orange'])
    ax3.set_title('Trav6-3 gene expression')
    ax3.set_ylabel('log2(FPKM + 1)')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(cells)


    ax4 = fig2.add_subplot(133)
    ab =up_xprR.loc['Trav1'][['RAG-mean','TLX3-mean','TAP-mean']]
    ax4.bar(x_pos, ab,align='center', color=['green','red','orange'])
    ax4.set_title('Trav1 gene expression')
    ax4.set_ylabel('log2(FPKM + 1)')
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(cells)



# Histogram padj
fig21 = plt.figure(figsize=(8, 6))
ax51 = fig21.add_subplot(111)
ax51.hist(tbl['padj'], bins=50)
ax51.set_title('Histogram of P-value')
ax51.set_ylabel('Frequency')
ax51.set_xlabel('P-value')

plt.show()




# ---------------------------------------------------------
# Re-align the ResultsTables to match the GTF file
#tbl = tbl.reindex_to(features, attribute='transcript_id')



# Get gene annotations
#~ mm9 = join('tracks', 'mm9.gtf')


#~ cmap='bwr',
#~ cmap='seismic'
#~ cmap='bwr'
#~ [u'seaborn-darkgrid', u'seaborn-notebook', u'classic', u'seaborn-ticks', 
#~ u'grayscale', u'bmh', u'seaborn-talk', u'dark_background', u'ggplot', 
#~ u'fivethirtyeight', u'seaborn-colorblind', u'seaborn-deep', 
#~ u'seaborn-whitegrid', u'seaborn-bright', u'seaborn-poster', 
#~ u'seaborn-muted', u'seaborn-paper', u'seaborn-white', 
#~ u'seaborn-pastel', u'seaborn-dark', u'seaborn-dark-palette']

























#~ key_single = ['tlx-H3K27me3']

#~ x3 = np.linspace(-3000, 3000, 100)

#~ old = [
    #~ 'rag-POLII',
    #~ 'rag-H3K27me3',
    #~ 'rag-H3K4me2',
    #~ 'rag-H3K4me3',
    #~ 'rag-H3K9ac',
    #~ 'rag-H3K36me3',
    #~ 'rag-H3K9me3']


#~ for chname in arrays3K.keys(): #key_single: #arrays3K.keys():
    #~ chnameR = 'rag'+chname[3:]
    
    #~ charr = arrays3K[chname]
    #~ charrR = arrays3K_rag[chnameR]
    
    #~ if any(chnameR in string for string in old):
        #~ k = (abs(np.max(charr)-np.min(charr)))/(abs(np.max(charrR)-np.min(charrR)))
        #~ charrR = charrR*k
        #~ print 'k_'+chnameR+' = ', k
    
    
    #~ up = charr[(tlx_vs_rag_table.log2FoldChange > 1).values, :]
    #~ down = charr[(tlx_vs_rag_table.log2FoldChange < -1).values, :]
    #~ unchn =  charr[((tlx_vs_rag_table.log2FoldChange >= -1) & (tlx_vs_rag_table.log2FoldChange <= 1)).values, :]
    
    #~ upR = charrR[(tlx_vs_rag_table.log2FoldChange > 1).values, :]
    #~ downR = charrR[(tlx_vs_rag_table.log2FoldChange < -1).values, :]
    #~ unchnR =  charrR[((tlx_vs_rag_table.log2FoldChange >= -1) & (tlx_vs_rag_table.log2FoldChange <= 1)).values, :]
    
    #~ xx = x3 
    
    #~ fig = metaseq.plotutils.imshow(
        #~ charr,
        #~ x=xx,
        #~ figsize=(6, 14),
        #~ vmin=5, vmax=99,  percentile=True,
        #~ line_kwargs=dict(color='k', label='All'),
        #~ fill_kwargs=dict(color='k', alpha=0.3),
        #~ sort_by=charr.mean(axis=1),
        #~ # Default was (3,1); here we add another number
        #~ height_ratios=(3, 1, 1)
    #~ )

    #~ # `fig.gs` contains the `matplotlib.gridspec.GridSpec` object,
    #~ # so we can now create the new axes.
    #~ bottom_axes = plt.subplot(fig.gs[2, 0])
    
    #~ # Signal over TSSs of transcripts that were activated upon knockdown.
    #~ metaseq.plotutils.ci_plot(
        #~ xx,
        #~ up,
        #~ line_kwargs=dict(color='#fe9829', label='up'),
        #~ fill_kwargs=dict(color='#fe9829', alpha=0.3),
        #~ ax=bottom_axes)

    #~ # Signal over TSSs of transcripts that were repressed upon knockdown
    #~ metaseq.plotutils.ci_plot(
        #~ xx,
        #~ down,
        #~ line_kwargs=dict(color='#8e3104', label='down'),
        #~ fill_kwargs=dict(color='#8e3104', alpha=0.3),
        #~ ax=bottom_axes)

    #~ # Signal over TSSs tof transcripts that did not change upon knockdown
    #~ metaseq.plotutils.ci_plot(
        #~ xx,
        #~ unchn,
        #~ line_kwargs=dict(color='.5', label='unchanged'),
        #~ fill_kwargs=dict(color='.5', alpha=0.3),
        #~ ax=bottom_axes);
    
    #~ # Clean up redundant x tick labels, and add axes labels
    #~ fig.line_axes.set_xticklabels([])
    #~ fig.array_axes.set_xticklabels([])
    #~ fig.line_axes.set_ylabel('Average\nenrichement')
    #~ fig.array_axes.set_ylabel('Transcripts on all chromosomes')
    #~ bottom_axes.set_ylabel('Average\nenrichment')
    #~ bottom_axes.set_xlabel('Distance from TSS (bp)')
    #~ fig.cax.set_ylabel('Enrichment')
    #~ fig.array_axes.set_title(chname)

    #~ # Add the vertical lines for TSS position to all axes
    #~ for ax in [fig.line_axes, fig.array_axes, bottom_axes]:
        #~ ax.axvline(0, linestyle=':', color='k')

    #~ # Nice legend
    #~ bottom_axes.legend(loc='best', frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)
    #~ fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)

    #~ #plt.savefig("results/"+chname+"-1K.svg", format="svg")
    #~ #plt.savefig("results/"+chname+"-1K.pdf", format="pdf")
    #~ plt.savefig("results/"+chname+"-3K.pdf", format="pdf")
    
    #~ fig = metaseq.plotutils.imshow(
        #~ # Same as before...
        #~ charrR,
        #~ x=xx,
        #~ figsize=(6, 14),
        #~ vmin=5, vmax=99,  percentile=True,
        #~ line_kwargs=dict(color='k', label='All'),
        #~ fill_kwargs=dict(color='k', alpha=0.3),
        #~ sort_by=charr.mean(axis=1),
        #~ # Default was (3,1); here we add another number
        #~ height_ratios=(3, 1, 1)
    #~ )

    #~ # `fig.gs` contains the `matplotlib.gridspec.GridSpec` object,
    #~ # so we can now create the new axes.
    #~ bottom_axes = plt.subplot(fig.gs[2, 0])
    #~ # Signal over TSSs of transcripts that were activated upon knockdown.
    #~ metaseq.plotutils.ci_plot(
        #~ xx,
        #~ upR,
        #~ line_kwargs=dict(color='#fe9829', label='up'),
        #~ fill_kwargs=dict(color='#fe9829', alpha=0.3),
        #~ ax=bottom_axes)

    #~ # Signal over TSSs of transcripts that were repressed upon knockdown
    #~ metaseq.plotutils.ci_plot(
        #~ xx,
        #~ downR,
        #~ line_kwargs=dict(color='#8e3104', label='down'),
        #~ fill_kwargs=dict(color='#8e3104', alpha=0.3),
        #~ ax=bottom_axes)

    #~ # Signal over TSSs tof transcripts that did not change upon knockdown
    #~ metaseq.plotutils.ci_plot(
        #~ xx,
        #~ unchnR,
        #~ line_kwargs=dict(color='.5', label='unchanged'),
        #~ fill_kwargs=dict(color='.5', alpha=0.3),
        #~ ax=bottom_axes);
    
    #~ # Clean up redundant x tick labels, and add axes labels
    #~ fig.line_axes.set_xticklabels([])
    #~ fig.array_axes.set_xticklabels([])
    #~ fig.line_axes.set_ylabel('Average\nenrichement')
    #~ fig.array_axes.set_ylabel('Transcripts on all chromosomes sorted by TLX3 signal')
    #~ bottom_axes.set_ylabel('Average\nenrichment')
    #~ bottom_axes.set_xlabel('Distance from TSS (bp)')
    #~ fig.cax.set_ylabel('Enrichment')
    #~ fig.array_axes.set_title(chnameR)

    #~ # Add the vertical lines for TSS position to all axes
    #~ for ax in [fig.line_axes, fig.array_axes, bottom_axes]:
        #~ ax.axvline(0, linestyle=':', color='k')

    #~ # Nice legend
    #~ bottom_axes.legend(loc='best', frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)
    #~ fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)

    #~ #plt.savefig("results/"+chnameR+"-1K.svg", format="svg")
    #~ #plt.savefig("results/"+chnameR+"-1K.pdf", format="pdf")
    #~ plt.savefig("results/"+chnameR+"-tlx-3K.pdf", format="pdf")

#~ plt.show()
