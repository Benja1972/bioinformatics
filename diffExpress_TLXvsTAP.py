


import numpy as np
import yaml
import os
from os.path import join 
import pybedtools as pb
import gffutils
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import matplotlib.pyplot as plt
import metaseq
import pandas as pd
from metaseq.results_table import ResultsTable, DESeq2Results


# Tweak some font settings so the results look nicer
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12


# RNA-seq expression tables


#tbl_f = 'TLX3vsRAG_DESeq2-results.txt'
# tbl_f = join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt')
tbl_f = join('tracks', 'TAPvsTLX3-results.txt')

tbl = ResultsTable(tbl_f, import_kwargs=dict(index_col=0))
tbl.data['TLX3-mean'] = (tbl.data['TLX3.1_1']+tbl.data['TLX3.1_5']+tbl.data['TLX3.1_P'])/3.
tbl.data['TAP-mean'] = (tbl.data['TAP']+tbl.data['TAP1B']+tbl.data['TAP2B'])/3.


# Get gene annotations
mm9 = join('tracks', 'mm9.gtf')

up = (tbl.log2FoldChange > 1.).values
dn = (tbl.log2FoldChange < -1.).values
unch =  ((tbl.log2FoldChange >= -1) & (tbl.log2FoldChange <= 1)).values

#fc =tbl.data.dropna().sort_values('log2FoldChange', axis=0, ascending=False)

upr = tbl[(tbl.log2FoldChange > 1.) & (tbl.padj < 0.05)]
dnr = tbl[(tbl.log2FoldChange < -1.) & (tbl.padj < 0.05)]



#~ print list(fc.head(20).index)

#~ print list(fc.tail(20).index)



## ---- Pictures --------

def log2p1(x):
    return np.log2(x + 1)
    
    
tbl.scatter(
    x='TLX3-mean',
    y='TAP-mean',
    #---------------- 
    xfunc=log2p1,
    yfunc=log2p1)

#~ # Scatterplot of control vs knockdown FPKM
#~ ax = tbl.scatter(
    #~ x='TLX3-mean',
    #~ y='TAP-mean',
    #~ #----------------
    #~ xfunc=log2p1,
    #~ yfunc=log2p1,
    #~ xlab='TLX3, log2(FPKM + 1)',
    #~ ylab='TAP, log2(FPKM + 1)',
    #~ one_to_one=dict(color='k', linestyle=':'),
    #~ marginal_histograms=False,
    #~ # add the "unchanged" label
    #~ general_kwargs=dict(marker='.', color='0.5', alpha=0.2, s=5, label='unchanged'),

    #~ genes_to_highlight=[
    #~ (
     #~ up,
     #~ dict(
          #~ color='#da3b3a', alpha=0.8,
          #~ marginal_histograms=True,
          #~ xhist_kwargs=dict(bins=50, linewidth=0),
          #~ yhist_kwargs=dict(bins=50, linewidth=0),

          #~ # add label
          #~ label='upregulated',
          #~ )
     #~ ),
    #~ (
     #~ dn,
     #~ dict(
          #~ color='#00748e', alpha=0.8,
          #~ marginal_histograms=True,
          #~ xhist_kwargs=dict(bins=50, linewidth=0),
          #~ yhist_kwargs=dict(bins=50, linewidth=0),

          #~ # add label
          #~ label='downregulated'
          #~ )
     #~ )
    #~ ],
#~ );


#~ # Get handles and labels, and then reverse their order
#~ handles, legend_labels = ax.get_legend_handles_labels()
#~ handles = handles[::-1]
#~ legend_labels = legend_labels[::-1]

#~ # Draw a legend using the flipped handles and labels.
#~ leg = ax.legend(handles,
          #~ legend_labels,

          #~ # These values may take some tweaking.
          #~ # By default they are in axes coordinates, so this means
          #~ # the legend is slightly outside the axes.
          #~ loc=(1.01, 1.05),

          #~ # Various style fixes to default legend.
          #~ fontsize=9,
          #~ scatterpoints=1,
          #~ borderpad=0.1,
          #~ handletextpad=0.05,
          #~ frameon=False,
          #~ title=' Transcripts',
          #~ );

#~ # Adjust the legend title after it's created
#~ leg.get_title().set_weight('bold')

#~ top_axes = tbl.marginal.top_hists[-1]
#~ top_axes.set_title('Differential expression, TAP vs TLX3');

#~ for ax in tbl.marginal.top_hists:
    #~ ax.set_ylabel('No.\ntranscripts', rotation=0, ha='right', va='center', size=8)

#~ for ax in tbl.marginal.right_hists:
    #~ ax.set_xlabel('No.\ntranscripts', rotation=-90, ha='left', va='top', size=8)





plt.show()






#~ print fc.head(20).index
#~ Index([u'ENSMUST00000037746', u'ENSMUST00000103276', u'ENSMUST00000119402',
       #~ u'ENSMUST00000103740', u'ENSMUST00000045583', u'ENSMUST00000071281',
       #~ u'ENSMUST00000053225', u'ENSMUST00000034751', u'ENSMUST00000087033',
       #~ u'ENSMUST00000058475', u'ENSMUST00000103676', u'ENSMUST00000103675',
       #~ u'ENSMUST00000057634', u'ENSMUST00000163134', u'ENSMUST00000028533',
       #~ u'ENSMUST00000085138', u'ENSMUST00000152439', u'ENSMUST00000095015',
       #~ u'ENSMUST00000147658', u'ENSMUST00000022186'],
      #~ dtype='object', name=u'Gene')


#~ print fc.tail(20).index
#~ Index([u'ENSMUST00000027377', u'ENSMUST00000076623', u'ENSMUST00000069476',
       #~ u'ENSMUST00000029076', u'ENSMUST00000081035', u'ENSMUST00000124830',
       #~ u'ENSMUST00000103284', u'ENSMUST00000065167', u'ENSMUST00000000314',
       #~ u'ENSMUST00000031273', u'ENSMUST00000053218', u'ENSMUST00000034026',
       #~ u'ENSMUST00000046332', u'ENSMUST00000048860', u'ENSMUST00000046384',
       #~ u'ENSMUST00000046285', u'ENSMUST00000092163', u'ENSMUST00000023952',
       #~ u'ENSMUST00000127786', u'ENSMUST00000023803'],
      #~ dtype='object', name=u'Gene')







# ---------------------------------------------------------
# Re-align the ResultsTables to match the GTF file
#tbl = tbl.reindex_to(features, attribute='transcript_id')


























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
    #~ fig.array_axes.set_ylabel('Transcripts on all chromosomes sorted by TAP signal')
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
