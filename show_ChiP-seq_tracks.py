#!/usr/bin/env python


# # RAG vs TLX3 combined ChiP-seq RNA-seq analysis

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


# Paths
#tracks/annot_tracks/references
ref_dir = 'tracks/annot_tracks/referencesmm9'
data_dir = 'tracks'
rel_path = "/home/sergio/media"


# Files list

file_t = open(join(data_dir,"TLX3_list.yaml"), "r")
file_r = open(join(data_dir,"RAG_list.yaml"), "r")

tlx_lst = yaml.load(file_t)
rag_lst = yaml.load(file_r)


# Martix of signals

features, arrays = metaseq.persistence.load_features_and_arrays(prefix='tlx-list1K')
features3K, arrays3K = metaseq.persistence.load_features_and_arrays(prefix='tlx-list3K')
features_rag, arrays_rag = metaseq.persistence.load_features_and_arrays(prefix='rag-list1K')
features3K_rag, arrays3K_rag = metaseq.persistence.load_features_and_arrays(prefix='rag-list3K')
#features5K, arrays5K = metaseq.persistence.load_features_and_arrays(prefix='tlx-vs-rag_tlx3_5Kb')





# Tweak some font settings so the results look nicer
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12




# RNA-seq expression tables

from metaseq.results_table import ResultsTable, DESeq2Results


tlx_tabf = rel_path+tlx_lst["table"][0]

tlx_vs_rag_table = ResultsTable(tlx_tabf, import_kwargs=dict(index_col=0))

# ---------------------------------------------------------
# Re-align the ResultsTables to match the GTF file
tlx_vs_rag_table = tlx_vs_rag_table.reindex_to(features, attribute='transcript_id')

key_single = ['tlx-H3K27me3']

x3 = np.linspace(-3000, 3000, 100)

old = [
    'rag-POLII',
    'rag-H3K27me3',
    'rag-H3K4me2',
    'rag-H3K4me3',
    'rag-H3K9ac',
    'rag-H3K36me3',
    'rag-H3K9me3']


for chname in key_single: #arrays3K.keys():
    chnameR = 'rag'+chname[3:]
    
    charr = arrays3K[chname]
    charrR = arrays3K_rag[chnameR]
    
    if any(chnameR in string for string in old):
        k = (abs(np.max(charr)-np.min(charr)))/(abs(np.max(charrR)-np.min(charrR)))
        charrR = charrR*k
        print 'k_'+chnameR+' = ', k
    
    
    up = charr[(tlx_vs_rag_table.log2FoldChange > 1).values, :]
    down = charr[(tlx_vs_rag_table.log2FoldChange < -1).values, :]
    unchn =  charr[((tlx_vs_rag_table.log2FoldChange >= -1) & (tlx_vs_rag_table.log2FoldChange <= 1)).values, :]
    
    upR = charrR[(tlx_vs_rag_table.log2FoldChange > 1).values, :]
    downR = charrR[(tlx_vs_rag_table.log2FoldChange < -1).values, :]
    unchnR =  charrR[((tlx_vs_rag_table.log2FoldChange >= -1) & (tlx_vs_rag_table.log2FoldChange <= 1)).values, :]
    
    xx = x3 
    
    fig = metaseq.plotutils.imshow(
        charr,
        x=xx,
        figsize=(6, 14),
        vmin=5, vmax=99,  percentile=True,
        line_kwargs=dict(color='k', label='All'),
        fill_kwargs=dict(color='k', alpha=0.3),
        sort_by=charr.mean(axis=1),
        # Default was (3,1); here we add another number
        height_ratios=(3, 1, 1)
    )

    # `fig.gs` contains the `matplotlib.gridspec.GridSpec` object,
    # so we can now create the new axes.
    bottom_axes = plt.subplot(fig.gs[2, 0])
    
    # Signal over TSSs of transcripts that were activated upon knockdown.
    metaseq.plotutils.ci_plot(
        xx,
        up,
        line_kwargs=dict(color='#fe9829', label='up'),
        fill_kwargs=dict(color='#fe9829', alpha=0.3),
        ax=bottom_axes)

    # Signal over TSSs of transcripts that were repressed upon knockdown
    metaseq.plotutils.ci_plot(
        xx,
        down,
        line_kwargs=dict(color='#8e3104', label='down'),
        fill_kwargs=dict(color='#8e3104', alpha=0.3),
        ax=bottom_axes)

    # Signal over TSSs tof transcripts that did not change upon knockdown
    metaseq.plotutils.ci_plot(
        xx,
        unchn,
        line_kwargs=dict(color='.5', label='unchanged'),
        fill_kwargs=dict(color='.5', alpha=0.3),
        ax=bottom_axes);
    
    # Clean up redundant x tick labels, and add axes labels
    fig.line_axes.set_xticklabels([])
    fig.array_axes.set_xticklabels([])
    fig.line_axes.set_ylabel('Average\nenrichement')
    fig.array_axes.set_ylabel('Transcripts on all chromosomes')
    bottom_axes.set_ylabel('Average\nenrichment')
    bottom_axes.set_xlabel('Distance from TSS (bp)')
    fig.cax.set_ylabel('Enrichment')
    fig.array_axes.set_title(chname)

    # Add the vertical lines for TSS position to all axes
    for ax in [fig.line_axes, fig.array_axes, bottom_axes]:
        ax.axvline(0, linestyle=':', color='k')

    # Nice legend
    bottom_axes.legend(loc='best', frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)
    fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)

    #plt.savefig("results/"+chname+"-1K.svg", format="svg")
    #plt.savefig("results/"+chname+"-1K.pdf", format="pdf")
    #plt.savefig("results/"+chname+"-3K.pdf", format="pdf")
    
    fig = metaseq.plotutils.imshow(
        # Same as before...
        charrR,
        x=xx,
        figsize=(6, 14),
        vmin=5, vmax=99,  percentile=True,
        line_kwargs=dict(color='k', label='All'),
        fill_kwargs=dict(color='k', alpha=0.3),
        sort_by=charr.mean(axis=1),
        # Default was (3,1); here we add another number
        height_ratios=(3, 1, 1)
    )

    # `fig.gs` contains the `matplotlib.gridspec.GridSpec` object,
    # so we can now create the new axes.
    bottom_axes = plt.subplot(fig.gs[2, 0])
    # Signal over TSSs of transcripts that were activated upon knockdown.
    metaseq.plotutils.ci_plot(
        xx,
        upR,
        line_kwargs=dict(color='#fe9829', label='up'),
        fill_kwargs=dict(color='#fe9829', alpha=0.3),
        ax=bottom_axes)

    # Signal over TSSs of transcripts that were repressed upon knockdown
    metaseq.plotutils.ci_plot(
        xx,
        downR,
        line_kwargs=dict(color='#8e3104', label='down'),
        fill_kwargs=dict(color='#8e3104', alpha=0.3),
        ax=bottom_axes)

    # Signal over TSSs tof transcripts that did not change upon knockdown
    metaseq.plotutils.ci_plot(
        xx,
        unchnR,
        line_kwargs=dict(color='.5', label='unchanged'),
        fill_kwargs=dict(color='.5', alpha=0.3),
        ax=bottom_axes);
    
    # Clean up redundant x tick labels, and add axes labels
    fig.line_axes.set_xticklabels([])
    fig.array_axes.set_xticklabels([])
    fig.line_axes.set_ylabel('Average\nenrichement')
    fig.array_axes.set_ylabel('Transcripts on all chromosomes sorted by TLX3 signal')
    bottom_axes.set_ylabel('Average\nenrichment')
    bottom_axes.set_xlabel('Distance from TSS (bp)')
    fig.cax.set_ylabel('Enrichment')
    fig.array_axes.set_title(chnameR)

    # Add the vertical lines for TSS position to all axes
    for ax in [fig.line_axes, fig.array_axes, bottom_axes]:
        ax.axvline(0, linestyle=':', color='k')

    # Nice legend
    bottom_axes.legend(loc='best', frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)
    fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)

    #plt.savefig("results/"+chnameR+"-1K.svg", format="svg")
    #plt.savefig("results/"+chnameR+"-1K.pdf", format="pdf")
    #plt.savefig("results/"+chnameR+"-tlx-3K.pdf", format="pdf")

plt.show()
