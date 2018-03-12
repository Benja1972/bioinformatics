#!/usr/bin/env python2

import numpy as np
import yaml
import os
from os.path import join 
import pybedtools as pb
import gffutils
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import matplotlib.pyplot as plt
import matplotlib
import metaseq
import pandas as pd
import gseapy as gp

from metaseq.results_table import ResultsTable, DESeq2Results


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
features3K, arrays3K = metaseq.persistence.load_features_and_arrays(prefix='tlx-list3K')

# Expression table
tbl_f = join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt')

tbl = ResultsTable(tbl_f, import_kwargs=dict(index_col=0))
tbl = tbl.reindex_to(features3K, attribute='transcript_id')

# Add ChiP-seq signals to table
tbl.data['tlx3_ChiP-seq_3Kb'] = arrays3K['tlx-TLX3'].mean(axis=1)

tbl_g = tbl[(tbl.padj < 0.05)].dropna()

# === Load gene names 
names = pd.read_table("tracks/UCSC_mm9_transcripID_to_geneSymbol.sort.txt", 
                           index_col=0,
                           names=["Geneid", "NAME"])


names = names.loc[tbl_g.index]
assert names.shape[0] == tbl_g.shape[0]

tbl_n=names.join(tbl_g, how ='right')


# === 500 gene set for GSEA in gsea_chip-seq_py3.py
tbl_s = tbl_n.sort_values('tlx3_ChiP-seq_3Kb', axis=0, ascending=False)

gm = ['ChiP_seq_TLX3_TLX3top500', 'ChiP_seq_TLX3_TLX3top500']  + list(tbl_s['NAME'].head(500).dropna().str.upper())

#~ with open('tracks/ChiP_seq_TLX3_TLX3top500.gmt', 'w') as fp:
    #~ fp.write("\t".join(gm))


# === Pictures

chname = 'tlx-TLX3'
xx = np.linspace(-3000, 3000, 100)
charr = arrays3K[chname]



up = charr[(tbl.log2FoldChange > 1).values, :]
down = charr[(tbl.log2FoldChange < -1).values, :]
unchn =  charr[((tbl.log2FoldChange >= -1) & (tbl.log2FoldChange <= 1)).values, :]




fig = metaseq.plotutils.imshow(
    charr,
    x=xx,
    figsize=(6, 14),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=charr.mean(axis=1),
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

ypos = len(charr) - 500

print ypos

fig.array_axes.axhline(ypos, color='r', linewidth=1)
#~ plt.savefig("results02/"+chname+"-3K.pdf", format="pdf")


plt.show()

#plt.savefig("results/"+chname+"-1K.svg", format="svg")
#plt.savefig("results/"+chname+"-1K.pdf", format="pdf")

