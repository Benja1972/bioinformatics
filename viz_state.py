#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import pandas
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import ImageGrid
import os
import yaml

wc = dict(
        state = 14,
        cell1 = 'TLX3',
        cell2 = 'RAG')

files = dict(
        emissions='../models/{state}/emissions_{state}.txt'.format(**wc),
        enrich1='../models/{state}/{cell1}_enrichment.txt'.format(**wc),
        enrich2='../models/{state}/{cell2}_enrichment.txt'.format(**wc))




config = yaml.load(open('../config.yaml'))
order = config["marks"]

#annot14 = [14,12,11,13,8,9,3,4,10,7,6,5,2,1]
#annot14 = [14,12,11,13,8,9,7,6,5,2,1,10,4,3]
annot14 = [14,12,13,11,8,9,7,6,5,2,1,10,4,3]

print files
print order


dfs = {}
for k, v in files.items():
    df = pandas.read_table(v, index_col=0)
    if 'enrich' in k:
        df = df.fillna(0)
    dfs[k] = df



#fig = plt.figure(figsize=(8, 8))
fig = plt.figure(figsize=(15, 6))
grid = ImageGrid(fig, 111, nrows_ncols=(1, 3), # nrows_ncols=(1, 5),
                 axes_pad=0.1,
                 add_all=True,
                 label_mode="R")


cmaps = {
    'emissions': cm.Blues,
    'enrich1': cm.Greys,
    'enrich2': cm.Greys,
}


# TLX3
ax0 = grid[0]
k = 'emissions'
em = dfs[k]
em = em.ix[:, order]
em = em.ix[annot14,:]
sns.heatmap(em, ax=ax0, cbar=False, cmap=cmaps[k], linewidths=0.5)
ax0.set_title(k)
ax0.patch.set_edgecolor('0.5')
ax0.patch.set_linewidth(2)
ax0.set_ylabel('State')
for txt in ax0.get_xticklabels():
    txt.set_rotation(90)
#~ for txt in ax0.get_xticklabels():
    #~ txt.set_visible(False)
#~ ax0.xaxis.get_label().set_visible(False)

ax1 = grid[1]
k = 'enrich1'
enr1 = dfs[k]
enr1.columns = [i.replace('.bed.gz', '').replace('.mm9', '') for i in enr1.columns]
enr1 = enr1.drop(['genes','tsses'], axis=1)
enr1 = enr1.drop('Base')
enr1 = enr1.ix[np.array(annot14)-1,:]
enr1 = enr1 / enr1.max()
sns.heatmap(enr1, ax=ax1, cbar=False, cmap=cmaps[k], linewidths=0.5)
ax1.set_title('enrichment')
ax1.patch.set_edgecolor('0.5')
ax1.patch.set_linewidth(2)
for txt in ax1.get_yticklabels():
    txt.set_visible(False)
ax1.yaxis.get_label().set_visible(False)
for txt in ax1.get_xticklabels():
    txt.set_rotation(90)


#~ for txt in ax1.get_yticklabels():
    #~ txt.set_visible(False)
#~ for txt in ax1.get_xticklabels():
    #~ txt.set_visible(False)
#~ ax1.yaxis.get_label().set_visible(False)
#~ for txt in ax1.get_xticklabels():
    #~ txt.set_rotation(90)
#~ ax1.xaxis.get_label().set_visible(False)


# RAG
#~ ax2 = grid[2]
#~ k = 'emissions'
#~ em = dfs[k]
#~ em = em.ix[:, order]
#~ em = em.ix[annot14,:]
#~ sns.heatmap(em, ax=ax2, cbar=False, cmap=cmaps[k], linewidths=0.5)
#~ ax2.set_title('')
#~ ax2.patch.set_edgecolor('0.5')
#~ ax2.patch.set_linewidth(2)
#~ ax2.set_ylabel('State')
#~ for txt in ax2.get_xticklabels():
    #~ txt.set_rotation(90)

ax3 = grid[2]
k = 'enrich2'
enr2 = dfs[k]
enr2.columns = [i.replace('.bed.gz', '').replace('.mm9', '') for i in enr2.columns]
enr2 = enr2.drop(['genes','tsses'], axis=1)
enr2 = enr2.drop('Base')
enr2 = enr2.ix[np.array(annot14)-1,:]
enr2 = enr2 / enr2.max()
sns.heatmap(enr2, ax=ax3, cbar=False, cmap=cmaps[k], linewidths=0.5)
ax3.set_title('enrichment')
ax3.patch.set_edgecolor('0.5')
ax3.patch.set_linewidth(2)
for txt in ax3.get_yticklabels():
    txt.set_visible(False)
ax3.yaxis.get_label().set_visible(False)
for txt in ax3.get_xticklabels():
    txt.set_rotation(90)

fig.tight_layout()
fig.savefig('TLX3-RAG_anno_{state}_ver4.eps'.format(**wc))

plt.show()

#~ for i, k in enumerate([
    #~ 'emissions',
    #~ 'enrichment',
    #~ 'uniform_enrichment',
    #~ 'transitions',
    #~ 'comparison'
#~ ]):
    #~ ax = grid[i]
    #~ v = dfs[k].copy()
    #~ v.columns = [i.replace('.bed.gz', '').replace('.bed', '').replace('.txt', '') for i in v.columns]
    #~ if k == 'uniform_enrichment':
        #~ v = v.drop('Base')
        #~ v = v.drop('Genome %', axis=1)
    #~ if k == 'enrichment':
        #~ v = v.drop('Base')
        #~ v = v / v.max()
    #~ if k == 'emissions':
        #~ v = v.ix[:, order]
    #~ #sns.heatmap(v, annot=True, ax=ax, cbar=False, cmap=cmaps[k], linewidths=0.5)
    #~ sns.heatmap(v, ax=ax, cbar=False, cmap=cmaps[k], linewidths=0.5)
    #~ ax.set_title(k)
    #~ ax.patch.set_edgecolor('0.5')
    #~ ax.patch.set_linewidth(2)
    #~ if i > 0:
        #~ for txt in ax.get_yticklabels():
            #~ txt.set_visible(False)
        #~ ax.yaxis.get_label().set_visible(False)
    #~ else:
        #~ ax.set_ylabel('State')
    #~ for txt in ax.get_xticklabels():
        #~ txt.set_rotation(90)

#~ fig.tight_layout()
#~ fig.suptitle('{cell}, {state}-states'.format(**wc), weight='bold', size=20)
#~ try:
    #~ #fig.savefig(snakemake.output.png)
    #~ fig.savefig(snakemake.output.pdf)
    #~ fig.savefig(snakemake.output.eps)
#~ except NameError:
    #~ plt.show()
