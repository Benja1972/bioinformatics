
import numpy as np
import pandas as pd
import sklearn.cluster as sc
import seaborn as sns
import matplotlib.pyplot as plt

## project settings ===
from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')
## =====================


def log2p1(x):
    return np.log2(x + 1)
    
# test

# === Load expression table 
tbl = pd.read_table(join(DATADIR,'tracks', 'TLX3vsRAG-results_genes.txt'), index_col=0)


# Filter genes (Note: this filter remove microRNA expression)
tbl = tbl[(tbl.padj < 0.05)].dropna()


# === Load gene names 
names = pd.read_table(join(DATADIR,"tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt"), 
                           index_col=0,
                           header=0,
                           names=['GeneID', 'TransID', 'Gene_name'])


names = names.drop('TransID', axis=1).drop_duplicates()
names = names.loc[tbl.index]
assert names.shape[0] == tbl.shape[0]

tbl=names.join(tbl, how ='right')
tbl=tbl.sort_values('log2FoldChange', axis=0, ascending=False)
tbn = tbl[['Gene_name', 'R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P']]
tbn.set_index(keys=tbn.columns[0], inplace=True)



df = tbn.head(50)
mat = df.as_matrix()

#===Test========
#f, ax2 = plt.subplots(figsize=(9, 6))
#fig=plt.figure()
#fig = plt.figure(figsize=(6, 11)) #11
#ax2 = fig.add_axes()
#sns.clustermap(log2p1(df), cmap='RdBu_r', col_cluster=False, ax = ax2)
sns.clustermap(log2p1(df), cmap='RdBu_r', col_cluster=False, figsize=(6, 12))
fig = plt.gcf()
ax2 = fig.axes[2]

for txt in ax2.get_yticklabels():
        txt.set_rotation(0)
for txt in ax2.get_xticklabels():
        txt.set_rotation(90)
#===============

cl = 6
km = sc.KMeans(n_clusters=cl)
km.fit(mat)
labels = km.labels_

results = pd.DataFrame(data=labels, columns=['cluster'], index=df.index)

dfc = pd.concat([log2p1(df), results], axis=1)

dfs = dfc.sort_values('cluster', axis=0, ascending=True)

k = log2p1(mat.max())/(cl-1)
dfs['cluster'] = k*dfs['cluster']


# ==== Figures
ttl = 'Cluster test'
fig1 = plt.figure(figsize=(6, 11)) #11
#ax1 = fig1.add_subplot(111)
ax1 = sns.heatmap(dfs, 
            cmap='RdBu_r',  
            linewidths=0.000)
ax1.set_title('diffExpressed')
for txt in ax1.get_yticklabels():
        txt.set_rotation(0)
for txt in ax1.get_xticklabels():
        txt.set_rotation(90)
    
fig1.suptitle(ttl, size='x-large')


#sns.heatmap(dfs, cmap='RdBu_r')





plt.show()
