import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns

import networkx as nx
import numpy as np


import EnrichRLib as erl


def cluser_enrich(enr,gl,pval=0.05, top_clusters=20):
    #2= Filter terms by p-Val
    enr = enr[enr['p-Val']<pval]

    #3= Make claster by kappa coeff
    enr = erl.cluster(set(gl), enr, deep=2)

    #3-1= Filter top clusters
    enr = enr[enr['cluster']<top_clusters]

    #4= Make clustered geneset 
    gs_clust,nt_cl = erl.cluster_genset(enr)

    #5= Enrich clustered geneset
    enr_clust = erl.enrich(gl,gs_clust)
    
    # deduplicate index -- TODO!!! in package
    nt_cl = nt_cl.loc[~nt_cl.index.duplicated(keep='first')]

    #6= Add cluster to table
    enr_clust = pd.concat([enr_clust,nt_cl.loc[enr_clust.index]['cluster']],axis=1, sort=False)

    #7= Make graphs
    G_gs = erl.make_graph_n(gl,enr, kappa=0.4)
    G_cl = erl.make_graph_n(gl,enr_clust, kappa=0.01)

    #8= Draw graphs
    erl.draw_graph(G_gs, spring=150)
    erl.draw_direct(G_cl)



    #9= Draw barplot for clustered terms
    enr.sort_values('cluster', axis=0, inplace = True)

    cm = ('tab20' if max(enr['cluster'])>10  else 'tab10')

    f, ax = plt.subplots(figsize=(8, 12))
    sns.barplot(y=enr.index,
                x='-log10(p-Val)',
                ax = ax, 
                hue ='cluster',
                dodge=False,
                data = enr,
                palette = cm)
    ax.set_title('Top terms in clusters ')



## =====================
## =====================


gl = [  'CD3E', 'BLK', 'PTPN22', 'PAG1', 'CTLA4', 'PIK3CD', 'LAT2', 'CSK', 
        'CD247', 'CD3G', 'THEMIS', 'PSMB8', 'LCP2', 'GATA3', 'LAT', 'SLA2', 
        'SKAP1', 'TRAT1', 'BCL2', 'CD3D', 'THY1','RUNX1', 'BLK', 'PTPN22']


lib_dir='_tmp'


gss = [ 
       'GO_Biological_Process_2018',
       #~ 'GO_Cellular_Component_2018',
       #~ 'GO_Molecular_Function_2018',
       'KEGG_2016',
       #'Reactome_2016'
       ]



#1= Enrich  all terms
enr = erl.enrich_gs(gl,gss, path_lib=lib_dir)

#2= Filter terms by p-Val
enr = enr[enr['p-Val']<0.005]

#3= Make claster by kappa coeff
enr = erl.cluster(set(gl), enr, deep=2)

#3-1= Filter top clusters
top_clusters = 18
enr = enr[enr['cluster']<top_clusters]

#4= Make clustered geneset 
gs_clust,nt_cl = erl.cluster_genset(enr)

#5= Enrich clustered geneset
enr_clust = erl.enrich(gl,gs_clust)

#6= Add cluster to table
enr_clust = pd.concat([enr_clust,nt_cl.loc[enr_clust.index]['cluster']],axis=1, sort=False)

#7= Make graphs
G_gs = erl.make_graph_n(gl,enr, kappa=0.4)
G_cl = erl.make_graph_n(gl,enr_clust, kappa=0.01)

#8= Draw graphs
erl.draw_graph(G_gs, spring=150)
erl.draw_direct(G_cl)



#9= Draw barplot for clustered terms
enr.sort_values('cluster', axis=0, inplace = True)

cm = ('tab20' if max(enr['cluster'])>10  else 'tab10')

f, ax = plt.subplots(figsize=(8, 12))
sns.barplot(y=enr.index,
            x='-log10(p-Val)',
            ax = ax, 
            hue ='cluster',
            dodge=False,
            data = enr,
            palette = cm)
ax.set_title('Top terms in clusters ')
plt.tight_layout()




#~ cluser_enrich(enr,gl)







plt.show()
