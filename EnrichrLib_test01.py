

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns

import networkx as nx
import numpy as np


import EnrichRLib as erl


## Functions ========
## ==================

def sigmoid(x):                                        
   return 1 / (1 + np.exp(-x))


def cluster_genset(enr,top_terms=5,top_keywords=5):
    #enr_grp = enr.groupby(by=['cluster','cluster_jacc'],axis=0)
    enr_grp = enr.groupby(by=['cluster'],axis=0)

    ng = len(enr_grp.groups.keys())
    cols = ['cloud_name', 'cluster','genes', 'top_terms', 'terms']
    df = pd.DataFrame(index=range(ng), columns=cols)

    for ic,k in enumerate(enr_grp.groups.keys()):
        lst_0 =  enr_grp.groups[k]
        clust = enr.loc[lst_0]
        clust.sort_values('p-Val', axis=0, inplace=True)
        clust_cut = clust.head(top_terms)
        
        gn_clust = set(clust_cut.iloc[0]['genes'])
        for i in range(len(clust_cut)):
            gn_clust = gn_clust & set(clust_cut.iloc[i]['genes'])

        lst = list(clust_cut.index)
        
        lsto = [x.split('_')[0].split('(')[0] for x in lst_0 ]
        fc = pd.DataFrame.from_dict(erl.word_count(lsto), orient='index', columns=['counts'])
        fc.sort_values('counts', axis=0, inplace=True, ascending=False)
        fc.index = [ind+' ('+str(fc.loc[ind]['counts'])+')' for ind in fc.index]
        fc_name = u'\u2295 '+u'\n\u2295 '.join(list(fc.index[:top_keywords]))
        
        df.loc[ic]['cloud_name'] = fc_name
        df.loc[ic]['top_terms'] = list(lst)
        df.loc[ic]['terms'] = list(lst_0)
        df.loc[ic]['genes'] = list(gn_clust)
        df.loc[ic]['cluster'] = k
        

    gs_cl = df[['cloud_name','genes']].set_index('cloud_name').T.to_dict('list')
    gs_clust = {k:v[0] for k,v in gs_cl.items()}
    df.set_index('cloud_name', inplace=True)
    
    return gs_clust, df


def categorical_cmap(nc, nsc, cmap="tab10", continuous=False):
    if nc > plt.get_cmap(cmap).N:
        raise ValueError("Too many categories for colormap.")
    if continuous:
        ccolors = plt.get_cmap(cmap)(np.linspace(0,1,nc))
    else:
        ccolors = plt.get_cmap(cmap)(np.arange(nc, dtype=int))
    cols = np.zeros((nc*nsc, 3))
    for i, c in enumerate(ccolors):
        chsv = matplotlib.colors.rgb_to_hsv(c[:3])
        arhsv = np.tile(chsv,nsc).reshape(nsc,3)
        arhsv[:,1] = np.linspace(chsv[1],0.25,nsc)
        arhsv[:,2] = np.linspace(chsv[2],1,nsc)
        rgb = matplotlib.colors.hsv_to_rgb(arhsv)
        cols[i*nsc:(i+1)*nsc,:] = rgb       
    cmap = matplotlib.colors.ListedColormap(cols)
    return cmap


def draw_graph(G, palette='gnuplot'):
    cl = [int(G.nodes[v]['cluster']) for v in G] #cl = range(len(G))
    nc = len(set(cl))
    #print(cl) 
    sz = [G.nodes[v]['mlog10pVal']*100 for v in G]
    ec = [G.edges[u, v]['weight'] for u,v in G.edges]
 
    #~ kappa = min(ec)
    
    f, ax = plt.subplots(figsize=(10, 8))
    pos=nx.spring_layout(G)#,k=0.1, weight = 'spring')     
    
    if nc >20:
        cmpl  = categorical_cmap(nc, 1, cmap=palette, continuous=True)
    elif nc>10:
        cmpl  = categorical_cmap(nc, 1, cmap="tab20")
    else:
        cmpl  = categorical_cmap(nc, 1, cmap="tab10")
    
    nx.draw(G, pos=pos,
            node_color = cl,
            node_size =  sz,
            edge_color =  ec,
            #~ edge_vmin = 0.6*kappa,
            #~ edge_vmax = 1.2+kappa,
            edge_cmap =plt.cm.Greys,
            cmap = cmpl) 
    
    nx.draw_networkx_labels(G, pos=pos,
                            font_size=12)
    plt.tight_layout()
    #~ return cmpl




def make_graph(setR, enr_in, kappa=0.001):
    setR = set(setR)
    terms = list(enr_in.index)

    G = nx.Graph()
    
    for tr in terms:
        G.add_node( tr,
                    genes = enr_in.loc[tr]['genes'],
                    num_list = enr_in.loc[tr]['num_list'],
                    num_term = enr_in.loc[tr]['num_term'],
                    ratio = enr_in.loc[tr]['num_list']/enr_in.loc[tr]['num_term'],
                    pVal = enr_in.loc[tr]['p-Val'],
                    padj = enr_in.loc[tr]['p-adj'],
                    mlog10pVal = enr_in.loc[tr]['-log10(p-Val)'],
                    cluster =  enr_in.loc[tr]['cluster'])

    for tri in terms:
        seti = set(list(enr_in['term_genes'].loc[tri]))
        for tro in terms:
            ed = tri,tro
            if not G.has_edge(*ed):
                seto = set(list(enr_in['term_genes'].loc[tro]))
                dist = erl.coeff_kappa(setR,seti,seto)
                #print(dist)
                if dist>kappa:
                    G.add_edge(*ed, weight=dist)
    

    return G




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
enr = enr[enr['p-Val']<0.007]

#3= Make claster by kappa coeff
enr = erl.cluster(set(gl), enr, deep=3)

#3-1= Filter top clusters
top_clusters = 20
enr = enr[enr['cluster']<top_clusters]

#4= Make clustered geneset 
gs_clust,nt_cl = cluster_genset(enr, top_terms=3)

#5= Enrich clustered geneset
enr_clust = erl.enrich(gl,gs_clust)

#6= Add cluster to table
#enr_clust = pd.concat([enr_clust,nt_cl[['cluster']]],axis=1, sort=False)
enr_clust = pd.concat([enr_clust,nt_cl.loc[enr_clust.index]['cluster']],axis=1, sort=False)

#7= Make graphs
G_gs = erl.make_graph_n(gl,enr, kappa=0.4)
G_cl = erl.make_graph_n(gl,enr_clust, kappa=0.01)

#8= Draw graphs
erl.draw_graph(G_gs, spring=150)
draw_graph(G_cl)



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




plt.show()




