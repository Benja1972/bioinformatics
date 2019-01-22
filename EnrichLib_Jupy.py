from os.path import join
import numpy as np
import pandas as pd
import EnrichRLib as erl
import seaborn as sns
import matplotlib.pyplot as plt



exp_name='GO_KEGG2'



#~ gt = pd.read_table('tracks/MARGE/relativeRP/DN_RP_enhancers_genes.txt')


gl = [  'CD3E', 'BLK', 'PTPN22', 'PAG1', 'CTLA4', 'PIK3CD', 'LAT2', 'CSK', 
        'CD247', 'CD3G', 'THEMIS', 'PSMB8', 'LCP2', 'GATA3', 'LAT', 'SLA2', 
        'SKAP1', 'TRAT1', 'BCL2', 'CD3D', 'THY1','RUNX1', 'BLK', 'PTPN22']

#setR = set(gl)
gs_fn ='EnrichrLibs/Reactome_2016.gmt'
gs = erl.read_gmt(gs_fn)
enr = erl.enrich(set(gl), gs)


#get_Enrichr(out_dir = 'EnrichrLibs')

# Several gene sets in one
gss = [ 
       'GO_Biological_Process_2018',
       'GO_Cellular_Component_2018',
       'GO_Molecular_Function_2018',
       'KEGG_2016',
       'Reactome_2016'
       ]


enrr = erl.enrich_gs(gl,gss, path_lib='EnrichrLibs')

## p-Value filter
enrr = enrr[enrr['p-Val']<0.001]

## Cluster
enrr = erl.cluster(gl,enrr)
enrr.sort_values('cluster', axis=0, inplace = True)
enrr.loc[:,'ass_genes_percnt'] = 100*enrr.loc[:,'num_list']/enrr.loc[:,'num_term']
enrr.sort_values('cluster', axis=0, inplace = True)

#cm = sns.color_palette("tab20", n_colors=enrr['cluster'].max())
cm = 'tab20'

ds = enrr.head(60)

f, ax = plt.subplots(figsize=(16.5, 6.5))
sns.barplot(y=ds.index,
            x='-log10(p-Val)',
            ax = ax, 
            hue ='cluster',#color="Red", 
            dodge=False,
            data = ds,
            palette = cm)
ax.set_title('All terms')








## Show cluster
dd = erl.dist_matrx(gl,enrr)
grid = sns.clustermap(dd, cmap='Blues', figsize=(16, 16))


## Graph and network
nt, nt_tb, G = erl.make_graph(gl, enrr, draw=True, palette=cm)

#~ cl = [G.nodes[v]['cluster'] for v in G] 
#~ sz = [G.nodes[v]['mlog10pVal']*100 for v in G] 

#~ f1, ax1 = plt.subplots(figsize=(10, 10))
     
#~ nx.draw(G, pos=nx.spring_layout(G),
        #~ with_labels = True, 
        #~ node_color = cl,
        #~ node_size =  sz,
        #~ font_size=8,
        #~ cmap = cm) 


plt.show()

