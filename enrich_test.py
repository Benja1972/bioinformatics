from os.path import join
import numpy as np
import pandas as pd
#import RNA_expression_processing as rn
import EnrichRLib as erl
import seaborn as sns
import matplotlib.pyplot as plt


SAVE = False

exp_name='GO_KEGG'


gt = pd.read_table('tracks/MARGE/relativeRP/DN_RP_enhancers_genes.txt')


setR = set(list(gt.iloc[:,0]))







gss = [ 'Pierre_gene_sets.gmt',
        #~ 'GO_Biological_Process_2018.gmt',
       #~ 'GO_Cellular_Component_2018.gmt',
       #~ 'GO_Molecular_Function_2018.gmt',
       'KEGG_2016.gmt',
       #~ 'Reactome_2016.gmt'
       ]



enrr = pd.DataFrame()



for gs in gss:
    pl = erl.read_gmt('EnrichrLibs/'+gs)
    enr = erl.enrich(setR, pl)
    print(len(enr))
    enrr = pd.concat([enrr, enr])
    ds = enr.head(20)
    f, ax = plt.subplots(figsize=(16.5, 6.5))
    sns.barplot(y=ds.index,
                x='-log10(p-Val)',
                ax = ax, 
                color="Red", 
                data = ds)
    ax.set_title(gs.split(".")[0])




# p-Value filter
enrfl = enrr[enrr['p-Val']<0.05]
#enrfl = enrr.copy()

dd = erl.dist_matrx(setR,enrfl)

grid = sns.clustermap(dd, cmap='Blues', figsize=(16, 16))

from scipy.cluster.hierarchy import fcluster
ass = fcluster(grid.dendrogram_col.linkage,2, 'distance')
enrfl['cluster'] = ass


plt.show()

if SAVE:
    enrfl.to_csv(exp_name+'-enrich.txt', sep="\t")



## Network
terms = list(enrfl.index)
n = len(terms)
col = ['source', 'target', 'dist', 'comm_genes']
nt = pd.DataFrame(columns=col)


for i in range(len(terms)):
    ter1 = terms[i]
    set1 = set(list(enrfl['genes'].loc[ter1]))
    for j in range(i):
        ter2 = terms[j]
        set2 = set(list(enrfl['genes'].loc[ter2]))
        cm = list(setR & set1 & set2)
        rw = pd.DataFrame([[ter1, ter2, erl.dist_kappa(setR,set1,set2), cm]], columns=col)
        nt = nt.append(rw,  ignore_index=True)


nt =  nt[nt['dist']>0.4]

ls = list(set(list(nt['source']))|set(list(nt['target'])))
nt_tb = enrfl[enrfl.index.isin(ls)]

if SAVE:
    nt.to_csv(exp_name+'-network.txt', sep="\t")
    nt_tb.to_csv(exp_name+'-network_table.txt', sep="\t")


## Graph 
import networkx as nx

G = nx.Graph()
nt_tb = nt_tb.sort_values('cluster', axis=0)

for tr in nt_tb.index:
    G.add_node( tr,
                mlog10pVal = nt_tb.loc[tr]['-log10(p-Val)'],
                cluster =  nt_tb.loc[tr]['cluster'])

for nm in nt.index:
    G.add_edge( nt.loc[nm]['source'],
                nt.loc[nm]['target'],
                weight = nt.loc[nm]['dist'])       


cl = nx.get_node_attributes(G, 'cluster')
sz = nx.get_node_attributes(G, 'mlog10pVal')

f, ax = plt.subplots()
     
nx.draw(G, pos=nx.spring_layout(G),
        with_labels = True, 
        node_color = list(cl.values()),
        node_size =  np.array(list((sz.values())))*100,
        ) #cmap = plt.cm.RdBu #font_weight='bold'

plt.show()
