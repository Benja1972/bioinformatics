import numpy as np
from os.path import join 
#import sys
import pandas as pd
import pybedtools as pb
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx

# Project settings
from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')





## FUNCTIONS
###==========================

def reg_network(gl,rn_df):
    
    rn = rn_df.copy()
    mr = list(rn[rn['regulator_symbol'].str.contains('mmu-')]['regulator_symbol'].unique())
    tf = list(rn[~rn['regulator_symbol'].str.contains('mmu-')]['regulator_symbol'].unique())

    gn = list(rn['target_symbol'].unique())
    gn = list(set(gn) - set(tf) -set(mr))


    rn_gl = rn[rn['target_symbol'].isin(gl)]


    terms = list(set(list(rn_gl['target_symbol'].unique())) | set(list(rn_gl['regulator_symbol'].unique())))



    G = nx.DiGraph()

    for tr in terms:
        if tr in tf:
            typ, cl = 'tf','c'
        elif tr in mr:
            typ, cl = 'mr','m'
        else:
            typ, cl = 'gn','g'
            
        G.add_node( tr,
                    typ = typ,
                    color=cl)


    for nm in rn_gl.index:
        G.add_edge( rn_gl.loc[nm]['regulator_symbol'],
                    rn_gl.loc[nm]['target_symbol'])   
    
    return G


def draw_reg_network(G):
    cl = [G.nodes[v]['color'] for v in G]

    #edges TF end miRNA
    ed_tf = [(u,v) for u,v in G.edges if G.nodes[u]['typ']=='tf']
    ed_mr = [(u,v) for u,v in G.edges if G.nodes[u]['typ']=='mr']



    f, ax = plt.subplots(figsize=(8, 8))
    pos=nx.spring_layout(G,k=0.4)


    nx.draw_networkx_nodes(G, pos=pos,
            node_color = cl,
            node_size = 400,
            alpha=0.85)

    nx.draw_networkx_edges(G, pos=pos, 
                            edgelist = ed_tf,
                            edge_color='#ED453F',#'r',
                            alpha = 0.95,
                            arrowstyle = '-|>',
                            arrowsize=12)

    nx.draw_networkx_edges(G, pos=pos, 
                            edgelist = ed_mr,
                            edge_color='#2588F3',# 'b',
                            alpha = 0.95,
                            arrowstyle = '-[',
                            arrowsize=12)




    nx.draw_networkx_labels(G, pos=pos,
                            font_size=12,
                            #font_color='grey',
                            alpha=0.7)



##  == RegNetworks analysis ==
### ==========================

gl  = ['Gucy1b2', 'Il5', 'Kras', 'Malt1']

rn_df = pd.read_csv(join(DATADIR,'RegNetworkDB/RegNetworkDB_Jun_2019.csv'))

G = reg_network(gl, rn_df)

draw_reg_network(G)

plt.show()

