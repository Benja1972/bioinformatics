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

def reg_network(gl,rn_df, type='target'):
    """
    Function construct regulatory network
    
    Parametes
    ---------
    rn_df : DataFrame based on csv file from 'http://www.regnetworkweb.org', ex rn_df = pd.read_csv('RegNetwork.csv')
    gl : gene list of interest
    
    Returns
    ------
    G :  networkx graph object
    """
    
    rn = rn_df.copy()
    mr = list(rn[rn['regulator_symbol'].str.contains('mmu-')]['regulator_symbol'].unique())
    tf = list(rn[~rn['regulator_symbol'].str.contains('mmu-')]['regulator_symbol'].unique())

    gn = list(rn['target_symbol'].unique())
    gn = list(set(gn) - set(tf) -set(mr))

    if type=='target':
        rn_gl = rn[rn['target_symbol'].isin(gl)]
    elif type=='regulator':
        rn_gl = rn[rn['regulator_symbol'].isin(gl)]


    terms = list(set(list(rn_gl['target_symbol'].unique())) | set(list(rn_gl['regulator_symbol'].unique())))



    G = nx.DiGraph()

    for tr in terms:
        if tr in tf:
            typ, cl = 'tf','c'
        elif tr in mr:
            typ, cl = 'mr','m'
        else:
            typ, cl = 'gn','g'
        
        if tr in gl:
            fl=1
        else:
            fl=0
            
            
        G.add_node( tr,
                    typ = typ,
                    color=cl,
                    from_list=fl)


    for nm in rn_gl.index:
        G.add_edge( rn_gl.loc[nm]['regulator_symbol'],
                    rn_gl.loc[nm]['target_symbol'])   
    
    return G


def draw_reg_network(G, spring=0.35):
    """
    Function draw regulatory network
    
    Parametes
    ---------
    G :  networkx graph object
    spring : spring  koefficient, smaller value -> tighter layout
    
    """
    cl =  [G.nodes[v]['color'] for v in G]
    edc = [G.nodes[v]['from_list'] for v in G]
    cm = {1:'#333333',0:'w'}
    lw = {1: 2,0:0.2}
    edm = [cm[i] for i in edc]
    elw = [lw[i] for i in edc]

    #edges TF end miRNA
    ed_tf = [(u,v) for u,v in G.edges if G.nodes[u]['typ']=='tf']
    ed_mr = [(u,v) for u,v in G.edges if G.nodes[u]['typ']=='mr']


    
    f, ax = plt.subplots(figsize=(8, 8))
    pos=nx.spring_layout(G,k=spring)


    nx.draw_networkx_nodes(G, pos=pos,
            node_color = cl,
            node_size = 400,
            edgecolors = edm,
            linewidths = elw,
            alpha=0.9)

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
                            font_size=10,
                            font_color='#414141',
                            alpha=1)



##  == RegNetworks analysis ==
### ==========================

gl  = ['Gucy1b2', 'Il5', 'Malt1', 'Tlx3', 'Foxn2'] # , 'Kras'

rn_df = pd.read_csv(join(DATADIR,'RegNetworkDB/RegNetworkDB_Jun_2019.csv'))

Gr = reg_network(gl, rn_df, type='regulator') #
Gt = reg_network(gl, rn_df, type='target') #

draw_reg_network(Gr)
draw_reg_network(Gt, spring=0.4)

plt.show()

