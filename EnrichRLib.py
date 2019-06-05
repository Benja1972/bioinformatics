from time import sleep
import requests
from os.path import join 
import os
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import scipy.stats as st
from scipy.cluster.hierarchy import fcluster, linkage
import networkx as nx

### == Math functions == 
def log2p1(x):
    return np.log2(x + 1)

def nlog10(x):
    return -np.log10(x)
    
# == Write gene set in gmt format
def write_gmt(st, name, fn):
    gmt = [name, name] + list(st)
    with open(fn, 'w') as fp:
        fp.write("\t".join(gmt))

# == Write gene set dic to gmt format
def write_dic2gmt(dic, fn):
    with open(fn, 'w') as fp:
        for db, term in dic.items():
            gmt = [db, db] + list(term)
            fp.write("\t".join(gmt))
            fp.write("\n")


# == Load gene set in gmt format
def read_gmt(gmt):
    with open(gmt) as f:
         gene_dict = { line.strip().split("\t")[0]: [ln.split(",")[0] for ln in line.strip().split("\t")[2:]]
                          for line in f.readlines()}
    return gene_dict

def coeff_kappa(setR, set1,set2):
    a = len(setR & set1 & set2)
    b = len((setR & set1) - set2)
    c = len((setR & set2) - set1)
    d = len((setR - set1) - set2)
    t = a + b + c + d
    p0 = np.float32(a + d)/t
    pY = np.float32((a + b)*(a + c))/(t*t)
    pN = np.float32((c + d)*(b + d))/(t*t)
    pe = pY + pN
    k = np.float32((p0 - pe))/np.float32(1 - pe)

    return np.float32(k)

def coeff_jaccard(set1,set2):
    a = len(set1 & set2)
    b = len(set1 | set2)
    jac = np.float32(a)/np.float32(b)

    return jac


def pValue_sets(setA, setB, tot=20000):
    """
    tot:    total number of genes 
    num_a:  total number of genes in the list with condition A
    num_b:  total number of genes in the list with condition B
    num_ab: number of genes with both condition A and B
    """
    num_a=len(setA)
    num_b=len(setB)
    num_ab=len(setA & setB)
    
    return st.hypergeom.sf(int(num_ab) - 1,int(tot),int(num_a),int(num_b))

def p_adjust(pvalue, method="fdr"):
    p = pvalue
    n = len(p)
    p0 = np.copy(p, order='K')
    nna = np.isnan(p)
    p = p[~nna]
    lp = len(p)
    if method == "bonferroni":
        p0[~nna] = np.fmin(1, lp * p)
    elif method == "fdr":
        i = np.arange(lp, 0, -1)
        o = (np.argsort(p))[::-1]
        ro = np.argsort(o)
        p0[~nna] = np.fmin(1, np.minimum.accumulate((p[o]/i*lp)))[ro]
    else:
        print("Method is not implemented")
        p0 = None
    return p0

# === Fetch all Enrichr libs
def get_Enrichr(out_dir = 'EnrichrLibs', libs='all'):
    # Create target directory & all intermediate directories if don't exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print("Directory " , out_dir ,  " created ")
    else:    
        print("Directory " , out_dir ,  " already exists")    

    if libs=='all':
        lib_url='http://amp.pharm.mssm.edu/Enrichr/datasetStatistics'
        libs_json = json.loads(requests.get(lib_url).text)
        gss = [lib['libraryName'] for lib in libs_json['statistics']]
    else:
        gss = libs

    link = 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName='

    for gs in gss:
        lk = link+gs
        fn = join(out_dir,gs+'.gmt')
        res = requests.get(lk, timeout=None)
        with open(fn, mode='w') as f:
            f.write(res.text)
        sleep(2)
        print(gs, " -- DONE")


# === Enrich single gene set
def enrich(setR,gs):
    setR = set(setR)
    enr = pd.DataFrame(index=gs.keys(), 
                        columns=['p-Val', 'num_list', 'num_term', 'genes', 'term_genes']) 

    for term, g in gs.items():
        enr.loc[term]['term_genes']  = list(g)
        set1 = set(g)
        enr.loc[term]['p-Val'] = pValue_sets(setR, set1)
        enr.loc[term]['num_list'] = int(len(setR & set1))
        enr.loc[term]['num_term'] = int(len(set1))
        enr.loc[term]['genes']  = list(setR & set1)

    enr = enr.astype({'p-Val':'float', 'num_list':'int', 'num_term':'int'})

    enr = enr[enr['num_list']>0 ]
    enr['p-adj'] = p_adjust(enr['p-Val'])

    enr['-log10(p-adj)'] =  - np.log10(enr.loc[:,'p-adj']).values
    enr['-log10(p-Val)'] =  - np.log10(enr.loc[:,'p-Val']).values
    enr.sort_values('p-Val', axis=0, inplace = True)
    
    return enr

# === Enrich many gene sets
def enrich_gs(setR,gss, path_lib='EnrichrLibs'):

    setR = set(setR)

    enrDf = pd.DataFrame()

    for gs in gss:
        pl = read_gmt(join(path_lib,gs)+'.gmt')
        enr = enrich(setR, pl)
        print(gs, len(enr))
        enrDf = pd.concat([enrDf, enr])
    enrDf = enrDf.drop(['p-adj', '-log10(p-adj)'], axis=1)
    enrDf['p-adj'] = p_adjust(enrDf['p-Val'])
    enrDf['-log10(p-adj)'] =  - np.log10(enrDf.loc[:,'p-adj']).values
    enrDf.sort_values('p-Val', axis=0, inplace = True)

    return enrDf


# ==== Distance matrix for clustering 
def kappa_matrx(setR,enr):
    setR = set(setR)
    terms = list(enr.index)

    # Similarity matrix
    dd = pd.DataFrame(index=terms, columns=terms)

    for ter1 in terms:
        set1 = set(list(enr['term_genes'].loc[ter1]))
        for ter2 in terms:
            set2 = set(list(enr['term_genes'].loc[ter2]))
            dd[ter1][ter2] = coeff_kappa(setR,set1,set2)

    dd = dd.astype(float)
    
    return dd


def jaccard_matrx(enr):
    terms = list(enr.index)

    # Similarity matrix
    dd = pd.DataFrame(index=terms, columns=terms)

    for ter1 in terms:
        set1 = set(list(enr['term_genes'].loc[ter1]))
        for ter2 in terms:
            set2 = set(list(enr['term_genes'].loc[ter2]))
            dd[ter1][ter2] = coeff_jaccard(set1,set2)

    dd = dd.astype(float)
    
    return dd



def cluster(setR, enr):
    enr = enr.copy()
    dd = kappa_matrx(setR,enr)

    ## Look to package fastclsuter !!!!!
    metric='euclidean'
    method='average'
    links = linkage(dd.values, method=method, metric=metric)
    cls = fcluster(links,2, 'distance')
    enr.loc[:,'cluster'] = cls
 
    return enr

def cluster_jacc(enr):
    enr = enr.copy()
    dd = jaccard_matrx(enr)

    ## Look to package fastclsuter !!!!!
    metric='euclidean'
    method='average'
    links = linkage(dd.values, method=method, metric=metric)
    cls = fcluster(links,2, 'distance')
    enr.loc[:,'cluster_jacc'] = cls
 
    return enr


def sigmoid(x):                                        
   return 1 / (1 + np.exp(-x))

def make_graph(setR, enr_in, kappa=0.4):
    setR = set(setR)

    ## Network
    terms = list(enr_in.index)
    n = len(terms)
    col = ['source', 'target', 'dist', 'comm_genes']
    nt = pd.DataFrame(columns=col)


    for i in range(len(terms)):
        ter1 = terms[i]
        set1 = set(list(enr_in['term_genes'].loc[ter1]))
        for j in range(i):
            ter2 = terms[j]
            set2 = set(list(enr_in['term_genes'].loc[ter2]))
            cm = list(setR & set1 & set2)
            rw = pd.DataFrame([[ter1, ter2, coeff_kappa(setR,set1,set2), cm]], columns=col)
            nt = nt.append(rw,  ignore_index=True)


    nt =  nt[nt['dist']>kappa]

    ls = list(set(list(nt['source']))|set(list(nt['target'])))
    enr_out = enr_in[enr_in.index.isin(ls)]
    
    ## Cluster
    enr_out = cluster(setR,enr_out)
    enr_out.sort_values('cluster', axis=0, inplace = True)
    
    ## Graph
    G = nx.Graph()
    
    for tr in enr_out.index:
        G.add_node( tr,
                    genes = enr_out.loc[tr]['genes'],
                    num_list = enr_out.loc[tr]['num_list'],
                    num_term = enr_out.loc[tr]['num_term'],
                    ratio = enr_out.loc[tr]['num_list']/enr_out.loc[tr]['num_term'],
                    pVal = enr_out.loc[tr]['p-Val'],
                    padj = enr_out.loc[tr]['p-adj'],
                    mlog10pVal = enr_out.loc[tr]['-log10(p-Val)'],
                    cluster =  enr_out.loc[tr]['cluster'])

    for nm in nt.index:
        G.add_edge( nt.loc[nm]['source'],
                    nt.loc[nm]['target'],
                    weight = nt.loc[nm]['dist'])       

    return G, enr_out , nt 



def draw_graph(G,spring = 1.0, pval_prcnt=0.8, palette='tab20'):
    cl = [float(G.nodes[v]['cluster']) for v in G] 
    sz = [G.nodes[v]['mlog10pVal']*100 for v in G]
    ec = [G.edges[u, v]['weight'] for u,v in G.edges]
    for u,v in G.edges:
        G.edges[u, v]['spring'] = sigmoid(G.edges[u, v]['weight'])*(spring/500)
    
    kappa = min(ec)
    mpValshr=pval_prcnt*max(sz)/100
    lb_M = {v:v for v in G if G.nodes[v]['mlog10pVal']>mpValshr}
    lb_m = {v:v for v in G if G.nodes[v]['mlog10pVal']<=mpValshr} 
    
    f, ax = plt.subplots(figsize=(15, 15))
    pos=nx.spring_layout(G,k=0.1, weight = 'spring')     
    
    
    #~ cmap = get_cmap(name='viridis')
    #~ colors = cmap(np.linspace(0, 1, max(cl)))
    nx.draw(G, pos=pos,
            node_color = cl,
            vmin = 1.,
            vmax = 20., ##FIX hardcoded number for tab20
            node_size =  sz,
            edge_color =  ec,
            edge_vmin = 0.6*kappa,
            edge_vmax = 1.2+kappa,
            edge_cmap =plt.cm.Greys,
            cmap = palette) 
    
    nx.draw_networkx_labels(G, pos=pos,
                            labels = lb_M,
                            font_size=14,
                            font_weight='bold')

    nx.draw_networkx_labels(G, pos=pos,
                            labels = lb_m,
                            font_size=6)
