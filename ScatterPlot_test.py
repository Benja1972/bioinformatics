import time
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from os.path import join 
import numpy as np
import scipy.stats as st
from scipy.stats.distributions import norm

def log2p1(x):
    return np.log2(x + 1)


#~ def p_value2(df, A, B, classes):
    #~ df = df.copy()
    #~ df_gr = df.groupby(by=classes, axis=1)

    #~ df['p-val']=np.nan

    #~ for i in range(len(df)):
        #~ df['p-val'][i] = st.ttest_ind(df[df_gr.groups[A]].iloc[i], df[df_gr.groups[B]].iloc[i])[1]
    #~ df['padj'] = df['p-val']*len(df)/len(classes)
    #~ return df



def p_value(df, A, B, classes):
    df = df.copy()
    la = classes.count(A)
    lb = classes.count(B)
    df_mean = df.groupby(by=classes, axis=1).mean()
    df_var = df.groupby(by=classes, axis=1).var(ddof=1)
    t = (df_mean[A] - df_mean[B])/ np.sqrt(df_var[A]/la+df_var[B]/lb)
    dgf = la+lb - 2
    pval =  2*st.t.cdf(-abs(t),df=dgf)
    return pval



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


def scatter(df, A, B, classes, n_top=0, geneList=[],  ttl=''):
    df = df.copy()
    df_mean= df.groupby(by=classes, axis=1).mean()


    df = pd.concat([df,df_mean], axis=1)
    df['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))

    if not geneList:
        geneList = list(df.index)
    
    gs = df.iloc[df.index.isin(geneList)]

    gs = gs.sort_values('lgFC', 
                        axis=0, 
                        ascending=False)
    
    up = gs[gs.lgFC > 1.]
    dn = gs[gs.lgFC < -1.]
    
    up_l = len(up)
    dn_l = len(dn)
    tt_l = len(gs)
    
    gs['lg'+A] = log2p1(gs[A])
    gs['lg'+B] = log2p1(gs[B])


    diff = "diff_"+A+"_vs_"+B

    gs[diff]='unchng'
    gs.loc[gs["lgFC"]>1,diff]='up'
    gs.loc[gs["lgFC"]<-1,diff]='down'

    cpp = [(0.86, 0.23, 0.22), (0.5,0.5,0.5),(0.03, 0.45, 0.56)]
    f, ax = plt.subplots(figsize=(6.5, 6.5))
    sns.scatterplot(x = 'lg'+B, 
                    y = 'lg'+A, 
                    hue=diff, 
                    data=gs, 
                    ax=ax, 
                    palette=cpp, 
                    linewidth=0, 
                    s=3.4)
    ax.set_xlabel(B +', log2(FPKM + 1)')
    ax.set_ylabel(A +', log2(FPKM + 1)')
    ax.set_title(ttl)
    ax.axis('equal')
    ax.axis('tight')
    
    # ==== Plot text of top genes
    for i in range(n_top):
        rx = .0003*np.random.randn()
        ry = .0003*np.random.randn()
        #~ print gss.index[i]
        ax.text(log2p1(gs[B][i]*(1.+rx)),
                log2p1(gs[A][i]*(1.+ry)),
                gs.index[i], 
                color=cpp[0])
        ax.text(log2p1(gs[B][-i-1]*(1.+rx)),
                log2p1(gs[A][-i-1]*(1.+ry)),
                gs.index[-i-1], 
                color=cpp[2])
                
    ax.text(0.37,0.95, 'UP: '+str(up_l),
            color=cpp[0],
            transform=ax.transAxes)
    ax.text(0.55,0.95, 'DOWN: '+str(dn_l),
            color=cpp[2],
            transform=ax.transAxes)
    ax.text(0.82,0.95, 'ALL: '+str(tt_l),
            transform=ax.transAxes)
    
    return up, dn, ax


def volcano(df, A, B, classes, n_top=0, geneList=[],  ttl=''):
    df = df.copy()
    df_mean= df.groupby(by=classes, axis=1).mean()
    df['p-val'] = p_value(df, A, B, classes)
    df['p-adj'] = p_adjust(df['p-val'])
    df['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))

    pv = 'p-adj'#'p-val'
    df['nlog10_pV'] = -np.log10(df[pv])
    
    if not geneList:
        geneList = list(df.index)
    
    gs = df.iloc[df.index.isin(geneList)]
    gs = gs.sort_values('lgFC', 
                        axis=0, 
                        ascending=False)
    
    up = gs[gs.lgFC > 1.]
    dn = gs[gs.lgFC < -1.]
    
    up_l = len(up)
    dn_l = len(dn)
    tt_l = len(gs)
    

    diff = "diff_"+A+"_vs_"+B

    gs[diff]='unchng'
    gs.loc[gs["lgFC"]>1,diff]='up'
    gs.loc[gs["lgFC"]<-1,diff]='down'

    cpp = [(0.86, 0.23, 0.22), (0.5,0.5,0.5),(0.03, 0.45, 0.56)]
    f, ax = plt.subplots(figsize=(6.5, 6.5))
    sns.scatterplot(x = 'lgFC', 
                    y = 'nlog10_pV', 
                    hue=diff, 
                    data=gs, 
                    ax=ax, 
                    palette=cpp, 
                    linewidth=0, 
                    s=3.4)
    ax.set_xlabel('log FC')
    ax.set_ylabel('-log10($p_{val}$)')
    ax.set_title(ttl)
    ax.axis('equal')
    ax.axis('tight')
    
    # ==== Plot text of top genes
    for i in range(n_top):
        rx = .0003*np.random.randn()
        ax.text(gs['lgFC'][i]*(1.+rx),
                gs['nlog10_pV'][i]+1e-300,
                gs.index[i], 
                color=cpp[0])
        ax.text(gs['lgFC'][-i-1]*(1.+rx),
                gs['nlog10_pV'][-i-1]+1e-300,
                gs.index[-i-1], 
                color=cpp[2])
                
    ax.text(0.037,0.95, 'UP: '+str(up_l),
            color=cpp[0],
            transform=ax.transAxes)
    ax.text(0.25,0.95, 'DOWN: '+str(dn_l),
            color=cpp[2],
            transform=ax.transAxes)
    ax.text(0.52,0.95, 'ALL: '+str(tt_l),
            transform=ax.transAxes)
    
    return up, dn, ax


def cluster_express(df, A, B, classes, n_top=0):
    df = df.copy()
    cols = df.columns
    df_mean= df.groupby(by=classes, axis=1).mean()

    df['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))
    
    gs = df[abs(df['lgFC'])>1].sort_values('lgFC', axis=0, ascending=False)
    
    if n_top != 0:
        n_top = int(min(n_top,len(gs)/2))
    
    gs = pd.concat([gs.head(n_top),gs.tail(n_top)], axis=0)
    
    dpic  = log2p1(gs[cols])
    
    # --Cluster
    grid = sns.clustermap(dpic, cmap='RdBu_r', 
                        linewidths=0.004,
                        #metric='correlation',
                        method='centroid',
                        col_cluster=False,
                        figsize=(8, int(n_top/2)))

    #grid.ax_heatmap.set_title(ttl+ ' diffExpressed' )
    #return gs
    # -----------
    






A = 'TLX3'
B = 'RAG'


tbl = pd.read_table(join('tracks', 'TLX3vsRAG-results_genesNames.txt'), index_col=0)

cols = ['Gene_name', 'TLX3.1_1','TLX3.1_5','TLX3.1_P','R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ']

df = tbl[cols]
df.set_index(keys=df.columns[0], inplace=True)


classes = [A]*3+[B]*3

dfs = df.copy()#sample(1000)
gs = cluster_express(dfs, A, B, classes, n_top=50)
#~ scatter(dfs, A, B, classes, n_top=2)
#~ volcano(dfs, A, B, classes, n_top=4)
plt.show()



#~ dfp = p_value(dfs, A, B, classes)   

t1 = time.time()

#~ dfs['p-val'] = p_value(dfs,A,B,classes)
#~ dfs['p-adj'] = p_adjust(dfs['p-val'])

t2 = time.time()


#print('t1= ', t1-t0)
print('t2= ', t2-t1)














