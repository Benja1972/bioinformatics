from os.path import join 

import numpy as np
import sys, logging
import multiprocessing
import itertools
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as st

# comment metaseq for Python3 (works only express)
if sys.version_info[0] < 3:
    from metaseq.results_table import ResultsTable

# RNA expression examples from gseapy package
# classes shoud be like this
# classes = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']
# expression DataFrame should have first column (e.g. 'NAME') for Gene names, 
# with one column per sample in class

### == Math functions == 
def log2p1(x):
    return np.log2(x + 1)

def nlog10(x):
    return -np.log10(x)
    
# == Write gene set in gmt format
def write_gmt(st, name, path=''):
    gmt = [name, name] + list(st)
    with open(join(path,name+'.gmt'), 'w') as fp:
        fp.write("\t".join(gmt))

# == Write gene set dic to gmt format
def write_dic2gmt(dic, name, path=''):
    with open(join(path,name+'.gmt'), 'w') as fp:
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


## == Proper read .narrowPeak files
def read_narrowPeak(name, path=''):
    pk = pd.read_table(join(path,name),sep='\t', names=['Chr', 
        'start', 'end', 'name', 'score', 'strand','FoldChange','-log10pvalue',
        '-log10qvalue', 'smmt_rel_to_start'])
    return pk


## == Write selected back to simple .bed file
def write_df2bed(df,name, path='' ):
    df[['Chr','start','end']].to_csv(join(path,name), 
                                        sep='\t', 
                                        index=False,
                                        header=False)


def dist_kappa(setR, set1,set2):
    a = len(setR & set1 & set2)
    b = len((setR & set1) - set2)
    c = len((setR & set2) - set1)
    d = len((setR - set1) - set2)
    t = a + b + c + d
    p0 = np.float64(a + d)/t
    pY = np.float64((a + b)*(a + c))/(t*t)
    pN = np.float64((c + d)*(b + d))/(t*t)
    pe = pY + pN
    k = np.float64((p0 - pe))/np.float64(1 - pe)

    return np.float64(k)

## == Enrichment analysis P-values
def pValue(tot, num_a, num_b, num_ab):
    """
    tot:    total number of genes 
    num_a:  total number of genes in the list with condition A
    num_b:  total number of genes in the list with condition B
    num_ab: number of genes with both condition A and B
    
    p_val_minus = st.hypergeom.cdf(int(num_ab),int(tot),int(num_a),int(num_b))
    p_val_plus  = st.hypergeom.sf(int(num_ab) - 1,int(tot),int(num_a),int(num_b)

    a p-value where by random chance number of genes with both condition A and B will be <= to your number with condition A and B
    a p-value where by random chance number of genes with both condition A and B will be >= to your number with condition A and B
    The second p-value is probably what you want
    """
    
    return st.hypergeom.sf(int(num_ab) - 1,int(tot),int(num_a),int(num_b))


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






def barplot(df, cutoff=0.05, figsize=(6.5,6), top_term=10, ttl=""):
    """ barplot for enrichr results"""

    # pvalue cut off
    d = df[df['Adjusted P-value'] <= cutoff]

    if len(d) < 1:
        return None
    d = d.assign(logAP = - np.log10(d.loc[:,'Adjusted P-value']).values )
    d = d.sort_values('logAP', ascending=False)
    dd = d.head(top_term).sort_values('logAP', ascending=False)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    sns.barplot(y='Term', x='logAP', color="Red", ax=ax, data = dd)
    ax.set_xlabel("-log$_{10}$ Adjust P-value")
    #ax.set_ylabel("")
    ax.set_title(ttl)
    ax.legend(loc=4)
    return fig




def preprocess(df):
    """pre-processed the data frame.new filtering methods will be implement here.
    """    
    
    df.drop_duplicates(subset=df.columns[0], inplace=True) #drop duplicate gene_names.    
    df.set_index(keys=df.columns[0], inplace=True)
    df.dropna(how='all', inplace=True)                     #drop rows with all NAs
    df2 = df.select_dtypes(include=['float64'])  + 0.001 #select numbers in DataFrame      
    
    return df2

def ranking_metric(df, method, phenoPos, phenoNeg, classes, ascending):
    """The main function to rank an expression table.
    
   :param df:      gene_expression DataFrame.    
   :param method:  The method used to calculate a correlation or ranking. Default: 'log2_ratio_of_classes'.
                   Others methods are:
   
                   1. 'signal_to_noise' 
      
                      You must have at least three samples for each phenotype to use this metric.
                      The larger the signal-to-noise ratio, the larger the differences of the means (scaled by the standard deviations);
                      that is, the more distinct the gene expression is in each phenotype and the more the gene acts as a "class marker." 
    
                   2. 't_test'
      
                      Uses the difference of means scaled by the standard deviation and number of samples. 
                      Note: You must have at least three samples for each phenotype to use this metric.
                      The larger the tTest ratio, the more distinct the gene expression is in each phenotype 
                      and the more the gene acts as a "class marker."
    
                   3. 'ratio_of_classes' (also referred to as fold change).
      
                      Uses the ratio of class means to calculate fold change for natural scale data.
    
                   4. 'diff_of_classes' 
      
                      Uses the difference of class means to calculate fold change for log scale data
    
                   5. 'log2_ratio_of_classes' 
      
                      Uses the log2 ratio of class means to calculate fold change for natural scale data.
                      This is the recommended statistic for calculating fold change for natural scale data.
   
      
   :param phenoPos: one of lables of phenotype's names.
   :param phenoNeg: one of lable of phenotype's names.   
   :param classes:  a list of phenotype labels, to specify which column of dataframe belongs to what catogry of phenotype.
   :param ascending:  bool or list of bool. Sort ascending vs. descending.
   :return: returns correlation to class of each variable. same format with .rnk file. gene_name in first coloum,
            correlation in second column.  
            
    visit here for more docs: http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html
    """ 
        
    A = phenoPos
    B = phenoNeg
    df2 = df.T   
    df2['class'] = classes
    df_mean= df2.groupby('class').mean().T
    df_std = df2.groupby('class').std().T  
    #exclude any zero stds.
    df_mean = df_mean[df_std.sum(axis=1) !=0]
    df_std = df_std[df_std.sum(axis=1) !=0]
    
    if method == 'signal_to_noise':
        sr = (df_mean[A] - df_mean[B])/(df_std[A] + df_std[B])
    elif method == 't_test':
        sr = (df_mean[A] - df_mean[B])/ np.sqrt(df_std[A]**2/len(df_std)+df_std[B]**2/len(df_std) )
    elif method == 'ratio_of_classes':
        sr = df_mean[A] / df_mean[B]
    elif method == 'diff_of_classes':
        sr  = df_mean[A] - df_mean[B]
    elif method == 'log2_ratio_of_classes':
        sr  =  np.log2(df_mean[A] / df_mean[B])
    else:
        logging.error("Please provide correct method name!!!")        
        sys.exit()
    sr.sort_values(ascending=ascending, inplace=True)
    df3 = sr.to_frame().reset_index()
    df3.columns = ['gene_name','rank']
    df3['rank2'] = df3['rank']

    return df3


def updn(df, phenoPos, phenoNeg, classes, geneList=[]):
    # ---
    A = phenoPos
    B = phenoNeg
    
    df = df.dropna()
    df = df.iloc[:,:len(classes)+1]
    
    df.set_index(keys=df.columns[0], inplace=True)
    df.index=df.index.str.upper()
    df2 = df.T
    df2['class'] = classes
    df_mean= df2.groupby('class').mean().T
    df_std = df2.groupby('class').std().T

    dfm = pd.concat([df,df_mean], axis=1)
    dfm['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))
    
    if not geneList:
        geneList = list(dfm.index)
    
    gs = dfm.iloc[dfm.index.isin(geneList)]
    
    up = gs[gs.lgFC > 1.]
    dn = gs[gs.lgFC < -1.]
    
    return up, dn



def express(df, phenoPos, phenoNeg, classes, n_top=0, geneList=[],  ttl='', sort=True, diffr=True):
    # ---
    A = phenoPos
    B = phenoNeg
    
    df = df.dropna()
    df = df.iloc[:,:len(classes)+1]


    df.set_index(keys=df.columns[0], inplace=True)
    df.index=df.index.str.upper()
    df2 = df.T
    df2['class'] = classes
    df_mean= df2.groupby('class').mean().T
    df_std = df2.groupby('class').std().T

    dfm = pd.concat([df,df_mean], axis=1)
    dfm['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))
    
    if not geneList:
        geneList = list(dfm.index)
    
    gs = dfm.iloc[dfm.index.isin(geneList)]
    gs = gs.reindex(geneList)
    gs['FCmod'] = abs(gs['lgFC'])
    
    up = gs[gs.lgFC > 1.].sort_values('lgFC', axis=0, ascending=False)
    dn = gs[gs.lgFC < -1.].sort_values('lgFC', axis=0, ascending=False)
    #diff = gs[abs(gs.lgFC) > 1.]
    
    
    up_l = len(up)
    dn_l = len(dn)
    tt_l = len(gs)
    
    #print tt_l, len(geneList)
    if diffr:
        #top = diff.copy()
        top = pd.concat([up[:n_top],dn[-n_top:]], axis=0)
    else:
        top = gs.copy()
        if sort:
            top = top.sort_values('FCmod', axis=0, ascending=False)
            #top = top.sort_values('lgFC', axis=0, ascending=False)

        if n_top==0:
            n_top=len(gs)

        top = top[:n_top]
        #top = top.sort_values(A, axis=0, ascending=False)
        #top = top.sort_values('lgFC', axis=0, ascending=False)



        

    
    fig = plt.figure(figsize=(6, max(11,int(n_top/9.)))) #11
    ax = fig.add_subplot(111)
    sns.heatmap(log2p1(top.iloc[:,:len(classes)]), 
                cmap='RdBu_r',  
                linewidths=0.004)
    ax.set_title('Top different')
    for txt in ax.get_yticklabels():
            txt.set_rotation(0)
    for txt in ax.get_xticklabels():
            txt.set_rotation(90)
    
    fig.suptitle(ttl, size='x-large')
    
    return top, up, dn, gs



def cluster_express(df, phenoPos, phenoNeg, classes, n_top=0, geneList=[],  ttl=''):
    # ---
    A = phenoPos
    B = phenoNeg
    
    df = df.dropna()
    df = df.iloc[:,:len(classes)+1]


    df.set_index(keys=df.columns[0], inplace=True)
    df.index=df.index.str.upper()
    df2 = df.T
    df2['class'] = classes
    df_mean= df2.groupby('class').mean().T
    df_std = df2.groupby('class').std().T

    dfm = pd.concat([df,df_mean], axis=1)
    dfm['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))
    
    if not geneList:
        geneList = list(dfm.index)
    
    gs = dfm.iloc[dfm.index.isin(geneList)]
    
    up = gs[gs.lgFC > 1.]
    dn = gs[gs.lgFC < -1.]
    gs = pd.concat([up,dn])
    
    
    up_l = len(up)
    dn_l = len(dn)
    tt_l = len(gs)
    
    
    gs['FCmod'] = abs(gs['lgFC'])
    top = gs.sort_values('FCmod', axis=0, ascending=False)
    if n_top != 0:
        ntp = min(2*n_top,len(gs))
        top = top.head(ntp)


    dpic  = log2p1(top.iloc[:,:len(classes)])
    
    # --Cluster
    grid = sns.clustermap(dpic, cmap='RdBu_r', 
                        linewidths=0.004,
                        #metric='cosine',
                        #method='complete',
                        #z_score= 1,# 0 - row, 1 - column
                        col_cluster=False,
                        figsize=(8, max(11,int(n_top/9.))))
    for txt in grid.ax_heatmap.get_yticklabels():
        txt.set_rotation(0)
    for txt in grid.ax_heatmap.get_xticklabels():
            txt.set_rotation(90)
    grid.ax_heatmap.set_title(ttl+ ' diffExpressed' )
    # -----------
    
    return top, up, dn, gs



def scatter(df, phenoPos, phenoNeg, classes, n_top=0, geneList=[],  ttl='', names_term=' Genes'):
    # ---
    A = phenoPos
    B = phenoNeg
    
    df = df.dropna()
    df = df.iloc[:,:len(classes)+1]


    df.set_index(keys=df.columns[0], inplace=True)
    df.index=df.index.str.upper()
    df2 = df.T
    df2['class'] = classes
    df_mean= df2.groupby('class').mean().T
    df_std = df2.groupby('class').std().T

    dfm = pd.concat([df,df_mean], axis=1)
    dfm['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))
    print(dfm.head())
    print(df_mean.head())
    if not geneList:
        geneList = list(dfm.index)
    
    gs = dfm.iloc[dfm.index.isin(geneList)]

    gss = gs.sort_values('lgFC', 
                        axis=0, 
                        ascending=False)
    
    up = gs[gs.lgFC > 1.]
    dn = gs[gs.lgFC < -1.]
    
    upv = (gs.lgFC > 1.).values
    dnv = (gs.lgFC < -1.).values
    
    up_l = len(up)
    dn_l = len(dn)
    tt_l = len(gs)
    
    gs = ResultsTable(gs)
    # ====== Scatterplot 
    ax = gs.scatter(
        x= B,
        y= A,
        xfunc= log2p1,
        yfunc= log2p1,
        #----------------
        xlab= B +', log2(FPKM + 1)',
        ylab= A +', log2(FPKM + 1)',
        one_to_one=dict(color='k', linestyle=':'),
        marginal_histograms=False,
        general_kwargs=dict(marker='.', color='0.5', alpha=1., s=10, label='unchanged'),

        genes_to_highlight=[
        (
         upv,
         dict(
              color='#da3b3a', alpha=1,
              marginal_histograms=True,
              xhist_kwargs=dict(bins=50, linewidth=0),
              yhist_kwargs=dict(bins=50, linewidth=0),
              label='upregulated',
              )
         ),
        (
         dnv,
         dict(
              color='#00748e', alpha=1,
              marginal_histograms=True,
              xhist_kwargs=dict(bins=50, linewidth=0),
              yhist_kwargs=dict(bins=50, linewidth=0),
              label='downregulated'
              )
         )
        ]
    )
    
    # ==== Plot text of top genes
    for i in range(n_top):
        rx = .0003*np.random.randn()
        ry = .0003*np.random.randn()
        #~ print gss.index[i]
        ax.text(log2p1(gss[B][i]*(1.+rx)),
                log2p1(gss[A][i]*(1.+ry)),
                gss.index[i], 
                color='#da3b3a')
        ax.text(log2p1(gss[B][-i-1]*(1.+rx)),
                log2p1(gss[A][-i-1]*(1.+ry)),
                gss.index[-i-1], 
                color='#00748e')
                
    ax.text(0.1,0.95, 'UP: '+str(up_l),
            color='#da3b3a',
            transform=ax.transAxes)
    ax.text(0.3,0.95, 'DOWN: '+str(dn_l),
            color='#00748e',
            transform=ax.transAxes)
    ax.text(0.5,0.95, 'ALL: '+str(tt_l),
            transform=ax.transAxes)
    # Get handles and labels, and then reverse their order
    handles, legend_labels = ax.get_legend_handles_labels()
    handles = handles[::-1]
    legend_labels = legend_labels[::-1]

    # Draw a legend using the flipped handles and labels.
    leg = ax.legend(handles,
              legend_labels,
              loc=(1.01, 1.05),
              # Various style fixes to default legend.
              fontsize=9,
              scatterpoints=1,
              borderpad=0.1,
              handletextpad=0.05,
              frameon=False,
              title=names_term,
              );

    # Adjust the legend title after it's created
    leg.get_title().set_weight('bold')

    top_axes = gs.marginal.top_hists[-1]
    top_axes.set_title(ttl)

    for ax in gs.marginal.top_hists:
        ax.set_ylabel('No.\n genes', rotation=0, ha='right', va='center', size=8)

    for ax in gs.marginal.right_hists:
        ax.set_xlabel('No.\n genes', rotation=-90, ha='left', va='top', size=8)
    
    return up, dn, ax


## === Volcano plot

def volcano(df, phenoPos, phenoNeg, classes, n_top, geneList=[],  ttl=''):
    # ---
    A = phenoPos
    B = phenoNeg
    
    df = df.dropna()
    dfo = df.copy()
    df = df.iloc[:,:len(classes)+1]


    df.set_index(keys=df.columns[0], inplace=True)
    dfo.set_index(keys=dfo.columns[0], inplace=True)
    df.index=df.index.str.upper()
    dfo.index=dfo.index.str.upper()
    df2 = df.T
    df2['class'] = classes
    df_mean= df2.groupby('class').mean().T
    df_std = df2.groupby('class').std().T

    dfm = pd.concat([df,df_mean, dfo['padj']], axis=1)
    dfm['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))
    
    if not geneList:
        geneList = list(dfm.index)
    
    gs = dfm.iloc[dfm.index.isin(geneList)]

    gss = gs.sort_values('lgFC', 
                        axis=0, 
                        ascending=False)
    
    up = (gs.lgFC > 1.).values
    dn = (gs.lgFC < -1.).values
    
    
    up_l = len(gs[gs.lgFC > 1.])
    dn_l = len(gs[gs.lgFC < -1.])
    tt_l = len(gs)
    
    gs = ResultsTable(gs)
    # ====== Volcano-plot
    ax = gs.scatter(
        x='lgFC', 
        y = 'padj', 
        yfunc=nlog10,
        xlab='log FC',
        ylab='-log10($p_{adj}$)',
        general_kwargs=dict(marker='.', color='0.5', alpha=1., s=10, label='unchanged'),
        genes_to_highlight=[
        (up,
         dict(color='#da3b3a', alpha=1,
              label='upregulated')
         ),
        (dn,
         dict(color='#00748e', alpha=1,
              label='downregulated')
         )
        ])
    # ==== Plot text of top genes
    for i in range(n_top):
        rx = .003*np.random.randn()
        ax.text(gss['lgFC'][i]*(1.+rx),
                nlog10(gss['padj'][i]+1e-300),
                gss.index[i], 
                color='#da3b3a')
        ax.text(gss['lgFC'][-i-1]*(1.+rx),
                nlog10(gss['padj'][-i-1]+1e-300),
                gss.index[-i-1], 
                color='#00748e')
                
    ax.text(0.1,0.95, 'UP: '+str(up_l),
            color='#da3b3a',
            transform=ax.transAxes)
    ax.text(0.3,0.95, 'DOWN: '+str(dn_l),
            color='#00748e',
            transform=ax.transAxes)
    ax.text(0.5,0.95, 'ALL: '+str(tt_l),
            transform=ax.transAxes)

    ax.set_title(ttl)

def scatter_n(df, A, B, classes, n_top=0, ttl=''):
    df = df.copy()
    df_mean= df.groupby(by=classes, axis=1).mean()

    df['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))
    df['lg'+A] = np.log2(df_mean[A]+1)
    df['lg'+B] = np.log2(df_mean[B]+1)
    df = df.sort_values('lgFC', axis=0, ascending=False)
    
    up = df[df.lgFC > 1.]
    dn = df[df.lgFC < -1.]
    
    diff = A+"_vs_"+B

    df[diff]='unchanged'
    df.loc[df["lgFC"]>1,diff]='up'
    df.loc[df["lgFC"]<-1,diff]='down'

    cpp = {'up':(0.86, 0.23, 0.22), 'unchanged':(0.5,0.5,0.5),'down':(0.03, 0.45, 0.56)}
    f, ax = plt.subplots(figsize=(6.5, 6.5))
    sns.scatterplot(x = 'lg'+B, 
                    y = 'lg'+A, 
                    hue=diff, 
                    data=df, 
                    ax=ax, 
                    palette=cpp, 
                    linewidth=0, 
                    s=3.4)
    ax.set_xlabel(B +', log2(FPKM + 1)')
    ax.set_ylabel(A +', log2(FPKM + 1)')
    ax.set_title(A+'/'+B+' '+ttl)
    ax.axis('equal')
    ax.axis('tight')
    
    # ==== Plot text of top genes
    for i in range(n_top):
        rx = np.log2(1+.0003*np.random.randn())
        ry = np.log2(1+.0003*np.random.randn())
        ax.text(df['lg'+B][i]+rx,
                df['lg'+A][i]+ry,
                df.index[i], 
                color=cpp['up'])
        ax.text(df['lg'+B][-i-1]+rx,
                df['lg'+A][-i-1]+ry,
                df.index[-i-1], 
                color=cpp['down'])
                
    ax.text(0.37,0.95, 'UP: '+str(len(up)),
            color=cpp['up'],
            transform=ax.transAxes)
    ax.text(0.55,0.95, 'DOWN: '+str(len(dn)),
            color=cpp['down'],
            transform=ax.transAxes)
    ax.text(0.82,0.95, 'ALL: '+str(len(df)),
            transform=ax.transAxes)
    
    return up, dn, ax


def volcano_n(df, A, B, classes, n_top=0, ttl=''):
    df = df.copy()
    df_mean= df.groupby(by=classes, axis=1).mean()
    df['p-val'] = p_value(df, A, B, classes)
    df['p-adj'] = p_adjust(df['p-val'])
    df['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))

    pv = 'p-adj'#'p-val'
    df['-log10_pV'] = -np.log10(df[pv])
    
    df = df.sort_values('lgFC', axis=0, ascending=False)
    
    up = df[df.lgFC > 1.]
    dn = df[df.lgFC < -1.]
    
    diff = A+"_vs_"+B

    df[diff]='unchanged'
    df.loc[df["lgFC"]>1,diff]='up'
    df.loc[df["lgFC"]<-1,diff]='down'

    cpp = {'up':(0.86, 0.23, 0.22), 'unchanged':(0.5,0.5,0.5),'down':(0.03, 0.45, 0.56)}
    f, ax = plt.subplots(figsize=(6.5, 6.5))
    sns.scatterplot(x = 'lgFC', 
                    y = '-log10_pV', 
                    hue=diff, 
                    data=df, 
                    ax=ax, 
                    palette=cpp, 
                    linewidth=0, 
                    s=3.4)
    ax.set_xlabel('log(FC)')
    ax.set_ylabel('$-log_{10}(p_{val}$)')
    ax.set_title(A+'/'+B+' '+ttl)
    ax.axis('equal')
    ax.axis('tight')
    
    # ==== Plot text of top genes
    for i in range(n_top):
        rx = 1+.0003*np.random.randn()
        ax.text(df['lgFC'][i]*rx,
                df['-log10_pV'][i]+1e-300,
                df.index[i], 
                color=cpp['up'])
        ax.text(df['lgFC'][-i-1]*rx,
                df['-log10_pV'][-i-1]+1e-300,
                df.index[-i-1], 
                color=cpp['down'])
                
    ax.text(0.037,0.95, 'UP: '+str(len(up)),
            color=cpp['up'],
            transform=ax.transAxes)
    ax.text(0.25,0.95, 'DOWN: '+str(len(dn)),
            color=cpp['down'],
            transform=ax.transAxes)
    ax.text(0.52,0.95, 'ALL: '+str(len(df)),
            transform=ax.transAxes)
    
    return up, dn, ax


def cluster(df, A, B, classes, n_top=0):
    df = df.copy()
    cols = df.columns
    df_mean= df.groupby(by=classes, axis=1).mean()

    df['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))
    
    df = df[abs(df['lgFC'])>1].sort_values('lgFC', axis=0, ascending=False)
    
    if n_top == 0:
        n_top = len(df)
    df = pd.concat([df.head(int(n_top/2)),df.tail(int(n_top/2))], axis=0)
    
    # --Cluster
    dpic  = log2p1(df[cols])
    grid = sns.clustermap(dpic, cmap='RdBu_r', 
                        linewidths=0.004,
                        #metric='correlation',
                        method='centroid',
                        col_cluster=False,
                        figsize=(8, int(n_top/2)))
    return grid
    # -----------



def exp_deg(df, A, B, classes):
    df = df.copy()
    df_mean= df.groupby(by=classes, axis=1).mean()
    df['p-val'] = p_value(df, A, B, classes)
    df['p-adj'] = p_adjust(df['p-val'])
    df['lgFC'] = np.log2((df_mean[A]+1.) / (df_mean[B]+1.))
    df = df.sort_values('lgFC', axis=0, ascending=False)
    return df


## END ====

# And if need remove of duplicated index to unique:
# df1.index = df1.index + df1.groupby(level=0).cumcount().astype(str).replace('0','')

#~ In [13]: dff = pd.DataFrame({'index':['b','a','a','c'], 'values':[10,20,30,40]})
#~ In [14]: dff.set_index(['index'], inplace=True)
#~ In [16]: dff.index
#~ Out[16]: Index(['b', 'a', 'a', 'c'], dtype='object', name='index')

#~ In [17]: dff.index.get_loc('a')
#~ Out[17]: array([False,  True,  True, False])

