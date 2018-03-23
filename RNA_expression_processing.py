from __future__ import  division
from functools import reduce
from os.path import join 

import numpy as np
import sys, logging
import multiprocessing
import itertools
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# comment metaseq for Python3 (works only express)
#from metaseq.results_table import ResultsTable

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

# == Load gene set in gmt format
def read_gmt(name, path=''):
    with open(join(path,name)) as f:
         gene_list = f.read().split()[2:]
    return list(gene_list)

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



def express(df, phenoPos, phenoNeg, classes, n_top=0, geneList=[],  ttl=''):
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
    
    
    up_l = len(up)
    dn_l = len(dn)
    tt_l = len(gs)
    
    #print tt_l, len(geneList) 
    
    gs['FCmod'] = abs(gs['lgFC'])
    #top = gs.copy()
    top = gs.sort_values('FCmod', axis=0, ascending=False)
    if n_top==0:
        n_top=len(gs)
        top = top[:n_top]
    else:
        top = top[:n_top].sort_values(A, axis=0, ascending=False)

    
    fig = plt.figure(figsize=(6, max(11,int(n_top/9.)))) #11
    ax = fig.add_subplot(111)
    sns.heatmap(log2p1(top.iloc[:,:len(classes)]), 
                cmap='RdBu_r',  
                linewidths=0.004)
    ax.set_title('diffExpressed')
    for txt in ax.get_yticklabels():
            txt.set_rotation(0)
    for txt in ax.get_xticklabels():
            txt.set_rotation(90)
    
    fig.suptitle(ttl, size='x-large')
    
    return top, up, dn, gs



def scatter(df, phenoPos, phenoNeg, classes, n_top, geneList=[],  ttl=''):
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
              title=' Genes',
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



## END ====
