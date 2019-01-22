import numpy as np
from os.path import join 
import pandas as pd
import pybedtools as pb
#import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt

# === User defined
import RNA_expression_processing as rn


def log2p1(x):
    return np.log2(x + 1)


SAVE = False

ntop=100


# === Load file
path = 'tracks/MARGE/relativeRP/DUXBL/'
#~ TLX3_regulated_TF2DNA.csv

df  =  pd.read_table(path+'TLX3_regulated_TF2DNA.csv')
#df  =  pd.read_table(path+'RAG_TLX_TAP_relativeRP_mm10mm9.txt')

# -- transform
#dfs = df.sort_values('TLX_rel_RP', axis=0, ascending=False)
#dfs.drop_duplicates(subset='gene_name', inplace=True)
#dfs.index=dfs['gene_name']




#~ A = 'TLX_rel_RP'
#~ if SAVE:
    #~ from matplotlib.backends.backend_pdf import PdfPages
    #~ pp = PdfPages(path+str(ntop)+'_RegPoten_TLX3.pdf')

#~ for B in ['RAG_rel_RP']: # ['RAG_rel_RP', 'TAP_rel_RP']:

    #~ cols = ['gene_name', A, B]

    #~ dfp = dfs[cols]

    #~ # classes
    #~ Ac = A+'c'
    #~ Bc = B+'c'
    #~ classe = [A+'c', B+'c']

    #~ rn.scatter(dfp, Ac, Bc, classes=classe, n_top=5, geneList=[],  ttl=A+'/'+B, names_term='Gene')
    #~ if SAVE:
        #~ plt.savefig(pp, format='pdf')
    #~ top, up, dn, gs = rn.express(dfp, Ac, Bc, classes=classe, n_top=ntop, geneList=[],  ttl=A+'/'+B)
    #~ if SAVE:
        #~ top.to_csv(path+'Top'+str(ntop)+'_'+A+'_vs_'+B+'.csv')
        #~ plt.savefig(pp, format='pdf')





## Expression analysis

# === Load expression table 
tbl = pd.read_table(join('tracks', 'TLX3vsRAG-results_genes.txt'), index_col=0)
tbl = tbl[(tbl.padj < 0.05)].dropna()
tbl = tbl.dropna()



# === Load gene names 
names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           index_col=0,
                           header=0,
                           names=['GeneID', 'TransID', 'Gene_name'])


names = names.drop('TransID', axis=1).drop_duplicates()
names = names.loc[tbl.index]
assert names.shape[0] == tbl.shape[0]

tbl=names.join(tbl, how ='right')


tbn = tbl[['Gene_name', 'TLX3.1_1','TLX3.1_5','TLX3.1_P','R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ', 'padj']]


## === Expresion analysis
classes = ['TLX3','TLX3','TLX3','RAG','RAG','RAG']


import RNA_expression_processing as rn

gl = list(df['gene_symbol'])
gl.insert(0,'TCF3')
#~ print(gl)


topN, upN, dnN, gsN = rn.express(tbn, 'TLX3', 'RAG', 
                    classes=classes, 
                    geneList=gl,  
                    ttl='ALL genes',
                    n_top=120,
                    sort=True)
if SAVE:
    plt.savefig(pp, format='pdf')

rn.scatter(tbn, 'TLX3', 'RAG', classes=classes, n_top=10, geneList=gl,  ttl='ALL genes')
#~ rn.scatter(tbn, 'TLX3', 'RAG', classes=classes, n_top=5, geneList=[],  ttl='ALL genes')
#~ rn.volcano(tbn, 'TLX3', 'RAG', classes=classes, n_top=5, geneList=[],  ttl='ALL genes')






#~ if SAVE:
    #~ dfs.head(500)['gene_name'].to_csv('topRPtlx3.txt')



if SAVE:
    pp.close()


plt.show()
