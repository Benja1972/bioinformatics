import numpy as np
from os.path import join 
import pandas as pd
#~ import pybedtools as pb
#import gseapy as gp
#~ import seaborn as sns
import matplotlib.pyplot as plt
import RNA_expression_processing as rn

#~ from pyliftover import LiftOver


def log2p1(x):
    return np.log2(x + 1)


SAVE = False

# Project settings
from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')

#RP = join(DATADIR,'tracks/MARGE/relativeRP/bam_input')


# === Load file

fn_rag = join(DATADIR,'tracks/TSS_RAG_H3K27ac_RP_3K.csv')
fn_tlx = join(DATADIR,'tracks/TSS_TLX3_H3K27ac_RP_3K.csv')




tlx_27ac = pd.read_csv(fn_tlx, index_col=0)
rag_27ac = pd.read_csv(fn_rag, index_col=0)

tlx_27ac.rename(columns={'potential':'RP_tlx'}, inplace=True)

rag_27ac.rename(columns={'potential':'RP_rag'}, inplace=True)

tlx_rag_27ac = pd.concat([tlx_27ac,rag_27ac.drop(['chrom','start','end', 'name'], axis=1)], axis=1)

tlx_rag_27ac.drop_duplicates(subset='name', inplace=True)

tlx_rag_27ac = tlx_rag_27ac[(tlx_rag_27ac['RP_rag']>0) | (tlx_rag_27ac['RP_tlx']>0)]


## Scatter of 27ac potential 

Ap,Bp = 'RP_tlx','RP_rag'
cols = ['name', Ap, Bp]

tlx_rag_27ac_s = tlx_rag_27ac[cols]

tlx_rag_27ac_s = tlx_rag_27ac_s.set_index(keys=tlx_rag_27ac_s.columns[0])

tlx_rag_27ac_s = tlx_rag_27ac_s.apply(pd.to_numeric)

Ac= Ap+'c'
Bc= Bp+'c'

classes=[Ac,Bc]
#df_mean= tlx_rag_27ac_s.groupby(by=classes, axis=1).mean()


df=rn.scatter_n(tlx_rag_27ac_s, Ac, Bc,classes=classes, n_top=16)





fn1_rag = 'tracks/MARGE/relativeRP/bam_input/RAG_H3K27ac_all_relativeRP.txt'
fn1_tlx = 'tracks/MARGE/relativeRP/bam_input/TLX3_H3K27ac_all_relativeRP.txt'


#~ fn1_rag = 'tracks/MARGE/RP/RAG_H3K27ac_repl1_all_RP.txt'
#~ fn1_tlx = 'tracks/MARGE/RP/TLX3_H3K27ac_repl1_all_RP.txt'


nm = ['chr','start','end', 'gene_id', 'raw_RP', 'rel_RP', 'gene_name', 'strand']
#nm = ['chr','start','end', 'gene_id', 'raw_RP', 'gene_name', 'strand']

df1_rag = pd.read_csv(join(DATADIR,fn1_rag), names=nm, sep='\t')
df1_tlx = pd.read_csv(join(DATADIR,fn1_tlx), names=nm, sep='\t')

df1_rag.rename(columns={'raw_RP':'RP_rag'}, inplace=True)
df1_tlx.rename(columns={'raw_RP':'RP_tlx'}, inplace=True)

df1_tlx_rag = pd.concat([df1_tlx,df1_rag.drop(['chr','start',
                                                'end', 'gene_id', 
                                                'rel_RP',
                                                'gene_name', 'strand'], axis=1)], axis=1)



df1_tlx_rag.drop_duplicates(subset='gene_name', inplace=True)

Ap,Bp = 'RP_tlx','RP_rag'
cols = ['gene_name', Ap, Bp]

df1_tlx_rag_s = df1_tlx_rag[cols]

df1_tlx_rag_s = df1_tlx_rag_s[(df1_tlx_rag_s['RP_rag']>0) | (df1_tlx_rag_s['RP_tlx']>0)]

df1_tlx_rag_s = df1_tlx_rag_s.set_index(keys=df1_tlx_rag_s.columns[0])

df1_tlx_rag_s = df1_tlx_rag_s.apply(pd.to_numeric)

Ac= Ap+'c'
Bc= Bp+'c'

classes=[Ac,Bc]

df1=rn.scatter_n(df1_tlx_rag_s, Ac, Bc,classes=classes, n_top=16)


#~ df_rag.rename(columns={'raw_RP':'RAG_raw_RP', 'rel_RP':'RAG_rel_RP'}, inplace=True)


#~ df_tlx = pd.read_table(path+fn_tlx, names=nm)

#~ df_tlx.rename(columns={'raw_RP':'TLX_raw_RP', 'rel_RP':'TLX_rel_RP'}, inplace=True)


#~ df_tap = pd.read_table(path+fn_tap, names=nm)

#~ df_tap.rename(columns={'raw_RP':'TAP_raw_RP', 'rel_RP':'TAP_rel_RP'}, inplace=True)

#~ df = pd.concat([df_rag[['chr','start','end', 'gene_id',  'gene_name', 'strand']],
                #~ df_rag[['RAG_raw_RP', 'RAG_rel_RP']], 
                #~ df_tlx[['TLX_raw_RP', 'TLX_rel_RP']], 
                #~ df_tap[['TAP_raw_RP', 'TAP_rel_RP']]], 
                #~ axis=1)



#~ dfc = dfc.dropna()

#~ dfc = dfc.astype({'start_mm9':int, 'end_mm9':int})


#~ ## add fold changes
#~ dfc['lgFC_TLXvsRAG'] = np.log2(dfc['TLX_rel_RP']/dfc['RAG_rel_RP'])


#~ if SAVE:
    #~ dfc.to_csv(join(DATADIR,'TLX3-vs-RAG_TSS_3K_RP_mm9.txt', index=False, sep='\t'))



### BED manipulation

#colb = ['chr','start','end', 'gene_name']
#~ colmm9 = ['chr_mm9','start_mm9','end_mm9', 'gene_name']

#~ dfb = dfc[colb]
#~ dfm = dfc[colmm9]
#~ dfm = dfm[dfm['end_mm9']-dfm['start_mm9']>0]

#rp_mm10 = pb.BedTool.from_dataframe(dfb) 
#~ rp_mm9 = pb.BedTool.from_dataframe(dfm)  #.saveas('tmp.bed')




#~ tlx_peak = pb.BedTool('tracks/TLX3_TLX3_peaks.bed')
#~ sl = 100
#~ tlx_peak = tlx_peak.slop(b=sl, genome='mm9')

#~ rp_mm9 = rp_mm9.slop(b=2000, genome='mm9')


#~ rp_tlx3 = rp_mm9+tlx_peak



plt.show()



