import numpy as np
from os.path import join 
import pandas as pd
import pybedtools as pb
#import gseapy as gp
#~ import seaborn as sns
#~ import matplotlib.pyplot as plt


from pyliftover import LiftOver


def log2p1(x):
    return np.log2(x + 1)


SAVE = True





# === Load file


prfx = 'bam_input/'
fn_tlx = 'TLX3_H3K27me3_all_RP.txt'
fn_tap = 'TAP_H3K27me3_all_RP.txt'


path = join('tracks/MARGE/dRP_PRC2/',prfx)




nm = ['chr','start','end', 'gene_id', 'raw_RP',  'gene_name', 'strand']

#~ df_rag = pd.read_table(path+fn_rag, names=nm)

#~ df_rag.rename(columns={'raw_RP':'RAG_raw_RP', 'rel_RP':'RAG_rel_RP'}, inplace=True)


df_tlx = pd.read_table(path+fn_tlx, names=nm)

df_tlx.rename(columns={'raw_RP':'TLX_raw_RP'}, inplace=True)


df_tap = pd.read_table(path+fn_tap, names=nm)

df_tap.rename(columns={'raw_RP':'TAP_raw_RP'}, inplace=True)

df = pd.concat([df_tlx[['chr','start','end', 'gene_id',  'gene_name', 'strand']],
                df_tlx[['TLX_raw_RP']], 
                df_tap[['TAP_raw_RP']]], 
                axis=1)

# -- Convert
lo = LiftOver('mm10', 'mm9')

def liftS(rw, col):
    lf = lo.convert_coordinate(rw['chr'],rw[col],rw['strand'])
    if len(lf) ==0:
        return np.nan, np.nan, np.nan
    else:
        return str(lf[0][0]), int(lf[0][1]), str(lf[0][2])


dfc = df


#dfc['bb'] = dfc.apply(lambda row: liftS(row)[1],axis=1)

dfc['chr_mm9'] = dfc.apply(lambda row: liftS(row,'start')[0],axis=1)
dfc['strand_mm9'] = dfc.apply(lambda row: liftS(row,'start')[2],axis=1)
dfc['start_mm9'] = dfc.apply(lambda row: liftS(row,'start')[1],axis=1)
dfc['end_mm9'] = dfc.apply(lambda row: liftS(row,'end')[1],axis=1)

dfc = dfc.dropna()

dfc = dfc.astype({'start_mm9':int, 'end_mm9':int})


## add fold changes
#~ dfc['lgFC_TLXvsRAG'] = np.log2(dfc['TLX_rel_RP']/dfc['RAG_rel_RP'])
#~ dfc['lgFC_TAPvsRAG'] = np.log2(dfc['TAP_rel_RP']/dfc['RAG_rel_RP'])
dfc['lgFC_TAPvsTLX'] = np.log2(dfc['TAP_raw_RP']/dfc['TLX_raw_RP'])

if SAVE:
    dfc.to_csv(path+'TLX_TAP_dRP_mm10mm9.txt', index=False, sep='\t')



### BED manipulation

#colb = ['chr','start','end', 'gene_name']
#colmm9 = ['chr_mm9','start_mm9','end_mm9', 'gene_name']

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






