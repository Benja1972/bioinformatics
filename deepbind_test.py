import sys,os
from Bio.Seq import MutableSeq
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects
from pyfaidx import Fasta
import allel


import warnings
warnings.filterwarnings('ignore')

### == Math functions == 
def log2p1(x):
    return np.log2(x + 1)



def deepbind(fa,md):
    cmd = 'echo {} | deepbind --no-head {}'.format(fa,md)
    res = os.popen(cmd).readlines()
    score = float(res[0].replace('\n',''))
    
    return score



def deepbind_list(fa,md_id):
    #cmd = 'echo {} | deepbind {}'.format(fa,' '.join(md_id))
    #res = os.popen(cmd).readlines()
    #cols = res[0].replace('\n','').split('\t')
    #vals = [float(i) for i in res[1].replace('\n','').split('\t')]
    score = pd.DataFrame([], index=md_id, columns=['score'])
    for idm in md_id:
        score.loc[idm,'score'] = deepbind(fs,idm)
    
    return score




def snp(seq,pos,nt):
    seq = MutableSeq(seq)
    seq[pos]=nt
    return str(seq)


def mut(fa,i,ref,alt):
    seq = MutableSeq(fa)
    ref = MutableSeq(ref)
    alt =  MutableSeq(alt)
    
    ref_s = seq[i:i+len(ref)]
    if (ref==ref_s):
        seq_mut = seq[:i]+alt+seq[i+len(ref):]
    else:
        seq_mut = MutableSeq('')
        print('It is not correct variant')
        
    return str(seq_mut)

def mut_map(fa,md):
    p0 = deepbind(fa,md)

    df = pd.DataFrame(index=['A','C','G','T'],columns=list(fa))
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            nt = df.index[i]
            #fm = snp(fa,j,nt)
            p1 = deepbind(snp(fa,j,nt),md)
            df.iloc[i,j] = (p1 - p0)*max([0.,p0,p1])
            
    df = df.astype(float)
    
    return df



def draw_motif(sm,ax):
    class Scale(matplotlib.patheffects.RendererBase):
        def __init__(self, sx, sy=None):
            self._sx = sx
            self._sy = sy

        def draw_path(self, renderer, gc, tpath, affine, rgbFace):
            affine = affine.identity().scale(self._sx, self._sy)+affine
            renderer.draw_path(gc, tpath, affine, rgbFace)


    colors = {'G': 'orange', 'A': 'darkgreen', 'C': 'blue', 'T': 'red'}

    ax.set_xlim(0,len(sm))
    ax.set_ylim(0,max(sm))

    for i in range(sm.shape[0]):
        nt = sm.index[i]
        txt = ax.text(i,0, nt, fontsize=18, weight='bold', color=colors[nt])
        txt.set_path_effects([Scale(1, sm[i]+1)])
    return ax

def plot_mutmap(fs,mdl):
    df = mut_map(fs,mdl)
    sc = deepbind(fs,mdl)
    
    fig = plt.figure(figsize=(10, 4))
    
    ax1 = fig.add_subplot(211)
    sns.heatmap(df, cmap='RdBu_r', ax=ax1, cbar = False, square=True)
    
    sm = (abs(df[df<0].sum(axis=0))+abs(sc))/abs(sc)
    
    ax2 = fig.add_subplot(212, sharex=ax1)
    draw_motif(sm,ax2)
    return fig


### main ###

db_f = '/home/sergio/tools/deepbind/db/db.tsv'

db = pd.read_csv(db_f, sep='\t',comment="#", index_col=0)
db = db[db['Labels'].isnull()]


#~ fa = 'CGTAAAGCCCTTGATAAACCCCTTCCCTGGA'
#~ md = 'D00411.003'

#~ sc0 = deepbind(fa,md)


#~ md_l = list(db.index[:25])

#~ sc_l = deepbind_list(fa,md_l)

#~ df = mut_map(fa,md)


fn = '/home/sergio/media/NAS4/PFlab/TLX3_project/WES-seq/references/mouse_mm9_reference_genome.fa'
vn = '/home/sergio/Res_CIML/TLX3_project/data/tracks/WGS-WES/Exomes/WES_TLX3_TAP.vcf'


var = allel.vcf_to_dataframe(vn,fields='*', numbers={'ALT': 2})
vnt  = var.loc[1]

fa = Fasta(fn)
rg = fa['chr1'][6372606:6372646]

if rg.name== vnt['CHROM']:
    pos = vnt['POS'] -rg.start
else:
    print('It is not correct variant')
    pos = np.nan

ref = vnt['REF'].upper()
alt = vnt['ALT_1'].upper()


fs = rg.seq.upper()
fm = mut(fs,pos,ref,alt)


print('pos = ', pos, 'REF = ', ref, 'ALT =', alt) 
print('fs = ', fs)
print('fm = ', fm)


# LIST

#sca = deepbind_list(fs,list(db.index))
#db['score'] = sca['score']

db['model'] = db.index
db['score'] = db.apply(lambda row: deepbind(fs,row['model']), axis=1)

dbl_s = db.sort_values('score',axis=0, ascending=False)

# best model
mdl = dbl_s.index[0]

#~ sc = deepbind(fs,mdl)
#~ df = mut_map(fs,mdl)







### Figures ####
if True:
    #~ fig = plt.figure(figsize=(10, 4))

    #~ ax1 = fig.add_subplot(211)
    #~ sns.heatmap(df, cmap='RdBu_r', ax=ax1, cbar = False, square=True)

    #~ sm = (abs(df[df<0].sum(axis=0))+sc)/(sc)

    #~ ax2 = fig.add_subplot(212, sharex=ax1)
    #~ draw_motif(sm,ax2)
    
    plot_mutmap(fs,mdl)
    plot_mutmap(fm,mdl)
    plt.show()

