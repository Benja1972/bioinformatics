import sys,os
from Bio.Seq import MutableSeq
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects
from pyfaidx import Fasta
import allel

from matplotlib import transforms
from matplotlib.font_manager import FontProperties
import matplotlib as mpl



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



#~ def deepbind_list(fa,md_id):
    #~ #cmd = 'echo {} | deepbind {}'.format(fa,' '.join(md_id))
    #~ #res = os.popen(cmd).readlines()
    #~ #cols = res[0].replace('\n','').split('\t')
    #~ #vals = [float(i) for i in res[1].replace('\n','').split('\t')]
    #~ score = pd.DataFrame([], index=md_id, columns=['score'])
    #~ for idm in md_id:
        #~ score.loc[idm,'score'] = deepbind(fs,idm)
    
    #~ return score




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
            p1 = deepbind(snp(fa,j,nt),md)
            df.iloc[i,j] = (p1 - p0)*max([0.,p0,p1])
            
    df = df.astype(float)
    
    return df

class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def draw_motif(sm,ax):

    colors = {'G': 'orange', 'A': 'darkgreen', 'C': 'blue', 'T': 'red'}

    ax.set_xlim(0,len(sm))
    ax.set_ylim(0,max(sm))

    for i in range(sm.shape[0]):
        nt = sm.index[i]
        txt = ax.text(i,0, nt, fontsize=18, weight='bold', color=colors[nt])
        txt.set_path_effects([Scale(1, sm[i])])
    return ax


def draw_logo(all_scores, fontfamily='Arial', size=80):

    mpl.rcParams['font.family'] = fontfamily
    
    colors = {'G': 'orange', 'A': 'darkgreen', 'C': 'blue', 'T': 'red'}

    fig, ax = plt.subplots(figsize=(len(all_scores), 2.5))

    font = FontProperties()
    font.set_size(size)
    font.set_weight('bold')

    ax.set_xticks(range(1,len(all_scores)+1))    
    ax.set_yticks(range(0,3))
    ax.set_xticklabels(range(1,len(all_scores)+1), rotation=90)
    ax.set_yticklabels(np.arange(0,3,1))    
    sns.despine(ax=ax, trim=True)
    
    trans_offset = transforms.offset_copy(ax.transData, 
                                          fig=fig, 
                                          x=1, 
                                          y=0, 
                                          units='dots')
   
    for index, scores in enumerate(all_scores):
        yshift = 0
        for base, score in scores:
            txt = ax.text(index+1, 
                          0, 
                          base, 
                          transform=trans_offset,
                          fontsize=80, 
                          color=colors[base],
                          ha='center',
                          fontproperties=font,

                         )
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height*score
            trans_offset = transforms.offset_copy(txt._transform, 
                                                  fig=fig,
                                                  y=yshift,
                                                  units='points')
        trans_offset = transforms.offset_copy(ax.transData, 
                                              fig=fig, 
                                              x=1, 
                                              y=0, 
                                              units='points')
    return fig


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

## test model from paper
#~ fa = 'CGTAAAGCCCTTGATAAACCCCTTCCCTGGA'
#~ md = 'D00411.003'



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


### Find best model

db['model'] = db.index
db['score'] = db.apply(lambda row: deepbind(fs,row['model']), axis=1)

dbl_s = db.sort_values('score',axis=0, ascending=False)

bst=3
mdl = dbl_s.index[bst]



scol = [[('C', 0.02247014831444764),
          ('T', 0.057903843733384308),
          ('A', 0.10370837683591219),
          ('G', 0.24803586793255664)],
         [('T', 0.046608227674354567),
          ('G', 0.048827667087419063),
          ('A', 0.084338697696451109),
          ('C', 0.92994511407402669)],
         [('G', 0.0),
          ('T', 0.011098351287382456),
          ('A', 0.022196702574764911),
          ('C', 1.8164301607015951)],
         [('C', 0.020803153636453006),
          ('T', 0.078011826136698756),
          ('G', 0.11268374886412044),
          ('A', 0.65529933954826969)],
         [('T', 0.017393530660176126),
          ('A', 0.030438678655308221),
          ('G', 0.22611589858228964),
          ('C', 0.45078233627623127)],
         [('G', 0.022364103549245576),
          ('A', 0.043412671595594352),
          ('T', 0.097349627214363091),
          ('C', 0.1657574733649966)],
         [('C', 0.03264675899941203),
          ('T', 0.045203204768416654),
          ('G', 0.082872542075430544),
          ('A', 1.0949220710572034)],
         [('C', 0.0),
          ('T', 0.0076232429756614498),
          ('A', 0.011434864463492175),
          ('G', 1.8867526364762088)],
         [('C', 0.0018955903000026028),
          ('T', 0.0094779515000130137),
          ('A', 0.35637097640048931),
          ('G', 0.58005063180079641)],
         [('A', 0.01594690817903021),
          ('C', 0.017541598996933229),
          ('T', 0.2774762023151256),
          ('G', 0.48638069946042134)],
         [('A', 0.003770051401807444),
          ('C', 0.0075401028036148881),
          ('T', 0.011310154205422331),
          ('G', 1.8624053924928772)],
         [('C', 0.036479877757360731),
          ('A', 0.041691288865555121),
          ('T', 0.072959755514721461),
          ('G', 1.1517218549109602)],
         [('G', 0.011831087684038642),
          ('T', 0.068620308567424126),
          ('A', 0.10174735408273231),
          ('C', 1.0009100180696691)],
         [('C', 0.015871770937774379),
          ('T', 0.018757547471915176),
          ('A', 0.32176408355669878),
          ('G', 0.36505073156881074)],
         [('A', 0.022798100897300954),
          ('T', 0.024064662058262118),
          ('G', 0.24571286522646588),
          ('C', 0.34070495229855319)]]


### Figures ####
if True:
    
    plot_mutmap(fs,mdl)
    #plt.title(dbl_s.iloc[bst]['Protein'])
    plot_mutmap(fm,mdl)
    #plt.title(dbl_s.iloc[bst]['Protein'])
    draw_logo(scol)
    plt.show()

