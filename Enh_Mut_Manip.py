import numpy as np
import pandas as pd
import pybedtools as pb
import deeptools.getScorePerBigWigBin as gs


def enh_gene(genl,df):
    if isinstance(genl, str):
        genl = [genl]
    enhll = list()
    for gen in genl:
        reg = df[df['gene_name']==gen]
        if len(reg)>0:
            dr = reg.iloc[0]['enhancers'].split(' ')
            enhl = [x for x in dr if 'enh' in x]
        else:
            enhl=[]
        enhll = enhll+enhl
    return list(set(enhll))

def gene_enh(enhl,df):
    if isinstance(enhl, str):
        enhl = [enhl]
    
    genl = list()
    for enh in enhl:
        reg = df[df['enhancers']==enh]
        if len(reg)>0:
            ls = reg.iloc[0]['gene_name'].split(', ')
            gl=[tr.split(' (')[0] for tr in ls]
        else:
            gl=[]
        genl = genl+gl
    return list(set(genl))

def gene_variants(gl,var,genes):
    """
    Function counts variants in gene list 
    
    Parametrs
    ---------
    gl : genes list
    var : variant DataFrame; ex. var = allel.vcf_to_dataframe('fn.vcf')
    genes : bed with genes; ex. genes = pb.BedTool('genes.bed')
    """
    genes = genes.to_dataframe()
    genes.loc[:,'name'] = genes.loc[:,'name'].str.upper()
    genes_l = genes[genes['name'].isin(gl)]
    
    var_gl = pd.DataFrame()
    for idx, row in genes_l.iterrows():
        vg = var[(var['CHROM']==row['chrom']) & (var['POS']>row['start']) & (var['POS']<row['end'])]
        genes_l.loc[idx,'num_mut'] = int(len(vg))
        var_gl = pd.concat([var_gl, vg])
    return  var_gl, genes_l


#~ def bed_variants_counts(var,bed):
    #~ """
    #~ Function counts variants in region bed file
    
    #~ Parametrs
    #~ ---------
    #~ var : variant DataFrame; ex. var = allel.vcf_to_dataframe('fn.vcf')
    #~ bed : bed file with regions; ex. bed = pb.BedTool('genes.bed')
    #~ """
    #~ bed = bed.to_dataframe()
    
    #~ var_bd = pd.DataFrame()
    #~ for idx, row in bed.iterrows():
        #~ vg = var[(var['CHROM']==row['chrom']) & (var['POS']>row['start']) & (var['POS']<row['end'])]
        #~ bed.loc[idx,'counts'] = int(len(vg))
        #~ var_bd = pd.concat([var_bd, vg])
    #~ return  var_bd , bed



def bed_variants(var,bed):
    """
    Function colect variants in region bed file and out to variants file
    
    Parametrs
    ---------
    var : variant bed file; ex. var = pb.BedTool('fn.vcf')
    bed : bed file with regions; ex. bed = pb.BedTool('genes.bed')
    """
    var_out = var.intersect(bed, header=True)
    
    return  var_out


def variants_bed_counts(var,bed):
    """
    Function counts variants in region bed file
    
    Parametrs
    ---------
    var : variant bed file; ex. var = pb.BedTool('fn.vcf')
    bed : bed file with regions; ex. bed = pb.BedTool('genes.bed')

    Output
    ------
    bed_out : bed with last 'counts' column 
    """
    bed_out = bed.intersect(var, c=True)
    #~ bed_df = bed_out.to_dataframe()
    #~ bed_df = bed_df[['chrom', 'start', 'end','thickStart', 'name']]
    #~ bed_out = pb.BedTool.from_dataframe(bed_df)
    return bed_out



def freq_on_gene(var):
    """
    Parametrs
    ---------
    var : variant DataFrame with ANN filed from SnpEff; 
            ex. var = allel.vcf_to_dataframe('fn.vcf', transformers=allel.ANNTransformer())
    """
    fr = var[['ANN_Gene_Name']]
    fr['num_mut'] = 1
    frq = fr.groupby('ANN_Gene_Name').sum()
    frq.sort_values('num_mut', axis=0, inplace=True, ascending=False)
    
    return frq



def bigWig2bed_pot(bw,bed,genome,pad=1e5,alpha=10, step=1):
    """
    Function calculates "regulatory potential" of bigWig track around summit points of bed
    
    Parametes
    ---------
    bw : bigWig file, ex. 'TLX3_H327ac.bw'
    bed : bed file with regions; ex. bed = pb.BedTool('genes.bed')
    genome : genome; ex. 'mm10'
    pad : padding; -pad|------summit------|+pad
    alpha :  exponential parameter for decay fucntion of weights 
    
    Returns
    ------
    df :  dataframe with potentials for all regions in bed
    
    Notes
    -----
    Weights
    .. math :: w(x)=\frac{2e^{-\alpha|x|}}{1+e^{-\alpha|x|}}
    """
    
    bed_df = bed.to_dataframe()
    bed_df['mid'] = (bed_df['end'] + bed_df['start'])/2
    bed_df['mid'] = bed_df['mid'].astype('int')
    bed_df['mid_end'] = bed_df['mid'] + 1

    col = ['chrom','mid','mid_end', 'name']

    tss = pb.BedTool.from_dataframe(bed_df[col])

    pad = int(pad)
    z = np.arange(-pad,pad+1, step)
    wt = 2.0*np.exp(-alpha*np.fabs(z)/1e5)/(1.0+np.exp(-alpha*np.fabs(z)/1e5))

    df = tss.slop(b=pad, genome=genome).to_dataframe()

    for i, row in df.iterrows():

        cnt = tss[i].start
        st = df.loc[i,'start']
        end = df.loc[i,'end']
        chrom = df.loc[i,'chrom']

        vl = gs.countFragmentsInRegions_worker(chrom, int(st), int(end), [bw], step, step, False)
        vl = np.transpose(np.squeeze(vl[0]))

        vl  = np.hstack((np.zeros(st - cnt + pad),vl,np.zeros(cnt + pad - end +1 )))
        df.loc[i,'potential'] = np.dot( vl, wt)
    
    return df


def deepbind(fa,md):
    """
    Parametes
    ---------
    md : ID of models
    fs : sequence, e.g. 'AAGTAAGCTGAACC'
    """
    cmd = 'echo {} | deepbind --no-head {}'.format(fa,md)
    res = os.popen(cmd).readlines()
    score = float(res[0].replace('\n',''))
    
    return score



def deepbind_list(md_l,fs):
    """
    Parametes
    ---------
    md_l : list of IDs of models
    fs : sequence, e.g. 'AAGTAAGCTGAACC'
    """
    df_md = pd.DataFrame(md_l, columns=['model'])
    df_md['score'] = df_md.apply(lambda row: deepbind(fs,row['model']), axis=1)
    
    return df_md 




def snp(seq,pos,nt):
    """
    Blind SNP
    """
    seq = MutableSeq(seq)
    seq[pos]=nt
    return str(seq)


def mut(fa,i,ref,alt):
    """
    General type mutation
    """
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
                          #fontsize=80, 
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
