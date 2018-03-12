#!/usr/bin/env python3

import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import metaseq
from mpl_toolkits.axes_grid1 import ImageGrid

## === Bed manipulation
path = 'tracks/state14_RXT/14/'

rxt= pb.BedTool(path+'RXT_14_segments.bed')
rxx = rxt.slop(b=50, genome='mm9')

tlx_peak = pb.BedTool('tracks/TLX3_TLX3_peaks.bed')

rxt_df = rxx.to_dataframe()

# --- States extraction
trj = '9'
#rxt_st1 = rxt_df[rxt_df['name'].isin(['E1','E2'])]
rxt_st1 = rxt_df[rxt_df['name']=='E'+trj]

#~ fig =plt.figure(10)
#~ lend =pd.DataFrame()
#~ lend['len'] = rxt_df['end'] - rxt_df['start']

#~ lend = lend[lend['len']<4000]
#~ lend = lend[lend['len']>1000]  
#~ print lend.head(100)
#plt.hist(np.log(lend['len']), bins=150)
#~ plt.hist(lend['len'], bins=150)


rxt_st1 = rxt_st1[rxt_st1['end']-rxt_st1['start']<4000]
rxt_st1 = rxt_st1[rxt_st1['end']-rxt_st1['start']>2000]

#~ st1 = pb.BedTool.from_dataframe(rxt_st1)
#~ reg = (tlx_peak+st1).sort()
#~ #reg = (st1-tlx_peak).sort()
#~ reg_dt = reg.to_dataframe()

reg_dt = rxt_st1.copy()

reg_dt['mid'] = ((reg_dt['end'] + reg_dt['start'])/2).astype('int')

reg_sm = reg_dt.copy()

reg_sm['start'],reg_sm['end'] = reg_sm['mid'],reg_sm['mid']

reg_sm = reg_sm.drop('mid', axis=1)


sm = pb.BedTool.from_dataframe(reg_sm)

shr = 2000
sm_Nkb = sm.slop(b=shr, genome='mm9')#, output='tsses-1kb.gtf')


## ==== Precallculation ====
import multiprocessing
processes = multiprocessing.cpu_count()

bn= 100
x = np.linspace(-shr, shr, bn)









tracks = [  'tracks/ChiP-seq_tracks/TLX3_TLX3_FE.bw',
            'tracks/ChiP-seq_tracks/RAG_TLX3_repl1_FE.bw',
            'tracks/ChiP-seq_tracks/TAP_TLX3_repl2_FE.bw',
            'tracks/ChiP-seq_tracks/TLX3_POLII_FE.bw',
            'tracks/ChiP-seq_tracks/RAG_PolII_FK97-98_sort.bw',
            'tracks/ChiP-seq_tracks/TAP_POLII_repl2_FE.bw',
            'tracks/ChiP-seq_tracks/TLX3_H3K27ac_FE.bw',
            'tracks/ChiP-seq_tracks/RAG_H3K27ac_repl1_FE.bw',
            'tracks/ChiP-seq_tracks/TAP_H3K27ac_repl2_FE.bw',
            'tracks/ChiP-seq_tracks/TLX3_H3K4me1_FE.bw',
            'tracks/ChiP-seq_tracks/RAG_H3K4me1_repl1_FE.bw',
            'tracks/ChiP-seq_tracks/TAP_H3K4me1_repl2_FE.bw',
            'tracks/ChiP-seq_tracks/TLX3_H3K4me2_FE.bw',
            'tracks/ChiP-seq_tracks/RAG_H3K4me2_FK89_sort.bw',
            'tracks/ChiP-seq_tracks/TAP_H3K4me2_repl2_FE.bw'
            ]


fig = plt.figure(1, (26., 14.))
grid = ImageGrid(fig, 111, nrows_ncols=(1, len(tracks)), # nrows_ncols=(1, 5),
                 axes_pad=0.1,
                 add_all=True,
                 label_mode="R",
                 aspect = False)


for i,tr in enumerate(tracks):
    tr_sig = metaseq.genomic_signal(tr,'bigwig')
    arr = tr_sig.array(sm_Nkb, bins=bn, processes=processes)
    if i==0:
        #i_a = arr.max(axis=1).argsort()
        i_a = arr.mean(axis=1).argsort()
        arr = arr[i_a,:]
    else:
        arr = arr[i_a,:]
    ax=grid[i]
    fig = metaseq.plotutils.imshow(
            arr,
            x=x, 
            figsize=(4, 10),
            ax=ax,
            vmin=5, vmax=99,  percentile=True,
            line_kwargs=dict(color='k', label='All'),
            fill_kwargs=dict(color='k', alpha=0.3))
            #sort_by=arr.max(axis=1)) # mean(axis=1))
    tit = '_'.join(tr.split('/')[-1].split('_')[:2])
    ax.set_title(tit)
    if i==0:
        ax.set_ylabel('trajectory '+trj+' ('+str(len(reg_dt))+')')
    #~ fig.array_axes.set_title(tr.split('/')[-1])
    #~ fig.array_axes.set_ylabel('trajectory '+trj+' ('+str(len(reg))+')')

plt.show()

#~ fig2 = metaseq.plotutils.imshow(
    #~ enh_tlx_arr,
    #~ x=x,
    #~ figsize=(5, 8),
    #~ vmin=5, vmax=99,  percentile=True,
    #~ line_kwargs=dict(color='k', label='All'),
    #~ fill_kwargs=dict(color='k', alpha=0.3),
    #~ sort_by=enh_tlx_arr.mean(axis=1)
#~ )




#~ rxt_gr = rxt_df.groupby(['name'])
#~ rxt_st1 = rxt_gr.get_group('E1')

#~ nf = 'RXT_14_st1.bed'
#~ rxt_st1[['chrom','start','end', 'name']].to_csv(join(path,nf), 
                                    #~ sep='\t', 
                                    #~ index=False,
                                    #~ header=False)



#~ st1 = pb.BedTool(join(path,nf))

#~ rxt_df['len'] = rxt_df['end'] - rxt_df['start']

#~ rxt_df['mid'] = (rxt_df['end'] + rxt_df['start'])/2
#~ rxt_df['mid'] = rxt_df['mid'].astype('int')

#~ rxt_df['mide'] = rxt_df['mid']


#~ plt.hist(rxt_df['len'], bins=50)

#~ plt.show()



#~ enh_df[['chrom','mid','mide']].to_csv('tracks/TLX3_6_FE_E4_mid.bed', 
                                    #~ sep='\t', 
                                    #~ index=False,
                                    #~ header=False)


#~ enhm= pb.BedTool('tracks/TLX3_6_FE_E4_mid.bed')

#~ shr = 2000
#~ enhm_Nkb = enhm.slop(b=shr, genome='mm9')#, output='tsses-1kb.gtf')



## ==== Precallculation ====
#~ import multiprocessing
#~ processes = multiprocessing.cpu_count()

#~ tlx_TLX3 = 'tracks/ChiP-seq_tracks/TLX3_TLX3_FE.bw'
#~ tlx_TLX3_sig = metaseq.genomic_signal(tlx_TLX3,'bigwig')

#~ bn= 100
#~ enh_tlx_arr = tlx_TLX3_sig.array(enhm_Nkb, bins=bn, processes=processes)

#~ metaseq.persistence.save_features_and_arrays(
        #~ features=enhm,
        #~ arrays={'enh_tlx': enh_tlx_arr},
        #~ prefix='enh_tlx1K',
        #~ link_features=True,
        #~ overwrite=True)

#~ features, arrays = metaseq.persistence.load_features_and_arrays(prefix='enh_tlx1K')


#~ enh_tlx_arr = arrays['enh_tlx']
#~ bn = enh_tlx_arr.shape[1]


#~ enh_df['tlx_sign2K']=enh_tlx_arr.mean(axis=1)

#~ enh_sr = enh_df.sort_values('tlx_sign2K', axis=0, ascending=False)

#~ top1K = enh_sr.head(1000)

#~ top500 = enh_sr.head(500)

#~ top1K[['chrom','start','end']].to_csv('tracks/enh_chromm_top1K_tlx_sign.bed', 
                                    #~ sep='\t', 
                                    #~ index=False,
                                    #~ header=False)

#~ top500[['chrom','start','end']].to_csv('tracks/enh_chromm_top500_tlx_sign.bed', 
                                    #~ sep='\t', 
                                    #~ index=False,
                                    #~ header=False)

 
#~ plt.figure(11)
#~ plt.hist(enh_sr['tlx_sign2K'], bins=100)


#~ x = np.linspace(-shr, shr, bn)

#~ fig2 = metaseq.plotutils.imshow(
    #~ enh_tlx_arr,
    #~ x=x,
    #~ figsize=(5, 8),
    #~ vmin=5, vmax=99,  percentile=True,
    #~ line_kwargs=dict(color='k', label='All'),
    #~ fill_kwargs=dict(color='k', alpha=0.3),
    #~ sort_by=enh_tlx_arr.mean(axis=1)
#~ )

#~ bottom_axes = plt.subplot(fig2.gs[1, 0])
#~ fig2.line_axes.set_xticklabels([])
#~ fig2.array_axes.set_xticklabels([])
#~ fig2.line_axes.set_ylabel('Average\nenrichement')
#~ fig2.array_axes.set_ylabel('Transcripts on all chromosomes')
#~ fig2.cax.set_ylabel('Enrichment')


## === Gene lists manipulation (tables comes from GREAT)
#~ enh_tlx_gt = pd.read_table('gene_lists/enh_chromm_top500_tlx_sign-genes.txt', header=1, names=['Genes','Regions'])

#~ enh_tlx_gn = list(enh_tlx_gt['Genes'].str.upper())

#~ # === Load expression table 
#~ tbl = pd.read_table(join('tracks', 'TLX3vsRAG-results_genes.txt'), index_col=0)


# Filter genes (Note: this filter remove microRNA expression)
#~ tbl = tbl[(tbl.padj < 0.05)].dropna()


# === Load gene names 
#~ names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           #~ index_col=0,
                           #~ header=0,
                           #~ names=['GeneID', 'TransID', 'Gene_name'])


#~ names = names.drop('TransID', axis=1).drop_duplicates()
#~ names = names.loc[tbl.index]
#~ assert names.shape[0] == tbl.shape[0]

#~ tbl=names.join(tbl, how ='right')


#~ tbn = tbl[['Gene_name', 'R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P', 'padj']]


## === Expresion analysis
#~ classes = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']

#~ import RNA_expression_processing as rn

#~ rn.rna_volcano(tbl, 8, ttl='ALL genes')

#~ rn.rna_scatter(tbl, 8, ttl='ALL genes')
#~ rn.rna_expression(tbl, 50, ttl='ALL genes')


#gn = rn.read_gmt(name = 'TLX3_peaks_RUNX_ETS_genes.gmt', path='tracks') 

#~ topN = rn.express(tbn, 'TLX3', 'RAG', 
                    #~ classes=classes, 
                    #~ n_top=90, 
                    #~ geneList=enh_tlx_gn,  
                    #~ ttl='Genes assoc. top TLX3 activated enhancers')


#~ rn.scatter(tbn, 'TLX3', 'RAG', 
            #~ classes=classes, 
            #~ n_top=5, 
            #~ geneList=enh_tlx_gn,  
            #~ ttl='Genes assoc. top TLX3 activated enhancers')

#~ plt.show()







#~ # === Load table 
#~ tbl = pd.read_table(join('tracks', 'TLX3vsRAGvsTAP_DESeq2-results.txt'), index_col=0)

#~ tbl = tbl[(tbl.padj < 0.05)].dropna()

#~ # === Load gene names 
#~ names = pd.read_table("tracks/UCSC_mm9_transcripID_to_geneSymbol.sort.txt", 
                           #~ index_col=0,
                           #~ names=["Geneid", "NAME"])


#~ names = names.loc[tbl.index]
#~ assert names.shape[0] == tbl.shape[0]


#~ tbl_raw = tbl[['R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ',
               #~ 'TLX3.1_1','TLX3.1_5','TLX3.1_P',
               #~ 'TAP','TAP1B','TAP2B']]


#~ tbl_n=names.join(tbl_raw, how ='right')

#~ tbl_n['NAME']=tbl_n['NAME'].str.upper()







# === Run GSEA
#~ tbl_c = tbl_n.copy() 

#~ tbl_c.index=tbl_n['NAME']

#gnc = list(gen_tb['Genes'].str.upper())
#tbl_c = tbl_c.loc[gen_tb['Genes'].str.upper()].dropna()
#tbl_cc=tbl_cc.dropna()

#~ tbl_c = tbl_c.groupby(tbl_c.index).agg({'NAME': 'first',
                                #~ 'R2.RAG1W.RAG1':sum,
                                #~ 'RAGS.RAGZ':sum,
                                #~ 'RAGZ':sum,
                                #~ 'TLX3.1_1':sum,
                                #~ 'TLX3.1_5':sum,
                                #~ 'TLX3.1_P':sum,
                                #~ 'TAP':sum,
                                #~ 'TAP1B':sum,
                                #~ 'TAP2B':sum})

#~ tbl_c = tbl_c[['NAME',
            #~ 'R2.RAG1W.RAG1',
            #~ 'RAGS.RAGZ',
            #~ 'RAGZ',
            #~ 'TLX3.1_1',
            #~ 'TLX3.1_5',
            #~ 'TLX3.1_P']]
            #~ 'TAP',
            #~ 'TAP1B',
            #~ 'TAP2B']]

#~ classi = ['RAG','RAG','RAG','TLX3','TLX3','TLX3'] #,'TLX3','TLX3','TLX3']

#~ gs_res = gp.gsea.call(data=tbl_c, 
                        #~ gene_sets= 'tracks/TLX_TSS3Kb_ChHMM.gmt', #gene_sets='KEGG_2016', 
                        #~ cls=classi,
                        #~ max_size = 2000,
                        #~ permutation_type='gene_set',#~ permutation_type='phenotype',
                        #~ outdir='gsea_TLX_TSS3Kb_ChHMM')

# === Pictures

#~ gsea_results = gs_res.reset_index().sort_values('fdr',axis=0,ascending=True)

#~ with plt.style.context('ggplot'):
    #~ gsea_results.head(40).plot.bar(y='fdr',x='Term', figsize=(12, 6),fontsize=12)

#~ plt.show() 




#------------------------------------------------------------------------

#~ genes = test1.closest('tracks/genes.bed')
#~ genes.saveas('tracks/TLX3_peaks_RUNX_ETS_genes.bed')



#~ gs = genes.to_dataframe()
#~ nm = list(gs['thickStart'].str.upper())
#~ nm = list(set(nm))
#~ gmt = {'TLX-RUNX-ETS':nm}


# -----------------------------------------------------------------
#~ b1 = [1,2,3,4,5,9,11,15]
#~ b2 = [4,5,6,7,8]
#~ b3 = [val for val in b1 if val in b2]
#~ or

#~ def intersect(a, b):
    #~ return list(set(a) & set(b))

#~ print intersect(b1, b2)

