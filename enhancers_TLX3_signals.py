#!/usr/bin/env python3

import pybedtools as pb
import numpy as np
import yaml
from os.path import join 
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import metaseq


## === Bed manipulation
enh= pb.BedTool('tracks/TLX3_6_FE_E4_sorted.bed')
enh_df = enh.to_dataframe()
enh_df['len'] = enh_df['end'] - enh_df['start']

enh_df['mid'] = (enh_df['end'] + enh_df['start'])/2
enh_df['mid'] = enh_df['mid'].astype('int')

enh_df['mide'] = enh_df['mid']


plt.hist(enh_df['len'], bins=50)





enh_df[['chrom','mid','mide']].to_csv('tracks/TLX3_6_FE_E4_mid.bed', 
                                    sep='\t', 
                                    index=False,
                                    header=False)


enhm= pb.BedTool('tracks/TLX3_6_FE_E4_mid.bed')

shr = 2000
enhm_Nkb = enhm.slop(b=shr, genome='mm9')#, output='tsses-1kb.gtf')



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

features, arrays = metaseq.persistence.load_features_and_arrays(prefix='enh_tlx1K')


enh_tlx_arr = arrays['enh_tlx']
bn = enh_tlx_arr.shape[1]


enh_df['tlx_sign2K']=enh_tlx_arr.mean(axis=1)

enh_sr = enh_df.sort_values('tlx_sign2K', axis=0, ascending=False)

top1K = enh_sr.head(1000)

top500 = enh_sr.head(500)

top1K[['chrom','start','end']].to_csv('tracks/enh_chromm_top1K_tlx_sign.bed', 
                                    sep='\t', 
                                    index=False,
                                    header=False)

top500[['chrom','start','end']].to_csv('tracks/enh_chromm_top500_tlx_sign.bed', 
                                    sep='\t', 
                                    index=False,
                                    header=False)

 
plt.figure(11)
plt.hist(enh_sr['tlx_sign2K'], bins=100)


x = np.linspace(-shr, shr, bn)

fig2 = metaseq.plotutils.imshow(
    enh_tlx_arr,
    x=x,
    figsize=(5, 8),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=enh_tlx_arr.mean(axis=1)
)

#~ bottom_axes = plt.subplot(fig2.gs[1, 0])
#~ fig2.line_axes.set_xticklabels([])
#~ fig2.array_axes.set_xticklabels([])
fig2.line_axes.set_ylabel('Average\nenrichement')
fig2.array_axes.set_ylabel('Transcripts on all chromosomes')
fig2.cax.set_ylabel('Enrichment')


## === Gene lists manipulation (tables comes from GREAT)
enh_tlx_gt = pd.read_table('gene_lists/enh_chromm_top500_tlx_sign-genes.txt', header=1, names=['Genes','Regions'])

enh_tlx_gn = list(enh_tlx_gt['Genes'].str.upper())

# === Load expression table 
tbl = pd.read_table(join('tracks', 'TLX3vsRAG-results_genes.txt'), index_col=0)


# Filter genes (Note: this filter remove microRNA expression)
tbl = tbl[(tbl.padj < 0.05)].dropna()


# === Load gene names 
names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           index_col=0,
                           header=0,
                           names=['GeneID', 'TransID', 'Gene_name'])


names = names.drop('TransID', axis=1).drop_duplicates()
names = names.loc[tbl.index]
assert names.shape[0] == tbl.shape[0]

tbl=names.join(tbl, how ='right')


tbn = tbl[['Gene_name', 'R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ','TLX3.1_1','TLX3.1_5','TLX3.1_P', 'padj']]


## === Expresion analysis
classes = ['RAG','RAG','RAG','TLX3','TLX3','TLX3']

import RNA_expression_processing as rn

#~ rn.rna_volcano(tbl, 8, ttl='ALL genes')

#~ rn.rna_scatter(tbl, 8, ttl='ALL genes')
#~ rn.rna_expression(tbl, 50, ttl='ALL genes')


#gn = rn.read_gmt(name = 'TLX3_peaks_RUNX_ETS_genes.gmt', path='tracks') 

topN = rn.express(tbn, 'TLX3', 'RAG', 
                    classes=classes, 
                    n_top=90, 
                    geneList=enh_tlx_gn,  
                    ttl='Genes assoc. top TLX3 activated enhancers')


rn.scatter(tbn, 'TLX3', 'RAG', 
            classes=classes, 
            n_top=5, 
            geneList=enh_tlx_gn,  
            ttl='Genes assoc. top TLX3 activated enhancers')

plt.show()







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

