import numpy as np
from os.path import join 
import sys
import pandas as pd
import pybedtools as pb
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt


# === User defined
import RNA_expression_processing as rn


def log2p1(x):
    return np.log2(x + 1)



SAVE = True
DO_EXPR = False
DO_ENH = False

DO_GSEA = False

ntop=40



# === Load file
prfx = 'bam_input/'

path = 'tracks/MARGE/dRP_PRC2/'

path_out = 'tracks/MARGE/dRP_PRC2/result/'

path_tb = join(path,prfx)


df  =  pd.read_table(path_tb+'TLX_TAP_dRP_mm10mm9.txt')

# -- transform
dfs = df.sort_values('lgFC_TAPvsTLX', axis=0, ascending=False)
dfs.drop_duplicates(subset='gene_name', inplace=True)

# dRP Diffirential analysis

A = 'TLX_raw_RP'
if SAVE:
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(path_out+'dRegPoten_TLX3_vs_TAP_analysis.pdf')

for B in ['TAP_raw_RP']: # ['RAG_rel_RP', 'TAP_rel_RP']:

    cols = ['gene_name', A, B]

    dfp = dfs[cols]
    dfsc = dfp.set_index(keys=dfp.columns[0])

    # classes
    Ac = A+'c'
    Bc = B+'c'
    classe = [A+'c', B+'c']

    up, dn, ax = rn.scatter_n(dfsc, Ac, Bc, 
                classes=classe, 
                n_top=4)
    
    if SAVE:
        plt.savefig(pp, format='pdf')
    
    #~ top, up, dn, gs = rn.express(dfp, Ac, Bc, classes=classe, 
                                    #~ n_top=ntop, geneList=[],  
                                    #~ ttl=A+'/'+B)
    #top = pd.concat([up,dn], axis=0)
    
    rn.cluster(dfsc, Ac, Bc, 
                classes=classe, 
                n_top=ntop)
    
    if SAVE:
        #top.to_csv(path+'TOP'+'_'+A+'_vs_'+B+'.csv')
        up.to_csv(path_out+'UP'+'_'+A+'_vs_'+B+'.csv')
        dn.to_csv(path_out+'DN'+'_'+A+'_vs_'+B+'.csv')
        plt.savefig(pp, format='pdf')



### STOP ###
#~ plt.show()
#~ sys.exit(0)


## === Expression analysis =====
## ============================= 

# === Load expression table 

if DO_EXPR:
    tb_in = pd.read_table(join('tracks', 'TAPvsTLX3-results_genes.txt'), index_col=0)
    names = pd.read_table("tracks/annot_tracks/references/mm9/mm9_EnsemblTransc_GeneNames.txt", 
                           index_col=0,
                           header=0,
                           names=['GeneID', 'TransID', 'Gene_name'])
    names = names.drop('TransID', axis=1).drop_duplicates()
    names = names.reindex(tb_in.index)
    tb_out = pd.concat([names,tb_in],axis=1)
    tb_out.to_csv(join('tracks', 'TAPvsTLX3-results_genesNames.txt'), sep='\t')


tbl = pd.read_table(join('tracks', 'TAPvsTLX3-results_genesNames.txt'), index_col=0)
#tbl = tbl[(tbl.padj < 0.05)].dropna()

# === Pheno ==
A,B = 'TLX3','TAP'
classes = [A]*3+[B]*3


cols = ['Gene_name', 'TLX3.1_1','TLX3.1_5','TLX3.1_P','TAP0','TAP1B','TAP2B']

tbn = tbl[cols]
tbv = tbn.set_index(keys=tbn.columns[0])
#tbv.index=tbv.index.str.upper()
# ---

### == UP analysis
gl_up = list(up.index)
tbu = tbv.iloc[tbv.index.isin(gl_up)]

# -- Scatter
rn.scatter_n(tbu, A, B, classes, n_top=3, ttl='UP_dRP') 
#rn.volcano_n(tbu, A, B, classes, n_top=3, ttl='UP_dRP') 

if SAVE:
    plt.savefig(pp, format='pdf')

# -- Cluster
gr = rn.cluster(tbu, A, B, classes, n_top=ntop)
gr.ax_heatmap.set_title('Cluster UP_dRP '+A+'/'+B)

if SAVE:
    plt.savefig(pp, format='pdf')
# ---


### == DN analysis
gl_dn = list(dn.index)
tbn = tbv.iloc[tbv.index.isin(gl_dn)]

# -- Scatter
rn.scatter_n(tbn, A, B, classes, n_top=4, ttl='DN_dRP') 
#rn.volcano_n(tbn, A, B, classes, n_top=3, ttl='DN_dRP')

if SAVE:
    plt.savefig(pp, format='pdf')

# -- Cluster
gr = rn.cluster(tbn, A, B, classes, n_top=ntop)
gr.ax_heatmap.set_title('Cluster DN_dRP '+A+'/'+B)

if SAVE:
    plt.savefig(pp, format='pdf')

### =================

#~ if SAVE:
    #~ pp.close()



### STOP ###
#~ plt.show()
#~ sys.exit(0)




## === TLX3 peaks analysis =====
## ============================= 

## Load TLX peaks file
colm9 = ['chr_mm9','start_mm9','end_mm9', 'gene_name']
tlx_peak = pb.BedTool('tracks/TLX3_TLX3_peaks.bed')
sl = 100
tlx_peak = tlx_peak.slop(b=sl, genome='mm9')

### == ALL genes with TLX3
tss =  3000

all_tss = dfs[dfs['end_mm9']-dfs['start_mm9']>0]
all_tss = pb.BedTool.from_dataframe(all_tss[colm9])
all_tss = all_tss.slop(b=tss, genome='mm9')

all_tss_tlx = all_tss+tlx_peak
all_tss_gene = all_tss_tlx.to_dataframe()

allTLX_list =list(all_tss_gene['name'])
allTLX_list.sort()


### == UP analysis
up_rp = dfs.loc[dfs['gene_name'].isin(list(up.index))]
up_rp = up_rp[up_rp['end_mm9']-up_rp['start_mm9']>0]
up_rp = pb.BedTool.from_dataframe(up_rp[colm9])
up_rp = up_rp.slop(b=tss, genome='mm9')


up_rp_tlx = up_rp+tlx_peak
up_rp_gene = up_rp_tlx.to_dataframe()

upTLX_list =list(up_rp_gene['name'])
upTLX_list.sort()

print('UP_dRP_with_TLX = ', upTLX_list)

if SAVE:
    up_rp_gene.to_csv(path_out+'UP_TLX3_peaks'+'_'+A+'_vs_'+B+'.csv')

### =================


### == DN analysis
dn_rp = dfs.loc[dfs['gene_name'].isin(list(dn.index))]
dn_rp = dn_rp[dn_rp['end_mm9']-dn_rp['start_mm9']>0]
dn_rp = pb.BedTool.from_dataframe(dn_rp[colm9])
dn_rp = dn_rp.slop(b=tss, genome='mm9')


dn_rp_tlx = dn_rp+tlx_peak
dn_rp_gene = dn_rp_tlx.to_dataframe()

dnTLX_list =list(dn_rp_gene['name'])
dnTLX_list.sort()

print('DN_dRP_with_TLX = ', dnTLX_list)

if SAVE:
    dn_rp_gene.to_csv(path_out+'DN_TLX3_peaks'+'_'+A+'_vs_'+B+'.csv')

### =================






## === Enhancers analysis ======
## ============================= 
flt=4000

if DO_ENH:
    enhm = pb.BedTool('tracks/TLX3_6_FE_E4_sorted.bed')
    enhm_df = enhm.to_dataframe()
    enhm_df.drop('name', axis=1, inplace=True)


    enhm_df['name'] = enhm_df.index
    enhm_df['name'] = 'enh_'+enhm_df['name'].astype(str)

    # Filter if needed
    plt.hist(enhm_df['end']-enhm_df['start'], bins=100)

    enhm_df = enhm_df[enhm_df['end']-enhm_df['start']<flt]

    enhm =  pb.BedTool.from_dataframe(enhm_df)
    if SAVE:
        enhm.saveas('tracks/Enhancers_ChromHMM_flt'+str(flt)+'.bed')

# Load enhancers
enhm = pb.BedTool('tracks/Enhancers_ChromHMM_flt'+str(flt)+'.bed')
enhm_df = enhm.to_dataframe()

# Load gene-enhancers table
en_g = pd.read_table('tracks/Enhancers_ChromHMM_genes2enh.txt', 
                    header=1, 
                    names=['gene_name','enhancers'])
#~ en_g['gene_name'] = en_g['gene_name'].str.upper()


### == UP analysis
up_gen_enh = list(set(list(en_g['gene_name']))&set(gl_up))

up_enh = pd.DataFrame(columns=['chrom', 'start', 'end', 'name'])

for gen in up_gen_enh:

    reg = en_g[en_g['gene_name']==gen]
    enhl = reg.iloc[0]['enhancers'].split(' ')[::2]

    ehGene = enhm_df.loc[enhm_df['name'].isin(enhl)]
    ehGene['name'] = ehGene['name']+'_'+gen
    up_enh = pd.concat([ehGene,up_enh], ignore_index=True)

up_enh_bed = pb.BedTool.from_dataframe(up_enh)
up_enh_bed_tlx = up_enh_bed+tlx_peak

df_up = up_enh_bed_tlx.to_dataframe()

up_enh_tlx_genes = [nm.split('_')[2] for nm in df_up['name']]
up_enh_tlx_genes = list(set(up_enh_tlx_genes))

if SAVE:
    up_enh_bed.saveas(path_out+'UP_dRP_enhancers.bed')
    up_enh_bed_tlx.saveas(path_out+'UP_dRP_enh_TLX3_pck.bed')
    with open(path_out+'UP_dRP_enhancers_genes.txt', 'w') as fp:
        fp.write("\n".join(up_gen_enh))
    with open(path_out+'UP_dRP_enh_TLX3_pck_genes.txt', 'w') as fp:
        fp.write("\n".join(up_enh_tlx_genes))


### =================


### == DN analysis
dn_gen_enh = list(set(list(en_g['gene_name']))&set(gl_dn))

dn_enh = pd.DataFrame(columns=['chrom', 'start', 'end', 'name'])

for gen in dn_gen_enh:

    reg = en_g[en_g['gene_name']==gen]
    enhl = reg.iloc[0]['enhancers'].split(' ')[::2]

    ehGene = enhm_df.loc[enhm_df['name'].isin(enhl)]
    ehGene['name'] = ehGene['name']+'_'+gen
    dn_enh = pd.concat([ehGene,dn_enh], ignore_index=True)

dn_enh_bed = pb.BedTool.from_dataframe(dn_enh)
dn_enh_bed_tlx = dn_enh_bed+tlx_peak

df_dn = dn_enh_bed_tlx.to_dataframe()

dn_enh_tlx_genes = [nm.split('_')[2] for nm in df_dn['name']]
dn_enh_tlx_genes = list(set(dn_enh_tlx_genes))


if SAVE:
    dn_enh_bed.saveas(path_out+'DN_dRP_enhancers.bed')
    dn_enh_bed_tlx.saveas(path_out+'DN_dRP_enh_TLX3_pck.bed')
    with open(path_out+'DN_dRP_enhancers_genes.txt', 'w') as fp:
        fp.write("\n".join(dn_gen_enh))
    with open(path_out+'DN_dRP_enh_TLX3_pck_genes.txt', 'w') as fp:
        fp.write("\n".join(dn_enh_tlx_genes))


### =================


### STOP ###
#~ plt.show()
#~ sys.exit(0)




##  == RegNetworks analysis ==
### ==========================


tf = pd.read_csv(path+'RegNetworkDB/RegNetworkDB_Aug_2018.csv')

# Drop micRNA
tf = tf[~tf['regulator_symbol'].str.contains('mmu-')]

tf['regulator_symbol'] = tf['regulator_symbol'].str.upper()
tf['target_symbol'] = tf['target_symbol'].str.upper()


allTLX_Lst = [x.upper() for x in allTLX_list]

# == UP analysis
tf_up = tf[tf['target_symbol'].isin(gl_up)]
tf_up['gene_name'] =  tf_up['regulator_symbol']

upTLX_Lst = [x.upper() for x in upTLX_list]

tf_up['TLX3-target'] = [1 if x in upTLX_Lst else 0 for x in tf_up['target_symbol']]
tf_up['TLX3-regulator'] = [1 if x in allTLX_Lst else 0 for x in tf_up['regulator_symbol']]

tf_up_tlx = tf_up[tf_up['target_symbol'].isin(upTLX_Lst)]
tf_up_not_tlx = tf_up[~tf_up['target_symbol'].isin(upTLX_Lst)]
tf_up_not_tlx = tf_up_not_tlx[tf_up_not_tlx['regulator_symbol'].isin(allTLX_Lst)]


if SAVE:
    tf_up.to_csv(path_out+'UP_RegNetwork.csv')
    tf_up_tlx.to_csv(path_out+'UP_RegNetwork_TLX3.csv')
    tf_up_not_tlx.to_csv(path_out+'UP_RegNetwork_notTLX3.csv')
### =================

# == DN analysis
tf_dn = tf[tf['target_symbol'].isin(gl_dn)]
tf_dn['gene_name'] =  tf_dn['regulator_symbol']

dnTLX_Lst = [x.upper() for x in dnTLX_list]

tf_dn['TLX3-target'] = [1 if x in dnTLX_Lst else 0 for x in tf_dn['target_symbol']]
tf_dn['TLX3-regulator'] = [1 if x in allTLX_Lst else 0 for x in tf_dn['regulator_symbol']]

tf_dn_tlx = tf_dn[tf_dn['target_symbol'].isin(dnTLX_Lst)]
tf_dn_not_tlx = tf_dn[~tf_dn['target_symbol'].isin(dnTLX_Lst)]
tf_dn_not_tlx = tf_dn_not_tlx[tf_dn_not_tlx['regulator_symbol'].isin(allTLX_Lst)]

if SAVE:
    tf_dn.to_csv(path_out+'DN_RegNetwork.csv')
    tf_dn_tlx.to_csv(path_out+'DN_RegNetwork_TLX3.csv')
    tf_dn_not_tlx.to_csv(path_out+'DN_RegNetwork_notTLX3.csv')
### =================



### ==== EnrichRlib =========
### =========================


import EnrichRLib as erl
gss = [
    #'BioCarta_2016',
    'Cancer_Cell_Line_Encyclopedia',
    #'ChEA_2016',
    'Disease_Perturbations_from_GEO_down',
    'Disease_Perturbations_from_GEO_up',
    #'Disease_Signatures_from_GEO_down_2014',
    #'Disease_Signatures_from_GEO_up_2014',
    'GO_Biological_Process_2018',
    'GO_Cellular_Component_2018',
    'GO_Molecular_Function_2018',
    #'Human_Phenotype_Ontology',
    'KEGG_2016',
    'MSigDB_Computational',
    'MSigDB_Oncogenic_Signatures',
    #'Mouse_Gene_Atlas',
    'NCI-60_Cancer_Cell_Lines',
    #'NCI-Nature_2016',
    #'OMIM_Disease',
    'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO',
    'Reactome_2016',
    #'WikiPathways_2016'
    ]


# == UP analysis

upTLX_list = [x.upper() for x in upTLX_list]

enrUP_tlx = erl.enrich_gs(upTLX_list,gss)

# --- Plot ---
enrUP_tlx.sort_values('p-Val', axis=0, inplace = True)
ds = enrUP_tlx.head(20)

f, ax = plt.subplots()
sns.barplot(y=ds.index,
            x='-log10(p-Val)',
            ax = ax, 
            color="Red", 
            data = ds)
ax.set_title('UP_dRP_Tlx_peaks')

if SAVE:
    plt.savefig(pp, format='pdf')


# == DN analysis

dnTLX_list = [x.upper() for x in dnTLX_list]

enrDN_tlx = erl.enrich_gs(dnTLX_list,gss)

# --- Plot ---
enrDN_tlx.sort_values('p-Val', axis=0, inplace = True)
ds = enrDN_tlx.head(20)

f, ax = plt.subplots()
sns.barplot(y=ds.index,
            x='-log10(p-Val)',
            ax = ax, 
            color="Red", 
            data = ds)
ax.set_title('DN_dRP_Tlx_peaks')





if SAVE:
    plt.savefig(pp, format='pdf')

### =================

if SAVE:
    pp.close()

plt.show()


### ==== GSEA Enrichr =================
### ===================================

if DO_GSEA:
    ### === Enrichr ======
    ### ==================

    gss = [
        'BioCarta_2016',
        'Cancer_Cell_Line_Encyclopedia',
        'ChEA_2016',
        'Disease_Perturbations_from_GEO_down',
        'Disease_Perturbations_from_GEO_up',
        'Disease_Signatures_from_GEO_down_2014',
        'Disease_Signatures_from_GEO_up_2014',
        'GO_Biological_Process_2018',
        'GO_Cellular_Component_2018',
        'GO_Molecular_Function_2018',
        'Human_Phenotype_Ontology',
        'KEGG_2016',
        'MSigDB_Computational',
        'MSigDB_Oncogenic_Signatures',
        'Mouse_Gene_Atlas',
        'NCI-60_Cancer_Cell_Lines',
        'NCI-Nature_2016',
        #'OMIM_Disease',
        'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO',
        #'Reactome_2016',
        'WikiPathways_2016'
        ]

    from matplotlib.backends.backend_pdf import PdfPages

    ### == UP analysis
    cutoff = 0.05
    
    pgs = PdfPages(path+'GSEA_UP_PR_analysis.pdf')

    for gs in gss:
        rs = gp.enrichr(gene_list=gl_up, 
                        description='up_'+gs,
                        outdir=path+'UP_Enrichr', 
                        gene_sets=gs, 
                        no_plot = True)
        if len(rs.res2d[rs.res2d['Adjusted P-value'] <= cutoff]) < 1:
            print("Warning: No enrich terms when cuttoff = %s"%cutoff)
        else:
            gf = rn.barplot(rs.res2d, figsize=(6, 6), ttl=gs)
            #gf.axes[0].set_title(gs)
            plt.savefig(pgs, format='pdf')


    pgs.close()

    ### =================


    ### == DN analysis
    pgs = PdfPages(path+'GSEA_DN_PR_analysis.pdf')

    for gs in gss:
        rs = gp.enrichr(gene_list=gl_dn, 
                        description='dn_'+gs,
                        outdir=path+'DN_Enrichr', 
                        gene_sets=gs, 
                        no_plot = True)
        if len(rs.res2d[rs.res2d['Adjusted P-value'] <= cutoff]) < 1:
            print("Warning: No enrich terms when cuttoff = %s"%cutoff)
        else:
            gf = rn.barplot(rs.res2d, figsize=(6, 6), ttl=gs)
            #gf.axes[0].set_title(gs)
            plt.savefig(pgs, format='pdf')


    pgs.close()










