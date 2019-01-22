#!/usr/bin/env python2

import numpy as np
from os.path import join 
import pandas as pd
import matplotlib.pyplot as plt
import RNA_expression_processing as rn

SAVE = True
path=''

if SAVE:
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(path+'DEG_TLX3_vs_RAG.pdf')


# === Load expression table 
tbl = pd.read_table(join('tracks', 'TLX3vsRAG-results_genesNames.txt'), index_col=0)
tbl = tbl[(tbl.padj < 0.05)].dropna()

# === Pheno ==
A,B = 'TLX3','RAG'
classes = [A]*3+[B]*3


cols = ['Gene_name', 'TLX3.1_1','TLX3.1_5','TLX3.1_P','R2.RAG1W.RAG1','RAGS.RAGZ','RAGZ']

tbn = tbl[cols]
tbn = tbn.set_index(keys=tbn.columns[0])


# === DEG and filter
#~ cls = tbn.columns
#~ tbn = rn.exp_deg(tbn, A, B, classes)
#~ tbn = tbn[(tbn['p-adj'] < 0.05)].dropna()
#~ tbn = tbn[cls]


# -- Scatter
rn.scatter_n(tbn, A, B, classes, n_top=3, ttl='') 

if SAVE:
    pp.savefig()

rn.volcano_n(tbn, A, B, classes, n_top=3, ttl='')

if SAVE:
    pp.savefig()

gr = rn.cluster(tbn, A, B, classes, n_top=20)
gr.ax_heatmap.set_title(A+'/'+B+ ' diff. expressed' )

if SAVE:
    pp.savefig()





plt.show() 






if SAVE:
    pp.close()















## === HTML parsing example -- TCRalpha locus
#~ df = pd.read_html('http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=GenePositions&species=mouse&group=TRA', header=2)

#~ tra_mm10 = df[0]

#~ tra_gn = list(tra_mm10['MGI symbol'].str.upper())

#~ rn.express(tbn, 'TLX3', 'RAG', 
            #~ classes=classes, 
            #~ #n_top=50, 
            #~ geneList=tra_gn,  
            #~ ttl=r'TCR$\alpha$ locus')



#~ df2 = pd.read_html('http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=GenePositions&species=mouse&group=TRB', header=2)


#~ trb_mm10 = df2[0]


#~ trb_gn = list(trb_mm10['MGI symbol'].str.upper())

#~ rn.express(tbn, 'TLX3', 'RAG', 
            #~ classes=classes, 
            #~ #n_top=50, 
            #~ geneList=trb_gn,  
            #~ ttl=r'TCR$\beta$ locus')

#~ plt.show()



### === Specific GENE ==
#~ cells = ('RAG','TLX3')
#~ cells2 = ('TLX3','TAP')
#~ x_pos = np.arange(len(cells))
#~ x_pos2 = np.arange(len(cells2))

#~ tbl['TLX3-mean'] = log2p1((tbl['TLX3.1_1']+tbl['TLX3.1_5']+tbl['TLX3.1_P'])/3.)
#~ #tbl.data['TAP-mean'] = log2p1((tbl.data['TAP']+tbl.data['TAP1B']+tbl.data['TAP2B'])/3.)
#~ tbl['RAG-mean'] = log2p1((tbl['R2.RAG1W.RAG1']+tbl['RAGS.RAGZ']+tbl['RAGZ'])/3.)

#~ tbl_p['TLX3-mean'] = log2p1((tbl_p['TLX3.1_1']+tbl_p['TLX3.1_5']+tbl_p['TLX3.1_P'])/3.)
#~ tbl_p['TAP-mean'] = log2p1((tbl_p['TAP0']+tbl_p['TAP1B']+tbl_p['TAP2B'])/3.)



#~ with plt.style.context('seaborn-talk'):
    #~ fig2 = plt.figure(figsize=(16, 6))
    #~ ax5 = fig2.add_subplot(131)
    #~ ac =tbl.loc[tbl['Gene_name']=='Arid5b'][['RAG-mean','TLX3-mean']]
    #~ ax5.bar(x_pos, np.array(ac).squeeze(),align='center', color=['green','red'])
    #~ ax5.set_title('Arid5b gene expression')
    #~ ax5.set_ylabel('log2(FPKM + 1)')
    #~ ax5.set_xticks(x_pos)
    #~ ax5.set_xticklabels(cells)
    
    
    #~ ax6 = fig2.add_subplot(132)
    #~ ac =tbl_p.loc[tbl_p['Gene_name']=='Arid5b'][['TLX3-mean','TAP-mean']]
    #~ ax6.bar(x_pos2, np.array(ac).squeeze(),align='center', color=['green','red'])
    #~ ax6.set_title('Arid5b gene expression')
    #~ ax6.set_ylabel('log2(FPKM + 1)')
    #~ ax6.set_xticks(x_pos2)
    #~ ax6.set_xticklabels(cells2)


#~ gt = pd.read_table('gene_lists/genes_list_Article.txt', names=['gene'])

#~ gn =list(gt['gene'])






## === GSEA analysis of all genes
## run in Python3
#tbn['Gene_name']=tbn.index


#~ tbn = tbn.dropna()
#~ tbn = tbn.iloc[:,:len(classes)+1]


#~ tbn.set_index(keys=tbn.columns[0], inplace=True)
#~ tbn.index=tbn.index.str.upper()





#~ tbn['Gene_name']=tbn.index

#~ tb_u = pd.concat([tbn['Gene_name'],tbn.iloc[:,:len(classes)]], axis=1)


#tbn = tbn.reset_index(drop=True)

#~ tbn['Gene_name']=tbn['Gene_name'].str.upper()

#~ tbn.set_index('Gene_name', drop=False, inplace=True)

#~ tb_u = tbn.drop(['padj'],axis=1)

#gs_res = gp.gsea.call(data=tbn,gene_sets = g_set, outdir='_test')



#~ tbn.drop(['padj'],axis=1, inplace=True)

#tb_u = pd.concat([tbn['Gene_name'],tbn.iloc[:,1:len(classes)+1]], axis=1)

#~ g_set = 'tracks/GSEA_gene_sets/c7.all.v6.0.symbols.gmt'
#~ out_dir = 'GSEA/gsea_TLXvsRAG_C7_Immuno'

#~ g_set = 'gene_lists/IMMPort/IMMPort.gmt'
#~ out_dir = 'GSEA/gsea_IMMPort'

#~ gs_res = gp.gsea.call(data=tb_u, 
                        #~ gene_sets = g_set, 
                        #~ cls=classes,
                        #~ max_size = 1000,
                        #~ top_enrich = 25,
                        #~ outdir=out_dir)

#~ gsea_results = gs_res.reset_index().sort_values('fdr',axis=0,ascending=True)

#~ with plt.style.context('ggplot'):
    #~ gsea_results.head(40).plot.bar(y='fdr',x='Term', figsize=(12, 6),fontsize=12)


#~ plt.savefig(out_dir+'/'+'TLX3vsRAG_IMMPort.pdf', format="pdf")









## === GSEA analysis of topN
## run in Python3
#~ upN['Gene_name']=upN.index

#~ tb_u = pd.concat([upN['Gene_name'],upN.iloc[:,:len(classes)]], axis=1)

#~ g_set = 'tracks/GSEA_gene_sets/c7.all.v6.0.symbols.gmt'
#~ out_dir = 'GSEA/gsea_UP500_C7_Immuno'
#~ gs_res = gp.gsea.call(data=tb_u, 
                        #~ gene_sets = g_set, 
                        #~ cls=classes,
                        #~ max_size = 2000,
                        #~ top_enrich = 50,
                        #~ outdir=out_dir)

#~ gsea_results = gs_res.reset_index().sort_values('fdr',axis=0,ascending=True)

#~ with plt.style.context('ggplot'):
    #~ gsea_results.head(40).plot.bar(y='fdr',x='Term', figsize=(12, 6),fontsize=12)


#~ plt.savefig(out_dir+'/'+'UP500_C7_Immuno.pdf', format="pdf")




#~ dnN['Gene_name']=dnN.index

#~ tb_d = pd.concat([dnN['Gene_name'],dnN.iloc[:,:len(classes)]], axis=1)

#~ g_set = 'tracks/GSEA_gene_sets/c7.all.v6.0.symbols.gmt'
#~ out_dir = 'GSEA/gsea_DN500_C7_Immuno'
#~ gs_res = gp.gsea.call(data=tb_d, 
                        #~ gene_sets = g_set, 
                        #~ cls=classes,
                        #~ max_size = 2000,
                        #~ top_enrich = 50,
                        #~ outdir=out_dir)

#~ gsea_results = gs_res.reset_index().sort_values('fdr',axis=0,ascending=True)

#~ with plt.style.context('ggplot'):
    #~ gsea_results.head(40).plot.bar(y='fdr',x='Term', figsize=(12, 6),fontsize=12)


#~ plt.savefig(out_dir+'/'+'DN500_C7_Immuno.pdf', format="pdf")
















#~ ## === Histogram padj
#~ fig21 = plt.figure(figsize=(8, 6))
#~ ax51 = fig21.add_subplot(111)
#~ ax51.hist(tbl['padj'], bins=50)
#~ ax51.set_title('Histogram of P-value')
#~ ax51.set_ylabel('Frequency')
#~ ax51.set_xlabel('P-value')



#~ ## === Gene specific
#~ tbl['TLX3-mean'] = log2p1((tbl['TLX3.1_1']+tbl['TLX3.1_5']+tbl['TLX3.1_P'])/3.)
#~ #tbl.data['TAP-mean'] = log2p1((tbl.data['TAP']+tbl.data['TAP1B']+tbl.data['TAP2B'])/3.)
#~ tbl['RAG-mean'] = log2p1((tbl['R2.RAG1W.RAG1']+tbl['RAGS.RAGZ']+tbl['RAGZ'])/3.)

#~ cells = ('RAG','TLX3')
#~ x_pos = np.arange(len(cells))

#~ with plt.style.context('seaborn-talk'):
    #~ fig2 = plt.figure(figsize=(16, 6))
    #~ ax5 = fig2.add_subplot(131)
    #~ ac =tbl.loc[tbl['Gene_name']=='Trac'][['RAG-mean','TLX3-mean']]
    #~ ax5.bar(x_pos, np.array(ac).squeeze(),align='center', color=['green','red'])
    #~ ax5.set_title('Trac gene expression')
    #~ ax5.set_ylabel('log2(FPKM + 1)')
    #~ ax5.set_xticks(x_pos)
    #~ ax5.set_xticklabels(cells)
    
    #~ ax3 = fig2.add_subplot(132)
    #~ aa =tbl.loc[tbl['Gene_name']=='Trav6-3'][['RAG-mean','TLX3-mean']]
    #~ ax3.bar(x_pos, np.array(aa).squeeze(),align='center', color=['green','red'])
    #~ ax3.set_title('Trav6-3 gene expression')
    #~ ax3.set_ylabel('log2(FPKM + 1)')
    #~ ax3.set_xticks(x_pos)
    #~ ax3.set_xticklabels(cells)


    #~ ax4 = fig2.add_subplot(133)
    #~ ab =tbl.loc[tbl['Gene_name']=='Trav1'][['RAG-mean','TLX3-mean']]
    #~ ax4.bar(x_pos, np.array(ab).squeeze(),align='center', color=['green','red'])
    #~ ax4.set_title('Trav1 gene expression')
    #~ ax4.set_ylabel('log2(FPKM + 1)')
    #~ ax4.set_xticks(x_pos)
    #~ ax4.set_xticklabels(cells)



#~ fig4 = plt.figure(figsize=(6, 6))
#~ ax51 = fig4.add_subplot(111)
#~ ac =tbl.loc[tbl['Gene_name']=='Tlx3'][['RAG-mean','TLX3-mean']]
#~ ax51.bar(x_pos, np.array(ac).squeeze(),align='center', color=['green','red'])
#~ ax51.set_title('TLX3 gene expression')
#~ ax51.set_ylabel('log2(FPKM + 1)')
#~ ax51.set_xticks(x_pos)
#~ ax51.set_xticklabels(cells)







