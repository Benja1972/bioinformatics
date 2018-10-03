from os.path import join
from gseapy.parser import gsea_gmt_parser
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt




SAVE = False

### ==== Gene sets formation

gs_dic = {  'Hallmark': 'h.all.v6.0.symbols.gmt',
            #~ 'BioCarta_2013': 'BioCarta_2013',
            #~ 'BioCarta_2015': 'BioCarta_2015',
            #~ 'BioCarta_2016': 'BioCarta_2016',
            #~ 'KEGG_2013': 'KEGG_2013',
            #~ 'KEGG_2015': 'KEGG_2015',
            #~ 'KEGG_2016': 'KEGG_2016',
            #~ 'Reactome_2013': 'Reactome_2013',
            #~ 'Reactome_2015': 'Reactome_2015',
            'Reactome_2016': 'Reactome_2016',
            'Canonical_pathways': 'c2.cp.v6.0.symbols.gmt',
            'Imunno':'c7.all.v6.0.symbols.gmt'}
            #~ 'Oncogenic_signatures': 'MSigDB_Oncogenic_Signatures',
            #~ 'Computational':'MSigDB_Computational',
            #~ 'NCI-60_Cancer_Cell_Lines': 'NCI-60_Cancer_Cell_Lines',
            #~ 'Transcription_factor_targets': 'tracks/GSEA_gene_sets/c3.tft.v6.0.symbols.gmt',
            #~ 'IMMPort': 'gene_lists/IMMPort/IMMPort_test.gmt'}



#~ pl = gsea_gmt_parser('Pierre_gene_set_list.gmt')


#~ pdic = dict()

#~ for db, term in pl.items():
    #~ pl_i = gsea_gmt_parser(gs_dic[db], max_size=5000)
    #~ ii = dict((k,pl_i[k]) for k in term)
    #~ pdic.update(ii)


#~ #/home/sergio/Res_CIML/TLX3_project/Down_stream/gene_lists/TLX3_TLX3_peaks-ALL-genesGREAT.txt
#~ tlx = pd.read_table('../../gene_lists/TLX3_TLX3_peaks-ALL-genesGREAT.txt', 
                    #~ sep='\t', header=1, 
                    #~ names=['Gene_name','peak'])

#~ gl = list(tlx['Gene_name'].str.upper())

#~ pdic.update({'TLX3_TLX3_peaks-ALL':gl})

#~ # == Write gene set dic to gmt format
#~ def write_dic2gmt(dic, name, path=''):
    #~ with open(join(path,name+'.gmt'), 'w') as fp:
        #~ for db, term in dic.items():
            #~ gmt = [db, db] + list(term)
            #~ fp.write("\t".join(gmt))
            #~ fp.write("\n")

#~ write_dic2gmt(pdic,'Pierre_gene_sets_v2')


### === Clustering

def dist_set(set1,set2):
    com = set1 & set2
    dist = 1-len(com)/min(len(set1),len(set2))
    return dist


def dist_kappa(setR, set1,set2):
    a = len(setR & set1 & set2)
    b = len((setR & set1) - set2)
    c = len((setR & set2) - set1)
    d = len((setR - set1) - set2)
    t = a + b + c + d
    p0 = (a + d)/t
    pY = (a + b)*(a + c)/(t*t)
    pN = (c + d)*(b + d)/(t*t)
    pe = pY + pN
    k = (p0 - pe)/(1 - pe)
    #dist = 1-len(com)/min(len(set1),len(set2))
    return k



path = 'GSEA/TLX3vsRAG_Pierre_sets_v2_classic_std/'


fn = 'TLX3vsRAG_Pierre_sets_v2.report.csv'

df = pd.read_csv(path+fn, index_col=0)
df = df.drop(['In Enh or TSS but not common', 'ALL:  TSS or Enh'], axis=0)


terms = list(df.index)

bb = pd.DataFrame(index=terms, columns=terms, dtype='float16')


for ter1 in terms:
    set1 = set(list(df['genes'].loc[ter1].split(',')))
    for ter2 in terms:
        set2 = set(list(df['genes'].loc[ter2].split(',')))
        bb[ter1][ter2] = dist_set(set1,set2)

bb = bb.astype(float)


#f, ax = plt.subplots(figsize=(11, 9))
grid = sns.clustermap(bb, cmap='BuPu_r', figsize=(16, 16)) #cmap='RdBu_r' 'Blues_r''gist_heat' 'ocean' Purples_r

#~ grid.ax_heatmap.set_xticklabels(rotation=30)

for txt in grid.ax_heatmap.get_yticklabels():
        txt.set_rotation(0)
for txt in grid.ax_heatmap.get_xticklabels():
            txt.set_rotation(90)


#~ To access the reordered row indices, use: clustergrid.dendrogram_row.reordered_ind
#~ Column indices, use: clustergrid.dendrogram_col.reordered_ind

#~ llt = list()
#~ for txt in grid.ax_heatmap.get_yticklabels():
        #~ llt.append(txt.get_text())

idd = grid.dendrogram_row.reordered_ind

llt = bb.iloc[idd].index



f, ax = plt.subplots(figsize=(6, 16))
dfff = df.reindex(llt)
dff = dfff.reset_index()
x = dff.index.astype(int)
sz = dff['matched_size'].astype(int)
cc = dff['nes'].astype(float)
ax.scatter(cc, x, s=sz*8, c=cc, 
            cmap="RdBu_r", vmin=-np.max(cc),
            alpha=0.85, edgecolors="grey", 
            linewidth=0.5)
ax.set_xlabel('NES')
ax.set_ylim(-3,95)
ax.set_yticks(np.arange(len(llt)))
ax.set_yticklabels(llt)


#~ plt.show()




import RNA_expression_processing as rn

term_t = 'TLX3_TLX3_peaks-ALL' # 'Tlx3_peaks_sc_gr1000-genes' # 'TSS and Enh'

dfin = df.drop(['TSS and Enh', 'Only TSS', 'Only Enh', term_t], axis=0)
terms = list(dfin.index)
totN = 21000



inter = pd.DataFrame(index=terms, columns=['genes_num', 'p-value', term_t], dtype='str')

set_t = set(list(dfff['genes'].loc[term_t].split(',')))
num_t = len(set_t)

for ter2 in terms:
    set_p = set(list(dfff['genes'].loc[ter2].split(',')))
    num_p = len(set_p)
    if ter2 != term_t:
        com = set_t & set_p
        num_com = len(com)
        inter[term_t][ter2] = list(com)
        inter['genes_num'][ter2] = len(list(com))
        inter['p-value'][ter2] = rn.pValue(totN,num_p,num_t,num_com)
        #~ print(term_t +' and '+ ter2+' = '+str(list(com)))
        #~ print(term_t +' and '+ ter2+' = '+str(len(list(com))))


#~ print(inter[['genes_num', 'p-value']].head(20))

if SAVE:
     inter.sort_values('p-value', ascending=True).to_csv(path+'Tlx3_pk_and_Gene_sets.csv')

inter['p-value'] = inter['p-value'].astype(float)
inter['-log10_pVal'] = -np.log10(inter['p-value'])
inter =inter.sort_values('-log10_pVal', ascending=True).reset_index()
inter.head(40).plot.barh(y='-log10_pVal',x='index', figsize=(18, 6), fontsize=10)
#plt.gca().invert_yaxis()
plt.title(term_t)




#~ inter = inter.sort_values('p-value', ascending=True).reset_index()
#~ inter['genes_num'] = inter['genes_num'].astype(int)
#~ inter['p-value'] = inter['p-value'].astype(float)
#~ inter.head(40).plot.scatter(y='p-value',x='genes_num', figsize=(18, 6), fontsize=10)


#~ f2, ax2 = plt.subplots(figsize=(6, 16))

#~ x = inter.index.astype(int)
#~ sz = inter['p-value'].astype(float)
#~ cc = inter['genes_num'].astype(int)


#~ ax2.scatter(x, sz, s=cc*10, c=sz, 
            #~ cmap="RdBu_r", vmin=-np.max(cc),
            #~ alpha=0.85, edgecolors="grey", 
            #~ linewidth=0.5)



plt.show()



#dff.iloc[::-1].to_csv('../../GSEA/TLX3vsRAG_Pierre_sets_classic_std/TLX3vsRAG_Pierre_cluster.csv')


