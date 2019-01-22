import gseapy as gp

import pandas as pd

import RNA_expression_processing as rn

import matplotlib.pyplot as plt

path = 'tracks/MARGE/relativeRP/'

dt = pd.read_csv(path+'DN_RegNetwork_TLX3.csv')

gl = list(dt['gene_name'])

gs = 'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO'

rs = gp.enrichr(gene_list=gl, gene_sets=gs)

rn.barplot(rs.res2d,ttl=gs)

plt.show()
