import pandas as pd

fn = 'gene_lists/IMMPort/Geneappend3.xls'

df = pd.read_excel(fn)

with open('gene_lists/IMMPort/IMMPort.gmt', 'w') as fp:
    for ct in df['Category'].unique():
        lg = [ct,ct]+list(df[df['Category']==ct]['Symbol'])
        fp.write("\t".join(lg))
        fp.write('\n')


