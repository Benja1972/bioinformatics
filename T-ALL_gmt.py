import pandas as pd
from os.path import join

## project settings ===
from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')
## =====================

def write_dic2gmt(dic, name, path=''):
    with open(join(path,name+'.gmt'), 'w') as fp:
        for db, term in dic.items():
            gmt = [db, db] + list(term)
            fp.write("\t".join(gmt))
            fp.write("\n")


df = pd.read_table(join(DATADIR,'gene_lists/Cancermine/cancermine_collated.tsv'))

tall = df[df['cancer_normalized']=='acute T cell leukemia']

tgmt = {'T-ALL Oncogenes':list(tall[tall['role']=='Oncogene']['gene_normalized'].unique())}
tgmt['T-ALL Tumor_Suppressor'] = list(tall[tall['role']=='Tumor_Suppressor']['gene_normalized'].unique())
tgmt['T-ALL Driver'] = list(tall[tall['role']=='Driver']['gene_normalized'].unique())
tgmt['T-ALL all'] = list(tall['gene_normalized'].unique())


print(tall.head())

write_dic2gmt(tgmt,'T-ALL',path=join(DATADIR,'gene_lists/Cancermine'))
