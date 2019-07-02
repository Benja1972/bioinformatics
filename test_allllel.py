import numpy as np
import scipy
import pandas as pd
import pybedtools as pb


#import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')
WGS = join(DATADIR,'tracks/WGS-WES/Germline')


# Load table
import allel


#vcf_tb = allel.vcf_to_dataframe(join(WGS,'Ehn_Active_TLX3_mut.vcf'),fields='*', numbers={'ALT': 1}, transformers=allel.ANNTransformer())
vcf_tb = allel.vcf_to_dataframe(join(WGS,'TLX3_WGS.vcf.gz'),fields='*', numbers={'ALT': 1})

vcf_tb = vcf_tb[vcf_tb['FILTER_PASS']==True]

cols = ['CHROM', 'POS', 'REF', 'ALT','is_snp']

vcf_tb = vcf_tb[cols]
