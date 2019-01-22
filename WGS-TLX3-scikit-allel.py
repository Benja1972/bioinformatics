
# coding: utf-8

# In[1]:

import numpy as np
import scipy
import pandas
import matplotlib as mpl
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')
import h5py
import allel; print('scikit-allel', allel.__version__)


# In[2]:

### Functions
def chrom2num(st):
    chrm = st.split(':')[0]
    pos = st.split(':')[1].split('-')

    pl = int(pos[0].replace(',',''))

    pr = int(pos[1].replace(',',''))
    
    return chrm, pl, pr

def plot_variant_density(pos, window_size, title=None):
    
    # setup windows 
    bins = np.arange(pos.min(), pos.max(), window_size)
    
    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2
    
    # compute variant density in each window
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size
    
    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)


# In[3]:

### VCF direct
ftlx = 'tracks/WGS/Germline/FERRIER_09_Germline.allchr.snpEff.p.SAL.SAL10_1.vcf'
#ftlxg = 'tracks/WGS/Germline/FERRIER_09_Germline.allchr.snpEff.p.SAL.SAL10_1.vcf.gz'

# read VCF file, transform SNPEFF to separated fields (optional)
cs = allel.read_vcf(ftlx,fields='*', numbers={'ALT': 4},transformers=allel.ANNTransformer())


# In[4]:

# variants data to DataFrame, transform SNPEFF to separated fields (optional)
var = allel.vcf_to_dataframe(ftlx,fields='*', numbers={'ALT': 4}, transformers=allel.ANNTransformer())


# In[5]:

# Genotype array to special class GenotypeArray
gt = allel.GenotypeArray(cs['calldata/GT'])
    ##- typical functions
    # gt.is_het()
    # gt.count_het(axis=1)
    # ac = gt.count_alleles()


# In[33]:

#gt
#sorted(var.columns)


# In[7]:

print(sorted(var['ANN_Feature_Type'].unique()))
print(sorted(var['ANN_Transcript_BioType'].unique()))


# In[8]:

cod_var = var[var['ANN_Feature_Type']=='transcript']


# In[9]:

print(len(var))
len(cod_var)


# In[10]:

cod_ind=cod_var.index
cod_gt=gt[cod_ind]


# ## Now we have pair {cod_var, cod_gt} for transcripts only

# In[11]:

#cod_var.head(12)
a,b,c = plt.hist(np.log(cod_var['QUAL']), bins=100)


# In[12]:

#cod_gt[2,:]


# ## Strip var data to region 

# In[14]:

st  ='chr12:77,033,211-78,041,433'
c,l,r = chrom2num(st)
print(c,l,r)


# In[15]:

var_reg = var[(var['CHROM']==c) & (var['POS']>l) & (var['POS']<r)]


# In[16]:

plot_variant_density(var_reg['POS'], window_size=35, title=c)


# In[ ]:

### Plot density for all chromosomes
# for ch in df['CHROM'].unique():
#     dfc = df[df['CHROM']==ch]
#    plot_windowed_variant_density(dfc['POS'], window_size=100000, title=ch+' , raw variant density')


# # Working with chunked table, we need HDF5 file

# In[17]:

### Save to hdf5
#import sys
#allel.vcf_to_hdf5(ftlx,'FERRIER_09_Germline.allchr.snpEff.p.SAL.SAL10_1_Shrt.h5', 
#                  fields='*', alt_number=4,transformers=allel.ANNTransformer(),log=sys.stdout, vlen=False)


# In[18]:

### HDF5 from VCF database
ftlxh5 ='tracks/WGS/Germline/FERRIER_09_Germline.allchr.snpEff.p.SAL.SAL10_1_Shrt.h5'

# read HDF5 file
csh = h5py.File(ftlxh5,mode='r')
var_tb = allel.VariantChunkedTable(csh['variants'], 
                                   names=['CHROM', 'POS', 'REF', 'ALT', 'DP', 'MQ', 'QD', 'ANN_AA_length',
                                             'ANN_Allele',
                                             'ANN_Annotation',
                                             'ANN_Annotation_Impact',
                                             'ANN_Feature_ID',
                                             'ANN_Feature_Type',
                                             'ANN_Gene_ID',
                                             'ANN_Gene_Name',
                                             'ANN_Rank',
                                             'ANN_Transcript_BioType','numalt'])


# In[19]:

#a,b,c=plt.hist(var_tb['DP'][:], bins=10)
#csh['variants/REF']


# ## Now we can work with filters

# In[20]:

#fltr_expr = '(QD > 5) & (MQ > 40) & (DP > 1500) & (DP < 3000)'
fltr_expr="ANN_Feature_Type==b'transcript'"

var_tb_fltr = var_tb.eval(fltr_expr)[:]

#var_tb
#var_tb_fltr
np.count_nonzero(var_tb_fltr)
#np.count_nonzero(~var_tb_fltr)

#list(csh['calldata'].keys())
#list(csh['variants'].keys())


# In[21]:

## apply filter
var_pass = var_tb.compress(var_tb_fltr)


# ## Genotype from HDF5 

# In[22]:

list(csh['calldata'].keys())


# In[23]:

gth = allel.GenotypeChunkedArray(csh['calldata/GT'])
gth


# In[24]:

list(csh['samples'])


# In[25]:

import pandas as pd
samples = pd.DataFrame({'sample':[b'AC3812', b'AC3813', b'AC3814', b'AC3815'], 'cell_type':['TAP','TAP','TLX3','TLX3']})
TLX = samples['cell_type'].isin(['TLX3'])
TAP = samples['cell_type'].isin(['TAP'])


# ## Subset genotype on transcrips and samples

# In[26]:

gth_tlx = gth.subset(var_tb_fltr, TLX)
gth_tap = gth.subset(var_tb_fltr, TAP)


# #### Now we have three tables: {var_pass, gth_tlx, gth_tap}  for transcripts only 

# In[27]:

n_variants = len(var_pass)
pc_missing_tlx = gth_tlx.count_missing(axis=0)[:] * 100 / n_variants
pc_het_tlx = gth_tlx.count_het(axis=0)[:] * 100 / n_variants

pc_missing_tap = gth_tap.count_missing(axis=0)[:] * 100 / n_variants
pc_het_tap = gth_tap.count_het(axis=0)[:] * 100 / n_variants



print('TLX3 missing = ', pc_missing_tlx)
print('TLX3 hetero = ', pc_het_tlx)

print('TAP missing = ', pc_missing_tap)
print('TAP hetero = ', pc_het_tap)


# In[31]:

tlx_seg = gth_tlx.count_alleles().count_segregating()
tap_seg = gth_tap.count_alleles().count_segregating()

print('TLX segregating = ', tlx_seg)
print('TAP segregating = ', tap_seg)


# In[ ]:

def plot_variant_hist_2d(f1, f2, variants, downsample):
    x = variants[f1][:][::downsample]
    y = variants[f2][:][::downsample]
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.despine(ax=ax, offset=10)
    ax.hexbin(x, y, gridsize=20)
    ax.set_xlabel(f1)
    ax.set_ylabel(f2)
    ax.set_title('Variant %s versus %s joint distribution' % (f1, f2))


# In[ ]:

#plot_variant_hist_2d('QD', 'MQ', var, downsample=500)


# In[ ]:

def plot_variant_hist(f, variants, bins=30, down=200):
    x = variants[f][:][::down]
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.despine(ax=ax, offset=10)
    ax.hist(x, bins=bins)
    ax.set_xlabel(f)
    ax.set_ylabel('No. variants')
    ax.set_title('Variant %s distribution' % f)


# In[ ]:

#“MQ” is average mapping quality across all samples.
plot_variant_hist('MQ', var, down=2)


# In[ ]:

#“QD” is a slightly odd statistic but turns out to be very useful 
# for finding poor quality SNPs. Roughly speaking, high numbers 
# mean that evidence for variation is strong (concentrated), 
# low numbers mean that evidence is weak (dilute).


x = var['QD'][:][::1000]
plot_variant_hist('QD', var, bins=30, down=500)


# In[ ]:

ac = gt.count_alleles()
ac


# In[ ]:



