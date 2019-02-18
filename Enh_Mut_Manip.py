import numpy as np
import pandas as pd
import pybedtools as pb


def enh_gene(genl,df):
    if isinstance(genl, str):
        genl = [genl]
    enhll = list()
    for gen in genl:
        reg = df[df['gene_name']==gen]
        if len(reg)>0:
            dr = reg.iloc[0]['enhancers'].split(' ')
            enhl = [x for x in dr if 'enh' in x]
        else:
            enhl=[]
        enhll = enhll+enhl
    return list(set(enhll))

def gene_enh(enhl,df):
    if isinstance(enhl, str):
        enhl = [enhl]
    
    genl = list()
    for enh in enhl:
        reg = df[df['enhancers']==enh]
        if len(reg)>0:
            ls = reg.iloc[0]['gene_name'].split(', ')
            gl=[tr.split(' (')[0] for tr in ls]
        else:
            gl=[]
        genl = genl+gl
    return list(set(genl))

def gene_variants(gl,var,genes):
    """
    Function counts variants in gene list 
    
    Parametrs
    ---------
    gl : genes list
    var : variant DataFrame; ex. var = allel.vcf_to_dataframe('fn.vcf')
    genes : bed with genes; ex. genes = pb.BedTool('genes.bed')
    """
    genes = genes.to_dataframe()
    genes.loc[:,'name'] = genes.loc[:,'name'].str.upper()
    genes_l = genes[genes['name'].isin(gl)]
    
    var_gl = pd.DataFrame()
    for idx, row in genes_l.iterrows():
        vg = var[(var['CHROM']==row['chrom']) & (var['POS']>row['start']) & (var['POS']<row['end'])]
        genes_l.loc[idx,'num_mut'] = int(len(vg))
        var_gl = pd.concat([var_gl, vg])
    return  var_gl, genes_l


#~ def bed_variants_counts(var,bed):
    #~ """
    #~ Function counts variants in region bed file
    
    #~ Parametrs
    #~ ---------
    #~ var : variant DataFrame; ex. var = allel.vcf_to_dataframe('fn.vcf')
    #~ bed : bed file with regions; ex. bed = pb.BedTool('genes.bed')
    #~ """
    #~ bed = bed.to_dataframe()
    
    #~ var_bd = pd.DataFrame()
    #~ for idx, row in bed.iterrows():
        #~ vg = var[(var['CHROM']==row['chrom']) & (var['POS']>row['start']) & (var['POS']<row['end'])]
        #~ bed.loc[idx,'counts'] = int(len(vg))
        #~ var_bd = pd.concat([var_bd, vg])
    #~ return  var_bd , bed



def bed_variants(var,bed):
    """
    Function colect variants in region bed file and out to variants file
    
    Parametrs
    ---------
    var : variant bed file; ex. var = pb.BedTool('fn.vcf')
    bed : bed file with regions; ex. bed = pb.BedTool('genes.bed')
    """
    var_out = var.intersect(bed, header=True)
    
    return  var_out


def variants_bed_counts(var,bed):
    """
    Function counts variants in region bed file
    
    Parametrs
    ---------
    var : variant bed file; ex. var = pb.BedTool('fn.vcf')
    bed : bed file with regions; ex. bed = pb.BedTool('genes.bed')

    Output
    ------
    bed_out : bed with last 'counts' column 
    """
    bed_out = bed.intersect(var, c=True)
    #~ bed_df = bed_out.to_dataframe()
    #~ bed_df = bed_df[['chrom', 'start', 'end','thickStart', 'name']]
    #~ bed_out = pb.BedTool.from_dataframe(bed_df)
    return bed_out



def freq_on_gene(var):
    """
    Parametrs
    ---------
    var : variant DataFrame with ANN filed from SnpEff; 
            ex. var = allel.vcf_to_dataframe('fn.vcf', transformers=allel.ANNTransformer())
    """
    fr = var[['ANN_Gene_Name']]
    fr['num_mut'] = 1
    frq = fr.groupby('ANN_Gene_Name').sum()
    frq.sort_values('num_mut', axis=0, inplace=True, ascending=False)
    
    return frq
