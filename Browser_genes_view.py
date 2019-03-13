import numpy as np
import pandas as pd
import pybedtools as pb

# my libs
import EnrichRLib as erl


# Project settings
from os.path import join 
WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')
WGS = join(DATADIR,'tracks/WGS-WES/Germline')
RP = join(DATADIR,'tracks/MARGE/relativeRP/bam_input')

# Load genes body regions
genes = pb.BedTool(join(DATADIR,'tracks/annot_tracks/references/mm9/mm9.refGene.bed'))
tall = erl.read_gmt(join(DATADIR,'gene_lists/Cancermine/T-ALL.gmt')) 

# Genes of interests

gl =     ['BCL11B',
        'EZH2',
        'RUNX1',
        'FBXW7',
        'FOS',
        'ETV6',
        'RPL5',
        'RB1',
        'GATA3',
        'LEF1',
        'TET1',
        'PTEN',
        'MSH2',
        'MSH6',
        'PMS2',
        'MLH1',
        'PHF6']
