#!/usr/bin/env python

import pandas as pd
import numpy as np
import shutil as sh

prefixRAW = '/home/sergio/media/NAS3B/PFlab/RawData/CHIPseq/'
prefixDATA = '/home/sergio/media/NAS3B/Bioinfoplatform/PFlab/ChiP-Seq/TAP/data/'


# Copy TAP data to data processing directory
df = pd.read_csv('TAP-data.csv')

for i in df.index:
	fileIN1 = df['FileName'][i].replace("$RAW/",prefixRAW)[:-11] + '_1.fastq.gz'
	fileIN2 = df['FileName'][i].replace("$RAW/",prefixRAW)[:-11] + '_2.fastq.gz'

	fileOUT1 = prefixDATA + df['SampleName'][i] + '_1.fastq.gz'
	fileOUT2 = prefixDATA + df['SampleName'][i] + '_2.fastq.gz'
	
	sh.copyfile(fileIN1, fileOUT1)
	sh.copyfile(fileIN2, fileOUT2)


# Copy inputs for TAP to data processing directory
df_inp = pd.read_csv('INPUTS-data.csv')

for i in [3,4]:
	fileIN1 = df_inp['FileName'][i].replace("$RAW/",prefixRAW)[:-11] + '_1.fastq.gz'
	fileIN2 = df_inp['FileName'][i].replace("$RAW/",prefixRAW)[:-11] + '_2.fastq.gz'

	fileOUT1 = prefixDATA + df_inp['SampleName'][i] + '_1.fastq.gz'
	fileOUT2 = prefixDATA + df_inp['SampleName'][i] + '_2.fastq.gz'
	
	sh.copyfile(fileIN1, fileOUT1)
	sh.copyfile(fileIN2, fileOUT2)


