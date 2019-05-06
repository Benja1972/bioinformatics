import pybedtools as pb
from matplotlib import pyplot as plt


from os.path import join 

WORKDIR = '/home/sergio/Res_CIML/TLX3_project'
SCRIPTS = join(WORKDIR,'scripts')
DATADIR = join(WORKDIR,'data')

bd = pb.BedTool(join(DATADIR,'tracks/ChromHMM_enhancers/Apply_HMM_FE_model_6_to_RAG/final_model/RAG_6_segments.bed'))
bd_dt = bd.to_dataframe()
bd_4 = bd_dt[bd_dt['name']=='E4']

st4 = pb.BedTool.from_dataframe(bd_4)
st4 = st4.merge(d=2000).sort()

bd_4 = st4.to_dataframe()

# Filter by lenght
bd_4['lenght'] = bd_4['end'] - bd_4['start']
bd_4 = bd_4[bd_4['lenght']>200]



# Add enhancer's name
bd_4['enh_name'] = bd_4.index
bd_4['enh_name'] = 'enh_'+bd_4['enh_name'].astype(str)

col2bed = ['chrom', 'start', 'end', 'enh_name']

# Final bed
st4 =  pb.BedTool.from_dataframe(bd_4[col2bed])
st4.saveas(join(DATADIR,'tracks/Enhancers_RAG_ChromHMM.bed'))




# fast pict
plt.hist(bd_4['lenght'], bins=80)

plt.show()
