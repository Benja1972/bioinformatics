import pybedtools as pb





tlx_st89 = pb.BedTool('TLX3_st8-9.bed')
tlx_peak = pb.BedTool('TLX3_TLX3_peaks.bed')
tlx_st3 = pb.BedTool('TLX3_st3.bed')

# in st89 AND in TLX-peaks
st89_peaks = tlx_st89 + tlx_peak
st3_peaks = tlx_st3 + tlx_peak
#st89_peaks2 =  tlx_peak + tlx_st89 

print(st89_peaks.count())
print(st3_peaks.count())


print(tlx_peak.count())

st89_peaks.saveas('st89_peaks.bed')
st3_peaks.saveas('st3_peaks.bed')






#~ st6_tlx =  pb.BedTool('TLX3_14_E6_sorted.bed')
#~ st6_rag =  pb.BedTool('RAG_14_E6_sorted.bed')

#~ # Regions in TLX3 but not in RAG
#~ st6_tlx_not_rag = st6_tlx - st6_rag
#~ st6_rag_not_tlx = st6_rag - st6_tlx

#~ #st6_tlx_not_rag.saveas('state6_tlx_not_rag.bed', trackline="track name='State 6 TLX3 -- not RAG' color=128,0,0")
#~ st6_tlx_not_rag.saveas('state6_tlx_not_rag.bed')

#~ print(st6_tlx_not_rag.count())
#~ print(st6_rag_not_tlx.count())
