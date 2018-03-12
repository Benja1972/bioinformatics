import pybedtools as pb
import os 
import subprocess

## ======== Extract specific states  from .bed
# E16 -- Promoter(A)

#~ bashTLX3 = "awk '($4==\"E16\")' tracks/state17_PolII_badRAG/TLX3_17_segments.bed | sort-bed - > tracks/state17_PolII_badRAG/TLX3_st16.bed"
#~ bashRAG = "awk '($4==\"E16\")' tracks/state17_PolII_badRAG/RAG_17_segments.bed | sort-bed - > tracks/state17_PolII_badRAG/RAG_st16.bed"
#~ bashTAP = "awk '($4==\"E16\")' tracks/state17_PolII_badRAG/TAP_17_segments.bed | sort-bed - > tracks/state17_PolII_badRAG/TAP_st16.bed"



#~ print(bashTLX3)
#~ print(bashRAG)
#~ print(bashTAP)


#~ subprocess.call(['bash','-c', bashTLX3])
#~ subprocess.call(['bash','-c', bashRAG])
#~ subprocess.call(['bash','-c', bashTAP])


## ======== Merge beds for 2 bins (200bp each)

#~ mTLX3 = "mergeBed -d 500 -i tracks/state17_PolII_badRAG/TLX3_st16.bed > tracks/state17_PolII_badRAG/TLX3_st16m.bed"
#~ mRAG = "mergeBed -d 500 -i tracks/state17_PolII_badRAG/RAG_st16.bed > tracks/state17_PolII_badRAG/RAG_st16m.bed"
#~ mTAP = "mergeBed -d 500 -i tracks/state17_PolII_badRAG/TAP_st16.bed > tracks/state17_PolII_badRAG/TAP_st16m.bed"

#~ print(mTLX3)
#~ subprocess.call(['bash','-c', mTLX3])

#~ print(mRAG)
#~ subprocess.call(['bash','-c', mRAG])

#~ print(mTAP)
#~ subprocess.call(['bash','-c', mTAP])


## ======= Tracks manipulations  

tlx = pb.BedTool('tracks/state17_PolII_badRAG/TLX3_st16m.bed')
rag = pb.BedTool('tracks/state17_PolII_badRAG/RAG_st16m.bed')
tap = pb.BedTool('tracks/state17_PolII_badRAG/TAP_st16m.bed')
tss2kb = pb.BedTool('tracks/RefSeqTSS2kb.mm9.bed')



# in st6 AND in TSS2Kb

tlx2Kb = (tss2kb + tlx).sort().merge() 
rag2Kb = (tss2kb + rag).sort().merge()  #rag.intersect(tss2kb, wa=True)
tap2Kb = (tss2kb + tap).sort().merge() 

#~ tlx2Kb.saveas('tracks/state17_PolII_badRAG/TLX_st16mTSS2k.bed')
#~ rag2Kb.saveas('tracks/state17_PolII_badRAG/RAG_st16mTSS2k.bed')
#~ tap2Kb.saveas('tracks/state17_PolII_badRAG/TAP_st16mTSS2k.bed')

#~ print(tlx2Kb.count())
#~ print(rag2Kb.count())
#~ print(tap2Kb.count())



tlx_ragtap = (tlx2Kb - rag2Kb - tap2Kb).sort().merge()
rag_tlxtap = (rag2Kb - tlx2Kb - tap2Kb).sort().merge()
tap_tlxrag = (tap2Kb - tlx2Kb - rag2Kb).sort().merge()


print(tlx_ragtap.count())
print(rag_tlxtap.count())
print(tap_tlxrag.count())

tlx_ragtap.saveas('tracks/state17_PolII_badRAG/TLX_st16only.bed')
rag_tlxtap.saveas('tracks/state17_PolII_badRAG/RAG_st16only.bed')
tap_tlxrag.saveas('tracks/state17_PolII_badRAG/TAP_st16only.bed')

#~ st89_peaks = tlx_st89 + tlx_peak
#~ st3_peaks = tlx_st3 + tlx_peak
#st89_peaks2 =  tlx_peak + tlx_st89 

#~ print(st89_peaks.count())
#~ print(st3_peaks.count())


#~ print(tlx_peak.count())

#~ st89_peaks.saveas('st89_peaks.bed')
#~ st3_peaks.saveas('st3_peaks.bed')






#~ st6_tlx =  pb.BedTool('TLX3_14_E6_sorted.bed')
#~ st6_rag =  pb.BedTool('RAG_14_E6_sorted.bed')

#~ # Regions in TLX3 but not in RAG
#~ st6_tlx_not_rag = st6_tlx - st6_rag
#~ st6_rag_not_tlx = st6_rag - st6_tlx

#~ #st6_tlx_not_rag.saveas('state6_tlx_not_rag.bed', trackline="track name='State 6 TLX3 -- not RAG' color=128,0,0")
#~ st6_tlx_not_rag.saveas('state6_tlx_not_rag.bed')

#~ print(st6_tlx_not_rag.count())
#~ print(st6_rag_not_tlx.count())
