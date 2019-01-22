#~ import vcf
#~ vcf_tlx = vcf.Reader(open('tracks/WGS/RAW_vcf/TLX3_AC3814_vs_AC3815.allchr.MuTect.snpEff.p.SAL.SAL10_2.vcf', 'r'))

#~ record = next(vcf_tlx)
#~ print(record)

#~ print(record.INFO)


st  ='chr12:77,033,211-78,041,433'

def chrom2num(st):
    chrm = st.split(':')[0]
    pos = st.split(':')[1].split('-')

    pl = int(pos[0].replace(',',''))

    pr = int(pos[1].replace(',',''))
    
    return chrm, pl, pr

print(chrom2num(st))


c,l,r = chrom2num(st)

print(c,l,r)
