from pyfaidx import Fasta, FastaVariant
import pybedtools as pb
from textwrap import TextWrapper
import subprocess


wr = TextWrapper()
wr.width=50


fn = '/home/sergio/media/NAS4/PFlab/TLX3_project/WES-seq/references/mouse_mm9_reference_genome.fa'
vn = 'WES_TLX3.vcf.gz'

fa = Fasta(fn)

bd  = pb.BedTool('test.bed')


inf = 'in_py.fa'
with open(inf,'w') as fp:
    for it in bd:
        rg = fa[it.chrom][it.start:it.end]
        fp.write('>'+rg.longname+'\n'+wr.fill(rg.seq)+'\n')






outf = 'out_py.fa'
cons_fa = "bcftools consensus -f {} {} -o {}".format(inf,vn,outf)

print('Running process ........ \n')
print(cons_fa)
subprocess.call(['bash','-c', cons_fa])

fv = Fasta(outf)



## Only SNP
#~ fv = FastaVariant(fn,vn)
#~ wr.width=60
#~ with open(outf,'w') as fp:
    #~ for it in bd:
        #~ rv = fv[it.chrom][it.start:it.end]
        #~ fp.write('>'+rv.longname+'\n'+wr.fill(rv.seq)+'\n')
