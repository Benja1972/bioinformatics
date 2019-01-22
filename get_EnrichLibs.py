from time import sleep
import requests
from os.path import join 
import json

out_dir = 'EnrichrLibs'

lib_url='http://amp.pharm.mssm.edu/Enrichr/datasetStatistics'
libs_json = json.loads(requests.get(lib_url).text)
gss = [lib['libraryName'] for lib in libs_json['statistics']]



link = 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName='


for gs in gss:
    lk = link+gs
    fn = join(out_dir,gs+'.gmt')
    res = requests.get(lk, timeout=None)
    with open(fn, mode='w') as f:
        f.write(res.text)
    sleep(2)
    print(gs)


